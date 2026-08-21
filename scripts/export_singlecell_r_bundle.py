#!/usr/bin/env python3
"""Export a compact, provenance-rich bundle for R-side single-cell plotting.

The exporter intentionally copies only metadata, embeddings, and selected marker
genes. It never materializes the full expression matrix as dense data.

Schema versions:
  v1 (default): CSV/CSV.gz for all files — backward-compatible, unchanged.
  v2: parquet for tables, optional mtx.gz for large marker matrices; richer manifest.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
from dataclasses import dataclass, field, replace as dataclass_replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import anndata as ad
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "anndata is required. Run inside the singlecell_factory sc_gpu environment."
    ) from exc

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    _HAVE_PYARROW = True
except ImportError:  # pragma: no cover
    _HAVE_PYARROW = False

try:
    import scipy.io as sio
    _HAVE_SCIPY_IO = True
except ImportError:  # pragma: no cover
    _HAVE_SCIPY_IO = False


SCHEMA_VERSION_V1 = "singlecell_r_bundle_v1"
BUNDLE_SCHEMA_V2 = "singlecell_r_bundle_v2"
BUNDLE_SCHEMA_V2_1 = "singlecell_r_bundle_v2.1"
BUNDLE_SCHEMA_V2_2 = "singlecell_r_bundle_v2.2"
# Backwards-compat alias (older callers / tests reference SCHEMA_VERSION_V2).
SCHEMA_VERSION_V2 = BUNDLE_SCHEMA_V2
SCHEMA_COMPATIBLE_WITH_V2 = ["singlecell_r_bundle_v1"]
SCHEMA_COMPATIBLE_WITH_V2_1 = ["singlecell_r_bundle_v1", "singlecell_r_bundle_v2"]
SCHEMA_COMPATIBLE_WITH_V2_2 = ["singlecell_r_bundle_v1", "singlecell_r_bundle_v2", "singlecell_r_bundle_v2.1"]
EXPRESSION_SOURCE_SLOT = "X"
EXPRESSION_VALUE_SCALE = "source_X_as_stored"
EXPRESSION_EXPORT_DTYPE = "float32"
INTENDED_USE = "plotting_and_visual_summary_only"
# H2 / Phase B: This claim_guard string MUST be surfaced by the R reader
# (see r_multiomics_factory/R_bundle/io_bundle.R::read_bundle_v2, which both
# emits a message() and attaches it as attr(expr_sparse, "claim_guard")).
# Do not change the literal value without updating the R-side validator and
# the cross-language parity tests.
#
# Phase B (v2.1) note: by default the SAME guard string is reused at the
# modality level for every optional bundle extension (protein, spatial,
# multimodal_obsm). Each extension is plotting-only until per-modality
# validation logic exists. Extension authors may override `claim_guard` when
# calling `add_extension(...)`, but they MUST NOT silently widen the
# semantics of the guard -- update both the R-side validator and the parity
# tests if the literal value or its meaning ever changes.
CLAIM_GUARD = "not_for_de_or_new_quantitative_claims_without_full_object_validation"

MTX_AUTO_THRESHOLD = 100  # switch to mtx when marker count >= this

DEFAULT_OBS_COLUMNS = (
    "cell_type",
    "cell_type_hint",
    "CellName",
    "leiden",
    "sample",
    "patient",
    "stage",
    "condition",
    "orig.ident",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "doublet_score",
    "predicted_doublet",
    "immune_subtype",
    "tme_subtype",
)
DEFAULT_MARKERS = (
    "CD3E",
    "LYZ",
    "MS4A1",
    "NKG7",
    "ELF3",
    "EPCAM",
    "KRT8",
    "KRT18",
    "PTPRC",
)
DEFAULT_OBSM = ("X_umap", "X_pca")


@dataclass(frozen=True)
class ExportConfig:
    input_h5ad: Path
    output_dir: Path
    source_run_dir: Path | None = None
    obs_columns: tuple[str, ...] = DEFAULT_OBS_COLUMNS
    markers: tuple[str, ...] = DEFAULT_MARKERS
    obsm_keys: tuple[str, ...] = DEFAULT_OBSM
    compression: str = "gzip"
    max_cells: int | None = None
    seed: int = 1
    marker_chunk_size: int = 50000
    schema_version: str = "v2.1"
    format: str = "auto"
    # Phase B (v2.1) extensions: opt-in modality exporters. Default off so existing
    # callers get bit-for-bit identical bundles. The protein exporter is the first
    # consumer of `add_extension(...)`.
    include_protein: bool = False
    protein_obsm_key: str = "protein_clr"
    protein_isotype_controls: tuple[str, ...] = ()
    # Phase B (v2.1) spatial extension: opt-in. Default off so existing
    # bundles are byte-for-byte unchanged.
    include_spatial: bool = False
    spatial_obsm_key: str = "spatial"
    spatial_include_image_paths: bool = True
    # Phase B (v2.1) multimodal_obsm extension (EXPERIMENTAL): opt-in.
    # Publishes whichever of ``multimodal_obsm_keys`` are actually present on
    # the AnnData; each one round-trips as a separate parquet file.
    include_multimodal_obsm: bool = False
    multimodal_obsm_keys: tuple[str, ...] = ("X_wnn", "X_mofa")
    # v2.2 marker_resolutions extension: opt-in. Exports the per-cell-type
    # marker resolution table from context_aware_annotation. No-op when
    # adata.uns["marker_db_index"] is absent or context_aware_annotation did
    # not run (obs["context_aware_celltype"] missing).
    include_marker_resolutions: bool = False
    # v2.2 ATAC extension: opt-in, active only when the source AnnData carries
    # the canonical ATAC contract written by ATACIngestModule/ATACLSIModule.
    include_atac: bool = False
    atac_lsi_obsm_key: str = "X_lsi"
    atac_peaks_uns_key: str = "atac_peaks"
    # v2.2 Hi-C / single-cell 3D genome extension: opt-in, active only when
    # hic_ingest/hic_tad wrote the canonical uns payloads.
    include_hic: bool = False
    hic_max_contacts: int = 1_000_000
    # v2.2 Ribo-seq extension: opt-in, active only when ribo_ingest wrote
    # adata.uns["ribo_translation_efficiency"].
    include_ribo: bool = False


# ---------------------------------------------------------------------------
# Phase B: bundle extension API (additive, opt-in)
# ---------------------------------------------------------------------------

# Known extension keys. Readers ignore unknown keys with a "skipping unknown
# extension" message; producers may freely register additional keys. Keep this
# list in sync with the R reader's known-extensions handling.
KNOWN_EXTENSION_KEYS = (
    "protein",
    "spatial",
    "multimodal_obsm",
    "marker_resolutions",
    "atac",
    "hic",
    "ribo",
)


def _unknown_hic_tad_metadata() -> dict[str, object]:
    return {
        "compartment_status": "unknown_or_unvalidated",
        "low_information_chromosomes": [],
        "compartment_status_by_chrom": {},
    }


def add_extension(manifest, name, *, version, files, claim_guard=CLAIM_GUARD, **fields):
    """Register an optional bundle extension on a manifest dict.

    This is the supported API for future module authors (proteomics / spatial
    / multimodal). It mutates ``manifest["extensions"]`` in place and is
    idempotent for the same ``name`` (a second call replaces the prior entry).

    Parameters
    ----------
    manifest : dict
        The bundle manifest (must already be a dict; raises TypeError otherwise).
        ``manifest["extensions"]`` is created on first use.
    name : str
        Extension key (e.g. ``"protein"``, ``"spatial"``, ``"multimodal_obsm"``).
        Unknown keys are accepted -- readers skip them gracefully -- but a
        debug-level message is implied for forward-compat.
    version : str
        Extension-internal schema version (e.g. ``"1.0"``).
    files : list[str]
        File stems within the bundle that belong to this extension. Their full
        records still live in ``manifest["files"]``; this list is just an
        index for readers.
    claim_guard : str
        By default the same module-level ``CLAIM_GUARD`` (plotting-only).
        Extensions MUST NOT silently weaken this guard.
    **fields
        Modality-specific small fields (e.g.
        ``protein.normalization``, ``spatial.coord_system``,
        ``multimodal_obsm.embeddings``).

    Returns
    -------
    dict
        The extension entry that was stored in ``manifest["extensions"][name]``.

    Examples
    --------
    Illustrative only -- no proteomics logic is wired up in this phase::

        # add_extension(
        #     manifest,
        #     "protein",
        #     version="1.0",
        #     files=["protein_expr"],
        #     normalization="CLR",
        # )
    """
    if not isinstance(manifest, dict):
        raise TypeError(
            f"add_extension expects manifest to be a dict, got {type(manifest).__name__}"
        )
    if not isinstance(name, str) or not name:
        raise TypeError("add_extension requires a non-empty string `name`")
    if not isinstance(version, str) or not version:
        raise TypeError("add_extension requires a non-empty string `version`")
    if not isinstance(files, (list, tuple)):
        raise TypeError(
            f"add_extension expects `files` to be a list/tuple, got {type(files).__name__}"
        )
    files_list = [str(f) for f in files]
    if not isinstance(claim_guard, str) or not claim_guard:
        raise TypeError("add_extension requires a non-empty string `claim_guard`")

    extensions = manifest.setdefault("extensions", {})
    if not isinstance(extensions, dict):
        raise TypeError(
            "manifest['extensions'] must be a dict; refusing to overwrite a "
            f"{type(extensions).__name__}"
        )

    entry: dict[str, object] = {
        "version": version,
        "files": files_list,
        "claim_guard": claim_guard,
    }
    # Modality-specific fields are flattened into the entry. We deliberately
    # let callers pass arbitrary keys here; the R reader does not validate
    # the inner shape (forward-compat).
    for key, value in fields.items():
        entry[key] = value

    extensions[name] = entry
    return entry


def maybe_export_protein(
    manifest: dict,
    output_dir: Path,
    adata,
    *,
    obsm_key: str = "protein_clr",
    isotype_controls: Sequence[str] = (),
) -> dict | None:
    """Write protein.parquet if adata has the modality; register via add_extension.

    Parameters
    ----------
    manifest : dict
        Bundle manifest (will be mutated: ``manifest["files"]["protein"]`` and
        ``manifest["extensions"]["protein"]`` are populated on success).
    output_dir : Path
        Bundle directory (typically the temp dir during a publish).
    adata : anndata.AnnData
        Source AnnData. Must expose ``adata.obsm[obsm_key]`` to be exported.
    obsm_key : str
        Key on ``adata.obsm`` holding the (cells x proteins) CLR matrix.
    isotype_controls : Sequence[str]
        Optional list of protein names that are isotype controls; surfaced
        in the extension entry for downstream R-side QC.

    Returns
    -------
    dict | None
        The extension entry that was registered, or ``None`` when the
        modality was not present (no-op, no file written).
    """
    if obsm_key not in getattr(adata, "obsm", {}):
        return None
    arr = np.asarray(adata.obsm[obsm_key])
    if arr.ndim != 2:
        raise ValueError(
            f"adata.obsm['{obsm_key}'] must be 2D (cells x proteins); got shape {arr.shape}"
        )
    n_cells, n_proteins = arr.shape

    cells = adata.obs_names[:n_cells].astype(str)
    protein_names = adata.uns.get("protein_names")
    if protein_names is None or len(protein_names) != n_proteins:
        protein_names = [f"ADT_{i + 1}" for i in range(n_proteins)]
    protein_names = [str(name) for name in protein_names]

    df = pd.DataFrame(arr.astype(np.float32, copy=False), index=cells, columns=protein_names)
    df.index.name = "cell"
    path = output_dir / "protein.parquet"
    _write_parquet(df, path)

    manifest_files = manifest.setdefault("files", {})
    manifest_files["protein"] = file_record(
        path, output_dir, n_cells, n_proteins, {"format": "parquet"}
    )

    return add_extension(
        manifest,
        "protein",
        version="1.0",
        files=["protein"],
        normalization="CLR",
        isotype_controls=list(isotype_controls),
        n_proteins=int(n_proteins),
    )


def maybe_export_spatial(
    manifest: dict,
    output_dir: Path,
    adata,
    *,
    obsm_key: str = "spatial",
    include_image_paths: bool = True,
) -> dict | None:
    """Write spatial.parquet if adata has the modality; register via add_extension.

    Parameters
    ----------
    manifest : dict
        Bundle manifest (mutated: ``manifest["files"]["spatial"]`` and
        ``manifest["extensions"]["spatial"]`` are populated on success).
    output_dir : Path
        Bundle directory.
    adata : anndata.AnnData
        Source AnnData; must expose ``adata.obsm[obsm_key]`` (cells x 2 coords).
    obsm_key : str
        Key on ``adata.obsm`` holding the (cells x 2) coordinate matrix.
    include_image_paths : bool
        When True (default), copy ``adata.uns['spatial']['library_id']`` path
        strings into the extension entry. We NEVER serialize image bytes; only
        path strings round-trip into the manifest.

    Returns
    -------
    dict | None
        The extension entry that was registered, or ``None`` if the modality
        was not present (no-op, no file written).
    """
    obsm = getattr(adata, "obsm", {})
    if obsm_key not in obsm:
        return None
    arr = np.asarray(obsm[obsm_key])
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(
            f"adata.obsm['{obsm_key}'] must be (cells x >=2); got shape {arr.shape}"
        )
    coords = arr[:, :2].astype(np.float32, copy=False)
    n_cells = coords.shape[0]
    cells = adata.obs_names[:n_cells].astype(str)

    df = pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1]}, index=cells)
    # Carry over optional sample / library_id obs columns when present so the
    # R loader can reconstruct per-library structure without separate files.
    if "spatial_sample" in adata.obs.columns:
        df["sample"] = adata.obs["spatial_sample"].astype(str).values[:n_cells]
    if "spatial_library_id" in adata.obs.columns:
        df["library_id"] = adata.obs["spatial_library_id"].astype(str).values[:n_cells]
    df.index.name = "cell"

    path = output_dir / "spatial.parquet"
    _write_parquet(df, path)

    manifest_files = manifest.setdefault("files", {})
    manifest_files["spatial"] = file_record(
        path, output_dir, n_cells, df.shape[1], {"format": "parquet"}
    )

    coord_system = {}
    spatial_uns = adata.uns.get("spatial_coord_system")
    if isinstance(spatial_uns, dict):
        coord_system = {str(k): _coerce_jsonable(v) for k, v in spatial_uns.items()}

    library_image_paths: dict[str, str] = {}
    if include_image_paths:
        spatial_pointer = adata.uns.get("spatial")
        if isinstance(spatial_pointer, dict):
            lib = spatial_pointer.get("library_id")
            if isinstance(lib, dict):
                # Force every value to a string. Bytes/Path/None are explicitly
                # rejected so an image blob can never sneak into the manifest.
                for key, value in lib.items():
                    if isinstance(value, (bytes, bytearray)):
                        raise ValueError(
                            "spatial library_image_paths must be string paths, "
                            f"not bytes (key={key!r})"
                        )
                    library_image_paths[str(key)] = str(value)

    return add_extension(
        manifest,
        "spatial",
        version="1.0",
        files=["spatial"],
        coord_system=coord_system,
        library_image_paths=library_image_paths,
        n_spots=int(n_cells),
    )


def maybe_export_multimodal_obsm(
    manifest: dict,
    output_dir: Path,
    adata,
    *,
    obsm_keys: Sequence[str] = ("X_wnn", "X_mofa"),
) -> dict | None:
    """Write multimodal_obsm parquets and register the ``multimodal_obsm`` extension.

    EXPERIMENTAL: this extension carries the same plotting-only ``CLAIM_GUARD``
    as the rest of the bundle. Embeddings are visualization aids; the
    cross-modality validation logic that would let a user make new
    quantitative claims is not yet wired up.

    Parameters
    ----------
    manifest : dict
        Bundle manifest (mutated: per-key entries land in
        ``manifest["files"]["multimodal_obsm_<key>"]`` and the extension
        record lands in ``manifest["extensions"]["multimodal_obsm"]``).
    output_dir : Path
        Bundle directory.
    adata : anndata.AnnData
        Source AnnData. For each key in ``obsm_keys`` that is actually
        present on ``adata.obsm``, a ``multimodal_obsm_<key>.parquet`` is
        written.
    obsm_keys : Sequence[str]
        Candidate obsm keys to publish. Keys not present on the AnnData are
        silently skipped (so the same call works on a partial run).

    Returns
    -------
    dict | None
        The extension entry that was registered, or ``None`` if none of the
        candidate keys were present (no-op, no files written).
    """
    obsm = getattr(adata, "obsm", {})
    present: list[tuple[str, np.ndarray]] = []
    for key in obsm_keys:
        if key not in obsm:
            continue
        arr = np.asarray(obsm[key])
        if arr.ndim != 2:
            raise ValueError(
                f"adata.obsm['{key}'] must be 2D (cells x dims); got shape {arr.shape}"
            )
        present.append((str(key), arr))

    if not present:
        return None

    file_stems: list[str] = []
    embeddings_meta: list[dict[str, object]] = []
    manifest_files = manifest.setdefault("files", {})

    for key, arr in present:
        n_cells, n_dims = arr.shape
        cells = adata.obs_names[:n_cells].astype(str)
        cols = [f"{key}_{i + 1}" for i in range(n_dims)]
        df = pd.DataFrame(arr.astype(np.float32, copy=False), index=cells, columns=cols)
        df.index.name = "cell"

        stem = f"multimodal_obsm_{key}"
        path = output_dir / f"{stem}.parquet"
        _write_parquet(df, path)

        manifest_files[stem] = file_record(
            path, output_dir, n_cells, n_dims, {"format": "parquet"}
        )
        file_stems.append(stem)
        embeddings_meta.append({
            "key": key,
            "n_dims": int(n_dims),
            "stem": stem,
        })

    engine_used = "unknown"
    multimodal_status = adata.uns.get("multimodal_status")
    if isinstance(multimodal_status, dict):
        engine_used = str(multimodal_status.get("engine", "unknown"))

    return add_extension(
        manifest,
        "multimodal_obsm",
        version="1.0",
        files=file_stems,
        embeddings=embeddings_meta,
        engine_used=engine_used,
        experimental=True,
    )


def maybe_export_marker_resolutions(
    manifest: dict,
    output_dir: Path,
    adata,
) -> dict | None:
    """Write extensions/marker_resolutions/markers.parquet and register the extension.

    Reads ``adata.uns["marker_db_index"]`` (written by marker_db_loader) and
    flattens all marker rows into the required schema:
    ``[cell_type, marker_set, source_db, source_version, score]``.

    No-op (returns None) when:
    - ``adata.uns["marker_db_index"]`` is absent or empty, OR
    - ``adata.obs["context_aware_celltype"]`` is absent (module did not run).

    Parameters
    ----------
    manifest : dict
        Bundle manifest (mutated on success).
    output_dir : Path
        Bundle directory (temp dir during publish).
    adata : anndata.AnnData
        Source AnnData (may be a backed view).

    Returns
    -------
    dict | None
        The extension entry registered in the manifest, or None on no-op.
    """
    marker_db_index = adata.uns.get("marker_db_index")
    if not marker_db_index:
        return None
    if "context_aware_celltype" not in getattr(adata, "obs", {}).columns if hasattr(adata, "obs") else True:
        return None

    rows = []
    for db_key, db_value in marker_db_index.items():
        # db_key is the source_db name; db_value is a dict or DataFrame-like
        # stored by marker_db_loader. Normalise to a list of dicts.
        if hasattr(db_value, "itertuples"):
            # It's a DataFrame (restored from checkpoint via anndata uns).
            for row in db_value.itertuples(index=False):
                rows.append({
                    "cell_type": str(getattr(row, "cell_type", "")),
                    "marker_set": str(getattr(row, "marker_set", db_key)),
                    "source_db": str(db_key),
                    "source_version": str(getattr(row, "source_version", "")),
                    "score": float(getattr(row, "score", 0.0)) if hasattr(row, "score") else None,
                })
        elif isinstance(db_value, dict):
            # Stored as a nested dict from anndata uns serialisation.
            for ct, markers in db_value.items():
                if isinstance(markers, (list, tuple)):
                    for gene in markers:
                        rows.append({
                            "cell_type": str(ct),
                            "marker_set": str(db_key),
                            "source_db": str(db_key),
                            "source_version": "",
                            "score": None,
                        })
                else:
                    rows.append({
                        "cell_type": str(ct),
                        "marker_set": str(db_key),
                        "source_db": str(db_key),
                        "source_version": "",
                        "score": None,
                    })

    if not rows:
        return None

    df = pd.DataFrame(rows, columns=["cell_type", "marker_set", "source_db", "source_version", "score"])
    df["score"] = pd.to_numeric(df["score"], errors="coerce").astype("float32")

    ext_dir = output_dir / "extensions" / "marker_resolutions"
    ext_dir.mkdir(parents=True, exist_ok=True)
    path = ext_dir / "markers.parquet"
    # marker_resolutions has no meaningful cell-index, so we don't reset index
    # via _write_parquet (which adds a "cell" column). Write directly.
    if not _HAVE_PYARROW:
        df.to_parquet(path, engine="fastparquet", index=False)
    else:
        table = pa.Table.from_pandas(df, preserve_index=False)
        pq.write_table(table, path, compression="snappy")

    manifest_files = manifest.setdefault("files", {})
    manifest_files["marker_resolutions"] = file_record(
        path, output_dir, df.shape[0], df.shape[1], {"format": "parquet"}
    )

    return add_extension(
        manifest,
        "marker_resolutions",
        version="1.0",
        files=["marker_resolutions"],
        columns=["cell_type", "marker_set", "source_db", "source_version", "score"],
        n_rows=int(df.shape[0]),
    )


def _write_table_parquet(df: pd.DataFrame, path: Path) -> None:
    """Write a non-cell-indexed extension table without adding a `cell` column."""
    if not _HAVE_PYARROW:
        df.to_parquet(path, engine="fastparquet", index=False)
        return
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(table, path, compression="snappy")


def maybe_export_atac(
    manifest: dict,
    output_dir: Path,
    adata,
    *,
    lsi_obsm_key: str = "X_lsi",
    peaks_uns_key: str = "atac_peaks",
) -> dict | None:
    """Write the v2.2 ATAC extension when LSI and peak metadata are present.

    The exported payload is intentionally plotting/reporting scale: dense LSI
    coordinates plus peak metadata. The sparse peak-count matrix remains in the
    AnnData and is not serialized into the compact R plotting bundle.
    """
    obsm = getattr(adata, "obsm", {})
    if lsi_obsm_key not in obsm:
        return None
    peaks = getattr(adata, "uns", {}).get(peaks_uns_key)
    if peaks is None:
        peaks = getattr(adata, "uns", {}).get("atac_var")
    if peaks is None:
        return None

    lsi = np.asarray(obsm[lsi_obsm_key])
    if lsi.ndim != 2:
        raise ValueError(
            f"adata.obsm['{lsi_obsm_key}'] must be 2D (cells x LSI dims); got shape {lsi.shape}"
        )
    n_cells, n_components = lsi.shape
    cells = adata.obs_names[:n_cells].astype(str)
    lsi_df = pd.DataFrame(
        lsi.astype(np.float32, copy=False),
        index=cells,
        columns=[f"lsi_{i + 1}" for i in range(n_components)],
    )
    lsi_df.index.name = "cell"

    if isinstance(peaks, pd.DataFrame):
        peaks_df = peaks.copy()
    else:
        peaks_df = pd.DataFrame(peaks)
    required_cols = {"chrom", "start", "end"}
    missing_cols = required_cols - set(peaks_df.columns)
    if missing_cols:
        raise ValueError(
            "ATAC peaks metadata must contain chrom/start/end columns; "
            f"missing {sorted(missing_cols)}"
        )
    peaks_df = peaks_df.copy()
    peaks_df["chrom"] = peaks_df["chrom"].astype(str)
    peaks_df["start"] = pd.to_numeric(peaks_df["start"], errors="raise").astype("int64")
    peaks_df["end"] = pd.to_numeric(peaks_df["end"], errors="raise").astype("int64")
    if "peak_id" not in peaks_df.columns:
        peaks_df["peak_id"] = (
            peaks_df["chrom"].astype(str) + ":" +
            peaks_df["start"].astype(str) + "-" +
            peaks_df["end"].astype(str)
        )
    peaks_df = peaks_df[["peak_id", "chrom", "start", "end"]]

    ext_dir = output_dir / "extensions" / "atac"
    ext_dir.mkdir(parents=True, exist_ok=True)
    lsi_path = ext_dir / "lsi.parquet"
    peaks_path = ext_dir / "peaks.parquet"
    _write_parquet(lsi_df, lsi_path)
    _write_table_parquet(peaks_df, peaks_path)

    manifest_files = manifest.setdefault("files", {})
    manifest_files["atac_lsi"] = file_record(
        lsi_path, output_dir, lsi_df.shape[0], lsi_df.shape[1], {"format": "parquet"}
    )
    manifest_files["atac_peaks"] = file_record(
        peaks_path, output_dir, peaks_df.shape[0], peaks_df.shape[1], {"format": "parquet"}
    )

    return add_extension(
        manifest,
        "atac",
        version="1.0",
        files=["atac_lsi", "atac_peaks"],
        status="active",
        table={
            "lsi_path": "extensions/atac/lsi.parquet",
            "peaks_path": "extensions/atac/peaks.parquet",
            "index_column": "cell",
            "lsi_index_column": "cell",
            "peak_id_column": "peak_id",
        },
        n_cells=int(n_cells),
        n_components=int(n_components),
        n_peaks=int(peaks_df.shape[0]),
        method="tfidf_lsi",
    )


def _normalize_hic_bins(raw_bins) -> pd.DataFrame:
    if raw_bins is None:
        raise ValueError("Hi-C export requires adata.uns['hic_bins']")
    bins_df = raw_bins.copy() if isinstance(raw_bins, pd.DataFrame) else pd.DataFrame(raw_bins)
    if "bin_id" not in bins_df.columns:
        bins_df = bins_df.reset_index(drop=True)
        bins_df.insert(0, "bin_id", np.arange(bins_df.shape[0], dtype=np.int64))
    required = {"bin_id", "chrom", "start", "end"}
    missing = required - set(bins_df.columns)
    if missing:
        raise ValueError(f"Hi-C bins table missing required columns: {sorted(missing)}")
    bins_df = bins_df[["bin_id", "chrom", "start", "end"]].copy()
    bins_df["bin_id"] = pd.to_numeric(bins_df["bin_id"], errors="raise").astype("int64")
    bins_df["chrom"] = bins_df["chrom"].astype(str)
    bins_df["start"] = pd.to_numeric(bins_df["start"], errors="raise").astype("int64")
    bins_df["end"] = pd.to_numeric(bins_df["end"], errors="raise").astype("int64")
    if bins_df["bin_id"].duplicated().any():
        raise ValueError("Hi-C bins table contains duplicate bin_id values")
    if not ((bins_df["end"] > bins_df["start"]).all()):
        raise ValueError("Hi-C bins table requires end > start for every bin")
    return bins_df.sort_values("bin_id").reset_index(drop=True)


def _hic_position_from_bins(df: pd.DataFrame, bins_df: pd.DataFrame) -> pd.Series:
    if "start" in df.columns and "end" in df.columns:
        start = pd.to_numeric(df["start"], errors="raise").astype("int64")
        end = pd.to_numeric(df["end"], errors="raise").astype("int64")
        return ((start + end) // 2).astype("int64")
    joined = df[["bin_id"]].merge(
        bins_df[["bin_id", "start", "end"]],
        on="bin_id",
        how="left",
        validate="many_to_one",
    )
    if joined[["start", "end"]].isna().any().any():
        raise ValueError("Hi-C table has bin_id values absent from hic_bins")
    return ((joined["start"].astype("int64") + joined["end"].astype("int64")) // 2).astype("int64")


def _normalize_hic_optional_table(raw_table, bins_df: pd.DataFrame, kind: str) -> pd.DataFrame:
    if kind == "boundaries":
        columns = ["bin_id", "chrom", "position", "insulation", "is_boundary"]
        empty = pd.DataFrame({
            "bin_id": pd.Series(dtype="int64"),
            "chrom": pd.Series(dtype="object"),
            "position": pd.Series(dtype="int64"),
            "insulation": pd.Series(dtype="float32"),
            "is_boundary": pd.Series(dtype="bool"),
        })
    elif kind == "compartments":
        columns = ["bin_id", "chrom", "position", "eigenvector_1", "compartment"]
        empty = pd.DataFrame({
            "bin_id": pd.Series(dtype="int64"),
            "chrom": pd.Series(dtype="object"),
            "position": pd.Series(dtype="int64"),
            "eigenvector_1": pd.Series(dtype="float32"),
            "compartment": pd.Series(dtype="object"),
        })
    else:  # pragma: no cover
        raise ValueError(f"Unknown Hi-C optional table kind: {kind}")

    if raw_table is None:
        return empty

    df = raw_table.copy() if isinstance(raw_table, pd.DataFrame) else pd.DataFrame(raw_table)
    if df.empty:
        return empty
    if "bin_id" not in df.columns:
        df = df.reset_index(drop=True)
        if df.shape[0] != bins_df.shape[0]:
            raise ValueError(f"Hi-C {kind} table lacks bin_id and row count does not match hic_bins")
        df.insert(0, "bin_id", bins_df["bin_id"].to_numpy())

    df = df.copy()
    df["bin_id"] = pd.to_numeric(df["bin_id"], errors="raise").astype("int64")
    unknown_bins = set(df["bin_id"]) - set(bins_df["bin_id"])
    if unknown_bins:
        raise ValueError(f"Hi-C {kind} table has bin_id values absent from hic_bins: {sorted(unknown_bins)[:5]}")
    if "chrom" not in df.columns:
        df = df.merge(bins_df[["bin_id", "chrom"]], on="bin_id", how="left", validate="many_to_one")
    if df["chrom"].isna().any():
        raise ValueError(f"Hi-C {kind} table has bin_id values absent from hic_bins")
    df["chrom"] = df["chrom"].astype(str)
    if "position" not in df.columns:
        df["position"] = _hic_position_from_bins(df, bins_df)
    df["position"] = pd.to_numeric(df["position"], errors="raise").astype("int64")

    if kind == "boundaries":
        required = {"insulation", "is_boundary"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Hi-C boundaries table missing required columns: {sorted(missing)}")
        df["insulation"] = pd.to_numeric(df["insulation"], errors="coerce").astype("float32")
        df["is_boundary"] = df["is_boundary"].astype(bool)
    else:
        if "eigenvector_1" not in df.columns:
            raise ValueError("Hi-C compartments table missing required column: eigenvector_1")
        df["eigenvector_1"] = pd.to_numeric(df["eigenvector_1"], errors="coerce").astype("float32")
        if "compartment" not in df.columns:
            df["compartment"] = np.where(df["eigenvector_1"] >= 0, "A", "B")
        df["compartment"] = df["compartment"].astype(str)
        allowed_compartments = {"A", "B", "low_information"}
        bad = set(df["compartment"].dropna().unique()) - allowed_compartments
        if bad:
            raise ValueError(
                "Hi-C compartments must be A/B/low_information labels; "
                f"observed {sorted(bad)}"
            )

    return df[columns].sort_values("bin_id").reset_index(drop=True)


def maybe_export_hic(
    manifest: dict,
    output_dir: Path,
    adata,
    *,
    max_contacts: int = 1_000_000,
) -> dict | None:
    """Write the v2.2 Hi-C/scHi-C extension from canonical AnnData uns payloads."""
    uns = getattr(adata, "uns", {})
    if "hic_contact_matrix" not in uns and "hic_bins" not in uns:
        return None
    if "hic_contact_matrix" not in uns or "hic_bins" not in uns:
        raise ValueError("Hi-C export requires both hic_contact_matrix and hic_bins")
    if max_contacts <= 0:
        raise ValueError("--hic-max-contacts must be positive")

    bins_df = _normalize_hic_bins(uns["hic_bins"])
    mat = uns["hic_contact_matrix"]
    if hasattr(mat, "to_memory"):
        mat = mat.to_memory()
    if sparse.issparse(mat):
        coo = mat.tocoo()
    else:
        arr = np.asarray(mat)
        if arr.ndim != 2:
            raise ValueError(f"Hi-C contact matrix must be 2D, got shape {arr.shape}")
        coo = sparse.coo_matrix(arr)
    if coo.shape[0] != bins_df.shape[0] or coo.shape[1] != bins_df.shape[0]:
        raise ValueError(
            "Hi-C contact matrix shape does not match hic_bins: "
            f"matrix={coo.shape}, bins={bins_df.shape[0]}"
        )
    if int(coo.nnz) > int(max_contacts):
        raise ValueError(
            f"Hi-C contact matrix has {coo.nnz} non-zero entries, above "
            f"--hic-max-contacts={max_contacts}. Export a coarser/bin-filtered matrix."
        )
    contacts_df = pd.DataFrame({
        "row": coo.row.astype(np.int64, copy=False),
        "col": coo.col.astype(np.int64, copy=False),
        "count": coo.data.astype(np.float32, copy=False),
    })

    boundaries_df = _normalize_hic_optional_table(
        uns.get("hic_tad_boundaries"), bins_df, "boundaries"
    )
    compartments_df = _normalize_hic_optional_table(
        uns.get("hic_compartments"), bins_df, "compartments"
    )

    ext_dir = output_dir / "extensions" / "hic"
    ext_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "hic_bins": ext_dir / "bins.parquet",
        "hic_contacts": ext_dir / "contacts.parquet",
        "hic_boundaries": ext_dir / "boundaries.parquet",
        "hic_compartments": ext_dir / "compartments.parquet",
    }
    _write_table_parquet(bins_df, paths["hic_bins"])
    _write_table_parquet(contacts_df, paths["hic_contacts"])
    _write_table_parquet(boundaries_df, paths["hic_boundaries"])
    _write_table_parquet(compartments_df, paths["hic_compartments"])

    manifest_files = manifest.setdefault("files", {})
    tables = {
        "hic_bins": bins_df,
        "hic_contacts": contacts_df,
        "hic_boundaries": boundaries_df,
        "hic_compartments": compartments_df,
    }
    for stem, df in tables.items():
        manifest_files[stem] = file_record(
            paths[stem],
            output_dir,
            df.shape[0],
            df.shape[1],
            {"format": "parquet"},
        )

    ingest_meta = uns.get("hic_ingest_metadata")
    if not isinstance(ingest_meta, dict):
        ingest_meta = {}
    hic_tad_metadata = uns.get("hic_tad_metadata")
    if isinstance(hic_tad_metadata, dict):
        hic_tad_metadata = _coerce_jsonable(hic_tad_metadata)
    else:
        hic_tad_metadata = _unknown_hic_tad_metadata()

    extension_fields = {
        "status": "active",
        "table": {
            "bins_path": "extensions/hic/bins.parquet",
            "contacts_path": "extensions/hic/contacts.parquet",
            "boundaries_path": "extensions/hic/boundaries.parquet",
            "compartments_path": "extensions/hic/compartments.parquet",
            "index_column": "bin_id",
            "contact_row_column": "row",
            "contact_col_column": "col",
            "contact_count_column": "count",
        },
        "n_bins": int(bins_df.shape[0]),
        "n_contacts": int(contacts_df.shape[0]),
        "n_boundaries": int(boundaries_df.shape[0]),
        "n_boundary_bins": int(boundaries_df["is_boundary"].sum()) if "is_boundary" in boundaries_df else 0,
        "n_compartments": int(compartments_df.shape[0]),
        "resolution_bp": int(ingest_meta["resolution_bp"]) if "resolution_bp" in ingest_meta else None,
        "method": "sparse_contact_tad_compartment",
        "normalization": str(ingest_meta.get("normalization", "raw_counts_or_unbalanced")),
        "matrix_format": str(ingest_meta.get("matrix_format", "csr_sparse")),
        "max_contacts": int(max_contacts),
    }
    extension_fields["hic_tad_metadata"] = hic_tad_metadata

    return add_extension(
        manifest,
        "hic",
        version="1.0",
        files=["hic_bins", "hic_contacts", "hic_boundaries", "hic_compartments"],
        **extension_fields,
    )


def maybe_export_ribo(
    manifest: dict,
    output_dir: Path,
    adata,
) -> dict | None:
    """Write the v2.2 Ribo-seq extension from ribo_ingest AnnData uns payloads."""
    uns = getattr(adata, "uns", {})
    if "ribo_translation_efficiency" not in uns:
        return None

    te = uns["ribo_translation_efficiency"]
    if not isinstance(te, pd.DataFrame):
        te = pd.DataFrame(te)
    required = {"gene", "sample", "footprint_count", "rna_count", "te"}
    missing = required - set(te.columns)
    if missing:
        raise ValueError(
            "Ribo-seq export requires columns "
            f"{sorted(required)}; missing {sorted(missing)}"
        )
    te_df = te.loc[:, ["gene", "sample", "footprint_count", "rna_count", "te"]].copy()
    if te_df.empty:
        raise ValueError("Ribo-seq export requires at least one row")

    for col in ("gene", "sample"):
        if te_df[col].isna().any():
            raise ValueError("Ribo-seq export requires non-empty gene and sample identifiers")
        te_df[col] = te_df[col].astype(str).str.strip()
    if (te_df[["gene", "sample"]] == "").any(axis=None):
        raise ValueError("Ribo-seq export requires non-empty gene and sample identifiers")
    if te_df.duplicated(["gene", "sample"]).any():
        raise ValueError("Ribo-seq export found duplicate gene/sample rows")

    for col in ("footprint_count", "rna_count", "te"):
        values = pd.to_numeric(te_df[col], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(
                "Ribo-seq export requires finite numeric values in "
                "footprint_count, rna_count, and te"
            )
        te_df[col] = values

    count_values = te_df[["footprint_count", "rna_count"]].to_numpy(dtype=float)
    if (count_values < 0).any() or (te_df["te"].to_numpy(dtype=float) < 0).any():
        raise ValueError("Ribo-seq export requires non-negative counts and te values")
    if not np.equal(count_values, np.floor(count_values)).all():
        raise ValueError("Ribo-seq export requires integer-valued footprint_count and rna_count")

    pseudocount = 1.0
    expected_te = te_df["footprint_count"] / (te_df["rna_count"] + pseudocount)
    if not np.allclose(
        te_df["te"].to_numpy(dtype=float),
        expected_te.to_numpy(dtype=float),
        rtol=1e-9,
        atol=1e-12,
    ):
        raise ValueError(
            "Ribo-seq export te does not match footprint_count / (rna_count + 1.0)"
        )

    ext_dir = output_dir / "extensions" / "ribo"
    ext_dir.mkdir(parents=True, exist_ok=True)
    te_path = ext_dir / "translation_efficiency.parquet"
    _write_table_parquet(te_df, te_path)

    manifest_files = manifest.setdefault("files", {})
    manifest_files["ribo_te"] = file_record(
        te_path,
        output_dir,
        te_df.shape[0],
        te_df.shape[1],
        {"format": "parquet"},
    )

    extension_fields = {
        "status": "active",
        "table": {
            "te_path": "extensions/ribo/translation_efficiency.parquet",
            "index_column": "gene",
        },
        "n_genes": int(te_df["gene"].nunique()) if not te_df.empty else 0,
        "n_samples": int(te_df["sample"].nunique()) if not te_df.empty else 0,
        "n_rows": int(te_df.shape[0]),
        "te_mean": float(te_df["te"].mean()),
        "te_median": float(te_df["te"].median()),
        "method": "footprint_over_rna_plus_pseudocount",
        "pseudocount": pseudocount,
    }
    return add_extension(
        manifest,
        "ribo",
        version="1.0",
        files=["ribo_te"],
        **extension_fields,
    )


def _coerce_jsonable(value):
    """Best-effort coercion of small uns metadata values into JSON-friendly types."""
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_coerce_jsonable(v) for v in value.tolist()]
    if isinstance(value, (bytes, bytearray)):
        # Reject silently dropping binary metadata into the manifest.
        raise ValueError("bundle extension metadata must not contain bytes.")
    if isinstance(value, dict):
        return {str(k): _coerce_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_coerce_jsonable(v) for v in value]
    return str(value)


def parse_csv_list(value: str | None, default: Sequence[str]) -> tuple[str, ...]:
    if value is None:
        return tuple(default)
    items = tuple(item.strip() for item in value.split(",") if item.strip())
    return items or tuple(default)


def output_path(output_dir: Path, stem: str, compression: str) -> Path:
    suffix = ".csv.gz" if compression == "gzip" else ".csv"
    return output_dir / f"{stem}{suffix}"


def open_text(path: Path, compression: str):
    if compression == "gzip":
        return gzip.open(path, "wt", newline="")
    return path.open("w", newline="")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, output_dir: Path, n_rows: int, n_cols: int,
                extra: dict | None = None) -> dict[str, object]:
    rec: dict[str, object] = {
        "path": path.relative_to(output_dir).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "n_rows": n_rows,
        "n_cols": n_cols,
    }
    if extra:
        rec.update(extra)
    return rec


def choose_cells(n_obs: int, max_cells: int | None, seed: int) -> np.ndarray:
    if max_cells is None or max_cells >= n_obs:
        return np.arange(n_obs)
    if max_cells <= 0:
        raise ValueError("--max-cells must be positive when provided.")
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_obs, size=max_cells, replace=False))


def materialize_small_matrix(matrix) -> np.ndarray:
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def write_dataframe_csv(df: pd.DataFrame, path: Path, compression: str) -> None:
    with open_text(path, compression) as handle:
        df.to_csv(handle, index=True, index_label="cell")


# ---------------------------------------------------------------------------
# v1 writers (unchanged)
# ---------------------------------------------------------------------------

def write_obs(adata, cell_idx: np.ndarray, config: ExportConfig) -> tuple[Path, list[str]]:
    present = [col for col in config.obs_columns if col in adata.obs.columns]
    obs = adata.obs.iloc[cell_idx][present].copy()
    path = output_path(config.output_dir, "obs", config.compression)
    write_dataframe_csv(obs, path, config.compression)
    return path, present


def embedding_columns(key: str, n_cols: int) -> list[str]:
    if key == "X_umap":
        return [f"UMAP_{i + 1}" for i in range(n_cols)]
    if key == "X_pca":
        return [f"PC_{i + 1}" for i in range(n_cols)]
    clean = key[2:] if key.startswith("X_") else key
    return [f"{clean}_{i + 1}" for i in range(n_cols)]


def write_obsm(adata, cell_idx: np.ndarray, config: ExportConfig) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    cells = adata.obs_names[cell_idx].astype(str)
    for key in config.obsm_keys:
        if key not in adata.obsm:
            raise ValueError(f"Missing required embedding adata.obsm['{key}'].")
        arr = np.asarray(adata.obsm[key][cell_idx])
        if arr.ndim != 2:
            raise ValueError(f"Embedding {key} must be 2D, got shape {arr.shape}.")
        stem = key
        path = output_path(config.output_dir, stem, config.compression)
        df = pd.DataFrame(arr, index=cells, columns=embedding_columns(key, arr.shape[1]))
        write_dataframe_csv(df, path, config.compression)
        records[stem] = file_record(path, config.output_dir, df.shape[0], df.shape[1])
    return records


def write_marker_expr(adata, cell_idx: np.ndarray, config: ExportConfig) -> tuple[Path, list[str]]:
    present = [gene for gene in config.markers if gene in adata.var_names]
    if not present:
        raise ValueError("None of the requested marker genes were found in adata.var_names.")

    path = output_path(config.output_dir, "marker_expr", config.compression)
    cells = adata.obs_names[cell_idx].astype(str)
    first = True
    with open_text(path, config.compression) as handle:
        for start in range(0, len(cell_idx), config.marker_chunk_size):
            stop = min(start + config.marker_chunk_size, len(cell_idx))
            rows = cell_idx[start:stop]
            block = materialize_small_matrix(adata[rows, present].X)
            block = block.astype(np.float32, copy=False)
            df = pd.DataFrame(block, index=cells[start:stop], columns=present)
            df.to_csv(handle, index=True, index_label="cell", header=first)
            first = False
    return path, present


def write_manifest_tsv(manifest: dict[str, object], path: Path) -> None:
    summary_rows = {
        "schema_version": manifest["schema_version"],
        "generated_at_utc": manifest["generated_at_utc"],
        "input_h5ad": manifest["source"]["input_h5ad"],
        "source_run_dir": manifest["source"].get("source_run_dir") or "",
        "source_n_obs": manifest["source"]["n_obs"],
        "source_n_vars": manifest["source"]["n_vars"],
        "input_h5ad_bytes": manifest["source"]["input_h5ad_bytes"],
        "input_h5ad_mtime_epoch": manifest["source"]["input_h5ad_mtime_epoch"],
        "n_cells_exported": manifest["bundle"]["n_cells_exported"],
        "max_cells": manifest["bundle"].get("max_cells") or "",
        "seed": manifest["bundle"]["seed"],
        "obs_columns": ",".join(manifest["bundle"]["obs_columns_present"]),
        "markers_requested": ",".join(manifest["bundle"]["markers_requested"]),
        "markers_present": ",".join(manifest["bundle"]["markers_present"]),
        "obsm_keys": ",".join(manifest["bundle"]["obsm_keys"]),
        "required_files": ",".join(manifest["bundle"]["required_files"]),
        "expression_source_slot": manifest["expression"]["source_slot"],
        "expression_value_scale": manifest["expression"]["value_scale"],
        "expression_export_dtype": manifest["expression"]["export_dtype"],
        "intended_use": manifest["expression"]["intended_use"],
        "claim_guard": manifest["expression"]["claim_guard"],
        "precision_policy": manifest["precision_policy"],
    }
    for stem, record in sorted(manifest["files"].items()):
        for key in ("path", "bytes", "sha256", "n_rows", "n_cols"):
            summary_rows[f"file_{stem}_{key}"] = record.get(key, "")

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["key", "value"])
        for key, value in summary_rows.items():
            writer.writerow([key, value])


# ---------------------------------------------------------------------------
# v2 writers
# ---------------------------------------------------------------------------

def _write_parquet(df: pd.DataFrame, path: Path) -> None:
    # M1: ensure the row-index column is literally named `cell` in the written
    # parquet, regardless of pyarrow vs fastparquet (pyarrow's preserve_index
    # path historically emits `__index_level_0__` when the index is unnamed).
    # Reset_index after naming so `cell` becomes a real column the R reader
    # can match by name without relying on engine-specific index conventions.
    df = df.copy()
    df.index.name = "cell"
    df = df.reset_index()
    if not _HAVE_PYARROW:
        df.to_parquet(path, engine="fastparquet", index=False)
        return
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(table, path, compression="snappy")


def write_obs_parquet(adata, cell_idx: np.ndarray, config: ExportConfig) -> tuple[Path, list[str]]:
    present = [col for col in config.obs_columns if col in adata.obs.columns]
    obs = adata.obs.iloc[cell_idx][present].copy()
    # M2: do NOT round-trip categoricals through astype(str).astype("category").
    # That destroys ordered/categories metadata. Pass categoricals through as-is;
    # parquet preserves dtype, ordered flag, and category levels natively.
    path = config.output_dir / "obs.parquet"
    _write_parquet(obs, path)
    return path, present


def write_obsm_parquet(adata, cell_idx: np.ndarray, config: ExportConfig) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    cells = adata.obs_names[cell_idx].astype(str)
    obsm_dir = config.output_dir / "obsm"
    obsm_dir.mkdir(parents=True, exist_ok=True)
    for key in config.obsm_keys:
        if key not in adata.obsm:
            raise ValueError(f"Missing required embedding adata.obsm['{key}'].")
        arr = np.asarray(adata.obsm[key][cell_idx])
        if arr.ndim != 2:
            raise ValueError(f"Embedding {key} must be 2D, got shape {arr.shape}.")
        df = pd.DataFrame(arr, index=cells, columns=embedding_columns(key, arr.shape[1]))
        df.index.name = "cell"
        path = obsm_dir / f"{key}.parquet"
        _write_parquet(df, path)
        records[key] = file_record(path, config.output_dir, df.shape[0], df.shape[1],
                                   {"format": "parquet"})
    return records


def _use_mtx_format(markers_present: list[str], fmt: str) -> bool:
    # M5: scipy.io.mmwrite is the only writer for mtx_gz. If scipy.io is not
    # importable, refuse to silently fall through to a non-mtx format -- callers
    # explicitly asking for "mtx" deserve a hard error pointing at the fix.
    if fmt == "mtx":
        if not _HAVE_SCIPY_IO:
            raise ImportError(
                "scipy.io is required for --format mtx (mtx_gz output). "
                "Install scipy >= 1.10 in the export environment, or rerun "
                "with --format parquet / --format csv."
            )
        return True
    if fmt == "csv" or fmt == "parquet":
        return False
    # auto: use mtx when marker count >= threshold AND scipy.io is available;
    # otherwise stay on parquet so we never produce a half-baked mtx bundle.
    if not _HAVE_SCIPY_IO:
        return False
    return len(markers_present) >= MTX_AUTO_THRESHOLD


def write_marker_expr_mtx(adata, cell_idx: np.ndarray, markers_present: list[str],
                          config: ExportConfig) -> dict[str, dict[str, object]]:
    """Write marker expression as mtx.gz + barcodes + genes sidecar files."""
    cells = adata.obs_names[cell_idx].astype(str)
    X_sub = adata[cell_idx, markers_present].X
    if hasattr(X_sub, "to_memory"):
        X_sub = X_sub.to_memory()
    if not sparse.issparse(X_sub):
        X_sub = sparse.csr_matrix(X_sub.astype(np.float32))
    else:
        X_sub = X_sub.astype(np.float32).tocsr()

    mtx_path = config.output_dir / "marker_expr.mtx.gz"
    barcodes_path = config.output_dir / "marker_expr.barcodes.tsv.gz"
    genes_path = config.output_dir / "marker_expr.genes.tsv.gz"

    with gzip.open(mtx_path, "wb") as fh:
        sio.mmwrite(fh, X_sub)

    with gzip.open(barcodes_path, "wt") as fh:
        fh.write("\n".join(cells) + "\n")

    with gzip.open(genes_path, "wt") as fh:
        fh.write("\n".join(markers_present) + "\n")

    n_rows, n_cols = X_sub.shape
    nnz = X_sub.nnz
    return {
        "marker_expr": file_record(
            mtx_path,
            config.output_dir,
            n_rows,
            n_cols,
            {"format": "mtx_gz", "sparse": True, "nnz": nnz},
        ),
        "marker_expr_barcodes": file_record(
            barcodes_path,
            config.output_dir,
            n_rows,
            1,
            {"format": "tsv_gz", "role": "marker_expr_row_index"},
        ),
        "marker_expr_genes": file_record(
            genes_path,
            config.output_dir,
            n_cols,
            1,
            {"format": "tsv_gz", "role": "marker_expr_col_index"},
        ),
    }


def write_marker_expr_parquet(adata, cell_idx: np.ndarray, markers_present: list[str],
                              config: ExportConfig) -> tuple[Path, dict[str, object]]:
    """Write marker expression as parquet (dense float32).

    Single-slice materialization: the marker matrix is the cell subset × the
    (small) marker gene set, so it fits in memory in one shot. Chunk-then-concat
    gave no memory benefit while issuing repeated backed-slice reads, so we read
    the full slice once. Output is byte/content-identical to the chunked path:
    same column order (``markers_present``), float32 dtype, row order, and the
    literal ``cell`` index name handled by ``_write_parquet``.
    """
    cells = adata.obs_names[cell_idx].astype(str)
    block = materialize_small_matrix(adata[cell_idx, markers_present].X)
    block = block.astype(np.float32, copy=False)
    df = pd.DataFrame(block, index=cells, columns=markers_present)
    df.index.name = "cell"
    path = config.output_dir / "marker_expr.parquet"
    _write_parquet(df, path)
    rec = file_record(path, config.output_dir, df.shape[0], df.shape[1],
                      {"format": "parquet"})
    return path, rec


def _extract_rank_genes_df(adata) -> pd.DataFrame | None:
    rgg = adata.uns.get("rank_genes_groups")
    if rgg is None:
        return None
    try:
        groups = list(rgg.get("names", {}).dtype.names or [])
        if not groups:
            return None
        rows = []
        for group in groups:
            names = rgg["names"][group]
            scores = rgg.get("scores", {}).get(group, [None] * len(names))
            logfcs = rgg.get("logfoldchanges", {}).get(group, [None] * len(names))
            pvals = rgg.get("pvals_adj", {}).get(group, [None] * len(names))
            for i, gene in enumerate(names):
                rows.append({
                    "group": group,
                    "gene": gene,
                    "score": float(scores[i]) if scores[i] is not None else None,
                    "logfoldchange": float(logfcs[i]) if logfcs[i] is not None else None,
                    "pval_adj": float(pvals[i]) if pvals[i] is not None else None,
                })
        return pd.DataFrame(rows)
    except Exception:
        return None


def _build_uns_summary(adata) -> dict:
    summary: dict = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
    }
    uns = adata.uns
    # module_status keys
    module_status = uns.get("module_status")
    if isinstance(module_status, dict):
        summary["module_status"] = {str(k): str(v) for k, v in module_status.items()}
    # gpu_mode / scale_mode from uns or obs
    for key in ("gpu_mode", "scale_mode"):
        if key in uns:
            summary[key] = str(uns[key])
        elif key in adata.obs.columns:
            summary[key] = str(adata.obs[key].iloc[0])
    return summary


def write_v2_files(adata, cell_idx: np.ndarray, config: ExportConfig) -> tuple[dict[str, dict], list[str], bool]:
    """Write v2-format files. Returns (files_dict, markers_present, used_mtx)."""
    files: dict[str, dict] = {}

    # obs.parquet
    obs_path, obs_present = write_obs_parquet(adata, cell_idx, config)
    files["obs"] = file_record(obs_path, config.output_dir, len(cell_idx), len(obs_present),
                               {"format": "parquet"})

    # obsm parquets
    obsm_records = write_obsm_parquet(adata, cell_idx, config)
    files.update(obsm_records)

    # marker expression
    markers_present = [gene for gene in config.markers if gene in adata.var_names]
    if not markers_present:
        raise ValueError("None of the requested marker genes were found in adata.var_names.")

    use_mtx = _use_mtx_format(markers_present, config.format)
    if use_mtx:
        files.update(write_marker_expr_mtx(adata, cell_idx, markers_present, config))
    else:
        _, marker_rec = write_marker_expr_parquet(adata, cell_idx, markers_present, config)
        files["marker_expr"] = marker_rec

    # rank_genes.parquet (optional)
    rg_df = _extract_rank_genes_df(adata)
    if rg_df is not None:
        rg_path = config.output_dir / "rank_genes.parquet"
        _write_parquet(rg_df, rg_path)
        files["rank_genes"] = file_record(rg_path, config.output_dir,
                                          rg_df.shape[0], rg_df.shape[1],
                                          {"format": "parquet"})

    # uns_summary.json
    uns_summary = _build_uns_summary(adata)
    uns_path = config.output_dir / "uns_summary.json"
    uns_path.write_text(json.dumps(uns_summary, indent=2) + "\n")
    files["uns_summary"] = file_record(uns_path, config.output_dir, 0, 0,
                                       {"format": "json"})

    return files, markers_present, use_mtx


# ---------------------------------------------------------------------------
# Bundle provenance helper (PREC-1)
# ---------------------------------------------------------------------------

def _r_factory_sha_at_export(r_factory_path: str | None = None) -> str:
    """Return the short HEAD SHA of the R factory repo at export time, or '' on failure."""
    if r_factory_path is None:
        r_factory_path = str(Path(__file__).resolve().parents[1].parent / "r_multiomics_factory")
    try:
        return subprocess.check_output(
            ["git", "-C", r_factory_path, "rev-parse", "--short=7", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return ""


def _bundle_sha256_concat(output_dir: Path) -> str:
    """SHA256 of canonical parquet files concatenated in order (obs, obsm/*, marker_expr, rank_genes)."""
    stems = ["obs.parquet"]
    obsm_dir = output_dir / "obsm"
    if obsm_dir.is_dir():
        stems += sorted(str(p.relative_to(output_dir).as_posix()) for p in obsm_dir.glob("*.parquet"))
    for optional in ("marker_expr.parquet", "marker_expr.mtx.gz", "rank_genes.parquet"):
        if (output_dir / optional).exists():
            stems.append(optional)
    digest = hashlib.sha256()
    for stem in stems:
        p = output_dir / stem
        if p.exists():
            with p.open("rb") as fh:
                for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                    digest.update(chunk)
    return digest.hexdigest()


def _backfill_project_manifest_bundle_sha256(run_dir: Path, bundle_output_dir: Path) -> str:
    """Write the real bundle digest into the project-root cross-factory manifest.

    The pipeline writes ``manifest.json`` before any bundle exists, so
    ``pipeline.py`` can only emit ``ctx.metadata.get("bundle_sha256", "")`` —
    and nothing ever sets that key. The field was therefore structurally empty
    in every run, including runs whose ``overall_status`` was ``complete``:
    the cross-factory manifest advertised a bundle checksum it never carried,
    so a consumer could not detect a truncated or substituted bundle.

    The authoritative digest is computed here at export time (see
    ``_bundle_sha256_concat``) and recorded in ``bundle/provenance.json``. This
    copies it into the manifest so producer and consumer agree on one value.
    Returns the digest, or "" when there is nothing to back-fill.
    """
    manifest_json = run_dir / "manifest.json"
    provenance_json = bundle_output_dir / "provenance.json"
    if not manifest_json.is_file() or not provenance_json.is_file():
        return ""
    try:
        digest = str(json.loads(provenance_json.read_text(encoding="utf-8")).get("bundle_sha256", ""))
        if not digest:
            return ""
        manifest = json.loads(manifest_json.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        warnings.warn(f"could not back-fill bundle_sha256 into {manifest_json}: {exc}", stacklevel=2)
        return ""
    manifest["bundle_sha256"] = digest
    manifest["bundle_sha256_source"] = str(provenance_json)
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return digest


def _write_bundle_provenance(output_dir: Path) -> None:
    """Write bundle/provenance.json with r_factory_sha_at_export and audit fields."""
    provenance = {
        "r_factory_sha_at_export": _r_factory_sha_at_export(),
        "exported_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "exported_by": socket.gethostname(),
        "bundle_sha256": _bundle_sha256_concat(output_dir),
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )


# ---------------------------------------------------------------------------
# Main export entry point
# ---------------------------------------------------------------------------

def export_bundle(config: ExportConfig) -> dict[str, object]:
    if config.compression not in {"gzip", "none"}:
        raise ValueError("--compression must be 'gzip' or 'none'.")
    if config.marker_chunk_size <= 0:
        raise ValueError("--marker-chunk-size must be positive.")
    if not config.input_h5ad.exists():
        raise FileNotFoundError(f"Input .h5ad not found: {config.input_h5ad}")
    if config.schema_version not in ("v1", "v2", "v2.1", "v2.2"):
        raise ValueError("--schema-version must be 'v1', 'v2', 'v2.1', or 'v2.2'.")
    if config.include_atac and config.schema_version != "v2.2":
        raise ValueError("--include-atac requires --schema-version v2.2")
    if config.include_hic and config.schema_version != "v2.2":
        raise ValueError("--include-hic requires --schema-version v2.2")
    if config.include_ribo and config.schema_version != "v2.2":
        raise ValueError("--include-ribo requires --schema-version v2.2")
    if config.format not in ("auto", "csv", "parquet", "mtx"):
        raise ValueError("--format must be 'auto', 'csv', 'parquet', or 'mtx'.")

    final_output_dir = config.output_dir
    parent_dir = final_output_dir.parent
    parent_dir.mkdir(parents=True, exist_ok=True)

    # C2: materialize the entire bundle into a sibling temp directory and only
    # publish it via os.replace at the very end (after the manifest is written).
    # This means readers either see no bundle dir at all, or a fully-populated
    # one with a final manifest -- never a partially written intermediate.
    temp_dir = parent_dir / f"{final_output_dir.name}.tmp.{os.getpid()}"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=False)

    temp_config = dataclass_replace(config, output_dir=temp_dir)

    backup_dir: Path | None = None
    try:
        adata = ad.read_h5ad(config.input_h5ad, backed="r")
        try:
            cell_idx = choose_cells(adata.n_obs, config.max_cells, config.seed)

            if temp_config.schema_version in ("v2", "v2.1", "v2.2"):
                manifest = _export_bundle_v2(adata, cell_idx, temp_config)
            else:
                manifest = _export_bundle_v1(adata, cell_idx, temp_config)
        finally:
            adata.file.close()

        # Atomic publish step. If the final dir already exists, move it aside
        # first so we can either restore it on a failed rename or remove it on
        # success -- this gives readers a single rename event to observe.
        if final_output_dir.exists():
            backup_dir = parent_dir / f"{final_output_dir.name}.bak.{os.getpid()}"
            if backup_dir.exists():
                shutil.rmtree(backup_dir)
            os.replace(final_output_dir, backup_dir)

        os.replace(temp_dir, final_output_dir)

        if backup_dir is not None and backup_dir.exists():
            shutil.rmtree(backup_dir, ignore_errors=True)

        return manifest
    except Exception:
        # Best-effort cleanup of the temp dir; restore the previous output_dir
        # from backup if we managed to move it aside but failed before rename.
        if temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)
        if backup_dir is not None and backup_dir.exists() and not final_output_dir.exists():
            try:
                os.replace(backup_dir, final_output_dir)
            except OSError:
                pass
        raise


def _export_bundle_v1(adata, cell_idx: np.ndarray, config: ExportConfig) -> dict[str, object]:
    obs_path, obs_present = write_obs(adata, cell_idx, config)
    obsm_records = write_obsm(adata, cell_idx, config)
    marker_path, markers_present = write_marker_expr(adata, cell_idx, config)

    files = {"obs": file_record(obs_path, config.output_dir, len(cell_idx), len(obs_present))}
    files.update(obsm_records)
    files["marker_expr"] = file_record(
        marker_path, config.output_dir, len(cell_idx), len(markers_present)
    )

    source_stat = config.input_h5ad.stat()
    manifest = {
        "schema_version": SCHEMA_VERSION_V1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "input_h5ad": str(config.input_h5ad),
            "source_run_dir": str(config.source_run_dir) if config.source_run_dir else None,
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "input_h5ad_bytes": source_stat.st_size,
            "input_h5ad_mtime": source_stat.st_mtime,
            "input_h5ad_mtime_epoch": int(source_stat.st_mtime),
            "x_storage": type(adata.X).__name__,
        },
        "bundle": {
            "schema_name": "compact_metadata_embedding_marker_bundle",
            "n_cells_exported": int(len(cell_idx)),
            "max_cells": config.max_cells,
            "seed": config.seed,
            "obs_columns_requested": list(config.obs_columns),
            "obs_columns_present": obs_present,
            "markers_requested": list(config.markers),
            "markers_present": markers_present,
            "obsm_keys": list(config.obsm_keys),
            "compression": config.compression,
            "required_files": sorted(files),
        },
        "expression": {
            "source_slot": EXPRESSION_SOURCE_SLOT,
            "value_scale": EXPRESSION_VALUE_SCALE,
            "export_dtype": EXPRESSION_EXPORT_DTYPE,
            "intended_use": INTENDED_USE,
            "claim_guard": CLAIM_GUARD,
        },
        "files": files,
        "precision_policy": (
            "Metadata, embeddings, and selected marker expression are copied from the "
            "source AnnData. The full expression matrix is not converted to dense data. "
            "This bundle is for plotting/reporting handoff, not differential expression "
            "or new quantitative biological claims without full-object validation."
        ),
    }
    manifest_json = config.output_dir / "bundle_manifest.json"
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    write_manifest_tsv(manifest, config.output_dir / "bundle_manifest.tsv")
    (config.output_dir / "README.md").write_text(
        "# singlecell_factory R bundle\n\n"
        "This compact bundle is intended for R-side plotting/reporting without copying "
        "or densifying the full AnnData expression matrix.\n\n"
        f"- Schema: `{SCHEMA_VERSION_V1}`\n"
        f"- Source: `{config.input_h5ad}`\n"
        f"- Cells exported: `{len(cell_idx)}`\n"
        f"- Marker genes present: `{', '.join(markers_present)}`\n"
        f"- Expression values: `{EXPRESSION_SOURCE_SLOT}` as `{EXPRESSION_VALUE_SCALE}` "
        f"exported as `{EXPRESSION_EXPORT_DTYPE}`\n"
        f"- Claim guard: `{CLAIM_GUARD}`\n"
    )
    return manifest


def _export_bundle_v2(adata, cell_idx: np.ndarray, config: ExportConfig) -> dict[str, object]:
    files, markers_present, used_mtx = write_v2_files(adata, cell_idx, config)

    # Phase B: schema-version negotiation. The default emitter writes v2.1
    # (additive, backwards-compatible). Callers passing --schema-version v2
    # get the legacy literal so older readers do not see an unfamiliar
    # schema string. The on-disk file layout is identical between v2 and
    # v2.1; only the manifest schema string and the optional `extensions`
    # field differ. v2.2 adds the marker_resolutions extension slot.
    is_v2_2 = config.schema_version == "v2.2"
    is_v2_1 = config.schema_version == "v2.1" or is_v2_2
    if is_v2_2:
        schema_literal = BUNDLE_SCHEMA_V2_2
        compatible_with = SCHEMA_COMPATIBLE_WITH_V2_2
    elif config.schema_version == "v2.1":
        schema_literal = BUNDLE_SCHEMA_V2_1
        compatible_with = SCHEMA_COMPATIBLE_WITH_V2_1
    else:
        schema_literal = BUNDLE_SCHEMA_V2
        compatible_with = SCHEMA_COMPATIBLE_WITH_V2

    source_stat = config.input_h5ad.stat()
    manifest = {
        "schema_version": schema_literal,
        "compatible_with": compatible_with,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "input_h5ad": str(config.input_h5ad),
            "source_run_dir": str(config.source_run_dir) if config.source_run_dir else None,
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "input_h5ad_bytes": source_stat.st_size,
            "input_h5ad_mtime": source_stat.st_mtime,
            "input_h5ad_mtime_epoch": int(source_stat.st_mtime),
            "x_storage": type(adata.X).__name__,
        },
        "bundle": {
            "schema_name": "compact_metadata_embedding_marker_bundle_v2",
            "n_cells_exported": int(len(cell_idx)),
            "max_cells": config.max_cells,
            "seed": config.seed,
            "obs_columns_requested": list(config.obs_columns),
            "obs_columns_present": [
                col for col in config.obs_columns if col in adata.obs.columns
            ],
            "markers_requested": list(config.markers),
            "markers_present": markers_present,
            "obsm_keys": list(config.obsm_keys),
            "compression": config.compression,
            "required_files": sorted(files),
            "marker_format": "mtx_gz" if used_mtx else "parquet",
        },
        "cell_alignment": {
            "primary_index": "barcode",
            "all_files_share_primary": True,
        },
        "expression": {
            "source_slot": EXPRESSION_SOURCE_SLOT,
            "value_scale": EXPRESSION_VALUE_SCALE,
            "export_dtype": EXPRESSION_EXPORT_DTYPE,
            "intended_use": INTENDED_USE,
            "claim_guard": CLAIM_GUARD,
        },
        "files": files,
        "precision_policy": (
            "Metadata, embeddings, and selected marker expression are copied from the "
            "source AnnData. The full expression matrix is not converted to dense data. "
            "This bundle is for plotting/reporting handoff, not differential expression "
            "or new quantitative biological claims without full-object validation."
        ),
    }
    # Phase B: only the v2.1 emitter advertises an `extensions` field. v2
    # bundles must remain bit-for-bit identical in shape to pre-Phase-B
    # output, so we deliberately omit the key when emitting v2.
    if is_v2_1:
        manifest["extensions"] = {}
        # Optional protein/ADT extension. No-op unless ``include_protein`` was
        # set on the config AND the modality is actually present on the adata.
        if getattr(config, "include_protein", False):
            # Subset the (possibly backed) adata by cell_idx so the parquet
            # rows align with obs.parquet. ad.read_h5ad(..., backed='r') supports
            # boolean/integer fancy indexing returning a view, but materializing
            # via `to_memory()` is safer for the small obsm slice we need.
            sub = adata[cell_idx]
            try:
                if hasattr(sub, "to_memory"):
                    sub = sub.to_memory()
            except Exception:
                # Fall back to the view; obsm access still works.
                pass
            maybe_export_protein(
                manifest,
                config.output_dir,
                sub,
                obsm_key=config.protein_obsm_key,
                isotype_controls=config.protein_isotype_controls,
            )
            # Refresh the required_files / cell_alignment slots so manifest stays
            # internally consistent if a downstream consumer iterates files.
            manifest["bundle"]["required_files"] = sorted(manifest["files"])
        # Optional spatial extension. No-op unless ``include_spatial`` was set
        # AND the modality is actually present on the adata.
        if getattr(config, "include_spatial", False):
            sub = adata[cell_idx]
            try:
                if hasattr(sub, "to_memory"):
                    sub = sub.to_memory()
            except Exception:
                pass
            maybe_export_spatial(
                manifest,
                config.output_dir,
                sub,
                obsm_key=config.spatial_obsm_key,
                include_image_paths=config.spatial_include_image_paths,
            )
            manifest["bundle"]["required_files"] = sorted(manifest["files"])
        # Optional multimodal_obsm extension (EXPERIMENTAL). No-op unless
        # ``include_multimodal_obsm`` was set AND at least one of the named
        # obsm keys is actually present on the adata.
        if getattr(config, "include_multimodal_obsm", False):
            sub = adata[cell_idx]
            try:
                if hasattr(sub, "to_memory"):
                    sub = sub.to_memory()
            except Exception:
                pass
            maybe_export_multimodal_obsm(
                manifest,
                config.output_dir,
                sub,
                obsm_keys=config.multimodal_obsm_keys,
            )
            manifest["bundle"]["required_files"] = sorted(manifest["files"])
        # v2.2 marker_resolutions extension. Available in v2.1 bundles too when
        # explicitly requested, but the schema slot is only declared in v2.2.
        if getattr(config, "include_marker_resolutions", False):
            maybe_export_marker_resolutions(
                manifest,
                config.output_dir,
                adata,
            )
            manifest["bundle"]["required_files"] = sorted(manifest["files"])
        if getattr(config, "include_atac", False):
            sub = adata[cell_idx]
            try:
                if hasattr(sub, "to_memory"):
                    sub = sub.to_memory()
            except Exception:
                pass
            maybe_export_atac(
                manifest,
                config.output_dir,
                sub,
                lsi_obsm_key=config.atac_lsi_obsm_key,
                peaks_uns_key=config.atac_peaks_uns_key,
            )
            manifest["bundle"]["required_files"] = sorted(manifest["files"])
        if getattr(config, "include_hic", False):
            maybe_export_hic(
                manifest,
                config.output_dir,
                adata,
                max_contacts=config.hic_max_contacts,
            )
            manifest["bundle"]["required_files"] = sorted(manifest["files"])
        if getattr(config, "include_ribo", False):
            maybe_export_ribo(
                manifest,
                config.output_dir,
                adata,
            )
            manifest["bundle"]["required_files"] = sorted(manifest["files"])
    # PREC-1: write provenance.json after all parquet/mtx files are in place so
    # the bundle_sha256 it records covers the complete file set.
    _write_bundle_provenance(config.output_dir)

    manifest_json = config.output_dir / "bundle_manifest.json"
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    write_manifest_tsv(manifest, config.output_dir / "bundle_manifest.tsv")
    (config.output_dir / "README.md").write_text(
        "# singlecell_factory R bundle v2\n\n"
        "This compact bundle uses parquet + optional mtx.gz for efficient R-side "
        "plotting/reporting without densifying the full AnnData expression matrix.\n\n"
        f"- Schema: `{schema_literal}`\n"
        f"- Compatible with: `{', '.join(compatible_with)}`\n"
        f"- Source: `{config.input_h5ad}`\n"
        f"- Cells exported: `{len(cell_idx)}`\n"
        f"- Marker genes present: `{', '.join(markers_present)}`\n"
        f"- Marker format: `{'mtx_gz' if used_mtx else 'parquet'}`\n"
        f"- Expression values: `{EXPRESSION_SOURCE_SLOT}` as `{EXPRESSION_VALUE_SCALE}` "
        f"exported as `{EXPRESSION_EXPORT_DTYPE}`\n"
        f"- Claim guard: `{CLAIM_GUARD}`\n"
    )
    return manifest


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Source final_adata.h5ad")
    parser.add_argument("--output", default=None, type=Path, help="Output bundle directory (required only without --project-root)")
    parser.add_argument("--source-run-dir", type=Path, default=None, help="Source run directory")
    parser.add_argument(
        "--project-root",
        default=None,
        metavar="PATH",
        type=Path,
        help=(
            "Project directory root. When set, bundle is written to "
            "<project-root>/runs/<run-id>/python/bundle/ and --output is ignored."
        ),
    )
    parser.add_argument(
        "--run-id",
        default=None,
        metavar="STR",
        help=(
            "Run identifier (format: YYYY-MM-DDTHHMMZ-<7hex>). "
            "Auto-generated when absent. Only used when --project-root is set."
        ),
    )
    parser.add_argument("--obs-cols", default=None, help="Comma-separated obs columns to export")
    parser.add_argument("--markers", default=None, help="Comma-separated marker genes to export")
    parser.add_argument("--obsm", default=None, help="Comma-separated obsm keys to export")
    parser.add_argument("--compression", choices=("gzip", "none"), default="gzip")
    parser.add_argument("--max-cells", type=int, default=None, help="Optional deterministic subset size")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--marker-chunk-size", type=int, default=50000)
    parser.add_argument("--schema-version", choices=("v1", "v2", "v2.1", "v2.2"), default="v2.1",
                        help=("Bundle schema version. Default: v2.1 (additive "
                              "Phase B schema with an `extensions` field). "
                              "Pass 'v2.2' to enable marker_resolutions slot, "
                              "'v2' for the legacy non-extension v2 "
                              "layout, or 'v1' for the original CSV bundle."))
    parser.add_argument("--format", choices=("auto", "csv", "parquet", "mtx"), default="auto",
                        help="Output format for v2 (default: auto)")
    parser.add_argument("--include-protein", action="store_true",
                        help=("v2.1 only: also write protein.parquet and register a "
                              "`protein` extension in the manifest. No-op when the "
                              "source AnnData lacks the protein modality."))
    parser.add_argument("--protein-obsm-key", default="protein_clr",
                        help="adata.obsm key holding the (cells x proteins) CLR matrix.")
    parser.add_argument("--protein-isotype-controls", default=None,
                        help="Comma-separated isotype-control protein names (forwarded to extension entry).")
    parser.add_argument("--include-spatial", action="store_true",
                        help=("v2.1 only: also write spatial.parquet and register a "
                              "`spatial` extension in the manifest. No-op when the "
                              "source AnnData lacks adata.obsm['spatial']."))
    parser.add_argument("--spatial-obsm-key", default="spatial",
                        help="adata.obsm key holding the (cells x 2) spatial coords.")
    parser.add_argument("--no-spatial-image-paths", action="store_true",
                        help="Do not copy adata.uns['spatial']['library_id'] image-path strings into the manifest.")
    parser.add_argument("--include-multimodal-obsm", action="store_true",
                        help=("v2.1 only [EXPERIMENTAL]: also write "
                              "multimodal_obsm_<key>.parquet for each present "
                              "obsm key in --multimodal-obsm-keys and register "
                              "the `multimodal_obsm` extension."))
    parser.add_argument("--multimodal-obsm-keys", nargs="+", default=["X_wnn", "X_mofa"],
                        help=("List of adata.obsm keys to publish under the "
                              "multimodal_obsm extension. Keys not present on "
                              "the adata are silently skipped."))
    parser.add_argument("--include-marker-resolutions", action="store_true",
                        help=("v2.1/v2.2: write extensions/marker_resolutions/markers.parquet "
                              "from adata.uns['marker_db_index'] and register the "
                              "`marker_resolutions` extension. No-op when "
                              "marker_db_loader did not run or context_aware_annotation "
                              "is absent."))
    parser.add_argument("--include-atac", action="store_true",
                        help=("v2.2 only: write extensions/atac/lsi.parquet and "
                              "extensions/atac/peaks.parquet from the canonical "
                              "ATAC AnnData contract and register the `atac` extension."))
    parser.add_argument("--atac-lsi-obsm-key", default="X_lsi",
                        help="adata.obsm key holding the ATAC LSI embedding.")
    parser.add_argument("--atac-peaks-uns-key", default="atac_peaks",
                        help="adata.uns key holding ATAC peak metadata.")
    parser.add_argument("--include-hic", action="store_true",
                        help=("v2.2 only: write extensions/hic/{bins,contacts,boundaries,compartments}.parquet "
                              "from hic_ingest/hic_tad AnnData uns payloads and register the `hic` extension."))
    parser.add_argument("--hic-max-contacts", type=int, default=1_000_000,
                        help=("Maximum non-zero contacts to serialize into the compact HIC extension "
                              "(default: 1,000,000). Use coarser bins or filtering for larger matrices."))
    parser.add_argument("--include-ribo", action="store_true",
                        help=("v2.2 only: write extensions/ribo/translation_efficiency.parquet "
                              "from ribo_ingest AnnData uns payloads and register the `ribo` extension."))
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    import sys
    import warnings

    args = build_parser().parse_args(argv)

    # GOV-2 cutover semantics:
    # When `SC_REQUIRE_PROJECT_ROOT` is UNSET: a missing `--project-root` records
    # a contracted access event, emits `DeprecationWarning`, and uses explicit
    # legacy `--output`. When
    # `SC_REQUIRE_PROJECT_ROOT=1`: a missing `--project-root` is a hard error
    # (exit code 2). Warning and hard-error are mutually exclusive (no
    # double-fire). Retirement is governed by the committed deprecation contract.
    if args.project_root is None:
        if os.environ.get("SC_REQUIRE_PROJECT_ROOT") == "1":
            print(
                "ERROR: --project-root is required (SC_REQUIRE_PROJECT_ROOT=1 is set). "
                "Pass --project-root or unset the env var.",
                file=sys.stderr,
            )
            sys.exit(2)

    # Resolve effective output directory.
    _run_dir = None
    if args.project_root is not None:
        import sys as _sys
        _factory_root = Path(__file__).resolve().parent.parent
        _project_paths_mod = _factory_root / "workflow" / "modular"
        if str(_project_paths_mod) not in _sys.path:
            _sys.path.insert(0, str(_factory_root))
        from workflow.modular.project_paths import resolve_run_dir, bundle_dir
        _run_dir = resolve_run_dir(args.project_root, run_id=args.run_id)
        effective_output = bundle_dir(_run_dir)
    else:
        if args.output is None:
            print(
                "ERROR: --output is required when --project-root is not provided.",
                file=sys.stderr,
            )
            sys.exit(2)
        _factory_root = Path(__file__).resolve().parent.parent
        if str(_factory_root) not in sys.path:
            sys.path.insert(0, str(_factory_root))
        from workflow.modular.legacy_output import (
            legacy_output_warning_message,
            record_legacy_output_access,
        )
        try:
            telemetry_path = record_legacy_output_access(
                output_dir=Path(args.output),
                project="bundle_export",
            )
        except (OSError, ValueError) as exc:
            print(
                "ERROR: legacy output telemetry could not be recorded; the "
                "deprecation contract blocks this compatibility export: "
                f"{exc}",
                file=sys.stderr,
            )
            raise SystemExit(2) from exc
        warning_message = legacy_output_warning_message(telemetry_path)
        warnings.warn(
            warning_message,
            DeprecationWarning,
            stacklevel=2,
        )
        print(f"DeprecationWarning: {warning_message}", file=sys.stderr)
        effective_output = args.output

    config = ExportConfig(
        input_h5ad=args.input,
        output_dir=effective_output,
        source_run_dir=args.source_run_dir,
        obs_columns=parse_csv_list(args.obs_cols, DEFAULT_OBS_COLUMNS),
        markers=parse_csv_list(args.markers, DEFAULT_MARKERS),
        obsm_keys=parse_csv_list(args.obsm, DEFAULT_OBSM),
        compression=args.compression,
        max_cells=args.max_cells,
        seed=args.seed,
        marker_chunk_size=args.marker_chunk_size,
        schema_version=args.schema_version,
        format=args.format,
        include_protein=bool(args.include_protein),
        protein_obsm_key=args.protein_obsm_key,
        protein_isotype_controls=parse_csv_list(args.protein_isotype_controls, ()),
        include_spatial=bool(args.include_spatial),
        spatial_obsm_key=args.spatial_obsm_key,
        spatial_include_image_paths=not bool(args.no_spatial_image_paths),
        include_multimodal_obsm=bool(args.include_multimodal_obsm),
        multimodal_obsm_keys=tuple(args.multimodal_obsm_keys),
        include_marker_resolutions=bool(args.include_marker_resolutions),
        include_atac=bool(args.include_atac),
        atac_lsi_obsm_key=args.atac_lsi_obsm_key,
        atac_peaks_uns_key=args.atac_peaks_uns_key,
        include_hic=bool(args.include_hic),
        hic_max_contacts=int(args.hic_max_contacts),
        include_ribo=bool(args.include_ribo),
    )
    manifest = export_bundle(config)
    result = {"output_dir": str(config.output_dir), "manifest": manifest["schema_version"]}
    if _run_dir is not None:
        digest = _backfill_project_manifest_bundle_sha256(_run_dir, config.output_dir)
        if digest:
            result["bundle_sha256"] = digest
    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
