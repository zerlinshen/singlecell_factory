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
from dataclasses import dataclass, field
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
SCHEMA_VERSION_V2 = "singlecell_r_bundle_v2"
SCHEMA_COMPATIBLE_WITH_V2 = ["singlecell_r_bundle_v1"]
EXPRESSION_SOURCE_SLOT = "X"
EXPRESSION_VALUE_SCALE = "source_X_as_stored"
EXPRESSION_EXPORT_DTYPE = "float32"
INTENDED_USE = "plotting_and_visual_summary_only"
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
    schema_version: str = "v1"
    format: str = "auto"


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
    if not _HAVE_PYARROW:
        df.to_parquet(path, engine="fastparquet", index=True)
        return
    table = pa.Table.from_pandas(df, preserve_index=True)
    pq.write_table(table, path, compression="snappy")


def write_obs_parquet(adata, cell_idx: np.ndarray, config: ExportConfig) -> tuple[Path, list[str]]:
    present = [col for col in config.obs_columns if col in adata.obs.columns]
    obs = adata.obs.iloc[cell_idx][present].copy()
    # Reset string/object columns to preserve categorical info
    for col in obs.columns:
        if hasattr(obs[col], "cat"):
            obs[col] = obs[col].astype(str).astype("category")
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
    if fmt == "mtx":
        return True
    if fmt == "csv" or fmt == "parquet":
        return False
    # auto: use mtx when marker count >= threshold
    return len(markers_present) >= MTX_AUTO_THRESHOLD


def write_marker_expr_mtx(adata, cell_idx: np.ndarray, markers_present: list[str],
                          config: ExportConfig) -> tuple[Path, dict[str, object]]:
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
    rec = file_record(mtx_path, config.output_dir, n_rows, n_cols,
                      {"format": "mtx_gz", "sparse": True, "nnz": nnz})
    return mtx_path, rec


def write_marker_expr_parquet(adata, cell_idx: np.ndarray, markers_present: list[str],
                              config: ExportConfig) -> tuple[Path, dict[str, object]]:
    """Write marker expression as parquet (dense float32)."""
    cells = adata.obs_names[cell_idx].astype(str)
    first = True
    chunks = []
    for start in range(0, len(cell_idx), config.marker_chunk_size):
        stop = min(start + config.marker_chunk_size, len(cell_idx))
        rows = cell_idx[start:stop]
        block = materialize_small_matrix(adata[rows, markers_present].X)
        block = block.astype(np.float32, copy=False)
        chunks.append(pd.DataFrame(block, index=cells[start:stop], columns=markers_present))
    df = pd.concat(chunks, axis=0)
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
        _, marker_rec = write_marker_expr_mtx(adata, cell_idx, markers_present, config)
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
# Main export entry point
# ---------------------------------------------------------------------------

def export_bundle(config: ExportConfig) -> dict[str, object]:
    if config.compression not in {"gzip", "none"}:
        raise ValueError("--compression must be 'gzip' or 'none'.")
    if config.marker_chunk_size <= 0:
        raise ValueError("--marker-chunk-size must be positive.")
    if not config.input_h5ad.exists():
        raise FileNotFoundError(f"Input .h5ad not found: {config.input_h5ad}")
    if config.schema_version not in ("v1", "v2"):
        raise ValueError("--schema-version must be 'v1' or 'v2'.")
    if config.format not in ("auto", "csv", "parquet", "mtx"):
        raise ValueError("--format must be 'auto', 'csv', 'parquet', or 'mtx'.")

    config.output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(config.input_h5ad, backed="r")
    try:
        cell_idx = choose_cells(adata.n_obs, config.max_cells, config.seed)

        if config.schema_version == "v2":
            return _export_bundle_v2(adata, cell_idx, config)
        else:
            return _export_bundle_v1(adata, cell_idx, config)
    finally:
        adata.file.close()


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

    source_stat = config.input_h5ad.stat()
    manifest = {
        "schema_version": SCHEMA_VERSION_V2,
        "compatible_with": SCHEMA_COMPATIBLE_WITH_V2,
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
    manifest_json = config.output_dir / "bundle_manifest.json"
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    write_manifest_tsv(manifest, config.output_dir / "bundle_manifest.tsv")
    (config.output_dir / "README.md").write_text(
        "# singlecell_factory R bundle v2\n\n"
        "This compact bundle uses parquet + optional mtx.gz for efficient R-side "
        "plotting/reporting without densifying the full AnnData expression matrix.\n\n"
        f"- Schema: `{SCHEMA_VERSION_V2}`\n"
        f"- Compatible with: `{', '.join(SCHEMA_COMPATIBLE_WITH_V2)}`\n"
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
    parser.add_argument("--output", required=True, type=Path, help="Output bundle directory")
    parser.add_argument("--source-run-dir", type=Path, default=None, help="Source run directory")
    parser.add_argument("--obs-cols", default=None, help="Comma-separated obs columns to export")
    parser.add_argument("--markers", default=None, help="Comma-separated marker genes to export")
    parser.add_argument("--obsm", default=None, help="Comma-separated obsm keys to export")
    parser.add_argument("--compression", choices=("gzip", "none"), default="gzip")
    parser.add_argument("--max-cells", type=int, default=None, help="Optional deterministic subset size")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--marker-chunk-size", type=int, default=50000)
    parser.add_argument("--schema-version", choices=("v1", "v2"), default="v1",
                        help="Bundle schema version (default: v1 for backward compat)")
    parser.add_argument("--format", choices=("auto", "csv", "parquet", "mtx"), default="auto",
                        help="Output format for v2 (default: auto)")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = ExportConfig(
        input_h5ad=args.input,
        output_dir=args.output,
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
    )
    manifest = export_bundle(config)
    print(json.dumps({"output_dir": str(config.output_dir), "manifest": manifest["schema_version"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
