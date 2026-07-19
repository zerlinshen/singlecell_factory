"""Strict render-only artifact contract for the public Figure 3A velocity lane."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

SCHEMA_NAME = "singlecell_velocity_figure3a_render"
SCHEMA_VERSION = "1.0"
CLAIM_CLASS = "exploratory"
CELL_ARTIFACT_NAME = "figure3a_velocity_cells.csv"
MARKER_ARTIFACT_NAME = "figure3a_velocity_marker_trends.csv"
MANIFEST_NAME = "figure3a_velocity_render_manifest.json"

CELL_COLUMNS = (
    "cell_id",
    "umap1",
    "umap2",
    "velocity_umap1",
    "velocity_umap2",
    "velocity_length",
    "velocity_pseudotime",
    "cell_type",
)
MARKER_COLUMNS = (
    "marker",
    "pseudotime_bin_index",
    "pseudotime_bin_left",
    "pseudotime_bin_right",
    "pseudotime_bin_midpoint",
    "median_normalized_expression",
    "n_cells",
)
REQUIRED_PARAMETER_KEYS = {
    "max_cells",
    "top_genes",
    "min_shared_counts",
    "moments_n_pcs",
    "moments_n_neighbors",
    "n_pseudotime_bins",
    "velocity_mode",
    "basis",
    "input_overlap_cells",
    "computed_cells",
    "computed_genes",
}
REQUIRED_SOFTWARE_KEYS = {"python", "anndata", "numpy", "pandas", "scanpy", "scvelo"}
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def sha256_file(path: Path) -> str:
    """Return the SHA256 digest of a file without loading it into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def choose_cells(obs: pd.DataFrame, max_cells: int, seed: int) -> list[str]:
    """Choose a deterministic, rare-cell-aware subset by cell identifier."""
    if max_cells <= 0:
        raise ValueError("max_cells must be positive")
    if not obs.index.is_unique:
        raise ValueError("cell identifiers must be unique")

    cell_ids = obs.index.astype(str)
    if len(obs) <= max_cells:
        return sorted(cell_ids.tolist())

    rng = np.random.default_rng(seed)
    if "cell_type" not in obs.columns:
        picked = rng.choice(len(obs), size=max_cells, replace=False)
        return sorted(cell_ids[picked].tolist())

    chosen: list[str] = []
    groups = obs.groupby("cell_type", observed=True, sort=True).indices
    floor = min(250, max(30, max_cells // max(1, len(groups) * 6)))
    remaining = max_cells
    ordered_groups = sorted(groups.items(), key=lambda item: (len(item[1]), str(item[0])))
    for _, positions in ordered_groups:
        positions = np.asarray(positions)
        take = min(len(positions), floor, remaining)
        if take:
            selected = rng.choice(positions, size=take, replace=False)
            chosen.extend(cell_ids[selected].tolist())
            remaining -= take

    if remaining:
        chosen_set = set(chosen)
        pool = np.asarray([i for i, cell_id in enumerate(cell_ids) if cell_id not in chosen_set])
        if len(pool):
            selected = rng.choice(pool, size=min(remaining, len(pool)), replace=False)
            chosen.extend(cell_ids[selected].tolist())
    return sorted(chosen[:max_cells])


def _expression_vector(adata: Any, marker: str) -> np.ndarray:
    values = adata[:, marker].X
    if hasattr(values, "toarray"):
        values = values.toarray()
    return np.asarray(values).reshape(-1).astype(float, copy=False)


def build_render_tables(
    adata: Any,
    markers: Sequence[str],
    *,
    n_pseudotime_bins: int = 20,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Build strict cell-vector and marker-trend tables from computed scVelo state."""
    required_obs = {"velocity_length", "velocity_pseudotime", "cell_type"}
    missing_obs = sorted(required_obs.difference(adata.obs.columns))
    if missing_obs:
        raise ValueError(f"computed velocity object missing obs fields: {missing_obs}")
    for key in ("X_umap", "velocity_umap"):
        if key not in adata.obsm:
            raise ValueError(f"computed velocity object missing obsm[{key!r}]")
    if n_pseudotime_bins < 2:
        raise ValueError("n_pseudotime_bins must be at least 2")

    umap = np.asarray(adata.obsm["X_umap"], dtype=float)
    velocity_umap = np.asarray(adata.obsm["velocity_umap"], dtype=float)
    if umap.shape != (adata.n_obs, 2) or velocity_umap.shape != (adata.n_obs, 2):
        raise ValueError("X_umap and velocity_umap must both have shape (n_cells, 2)")

    cell_table = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "umap1": umap[:, 0],
            "umap2": umap[:, 1],
            "velocity_umap1": velocity_umap[:, 0],
            "velocity_umap2": velocity_umap[:, 1],
            "velocity_length": np.asarray(adata.obs["velocity_length"], dtype=float),
            "velocity_pseudotime": np.asarray(adata.obs["velocity_pseudotime"], dtype=float),
            "cell_type": adata.obs["cell_type"].astype(str).to_numpy(),
        },
        columns=CELL_COLUMNS,
    ).sort_values("cell_id", kind="stable", ignore_index=True)

    pseudotime = np.asarray(adata.obs["velocity_pseudotime"], dtype=float)
    if not np.isfinite(pseudotime).all():
        raise ValueError("velocity_pseudotime must be finite for every rendered cell")
    pt_min = float(pseudotime.min())
    pt_max = float(pseudotime.max())
    if not pt_max > pt_min:
        raise ValueError("velocity_pseudotime must span a non-zero range")
    bin_edges = np.linspace(pt_min, pt_max, n_pseudotime_bins + 1)

    marker_rows: list[dict[str, Any]] = []
    missing_markers: list[str] = []
    for marker in markers:
        if marker not in adata.var_names:
            missing_markers.append(marker)
            continue
        expression = _expression_vector(adata, marker)
        valid = np.isfinite(expression)
        for bin_index, (left, right) in enumerate(zip(bin_edges[:-1], bin_edges[1:], strict=True)):
            in_bin = valid & (pseudotime >= left)
            in_bin &= (
                pseudotime <= right if bin_index == n_pseudotime_bins - 1 else pseudotime < right
            )
            n_cells = int(in_bin.sum())
            if not n_cells:
                continue
            marker_rows.append(
                {
                    "marker": marker,
                    "pseudotime_bin_index": bin_index,
                    "pseudotime_bin_left": float(left),
                    "pseudotime_bin_right": float(right),
                    "pseudotime_bin_midpoint": float((left + right) / 2.0),
                    "median_normalized_expression": float(np.median(expression[in_bin])),
                    "n_cells": n_cells,
                }
            )

    marker_table = pd.DataFrame(marker_rows, columns=MARKER_COLUMNS)
    if marker_table.empty:
        raise ValueError("none of the requested marker genes produced trend rows")
    marker_table = marker_table.sort_values(
        ["marker", "pseudotime_bin_index"], kind="stable", ignore_index=True
    )
    validate_render_tables(cell_table, marker_table)
    return cell_table, marker_table, missing_markers


def validate_render_tables(cell_table: pd.DataFrame, marker_table: pd.DataFrame) -> None:
    """Fail loudly unless both render inputs match the named v1 schema."""
    if tuple(cell_table.columns) != CELL_COLUMNS:
        raise ValueError(f"cell artifact columns must be exactly {list(CELL_COLUMNS)}")
    if tuple(marker_table.columns) != MARKER_COLUMNS:
        raise ValueError(f"marker artifact columns must be exactly {list(MARKER_COLUMNS)}")
    if cell_table.empty or marker_table.empty:
        raise ValueError("render artifacts must not be empty")
    if cell_table["cell_id"].duplicated().any():
        raise ValueError("cell_id must be unique")
    if cell_table["cell_id"].isna().any() or (cell_table["cell_id"].astype(str) == "").any():
        raise ValueError("cell_id must be non-empty")
    if cell_table["cell_type"].isna().any() or (cell_table["cell_type"].astype(str) == "").any():
        raise ValueError("cell_type must be non-empty")

    cell_numeric = list(CELL_COLUMNS[1:7])
    if not np.isfinite(cell_table[cell_numeric].to_numpy(dtype=float)).all():
        raise ValueError("all cell render coordinates and statistics must be finite")
    if (cell_table["velocity_length"] < 0).any():
        raise ValueError("velocity_length must be non-negative")
    if not cell_table["velocity_pseudotime"].between(0.0, 1.0).all():
        raise ValueError("velocity_pseudotime must be within [0, 1]")

    marker_numeric = list(MARKER_COLUMNS[1:6])
    if not np.isfinite(marker_table[marker_numeric].to_numpy(dtype=float)).all():
        raise ValueError("all marker trend fields must be finite")
    if (marker_table["n_cells"] <= 0).any():
        raise ValueError("marker trend n_cells must be positive")
    if (marker_table["pseudotime_bin_index"] < 0).any():
        raise ValueError("pseudotime_bin_index must be non-negative")
    bin_indices = marker_table["pseudotime_bin_index"].to_numpy(dtype=float)
    if not np.equal(bin_indices, np.floor(bin_indices)).all():
        raise ValueError("pseudotime_bin_index must contain integers")
    for column in (
        "pseudotime_bin_left",
        "pseudotime_bin_right",
        "pseudotime_bin_midpoint",
    ):
        if not marker_table[column].between(0.0, 1.0).all():
            raise ValueError(f"{column} must be within [0, 1]")
    if (
        not (marker_table["pseudotime_bin_left"] <= marker_table["pseudotime_bin_midpoint"]).all()
        or not (
            marker_table["pseudotime_bin_midpoint"] <= marker_table["pseudotime_bin_right"]
        ).all()
    ):
        raise ValueError("marker pseudotime bin boundaries are inconsistent")


def _file_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _reproducibility_payload(manifest: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_name": manifest["schema_name"],
        "schema_version": manifest["schema_version"],
        "source_sha256": {name: record["sha256"] for name, record in manifest["sources"].items()},
        "producer_sha256": manifest["producer"]["sha256"],
        "parameters": manifest["parameters"],
        "seed": manifest["seed"],
        "software_versions": manifest["software_versions"],
        "artifact_sha256": {
            name: record["sha256"] for name, record in manifest["artifacts"].items()
        },
    }


def _reproducibility_key(manifest: Mapping[str, Any]) -> str:
    payload = _reproducibility_payload(manifest)
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _validate_file_record(record: Any, *, label: str) -> None:
    if not isinstance(record, Mapping):
        raise ValueError(f"{label} must be a mapping")
    if not isinstance(record.get("path"), str) or not record["path"].strip():
        raise ValueError(f"{label}.path must be nonblank")
    if not isinstance(record.get("sha256"), str) or not _SHA256_PATTERN.fullmatch(record["sha256"]):
        raise ValueError(f"{label}.sha256 must be lowercase 64-hex")
    if not isinstance(record.get("size_bytes"), int) or record["size_bytes"] <= 0:
        raise ValueError(f"{label}.size_bytes must be a positive integer")


def write_render_bundle(
    adata: Any,
    output_dir: Path,
    *,
    markers: Sequence[str],
    source_files: Mapping[str, Path],
    parameters: Mapping[str, Any],
    seed: int,
    software_versions: Mapping[str, str | None],
    producer_path: Path,
    truth_boundary: str,
    n_pseudotime_bins: int = 20,
) -> dict[str, Any]:
    """Write an immutable, hash-linked render bundle and return its manifest."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = [
        output_dir / CELL_ARTIFACT_NAME,
        output_dir / MARKER_ARTIFACT_NAME,
        output_dir / MANIFEST_NAME,
    ]
    existing = [str(path) for path in targets if path.exists()]
    if existing:
        raise FileExistsError(f"immutable render artifact target already exists: {existing}")

    cell_table, marker_table, missing_markers = build_render_tables(
        adata, markers, n_pseudotime_bins=n_pseudotime_bins
    )
    cell_path, marker_path, manifest_path = targets
    cell_table.to_csv(cell_path, index=False, float_format="%.10g", lineterminator="\n")
    marker_table.to_csv(marker_path, index=False, float_format="%.10g", lineterminator="\n")

    source_records = {name: _file_record(Path(path)) for name, path in sorted(source_files.items())}
    producer_record = _file_record(Path(producer_path))
    artifact_records = {
        "cells": {
            "filename": CELL_ARTIFACT_NAME,
            "sha256": sha256_file(cell_path),
            "size_bytes": cell_path.stat().st_size,
            "rows": len(cell_table),
            "columns": list(CELL_COLUMNS),
        },
        "marker_trends": {
            "filename": MARKER_ARTIFACT_NAME,
            "sha256": sha256_file(marker_path),
            "size_bytes": marker_path.stat().st_size,
            "rows": len(marker_table),
            "columns": list(MARKER_COLUMNS),
        },
    }
    manifest = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "status": "complete",
        "claim_class": CLAIM_CLASS,
        "truth_boundary": truth_boundary,
        "seed": int(seed),
        "parameters": dict(parameters),
        "sources": source_records,
        "producer": producer_record,
        "software_versions": dict(software_versions),
        "markers": {
            "requested": list(markers),
            "missing": missing_markers,
            "rendered": sorted(marker_table["marker"].unique().tolist()),
        },
        "artifacts": artifact_records,
    }
    manifest["reproducibility_key_sha256"] = _reproducibility_key(manifest)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    validate_render_bundle(output_dir)
    return manifest


def validate_render_bundle(bundle_dir: Path) -> dict[str, Any]:
    """Validate schema, exact columns, row counts, and SHA256 integrity."""
    bundle_dir = Path(bundle_dir)
    manifest_path = bundle_dir / MANIFEST_NAME
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_name") != SCHEMA_NAME:
        raise ValueError(
            f"unsupported velocity render schema name: {manifest.get('schema_name')!r}"
        )
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"unsupported velocity render schema version: {manifest.get('schema_version')!r}"
        )
    if manifest.get("status") != "complete" or manifest.get("claim_class") != CLAIM_CLASS:
        raise ValueError("velocity render manifest must be complete and exploratory")
    if (
        not isinstance(manifest.get("truth_boundary"), str)
        or not manifest["truth_boundary"].strip()
    ):
        raise ValueError("velocity render manifest truth_boundary must be nonblank")
    if not isinstance(manifest.get("seed"), int) or isinstance(manifest["seed"], bool):
        raise ValueError("velocity render manifest seed must be an integer")

    parameters = manifest.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ValueError("velocity render manifest parameters must be a mapping")
    missing_parameters = sorted(REQUIRED_PARAMETER_KEYS.difference(parameters))
    if missing_parameters:
        raise ValueError(f"velocity render manifest missing parameters: {missing_parameters}")

    software_versions = manifest.get("software_versions")
    if not isinstance(software_versions, Mapping):
        raise ValueError("velocity render manifest software_versions must be a mapping")
    missing_software = sorted(REQUIRED_SOFTWARE_KEYS.difference(software_versions))
    if missing_software:
        raise ValueError(f"velocity render manifest missing software versions: {missing_software}")
    if any(
        not isinstance(software_versions[key], str) or not software_versions[key].strip()
        for key in REQUIRED_SOFTWARE_KEYS
    ):
        raise ValueError("velocity render manifest software versions must be nonblank strings")

    sources = manifest.get("sources")
    if not isinstance(sources, Mapping) or not sources:
        raise ValueError("velocity render manifest sources must be a nonempty mapping")
    for name, record in sources.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("velocity render manifest source names must be nonblank")
        _validate_file_record(record, label=f"sources.{name}")
    _validate_file_record(manifest.get("producer"), label="producer")

    marker_record = manifest.get("markers")
    if not isinstance(marker_record, Mapping):
        raise ValueError("velocity render manifest markers must be a mapping")
    for key in ("requested", "missing", "rendered"):
        values = marker_record.get(key)
        if not isinstance(values, list) or any(
            not isinstance(value, str) or not value.strip() for value in values
        ):
            raise ValueError(f"velocity render manifest markers.{key} must be a string list")
    if not marker_record["requested"] or not marker_record["rendered"]:
        raise ValueError("velocity render manifest requested/rendered markers must be nonempty")
    if not set(marker_record["rendered"]).issubset(marker_record["requested"]):
        raise ValueError("rendered markers must be a subset of requested markers")

    artifact_manifest = manifest.get("artifacts")
    if not isinstance(artifact_manifest, Mapping):
        raise ValueError("velocity render manifest artifacts must be a mapping")
    if set(artifact_manifest) != {"cells", "marker_trends"}:
        raise ValueError(
            "velocity render manifest artifacts must contain exactly cells and marker_trends"
        )

    cells = pd.read_csv(bundle_dir / CELL_ARTIFACT_NAME)
    markers = pd.read_csv(bundle_dir / MARKER_ARTIFACT_NAME)
    validate_render_tables(cells, markers)
    for key, table, expected_filename in (
        ("cells", cells, CELL_ARTIFACT_NAME),
        ("marker_trends", markers, MARKER_ARTIFACT_NAME),
    ):
        record = artifact_manifest[key]
        if not isinstance(record, Mapping):
            raise ValueError(f"manifest {key} artifact record must be a mapping")
        if record.get("filename") != expected_filename:
            raise ValueError(f"manifest {key} filename must be {expected_filename}")
        path = bundle_dir / expected_filename
        if not isinstance(record.get("sha256"), str) or not _SHA256_PATTERN.fullmatch(
            record["sha256"]
        ):
            raise ValueError(f"manifest {key} sha256 must be lowercase 64-hex")
        if record.get("sha256") != sha256_file(path):
            raise ValueError(f"manifest {key} sha256 does not match on-disk artifact")
        if not isinstance(record.get("size_bytes"), int) or record["size_bytes"] <= 0:
            raise ValueError(f"manifest {key} size_bytes must be a positive integer")
        if record["size_bytes"] != path.stat().st_size:
            raise ValueError(f"manifest {key} size_bytes does not match artifact")
        if not isinstance(record.get("rows"), int) or record["rows"] <= 0:
            raise ValueError(f"manifest {key} rows must be a positive integer")
        if record.get("rows") != len(table) or record.get("columns") != list(table.columns):
            raise ValueError(f"manifest {key} shape/schema metadata does not match artifact")

    expected_key = _reproducibility_key(manifest)
    actual_key = manifest.get("reproducibility_key_sha256")
    if not isinstance(actual_key, str) or not _SHA256_PATTERN.fullmatch(actual_key):
        raise ValueError("reproducibility_key_sha256 must be lowercase 64-hex")
    if actual_key != expected_key:
        raise ValueError("reproducibility_key_sha256 does not match manifest fields")
    return manifest
