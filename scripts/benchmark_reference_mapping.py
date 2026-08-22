#!/usr/bin/env python3
"""Process-isolated CPU/GPU benchmark for reference KNN OOD mapping.

This is one bounded benchmark entrypoint, not a workflow framework.  Its
coordinator freezes every consumed input exactly once, then starts CPU and GPU
children with the same interpreter.  The children can only read hash-bound
serialized inputs, so their timing, RSS, output determinism, and parity claims
are attributable to a single backend lane.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import logging
import os
import resource
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics import f1_score

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from workflow.modular._reference_mapping import (
    ASSIGNMENT_STATUS_ACCEPTED,
    calibrate_reference_ood_threshold,
    compute_parity_metrics,
    map_knn_reference,
)
from workflow.modular.manifest_writer import factory_git_state, write_manifest
from workflow.modular.project_paths import resolve_run_dir

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("benchmark_reference_mapping")

DEFAULT_TREVINO_SOURCE = Path(
    "/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/"
    "python/wave5_trevino_public_rna_20260518_040458/final_adata.h5ad"
)
DEFAULT_TREVINO_SOURCE_SHA256 = "44b3cf6f98ef890ad9a2723f275f8ac63b51148f5a685e38cfe802d561282ad7"
_CONSUMED_OUTPUT_COLUMNS = (
    "reference_predicted_label",
    "reference_confidence",
    "reference_distance",
    "reference_assignment_status",
    "reference_cell_type",
    "reference_ood",
)
_WARMUP_ROWS = 256
_MEASURED_REPETITIONS = 3
_WITHIN_LANE_DISTANCE_RTOL = 1e-6
_WITHIN_LANE_DISTANCE_ATOL = 1e-6
_EXACT_SCIENTIFIC_OUTPUT_COLUMNS = tuple(
    column for column in _CONSUMED_OUTPUT_COLUMNS if column != "reference_distance"
)


def compute_sha256(file_path: Path) -> str:
    digest = hashlib.sha256()
    with Path(file_path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json_sha256(payload: Any) -> str:
    return _sha256_bytes(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8"))


def _atomic_write_json(file_path: Path, payload: dict[str, Any]) -> None:
    destination = Path(file_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, destination)


def _atomic_write_text(file_path: Path, text: str) -> None:
    destination = Path(file_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, destination)


def _atomic_write_dataframe(file_path: Path, frame: pd.DataFrame) -> None:
    destination = Path(file_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    frame.to_csv(temporary, lineterminator="\n")
    os.replace(temporary, destination)


def _write_lines(file_path: Path, values: list[str]) -> None:
    _atomic_write_text(file_path, "\n".join(values) + "\n")


def _read_lines(file_path: Path) -> list[str]:
    return Path(file_path).read_text(encoding="utf-8").splitlines()


def _file_record(file_path: Path) -> dict[str, Any]:
    return {"sha256": compute_sha256(file_path), "size_bytes": int(file_path.stat().st_size)}


def get_peak_rss_mb() -> float:
    """Return Linux process ru_maxrss in MiB (child-lane attributable)."""
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0)


def select_reference_hvg_genes(
    adata: ad.AnnData,
    reference_cells: list[str],
    max_genes: int,
) -> list[str]:
    """Freeze the upstream retained feature set without recomputing HVGs.

    ``highly_variable`` belongs to the source object, not this benchmark.  The
    selection remains deterministic by source feature order and is capped at
    3,000.  It is deliberately recorded as *not independently query-blind*:
    the upstream retained feature set may have seen cells outside the frozen
    benchmark reference split.
    """
    if max_genes < 1:
        raise ValueError(f"max_genes must be >= 1, got {max_genes}.")
    if not reference_cells:
        raise ValueError("Reference split must contain at least one cell.")
    if np.any(adata.obs_names.get_indexer(reference_cells) < 0):
        raise ValueError("Reference split contains cell IDs absent from the source AnnData.")
    if "highly_variable" not in adata.var.columns:
        raise ValueError(
            "Source AnnData must contain an upstream 'highly_variable' field; "
            "the benchmark must not recompute Seurat HVGs."
        )
    upstream_mask = adata.var["highly_variable"].fillna(False).astype(bool).to_numpy()
    selected = adata.var_names.astype(str)[upstream_mask].tolist()
    if not selected:
        raise ValueError("Source AnnData highly_variable field contains no selected features.")
    return selected[: min(int(max_genes), 3000)]


def prepare_trevino_split(
    source_h5ad: Path,
    out_dir: Path,
    *,
    expected_source_sha256: str = DEFAULT_TREVINO_SOURCE_SHA256,
    max_ref_cells: int = 15000,
    max_query_cells: int = 6000,
    max_genes: int = 3000,
    ood_label: str = "Microglia",
    seed: int = 42,
) -> dict[str, Any]:
    """Create a sample-disjoint Trevino split and serialize its consumed inputs."""
    source_h5ad = Path(source_h5ad)
    if not source_h5ad.is_file():
        raise FileNotFoundError(f"Trevino source H5AD is absent: {source_h5ad}")
    source_sha256 = compute_sha256(source_h5ad)
    if expected_source_sha256 and source_sha256 != expected_source_sha256:
        raise RuntimeError(
            "Trevino source H5AD SHA-256 changed from the frozen benchmark contract: "
            f"expected {expected_source_sha256}, found {source_sha256}."
        )

    out_dir = Path(out_dir)
    frozen_dir = out_dir / "frozen_inputs"
    frozen_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(source_h5ad)
    required_obs = {"Tissue.ID", "sample", "cell_type"}
    missing_obs = sorted(required_obs - set(adata.obs.columns))
    if missing_obs:
        raise ValueError(f"Trevino source is missing benchmark columns: {missing_obs}")

    query_mask = adata.obs["Tissue.ID"].astype(str) == "HFT3"
    reference_tissues = {"HFT7", "HFT5", "HFT6"}
    reference_mask = adata.obs["Tissue.ID"].astype(str).isin(reference_tissues)
    query_obs = adata.obs.loc[query_mask]
    reference_obs = adata.obs.loc[reference_mask]
    if query_obs.empty or reference_obs.empty:
        raise ValueError("Trevino source produced an empty HFT3 query or HFT7/HFT5/HFT6 reference split.")

    reference_obs = reference_obs.loc[reference_obs["cell_type"].astype(str) != ood_label]
    if reference_obs.empty:
        raise ValueError(f"Purging OOD label {ood_label!r} emptied the Trevino reference split.")
    rng = np.random.default_rng(seed)

    reference_cells: list[str] = []
    max_per_reference_label = int(np.ceil(max_ref_cells / max(reference_obs["cell_type"].nunique(), 1)))
    for _, group in reference_obs.groupby("cell_type", sort=True, observed=True):
        cell_ids = np.asarray(sorted(group.index.astype(str).tolist()), dtype=object)
        rng.shuffle(cell_ids)
        reference_cells.extend(cell_ids[:max_per_reference_label].tolist())
    if len(reference_cells) > max_ref_cells:
        reference_cells = rng.choice(np.asarray(reference_cells, dtype=object), size=max_ref_cells, replace=False).tolist()
    reference_cells = sorted(str(cell) for cell in reference_cells)

    ood_query_cells = sorted(
        query_obs.index[query_obs["cell_type"].astype(str) == ood_label].astype(str).tolist()
    )
    if not ood_query_cells:
        raise ValueError(f"Trevino HFT3 query has no held-out OOD label {ood_label!r}.")
    if len(ood_query_cells) > max_query_cells:
        raise ValueError(
            f"Held-out OOD label {ood_label!r} has {len(ood_query_cells)} cells, exceeding max_query_cells={max_query_cells}."
        )
    remaining_query_slots = max_query_cells - len(ood_query_cells)
    known_query = query_obs.loc[query_obs["cell_type"].astype(str) != ood_label]
    known_selected: list[str] = []
    max_per_known_label = int(np.ceil(remaining_query_slots / max(known_query["cell_type"].nunique(), 1)))
    for _, group in known_query.groupby("cell_type", sort=True, observed=True):
        cell_ids = np.asarray(sorted(group.index.astype(str).tolist()), dtype=object)
        rng.shuffle(cell_ids)
        known_selected.extend(cell_ids[:max_per_known_label].tolist())
    if len(known_selected) > remaining_query_slots:
        known_selected = rng.choice(np.asarray(known_selected, dtype=object), size=remaining_query_slots, replace=False).tolist()
    query_cells = sorted(str(cell) for cell in [*ood_query_cells, *known_selected])
    if len(query_cells) < _WARMUP_ROWS:
        raise ValueError(f"Benchmark requires at least {_WARMUP_ROWS} query cells, found {len(query_cells)}.")

    source_hvg_count = int(adata.var["highly_variable"].fillna(False).astype(bool).sum())
    selected_genes = select_reference_hvg_genes(adata, reference_cells, max_genes=max_genes)
    ref_adata = adata[reference_cells, selected_genes].copy()
    query_adata = adata[query_cells, selected_genes].copy()
    del adata

    if set(ref_adata.obs["sample"].astype(str)) & set(query_adata.obs["sample"].astype(str)):
        raise RuntimeError("Reference and query sample sets overlap; the sample-disjoint benchmark contract failed.")
    if ood_label in set(ref_adata.obs["cell_type"].astype(str)):
        raise RuntimeError(f"Held-out OOD label {ood_label!r} remains in the reference split.")
    ref_adata.X = ref_adata.X.tocsr() if sp.issparse(ref_adata.X) else sp.csr_matrix(ref_adata.X)
    query_adata.X = query_adata.X.tocsr() if sp.issparse(query_adata.X) else sp.csr_matrix(query_adata.X)

    file_paths = {
        "reference_matrix.npz": frozen_dir / "reference_matrix.npz",
        "query_matrix.npz": frozen_dir / "query_matrix.npz",
        "reference_cells.txt": frozen_dir / "reference_cells.txt",
        "query_cells.txt": frozen_dir / "query_cells.txt",
        "reference_labels.txt": frozen_dir / "reference_labels.txt",
        "query_proxy_labels.txt": frozen_dir / "query_proxy_labels.txt",
        "genes.txt": frozen_dir / "genes.txt",
    }
    sp.save_npz(file_paths["reference_matrix.npz"], ref_adata.X)
    sp.save_npz(file_paths["query_matrix.npz"], query_adata.X)
    _write_lines(file_paths["reference_cells.txt"], ref_adata.obs_names.astype(str).tolist())
    _write_lines(file_paths["query_cells.txt"], query_adata.obs_names.astype(str).tolist())
    _write_lines(file_paths["reference_labels.txt"], ref_adata.obs["cell_type"].astype(str).tolist())
    _write_lines(file_paths["query_proxy_labels.txt"], query_adata.obs["cell_type"].astype(str).tolist())
    _write_lines(file_paths["genes.txt"], ref_adata.var_names.astype(str).tolist())

    source_sha_after = compute_sha256(source_h5ad)
    if source_sha_after != source_sha256:
        raise RuntimeError("Trevino source H5AD changed while the frozen split was being prepared.")

    split_receipt = {
        "schema_version": "2.0",
        "source_h5ad": str(source_h5ad),
        "source_sha256": source_sha256,
        "expected_source_sha256": expected_source_sha256,
        "seed": int(seed),
        "reference_tissue_ids": sorted(reference_tissues),
        "query_tissue_id": "HFT3",
        "held_out_ood_label": ood_label,
        "n_reference_cells": int(ref_adata.n_obs),
        "n_query_cells": int(query_adata.n_obs),
        "n_query_ood_cells": int((query_adata.obs["cell_type"].astype(str) == ood_label).sum()),
        "n_genes": int(ref_adata.n_vars),
        "reference_samples": sorted(ref_adata.obs["sample"].astype(str).unique().tolist()),
        "query_samples": sorted(query_adata.obs["sample"].astype(str).unique().tolist()),
        "feature_selection": {
            "selection_route": "source_highly_variable_ordered_cap",
            "source_highly_variable_required": True,
            "source_highly_variable_count": source_hvg_count,
            "max_cap": min(int(max_genes), 3000),
            "selected_count": int(ref_adata.n_vars),
            "ordered_genes_sha256": compute_sha256(file_paths["genes.txt"]),
            "upstream_retained_feature_set_independently_query_blind": False,
            "note": "The benchmark did not recompute HVGs; source highly_variable provenance may predate this split.",
        },
        "files": {name: _file_record(path) for name, path in file_paths.items()},
    }
    _atomic_write_json(out_dir / "source_and_split_receipt.json", split_receipt)
    return {
        "reference_adata": ref_adata,
        "query_adata": query_adata,
        "genes": ref_adata.var_names.astype(str).tolist(),
        "split_receipt": split_receipt,
        "frozen_dir": frozen_dir,
    }


def _write_frozen_input_manifest(
    frozen_dir: Path,
    split_receipt: dict[str, Any],
    threshold_receipt: dict[str, Any],
    parameters: dict[str, Any],
) -> tuple[Path, str]:
    """Bind every serialized child input to hashes before either lane starts."""
    threshold_path = frozen_dir / "ood_threshold_receipt.json"
    parameters_path = frozen_dir / "parameters.json"
    _atomic_write_json(threshold_path, threshold_receipt)
    _atomic_write_json(parameters_path, parameters)
    names = [
        "reference_matrix.npz",
        "query_matrix.npz",
        "reference_cells.txt",
        "query_cells.txt",
        "reference_labels.txt",
        "query_proxy_labels.txt",
        "genes.txt",
        "ood_threshold_receipt.json",
        "parameters.json",
    ]
    manifest = {
        "schema_version": "2.0",
        "source_sha256": split_receipt["source_sha256"],
        "files": {name: _file_record(frozen_dir / name) for name in names},
        "matrix_shapes": {
            "reference": [int(split_receipt["n_reference_cells"]), int(split_receipt["n_genes"])],
            "query": [int(split_receipt["n_query_cells"]), int(split_receipt["n_genes"])],
        },
        "parameters_sha256": _json_sha256(parameters),
        "threshold_receipt_sha256": compute_sha256(threshold_path),
        "ordered_gene_sha256": compute_sha256(frozen_dir / "genes.txt"),
        "contract": "Both lane children must consume these exact hash-bound files before mapping.",
    }
    manifest_path = frozen_dir / "input_manifest.json"
    _atomic_write_json(manifest_path, manifest)
    return manifest_path, compute_sha256(manifest_path)


def verify_frozen_inputs(input_dir: Path, expected_manifest_sha256: Optional[str] = None) -> dict[str, Any]:
    """Fail closed on any matrix, order, threshold, or parameter drift."""
    input_dir = Path(input_dir)
    manifest_path = input_dir / "input_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Frozen input manifest is missing: {manifest_path}")
    actual_manifest_sha = compute_sha256(manifest_path)
    if expected_manifest_sha256 and actual_manifest_sha != expected_manifest_sha256:
        raise ValueError("Frozen input manifest SHA-256 differs from the coordinator contract.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for name, record in manifest.get("files", {}).items():
        file_path = input_dir / name
        if not file_path.is_file():
            raise FileNotFoundError(f"Frozen input artifact is missing: {name}")
        if record.get("sha256") != compute_sha256(file_path):
            raise ValueError(f"Frozen input SHA-256 mismatch: {name}")
        if int(record.get("size_bytes", -1)) != int(file_path.stat().st_size):
            raise ValueError(f"Frozen input size mismatch: {name}")

    reference_matrix = sp.load_npz(input_dir / "reference_matrix.npz").tocsr()
    query_matrix = sp.load_npz(input_dir / "query_matrix.npz").tocsr()
    reference_cells = _read_lines(input_dir / "reference_cells.txt")
    query_cells = _read_lines(input_dir / "query_cells.txt")
    reference_labels = _read_lines(input_dir / "reference_labels.txt")
    query_proxy_labels = _read_lines(input_dir / "query_proxy_labels.txt")
    genes = _read_lines(input_dir / "genes.txt")
    expected_shapes = manifest.get("matrix_shapes", {})
    if tuple(expected_shapes.get("reference", [])) != tuple(reference_matrix.shape):
        raise ValueError("Frozen reference matrix shape differs from the coordinator contract.")
    if tuple(expected_shapes.get("query", [])) != tuple(query_matrix.shape):
        raise ValueError("Frozen query matrix shape differs from the coordinator contract.")
    if not (len(reference_cells) == len(reference_labels) == reference_matrix.shape[0]):
        raise ValueError("Frozen reference cell-list/label/matrix lengths disagree.")
    if not (len(query_cells) == len(query_proxy_labels) == query_matrix.shape[0]):
        raise ValueError("Frozen query cell-list/label/matrix lengths disagree.")
    if len(genes) != reference_matrix.shape[1] or len(genes) != query_matrix.shape[1]:
        raise ValueError("Frozen ordered gene list does not match serialized matrix columns.")
    if len(set(reference_cells)) != len(reference_cells) or len(set(query_cells)) != len(query_cells):
        raise ValueError("Frozen cell IDs must be unique within each lane input.")

    parameters = json.loads((input_dir / "parameters.json").read_text(encoding="utf-8"))
    if _json_sha256(parameters) != manifest.get("parameters_sha256"):
        raise ValueError("Frozen parameter payload differs from the coordinator contract.")
    if int(parameters.get("warmup_rows", 0)) != _WARMUP_ROWS:
        raise ValueError(f"Benchmark warmup must be exactly {_WARMUP_ROWS} rows.")
    if int(parameters.get("measured_repetitions", 0)) != _MEASURED_REPETITIONS:
        raise ValueError(f"Benchmark requires exactly {_MEASURED_REPETITIONS} measured repetitions.")
    threshold_receipt = json.loads((input_dir / "ood_threshold_receipt.json").read_text(encoding="utf-8"))
    if compute_sha256(input_dir / "ood_threshold_receipt.json") != manifest.get("threshold_receipt_sha256"):
        raise ValueError("Frozen OOD threshold receipt differs from the coordinator contract.")
    threshold = float(threshold_receipt.get("distance_threshold"))
    if not 0.0 <= threshold <= 2.0:
        raise ValueError("Frozen OOD distance threshold is outside [0, 2].")

    return {
        "manifest": manifest,
        "manifest_sha256": actual_manifest_sha,
        "reference_matrix": reference_matrix,
        "query_matrix": query_matrix,
        "reference_cells": reference_cells,
        "query_cells": query_cells,
        "reference_labels": reference_labels,
        "query_proxy_labels": query_proxy_labels,
        "genes": genes,
        "parameters": parameters,
        "threshold_receipt": threshold_receipt,
    }


def _mapping_sha256(frame: pd.DataFrame) -> str:
    consumed = frame.loc[:, list(_CONSUMED_OUTPUT_COLUMNS)]
    return _sha256_bytes(consumed.to_csv(lineterminator="\n").encode("utf-8"))


def _scientific_output_sha256(frame: pd.DataFrame) -> str:
    consumed = frame.loc[:, list(_EXACT_SCIENTIFIC_OUTPUT_COLUMNS)]
    return _sha256_bytes(consumed.to_csv(lineterminator="\n").encode("utf-8"))


def _assess_within_lane_repeatability(frames: list[pd.DataFrame]) -> dict[str, Any]:
    """Require exact scientific decisions and tightly bounded float drift.

    Brute-force CUDA reductions are not guaranteed to be bitwise stable across
    repeated launches.  Scientific reproducibility is therefore evaluated at
    the consumed-output boundary: every discrete decision and confidence value
    must be exactly identical, while diagnostic cosine distances must remain
    finite and within a predeclared 1e-6 tolerance.  Raw hashes are retained so
    bitwise drift remains visible rather than being normalized away.
    """
    if len(frames) != _MEASURED_REPETITIONS:
        raise RuntimeError(
            f"Within-lane repeatability requires exactly {_MEASURED_REPETITIONS} measured outputs."
        )
    baseline = frames[0].loc[:, list(_CONSUMED_OUTPUT_COLUMNS)]
    baseline_distances = baseline["reference_distance"].to_numpy(dtype=float)
    if not np.isfinite(baseline_distances).all():
        raise RuntimeError("Within-lane reference distances contain non-finite values.")

    max_abs_diff = 0.0
    mean_abs_diffs: list[float] = []
    raw_hashes: list[str] = []
    scientific_hashes: list[str] = []
    for repetition, frame in enumerate(frames, start=1):
        consumed = frame.loc[:, list(_CONSUMED_OUTPUT_COLUMNS)]
        if not baseline.index.equals(consumed.index):
            raise RuntimeError(f"Within-lane cell order changed in repetition {repetition}.")
        for column in _EXACT_SCIENTIFIC_OUTPUT_COLUMNS:
            if not baseline[column].equals(consumed[column]):
                raise RuntimeError(
                    f"Within-lane scientific output {column!r} changed in repetition {repetition}."
                )
        distances = consumed["reference_distance"].to_numpy(dtype=float)
        if not np.isfinite(distances).all():
            raise RuntimeError(
                f"Within-lane reference distances contain non-finite values in repetition {repetition}."
            )
        absolute_diff = np.abs(baseline_distances - distances)
        max_abs_diff = max(max_abs_diff, float(np.max(absolute_diff, initial=0.0)))
        mean_abs_diffs.append(float(np.mean(absolute_diff)))
        if not np.allclose(
            baseline_distances,
            distances,
            rtol=_WITHIN_LANE_DISTANCE_RTOL,
            atol=_WITHIN_LANE_DISTANCE_ATOL,
        ):
            raise RuntimeError(
                "Within-lane reference distances exceed the predeclared "
                f"rtol={_WITHIN_LANE_DISTANCE_RTOL:g}, atol={_WITHIN_LANE_DISTANCE_ATOL:g} contract "
                f"in repetition {repetition} (max_abs_diff={float(np.max(absolute_diff)):.9g})."
            )
        raw_hashes.append(_mapping_sha256(consumed))
        scientific_hashes.append(_scientific_output_sha256(consumed))

    if len(set(scientific_hashes)) != 1:
        raise RuntimeError("Within-lane exact scientific-output hashes differ across repetitions.")
    return {
        "status": "pass",
        "contract": "exact_scientific_outputs_plus_bounded_float_distance",
        "exact_columns": list(_EXACT_SCIENTIFIC_OUTPUT_COLUMNS),
        "distance_rtol": _WITHIN_LANE_DISTANCE_RTOL,
        "distance_atol": _WITHIN_LANE_DISTANCE_ATOL,
        "distance_max_abs_diff": max_abs_diff,
        "distance_mean_abs_diff_max": max(mean_abs_diffs, default=0.0),
        "scientific_output_sha256": scientific_hashes[0],
        "raw_mapping_sha256s": raw_hashes,
        "raw_bitwise_identical": len(set(raw_hashes)) == 1,
    }


def run_lane(
    *,
    device: str,
    input_dir: Path,
    output_dir: Path,
    expected_input_manifest_sha256: str,
) -> dict[str, Any]:
    """Private child subcommand: load frozen inputs and measure one backend."""
    device = str(device).strip().lower()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    receipt_path = output_dir / "lane_receipt.json"
    started = time.perf_counter()
    try:
        inputs = verify_frozen_inputs(input_dir, expected_input_manifest_sha256)
        parameters = inputs["parameters"]
        ref = ad.AnnData(
            X=inputs["reference_matrix"],
            obs=pd.DataFrame({"cell_type": inputs["reference_labels"]}, index=inputs["reference_cells"]),
            var=pd.DataFrame(index=inputs["genes"]),
        )
        query = ad.AnnData(
            X=inputs["query_matrix"],
            obs=pd.DataFrame({"cell_type": inputs["query_proxy_labels"]}, index=inputs["query_cells"]),
            var=pd.DataFrame(index=inputs["genes"]),
        )
        if query.n_obs < _WARMUP_ROWS:
            raise ValueError(f"Frozen query contains fewer than {_WARMUP_ROWS} rows.")
        mapping_kwargs = {
            "query_adata": query,
            "ref_adata": ref,
            "label_key": "cell_type",
            "k": int(parameters["k"]),
            "min_confidence": float(parameters["min_confidence"]),
            "distance_threshold": float(inputs["threshold_receipt"]["distance_threshold"]),
            "device": device,
            "seed": int(parameters["seed"]),
            "min_shared_genes": 1,
        }
        # Exactly one 256-row warm-up.  It is deliberately excluded from the
        # three reported repetitions while remaining inside lane RSS scope.
        map_knn_reference(query_adata=query[:_WARMUP_ROWS].copy(), **{key: value for key, value in mapping_kwargs.items() if key != "query_adata"})

        repetitions: list[dict[str, Any]] = []
        repetition_mappings: list[pd.DataFrame] = []
        canonical_mapping: Optional[pd.DataFrame] = None
        canonical_metadata: Optional[dict[str, Any]] = None
        for repetition in range(_MEASURED_REPETITIONS):
            started_rep = time.perf_counter()
            mapping, metadata = map_knn_reference(**mapping_kwargs)
            elapsed = time.perf_counter() - started_rep
            mapping_hash = _mapping_sha256(mapping)
            repetition_path = output_dir / f"mapping_repetition_{repetition + 1}.csv"
            _atomic_write_dataframe(
                repetition_path,
                mapping.loc[:, list(_CONSUMED_OUTPUT_COLUMNS)],
            )
            repetitions.append(
                {
                    "repetition": repetition + 1,
                    "wall_seconds": elapsed,
                    "mapping_path": str(repetition_path),
                    "mapping_sha256": mapping_hash,
                    "mapping_file_sha256": compute_sha256(repetition_path),
                    "scientific_output_sha256": _scientific_output_sha256(mapping),
                }
            )
            repetition_mappings.append(mapping)
            if canonical_mapping is None:
                canonical_mapping = mapping
                canonical_metadata = metadata
        assert canonical_mapping is not None and canonical_metadata is not None
        repeatability = _assess_within_lane_repeatability(repetition_mappings)
        device_info = dict(canonical_metadata.get("device_info", {}))
        if device == "gpu" and (
            device_info.get("backend_residency") != "cuda"
            or not bool(device_info.get("cuda_residency_verified", False))
        ):
            raise RuntimeError("GPU lane did not provide direct CUDA residency evidence.")

        mapping_path = output_dir / "mapping.csv"
        _atomic_write_dataframe(mapping_path, canonical_mapping.loc[:, list(_CONSUMED_OUTPUT_COLUMNS)])
        metrics = {
            "schema_version": "2.0",
            "status": "completed",
            "device": device,
            "input_manifest_sha256": inputs["manifest_sha256"],
            "interpreter": sys.executable,
            "warmup_rows": _WARMUP_ROWS,
            "measured_repetitions": _MEASURED_REPETITIONS,
            "repetitions": repetitions,
            "median_wall_seconds": float(np.median([record["wall_seconds"] for record in repetitions])),
            "mean_wall_seconds": float(np.mean([record["wall_seconds"] for record in repetitions])),
            "lane_end_to_end_seconds": time.perf_counter() - started,
            "ru_maxrss_mb": get_peak_rss_mb(),
            "ru_maxrss_scope": "independent child process including serialized-input loading and warmup",
            "mapping_path": str(mapping_path),
            "mapping_sha256": compute_sha256(mapping_path),
            "within_lane_repeatability": repeatability,
            "metadata": canonical_metadata,
        }
        metrics_path = output_dir / "metrics.json"
        _atomic_write_json(metrics_path, metrics)
        receipt = {
            "schema_version": "2.0",
            "status": "completed",
            "device": device,
            "input_manifest_sha256": inputs["manifest_sha256"],
            "metrics_path": str(metrics_path),
            "metrics_sha256": compute_sha256(metrics_path),
            "mapping_path": str(mapping_path),
            "mapping_sha256": compute_sha256(mapping_path),
        }
        _atomic_write_json(receipt_path, receipt)
        return receipt
    except Exception as exc:
        failure = {
            "schema_version": "2.0",
            "status": "failed",
            "device": device,
            "input_manifest_sha256": expected_input_manifest_sha256,
            "error_type": exc.__class__.__name__,
            "error_message": str(exc),
            "lane_end_to_end_seconds": time.perf_counter() - started,
            "ru_maxrss_mb": get_peak_rss_mb(),
        }
        _atomic_write_json(receipt_path, failure)
        raise


def verify_lane_output(
    output_dir: Path,
    *,
    device: str,
    expected_input_manifest_sha256: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Read a child lane only after its receipt and artifact hashes validate."""
    output_dir = Path(output_dir)
    receipt = json.loads((output_dir / "lane_receipt.json").read_text(encoding="utf-8"))
    if receipt.get("status") != "completed":
        raise RuntimeError(f"{device} child did not complete: {receipt.get('error_message', receipt)}")
    if receipt.get("device") != device or receipt.get("input_manifest_sha256") != expected_input_manifest_sha256:
        raise ValueError(f"{device} child receipt does not match the frozen lane contract.")
    mapping_path = Path(receipt["mapping_path"])
    metrics_path = Path(receipt["metrics_path"])
    if compute_sha256(mapping_path) != receipt.get("mapping_sha256"):
        raise ValueError(f"{device} mapping artifact SHA-256 mismatch.")
    if compute_sha256(metrics_path) != receipt.get("metrics_sha256"):
        raise ValueError(f"{device} metrics artifact SHA-256 mismatch.")
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    if metrics.get("input_manifest_sha256") != expected_input_manifest_sha256:
        raise ValueError(f"{device} metrics did not consume the coordinator input manifest.")
    repetition_frames: list[pd.DataFrame] = []
    for record in metrics.get("repetitions", []):
        repetition_path = Path(record.get("mapping_path", ""))
        if not repetition_path.is_file():
            raise ValueError(f"{device} child repetition artifact is missing: {repetition_path}")
        repetition_file_sha256 = compute_sha256(repetition_path)
        if repetition_file_sha256 != record.get("mapping_file_sha256"):
            raise ValueError(f"{device} child repetition artifact SHA-256 mismatch.")
        repetition_frame = pd.read_csv(repetition_path, index_col=0)
        if _mapping_sha256(repetition_frame) != record.get("mapping_sha256"):
            raise ValueError(f"{device} child repetition mapping hash mismatch.")
        if _scientific_output_sha256(repetition_frame) != record.get("scientific_output_sha256"):
            raise ValueError(f"{device} child repetition scientific-output hash mismatch.")
        repetition_frames.append(repetition_frame)
    try:
        verified_repeatability = _assess_within_lane_repeatability(repetition_frames)
    except RuntimeError as exc:
        raise ValueError(f"{device} child output repeatability failed verification: {exc}") from exc
    recorded_repeatability = metrics.get("within_lane_repeatability", {})
    for field in (
        "status",
        "contract",
        "exact_columns",
        "distance_rtol",
        "distance_atol",
        "scientific_output_sha256",
        "raw_mapping_sha256s",
        "raw_bitwise_identical",
    ):
        if recorded_repeatability.get(field) != verified_repeatability.get(field):
            raise ValueError(f"{device} child repeatability receipt mismatch for {field}.")
    for field in ("distance_max_abs_diff", "distance_mean_abs_diff_max"):
        if not np.isclose(
            float(recorded_repeatability.get(field, np.nan)),
            float(verified_repeatability[field]),
            rtol=0.0,
            atol=1e-12,
        ):
            raise ValueError(f"{device} child repeatability receipt mismatch for {field}.")
    if device == "gpu":
        device_info = metrics.get("metadata", {}).get("device_info", {})
        if device_info.get("backend_residency") != "cuda" or not device_info.get("cuda_residency_verified"):
            raise ValueError("GPU child lacks CUDA residency evidence.")
    mapping = pd.read_csv(mapping_path, index_col=0)
    if list(mapping.columns) != list(_CONSUMED_OUTPUT_COLUMNS):
        raise ValueError(f"{device} mapping output columns differ from the consumed-output contract.")
    return mapping, metrics


def _poll_gpu_vram(pid: int) -> dict[str, Any]:
    executable = shutil.which("nvidia-smi")
    if not executable:
        return {"peak_vram_mb": None, "reason": "nvidia-smi is unavailable", "samples": []}
    try:
        completed = subprocess.run(
            [executable, "--query-compute-apps=pid,used_gpu_memory", "--format=csv,noheader,nounits"],
            text=True,
            capture_output=True,
            check=True,
            timeout=5,
        )
    except Exception as exc:
        return {"peak_vram_mb": None, "reason": f"nvidia-smi polling failed: {exc}", "samples": []}
    samples: list[int] = []
    for line in completed.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) != 2:
            continue
        try:
            if int(parts[0]) == int(pid):
                samples.append(int(parts[1]))
        except ValueError:
            continue
    return {
        "peak_vram_mb": max(samples) if samples else None,
        "reason": None if samples else "nvidia-smi returned no active allocation for the child PID",
        "samples": samples,
    }


def _run_lane_process(
    *,
    device: str,
    input_dir: Path,
    output_dir: Path,
    expected_input_manifest_sha256: str,
    logs_dir: Path,
) -> dict[str, Any]:
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "_lane",
        "--device",
        device,
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--expected-input-manifest-sha256",
        expected_input_manifest_sha256,
    ]
    child_environment = os.environ.copy()
    cuda_bootstrap: dict[str, Any] = {}
    if device == "gpu":
        interpreter_prefix = Path(sys.executable).resolve().parent.parent
        toolkit_root = interpreter_prefix / "targets" / "x86_64-linux"
        if (toolkit_root / "include" / "cuda_runtime.h").is_file():
            # AnnData imports CuPy before lane code can run. Set this in the
            # child process environment, rather than after import, so CuPy's
            # lazy CUB compiler sees the interpreter's tested CUDA headers.
            child_environment["CUDA_PATH"] = str(toolkit_root)
            toolkit_bin = toolkit_root / "bin"
            path_entries = [entry for entry in child_environment.get("PATH", "").split(os.pathsep) if entry]
            if toolkit_bin.is_dir() and str(toolkit_bin) not in path_entries:
                child_environment["PATH"] = os.pathsep.join([str(toolkit_bin), *path_entries])
            cuda_bootstrap = {
                "cuda_path": str(toolkit_root),
                "cuda_headers_verified": True,
                "source": "child_interpreter_target_toolkit",
            }
        else:
            cuda_bootstrap = {
                "cuda_path": None,
                "cuda_headers_verified": False,
                "source": "child_interpreter_target_toolkit_missing",
            }
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=child_environment,
    )
    vram_samples: list[int] = []
    vram_poll_reasons: list[str] = []
    while process.poll() is None:
        if device == "gpu":
            poll = _poll_gpu_vram(process.pid)
            vram_samples.extend(poll["samples"])
            if poll.get("reason"):
                vram_poll_reasons.append(str(poll["reason"]))
        time.sleep(0.1)
    stdout, stderr = process.communicate()
    log_path = logs_dir / f"reference_mapping_{device}_lane.log"
    _atomic_write_text(
        log_path,
        f"command={json.dumps(command)}\n"
        f"cuda_bootstrap={json.dumps(cuda_bootstrap, sort_keys=True)}\n"
        f"stdout:\n{stdout}\nstderr:\n{stderr}\n",
    )
    polling = {
        "peak_vram_mb": max(vram_samples) if vram_samples else None,
        "reason": None if vram_samples else (vram_poll_reasons[-1] if vram_poll_reasons else "CPU lane"),
        "sample_count": len(vram_samples),
    }
    if device == "gpu":
        _atomic_write_json(output_dir / "vram_poll.json", polling)
    return {
        "command": command,
        "exit_code": int(process.returncode),
        "log_path": str(log_path),
        "vram_poll": polling,
        "cuda_bootstrap": cuda_bootstrap,
    }


def _compute_ood_metrics(mapping: pd.DataFrame, proxy_labels: pd.Series, ood_label: str) -> dict[str, Any]:
    labels = proxy_labels.astype(str).reindex(mapping.index)
    if labels.isna().any():
        raise ValueError("Query proxy labels do not align with the CPU mapping output.")
    ood_mask = labels == ood_label
    known_mask = ~ood_mask
    if not ood_mask.any() or not known_mask.any():
        raise ValueError("Held-out OOD metric requires both OOD and known query cells.")
    rejected = mapping["reference_assignment_status"] != ASSIGNMENT_STATUS_ACCEPTED
    ood_rejection = float(rejected.loc[ood_mask].mean())
    known_rejection = float(rejected.loc[known_mask].mean())
    known_coverage = 1.0 - known_rejection
    accepted_known = known_mask & ~rejected
    accepted_f1 = (
        float(f1_score(labels.loc[accepted_known], mapping.loc[accepted_known, "reference_cell_type"], average="macro", zero_division=0))
        if accepted_known.any()
        else 0.0
    )
    all_known_f1 = float(
        f1_score(labels.loc[known_mask], mapping.loc[known_mask, "reference_cell_type"], average="macro", zero_division=0)
    )
    separation = ood_rejection - known_rejection
    gates = {
        "known_coverage_ge_80": known_coverage >= 0.80,
        "known_macro_f1_ge_70": all_known_f1 >= 0.70,
        "ood_rejection_recall_ge_80": ood_rejection >= 0.80,
        "delta_rejection_ge_30": separation >= 0.30,
    }
    passed = all(gates.values())
    return {
        "held_out_label": ood_label,
        "n_held_out_cells": int(ood_mask.sum()),
        "n_known_cells": int(known_mask.sum()),
        "held_out_ood_rejection_recall": ood_rejection,
        "known_cells_rejection_rate": known_rejection,
        "known_cells_acceptance_coverage": known_coverage,
        "delta_rejection_rate_ood_vs_known": separation,
        "known_accepted_macro_f1": accepted_f1,
        "known_all_cells_macro_f1_rejected_as_unknown": all_known_f1,
        "gates": gates,
        "all_falsifiability_gates_passed": passed,
        "verdict": "PASS" if passed else "FAIL_NOT_PROMOTED",
        "claim_class": "technical_ood_falsifiability_using_pipeline_derived_proxy_labels",
    }


def _write_validation_summary(
    path: Path,
    *,
    run_dir: Path,
    technical_verdict: str,
    split_receipt: dict[str, Any],
    threshold_receipt: dict[str, Any],
    cpu_metrics: Optional[dict[str, Any]],
    gpu_metrics: Optional[dict[str, Any]],
    parity: Optional[dict[str, Any]],
    ood_metrics: Optional[dict[str, Any]],
    failure_reason: Optional[str] = None,
) -> None:
    lines = [
        "# Reference Atlas OOD Technical Benchmark",
        "",
        f"- Run: `{run_dir}`",
        f"- Technical verdict: `{technical_verdict}`",
        f"- Source SHA-256: `{split_receipt['source_sha256']}`",
        f"- Feature route: `{split_receipt['feature_selection']['selection_route']}`",
        f"- Features: `{split_receipt['n_genes']}` (upstream retained set; not independently query-blind)",
        f"- Full-precision OOD threshold: `{threshold_receipt['distance_threshold']!r}`",
        "- Claim boundary: Trevino labels are pipeline-derived proxies; this is technical P0 evidence, not biological validation or production promotion.",
    ]
    if failure_reason:
        lines.extend(["", f"## Failure", "", failure_reason])
    if cpu_metrics:
        lines.extend(["", "## CPU lane", "", f"- Median wall seconds: `{cpu_metrics['median_wall_seconds']:.6f}`", f"- Child ru_maxrss MiB: `{cpu_metrics['ru_maxrss_mb']:.2f}`"])
    if gpu_metrics:
        device_info = gpu_metrics.get("metadata", {}).get("device_info", {})
        lines.extend(["", "## GPU lane", "", f"- Median wall seconds: `{gpu_metrics['median_wall_seconds']:.6f}`", f"- Child ru_maxrss MiB: `{gpu_metrics['ru_maxrss_mb']:.2f}`", f"- CUDA residency: `{device_info.get('cuda_residency_verified')}`"])
    if parity:
        lines.extend(["", "## Parity", "", f"- Verdict: `{parity.get('parity_verdict')}`", f"- Candidate-label agreement: `{parity.get('candidate_label_agreement')}`", f"- Status agreement: `{parity.get('assignment_status_agreement')}`", f"- Distance allclose: `{parity.get('distance_allclose_1e4')}`"])
    if ood_metrics:
        lines.extend(["", "## Held-out OOD", "", f"- Verdict: `{ood_metrics['verdict']}`", f"- Held-out rejection recall: `{ood_metrics['held_out_ood_rejection_recall']:.6f}`", f"- Known acceptance coverage: `{ood_metrics['known_cells_acceptance_coverage']:.6f}`"])
    _atomic_write_text(path, "\n".join(lines) + "\n")


def _write_run_manifest(
    *,
    project_root: Path,
    run_dir: Path,
    git_state: dict[str, Any],
    technical_verdict: str,
    validation_dir: Path,
    input_manifest_sha256: str,
) -> Path:
    artifact_paths = [
        validation_dir / "benchmark_contract.json",
        validation_dir / "source_and_split_receipt.json",
        validation_dir / "ood_threshold_receipt.json",
        validation_dir / "cpu" / "metrics.json",
        validation_dir / "gpu" / "metrics.json",
        validation_dir / "cpu_gpu_parity.json",
        validation_dir / "heldout_ood_metrics.json",
        validation_dir / "validation_summary.md",
    ]
    artifact_hashes = {
        str(path.relative_to(run_dir)): compute_sha256(path)
        for path in artifact_paths
        if path.is_file()
    }
    return write_manifest(
        run_dir,
        project_id=project_root.name,
        run_id=run_dir.name,
        modules_run=["reference_mapping_benchmark"],
        produced_on=os.uname().nodename,
        factory_python_path=_REPO_ROOT,
        factory_python_state=git_state,
        overall_status="complete" if technical_verdict == "PASS" else "failed",
        extra={
            "reference_mapping_benchmark": {
                "technical_verdict": technical_verdict,
                "input_manifest_sha256": input_manifest_sha256,
                "artifact_sha256": artifact_hashes,
                "claim_class": "technical_p0_validation_not_biological_validation_or_production_promotion",
            }
        },
    )


def run_benchmark(
    project_root: Path,
    run_id: Optional[str] = None,
    source_h5ad: Path = DEFAULT_TREVINO_SOURCE,
    expected_source_sha256: str = DEFAULT_TREVINO_SOURCE_SHA256,
    k: int = 15,
    min_confidence: float = 0.6,
    cal_quantile: float = 0.95,
    seed: int = 42,
    allow_dirty: bool = False,
) -> dict[str, Any]:
    """Coordinate isolated child lanes and publish PASS or FAIL_NOT_PROMOTED evidence."""
    if int(k) < 1:
        raise ValueError("k must be >= 1.")
    if not 0.0 <= float(min_confidence) <= 1.0:
        raise ValueError("min_confidence must be in [0, 1].")
    if not 0.0 < float(cal_quantile) < 1.0:
        raise ValueError("cal_quantile must be in (0, 1).")
    git_state = factory_git_state(_REPO_ROOT)
    if git_state["dirty"] and not allow_dirty:
        raise RuntimeError("Reference benchmark requires a clean implementation SHA; commit or pass --allow-dirty for non-promotable diagnostics.")
    project_root = Path(project_root)
    run_dir = resolve_run_dir(project_root, run_id=run_id, factory_sha=git_state["sha"] or None)
    validation_dir = run_dir / "python" / "validation"
    if validation_dir.exists() and any(validation_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite an existing benchmark validation directory: {validation_dir}")
    validation_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = run_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    prepared = prepare_trevino_split(
        source_h5ad=source_h5ad,
        out_dir=validation_dir,
        expected_source_sha256=expected_source_sha256,
        seed=seed,
    )
    ref_adata = prepared["reference_adata"]
    query_adata = prepared["query_adata"]
    threshold, threshold_receipt = calibrate_reference_ood_threshold(
        ref_adata=ref_adata,
        ref_genes=prepared["genes"],
        mode="reference_quantile",
        quantile=cal_quantile,
        group_key="sample",
        cal_fraction=0.2,
        seed=seed,
        k=k,
    )
    threshold_receipt["source_sha256"] = prepared["split_receipt"]["source_sha256"]
    threshold_receipt["feature_order_sha256"] = prepared["split_receipt"]["feature_selection"]["ordered_genes_sha256"]
    parameters = {
        "k": int(k),
        "min_confidence": float(min_confidence),
        "seed": int(seed),
        "warmup_rows": _WARMUP_ROWS,
        "measured_repetitions": _MEASURED_REPETITIONS,
        "consumed_output_columns": list(_CONSUMED_OUTPUT_COLUMNS),
    }
    input_manifest_path, input_manifest_sha256 = _write_frozen_input_manifest(
        prepared["frozen_dir"], prepared["split_receipt"], threshold_receipt, parameters
    )
    benchmark_contract = {
        "schema_version": "2.0",
        "factory_git_state": git_state,
        "implementation_sha256": {
            "benchmark_reference_mapping.py": compute_sha256(Path(__file__).resolve()),
            "workflow/modular/_reference_mapping.py": compute_sha256(_REPO_ROOT / "workflow" / "modular" / "_reference_mapping.py"),
        },
        "source_sha256": prepared["split_receipt"]["source_sha256"],
        "input_manifest_path": str(input_manifest_path),
        "input_manifest_sha256": input_manifest_sha256,
        "threshold_receipt_sha256": compute_sha256(prepared["frozen_dir"] / "ood_threshold_receipt.json"),
        "parameters": parameters,
        "execution": {
            "cpu_gpu_process_isolation": True,
            "interpreter_for_both_children": sys.executable,
            "gpu_failure_policy": "failure_receipt_and_nonzero_no_cpu_fallback",
            "rss_scope": "independent child ru_maxrss including frozen-input loading",
        },
        "claim_class": "technical_p0_validation_not_biological_validation_or_production_promotion",
    }
    _atomic_write_json(validation_dir / "benchmark_contract.json", benchmark_contract)

    cpu_dir = validation_dir / "cpu"
    gpu_dir = validation_dir / "gpu"
    cpu_launch = _run_lane_process(
        device="cpu", input_dir=prepared["frozen_dir"], output_dir=cpu_dir,
        expected_input_manifest_sha256=input_manifest_sha256, logs_dir=logs_dir,
    )
    if cpu_launch["exit_code"] != 0:
        reason = f"CPU child failed with exit {cpu_launch['exit_code']}; see {cpu_launch['log_path']}"
        _write_validation_summary(
            validation_dir / "validation_summary.md", run_dir=run_dir, technical_verdict="FAIL_NOT_PROMOTED",
            split_receipt=prepared["split_receipt"], threshold_receipt=threshold_receipt,
            cpu_metrics=None, gpu_metrics=None, parity=None, ood_metrics=None, failure_reason=reason,
        )
        _write_run_manifest(project_root=project_root, run_dir=run_dir, git_state=git_state, technical_verdict="FAIL_NOT_PROMOTED", validation_dir=validation_dir, input_manifest_sha256=input_manifest_sha256)
        return {"run_id": run_dir.name, "technical_verdict": "FAIL_NOT_PROMOTED", "failure_reason": reason}
    cpu_mapping, cpu_metrics = verify_lane_output(cpu_dir, device="cpu", expected_input_manifest_sha256=input_manifest_sha256)

    gpu_launch = _run_lane_process(
        device="gpu", input_dir=prepared["frozen_dir"], output_dir=gpu_dir,
        expected_input_manifest_sha256=input_manifest_sha256, logs_dir=logs_dir,
    )
    if gpu_launch["exit_code"] != 0:
        reason = f"GPU child failed with exit {gpu_launch['exit_code']}; no CPU fallback occurred; see {gpu_launch['log_path']}"
        ood_metrics = _compute_ood_metrics(cpu_mapping, query_adata.obs["cell_type"], prepared["split_receipt"]["held_out_ood_label"])
        _atomic_write_json(validation_dir / "heldout_ood_metrics.json", ood_metrics)
        _write_validation_summary(
            validation_dir / "validation_summary.md", run_dir=run_dir, technical_verdict="FAIL_NOT_PROMOTED",
            split_receipt=prepared["split_receipt"], threshold_receipt=threshold_receipt,
            cpu_metrics=cpu_metrics, gpu_metrics=None, parity=None, ood_metrics=ood_metrics, failure_reason=reason,
        )
        _write_run_manifest(project_root=project_root, run_dir=run_dir, git_state=git_state, technical_verdict="FAIL_NOT_PROMOTED", validation_dir=validation_dir, input_manifest_sha256=input_manifest_sha256)
        return {"run_id": run_dir.name, "technical_verdict": "FAIL_NOT_PROMOTED", "failure_reason": reason, "ood_summary": ood_metrics}
    gpu_mapping, gpu_metrics = verify_lane_output(gpu_dir, device="gpu", expected_input_manifest_sha256=input_manifest_sha256)

    parity = compute_parity_metrics(
        cpu_mapping,
        gpu_mapping,
        ground_truth=query_adata.obs["cell_type"],
        ood_label=prepared["split_receipt"]["held_out_ood_label"],
    )
    _atomic_write_json(validation_dir / "cpu_gpu_parity.json", parity)
    ood_metrics = _compute_ood_metrics(cpu_mapping, query_adata.obs["cell_type"], prepared["split_receipt"]["held_out_ood_label"])
    _atomic_write_json(validation_dir / "heldout_ood_metrics.json", ood_metrics)
    technical_verdict = "PASS" if parity.get("parity_gates_passed") and ood_metrics["all_falsifiability_gates_passed"] else "FAIL_NOT_PROMOTED"
    _write_validation_summary(
        validation_dir / "validation_summary.md", run_dir=run_dir, technical_verdict=technical_verdict,
        split_receipt=prepared["split_receipt"], threshold_receipt=threshold_receipt,
        cpu_metrics=cpu_metrics, gpu_metrics=gpu_metrics, parity=parity, ood_metrics=ood_metrics,
    )
    manifest_path = _write_run_manifest(
        project_root=project_root, run_dir=run_dir, git_state=git_state,
        technical_verdict=technical_verdict, validation_dir=validation_dir,
        input_manifest_sha256=input_manifest_sha256,
    )
    return {
        "run_id": run_dir.name,
        "technical_verdict": technical_verdict,
        "manifest_path": str(manifest_path),
        "cal_receipt": threshold_receipt,
        "cpu_metrics": cpu_metrics,
        "gpu_metrics": gpu_metrics,
        "parity_summary": parity,
        "ood_summary": ood_metrics,
    }


def _lane_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Private isolated reference-mapping benchmark lane")
    parser.add_argument("--device", required=True, choices=["cpu", "gpu"])
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--expected-input-manifest-sha256", required=True)
    args = parser.parse_args(argv)
    try:
        run_lane(
            device=args.device,
            input_dir=Path(args.input_dir),
            output_dir=Path(args.output_dir),
            expected_input_manifest_sha256=args.expected_input_manifest_sha256,
        )
        return 0
    except Exception as exc:
        logger.error("%s lane failed: %s", args.device, exc, exc_info=True)
        return 1


def main(argv: Optional[list[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "_lane":
        return _lane_main(argv[1:])
    parser = argparse.ArgumentParser(description="Run isolated Reference Atlas OOD CPU/GPU technical benchmark")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--source-h5ad", default=str(DEFAULT_TREVINO_SOURCE))
    parser.add_argument("--expected-source-sha256", default=DEFAULT_TREVINO_SOURCE_SHA256)
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument("--min-confidence", type=float, default=0.6)
    parser.add_argument("--cal-quantile", type=float, default=0.95)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--allow-dirty", action="store_true", help="Diagnostic-only; clean SHA is required for promotable evidence")
    args = parser.parse_args(argv)
    try:
        result = run_benchmark(
            project_root=Path(args.project_root),
            run_id=args.run_id,
            source_h5ad=Path(args.source_h5ad),
            expected_source_sha256=args.expected_source_sha256,
            k=args.k,
            min_confidence=args.min_confidence,
            cal_quantile=args.cal_quantile,
            seed=args.seed,
            allow_dirty=args.allow_dirty,
        )
        print(json.dumps(result, indent=2, default=str, sort_keys=True))
        return 0 if result["technical_verdict"] == "PASS" else 1
    except Exception as exc:
        logger.error("Benchmark execution failed: %s", exc, exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
