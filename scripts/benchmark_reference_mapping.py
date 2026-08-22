#!/usr/bin/env python3
"""Technical benchmark and parity evaluation for Reference Atlas OOD KNN mapping.

Evaluates CPU (scikit-learn) vs GPU (cuML RAPIDS 26.08) on exact-input Trevino splits:
- Held-out biological query: Tissue.ID == 'HFT3'
- Reference index: Tissue.ID in {'HFT7', 'HFT5', 'HFT6'}
- Held-out OOD negative control: 'Microglia' (completely purged from reference fit/cal)
- Validates 5 predeclared analytical parity gates + 4 biological falsifiability criteria.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import logging
import os
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# Add repo root to path
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from workflow.modular._reference_mapping import (
    ASSIGNMENT_STATUS_ACCEPTED,
    calibrate_reference_ood_threshold,
    compute_parity_metrics,
    map_knn_reference,
)
from workflow.modular.manifest_writer import factory_git_state
from workflow.modular.project_paths import generate_run_id, resolve_run_dir

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("benchmark_reference_mapping")

DEFAULT_TREVINO_SOURCE = Path(
    "/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/final_adata.h5ad"
)


def compute_sha256(filepath: Path) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def get_peak_rss_mb() -> float:
    """Return peak RSS in MB for current process."""
    usage = resource.getrusage(resource.RUSAGE_SELF)
    # ru_maxrss is in kilobytes on Linux
    return float(usage.ru_maxrss / 1024.0)


def select_reference_hvg_genes(
    adata: ad.AnnData,
    reference_cells: list[str],
    max_genes: int,
) -> list[str]:
    """Select Scanpy Seurat-style HVGs using the purged reference cells only."""
    if max_genes < 1:
        raise ValueError(f"max_genes must be >= 1, got {max_genes}")
    if not reference_cells:
        raise ValueError("Reference-only feature selection requires at least one reference cell.")

    if np.any(adata.obs_names.get_indexer(reference_cells) < 0):
        raise ValueError("Reference-only feature selection received cell IDs absent from the source AnnData.")

    reference = adata[reference_cells, :].copy()
    sc.pp.highly_variable_genes(
        reference,
        flavor="seurat",
        n_top_genes=min(max_genes, reference.n_vars),
        subset=False,
        inplace=True,
        check_values=False,
    )
    if "highly_variable" not in reference.var or not bool(reference.var["highly_variable"].any()):
        raise ValueError("Reference-only Scanpy HVG selection returned no genes.")
    selected_var = reference.var.loc[reference.var["highly_variable"]].copy()
    selected_var["_gene_name"] = selected_var.index.astype(str)
    selected_var = selected_var.sort_values(
        ["dispersions_norm", "_gene_name"],
        ascending=[False, True],
        kind="mergesort",
    )
    selected = selected_var["_gene_name"].head(max_genes).tolist()
    logger.info(
        "Selected %d Scanpy Seurat HVGs from %d reference cells only",
        len(selected),
        len(reference_cells),
    )
    return selected


def prepare_trevino_split(
    source_h5ad: Path,
    out_dir: Path,
    max_ref_cells: int = 15000,
    max_query_cells: int = 6000,
    max_genes: int = 3000,
    ood_label: str = "Microglia",
    seed: int = 42,
) -> dict[str, Any]:
    """Prepare and serialize exact frozen split from Trevino H5AD."""
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Loading Trevino source: %s", source_h5ad)
    source_sha = compute_sha256(source_h5ad)

    adata = ad.read_h5ad(source_h5ad)
    logger.info("Loaded AnnData shape: %s", adata.shape)

    # 1. Partition query vs reference by Tissue.ID
    # Query: Tissue.ID == 'HFT3' (GW20_3148, GW20_3408)
    # Reference: Tissue.ID in {'HFT7', 'HFT5', 'HFT6'}
    is_hft3 = adata.obs["Tissue.ID"].astype(str) == "HFT3"
    query_obs_all = adata.obs[is_hft3]
    reference_tissues = {"HFT7", "HFT5", "HFT6"}
    ref_obs_all = adata.obs[adata.obs["Tissue.ID"].astype(str).isin(reference_tissues)]

    if query_obs_all.empty or ref_obs_all.empty:
        raise ValueError(
            f"Trevino split is empty (query HFT3={len(query_obs_all)}, "
            f"reference HFT7/HFT5/HFT6={len(ref_obs_all)})."
        )

    logger.info("Initial split: Query (HFT3)=%d, Reference=%d", len(query_obs_all), len(ref_obs_all))

    # 2. Purge OOD label (Microglia) from Reference
    ref_obs_clean = ref_obs_all[ref_obs_all["cell_type"].astype(str) != ood_label]
    logger.info(
        "Reference after purging %s: %d cells (purged %d)",
        ood_label,
        len(ref_obs_clean),
        len(ref_obs_all) - len(ref_obs_clean),
    )

    rng = np.random.default_rng(seed)

    # 3. Stratified downsample Reference to max_ref_cells
    ref_grouped = ref_obs_clean.groupby("cell_type")
    ref_selected_cells: list[str] = []
    max_per_ref_lbl = int(np.ceil(max_ref_cells / max(len(ref_grouped), 1)))

    for lbl, group in sorted(ref_grouped.groups.items()):
        cell_ids = np.array(group)
        rng.shuffle(cell_ids)
        ref_selected_cells.extend(cell_ids[:max_per_ref_lbl].tolist())

    if len(ref_selected_cells) > max_ref_cells:
        rng.shuffle(ref_selected_cells)
        ref_selected_cells = sorted(ref_selected_cells[:max_ref_cells])
    else:
        ref_selected_cells = sorted(ref_selected_cells)

    # 4. Downsample Query to max_query_cells, preserving ALL Microglia cells
    query_ood_cells = query_obs_all[query_obs_all["cell_type"].astype(str) == ood_label].index.tolist()
    if len(query_ood_cells) > max_query_cells:
        raise ValueError(
            f"Held-out label {ood_label!r} has {len(query_ood_cells)} cells, exceeding max_query_cells={max_query_cells}; "
            "refusing to silently downsample the OOD control."
        )
    query_known_obs = query_obs_all[query_obs_all["cell_type"].astype(str) != ood_label]

    remaining_query_slots = max(0, max_query_cells - len(query_ood_cells))
    known_grouped = query_known_obs.groupby("cell_type")
    query_known_selected: list[str] = []
    max_per_known = int(np.ceil(remaining_query_slots / max(len(known_grouped), 1)))

    for lbl, group in sorted(known_grouped.groups.items()):
        cell_ids = np.array(group)
        rng.shuffle(cell_ids)
        query_known_selected.extend(cell_ids[:max_per_known].tolist())

    if len(query_known_selected) > remaining_query_slots:
        rng.shuffle(query_known_selected)
        query_known_selected = query_known_selected[:remaining_query_slots]

    query_selected_cells = sorted(query_ood_cells + query_known_selected)

    logger.info(
        "Final split cells: Reference=%d, Query=%d (of which %d are %s)",
        len(ref_selected_cells),
        len(query_selected_cells),
        len(query_ood_cells),
        ood_label,
    )

    # 5. Rank features using only the purged, selected reference cells. Query
    # cells and the held-out OOD class cannot influence the feature space.
    selected_genes = select_reference_hvg_genes(
        adata,
        reference_cells=ref_selected_cells,
        max_genes=max_genes,
    )

    # Extract sub-AnnData objects
    ref_adata = adata[ref_selected_cells, selected_genes].copy()
    query_adata = adata[query_selected_cells, selected_genes].copy()
    del adata

    source_sha_after = compute_sha256(source_h5ad)
    if source_sha_after != source_sha:
        raise RuntimeError("Trevino source H5AD changed while preparing the frozen benchmark split.")

    reference_sample_set = set(ref_adata.obs["sample"].astype(str))
    query_sample_set = set(query_adata.obs["sample"].astype(str))
    if reference_sample_set & query_sample_set:
        raise RuntimeError("Reference and query sample sets overlap; sample-disjoint benchmark contract violated.")
    if ood_label in set(ref_adata.obs["cell_type"].astype(str)):
        raise RuntimeError(f"Held-out OOD label {ood_label!r} remains in the reference after purging.")

    # Ensure sparse CSR matrices
    if not sp.issparse(ref_adata.X):
        ref_adata.X = sp.csr_matrix(ref_adata.X)
    if not sp.issparse(query_adata.X):
        query_adata.X = sp.csr_matrix(query_adata.X)

    # Serialize files
    ref_mat_path = out_dir / "reference_matrix.npz"
    query_mat_path = out_dir / "query_matrix.npz"
    ref_cells_path = out_dir / "reference_cells.txt"
    query_cells_path = out_dir / "query_cells.txt"
    genes_path = out_dir / "genes.txt"
    receipt_path = out_dir / "source_and_split_receipt.json"

    sp.save_npz(ref_mat_path, ref_adata.X)
    sp.save_npz(query_mat_path, query_adata.X)
    np.savetxt(ref_cells_path, ref_selected_cells, fmt="%s")
    np.savetxt(query_cells_path, query_selected_cells, fmt="%s")
    np.savetxt(genes_path, selected_genes, fmt="%s")

    split_receipt = {
        "source_h5ad": str(source_h5ad),
        "source_sha256": source_sha,
        "seed": seed,
        "n_genes": len(selected_genes),
        "genes_sha256": compute_sha256(genes_path),
        "feature_selection": {
            "method": "reference_only_scanpy_seurat_hvg",
            "scanpy_version": importlib.metadata.version("scanpy"),
            "reference_cells_used": len(ref_selected_cells),
            "query_cells_used": 0,
            "held_out_ood_cells_used": 0,
            "max_genes": int(max_genes),
        },
        "n_reference_cells": len(ref_selected_cells),
        "reference_tissue_ids": list(ref_adata.obs["Tissue.ID"].unique()),
        "reference_samples": list(ref_adata.obs["sample"].unique()),
        "reference_cell_types": list(ref_adata.obs["cell_type"].unique()),
        "n_query_cells": len(query_selected_cells),
        "query_tissue_ids": list(query_adata.obs["Tissue.ID"].unique()),
        "query_samples": list(query_adata.obs["sample"].unique()),
        "query_cell_types": list(query_adata.obs["cell_type"].unique()),
        "held_out_ood_label": ood_label,
        "n_query_ood_cells": len(query_ood_cells),
        "artifacts": {
            "ref_matrix_sha256": compute_sha256(ref_mat_path),
            "query_matrix_sha256": compute_sha256(query_mat_path),
            "ref_cells_sha256": compute_sha256(ref_cells_path),
            "query_cells_sha256": compute_sha256(query_cells_path),
            "reference_labels_sha256": hashlib.sha256(
                "\n".join(ref_adata.obs["cell_type"].astype(str)).encode("utf-8")
            ).hexdigest(),
            "query_labels_sha256": hashlib.sha256(
                "\n".join(query_adata.obs["cell_type"].astype(str)).encode("utf-8")
            ).hexdigest(),
        },
    }
    receipt_path.write_text(json.dumps(split_receipt, indent=2), encoding="utf-8")

    return {
        "ref_adata": ref_adata,
        "query_adata": query_adata,
        "genes": selected_genes,
        "receipt": split_receipt,
    }


def run_benchmark(
    project_root: Path,
    run_id: Optional[str] = None,
    source_h5ad: Path = DEFAULT_TREVINO_SOURCE,
    k: int = 15,
    min_confidence: float = 0.6,
    cal_quantile: float = 0.95,
    n_reps: int = 3,
    seed: int = 42,
) -> dict[str, Any]:
    """Execute complete technical benchmark across CPU and GPU lanes."""
    git_state = factory_git_state(_REPO_ROOT)
    run_dir = resolve_run_dir(project_root, run_id=run_id, factory_sha=git_state["sha"] or None)

    val_dir = run_dir / "python" / "validation"
    cpu_dir = val_dir / "cpu"
    gpu_dir = val_dir / "gpu"
    cpu_dir.mkdir(parents=True, exist_ok=True)
    gpu_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Prepare Split
    split_info = prepare_trevino_split(
        source_h5ad=source_h5ad,
        out_dir=val_dir,
        seed=seed,
    )
    ref_adata = split_info["ref_adata"]
    query_adata = split_info["query_adata"]
    genes = split_info["genes"]
    ood_label = split_info["receipt"]["held_out_ood_label"]

    # Step 2: Calibrate OOD distance threshold on Reference-Only partition
    logger.info("Calibrating OOD distance threshold on reference-only data...")
    distance_thresh, cal_receipt = calibrate_reference_ood_threshold(
        ref_adata=ref_adata,
        ref_genes=genes,
        mode="reference_quantile",
        quantile=cal_quantile,
        group_key="sample",
        cal_fraction=0.2,
        seed=seed,
        k=k,
    )
    cal_receipt_path = val_dir / "ood_threshold_receipt.json"
    cal_receipt_path.write_text(json.dumps(cal_receipt, indent=2), encoding="utf-8")
    logger.info("Calibrated distance threshold: %.6f (written to %s)", distance_thresh, cal_receipt_path)

    benchmark_contract = {
        "schema_version": "1.0",
        "factory_git_state": git_state,
        "implementation_sha256": {
            "benchmark_reference_mapping.py": compute_sha256(Path(__file__).resolve()),
            "workflow/modular/_reference_mapping.py": compute_sha256(
                _REPO_ROOT / "workflow" / "modular" / "_reference_mapping.py"
            ),
        },
        "source_sha256": split_info["receipt"]["source_sha256"],
        "matrix_and_order_hashes": split_info["receipt"]["artifacts"],
        "genes_sha256": split_info["receipt"]["genes_sha256"],
        "threshold_receipt_sha256": compute_sha256(cal_receipt_path),
        "parameters": {
            "k": int(k),
            "min_confidence": float(min_confidence),
            "cal_quantile": float(cal_quantile),
            "seed": int(seed),
            "warmup_cells": int(min(256, query_adata.n_obs)),
            "measured_repetitions": int(n_reps),
            "consumed_outputs": [
                "reference_predicted_label",
                "reference_confidence",
                "reference_distance",
                "reference_assignment_status",
                "reference_cell_type",
                "reference_ood",
            ],
        },
        "lane_contract": "CPU and GPU consume the same in-memory matrices, cell order, gene order, threshold, seed, and outputs.",
    }
    parameter_bytes = json.dumps(benchmark_contract["parameters"], sort_keys=True).encode("utf-8")
    benchmark_contract["parameters_sha256"] = hashlib.sha256(parameter_bytes).hexdigest()
    benchmark_contract_path = val_dir / "benchmark_contract.json"
    benchmark_contract_path.write_text(json.dumps(benchmark_contract, indent=2), encoding="utf-8")

    # Prepare warm-up subset (first 256 rows)
    query_warmup = query_adata[: min(256, query_adata.n_obs)].copy()

    # Step 3: Run CPU Benchmark
    logger.info("Starting CPU Benchmark (repetitions=%d)...", n_reps)
    # Warm-up
    map_knn_reference(
        query_adata=query_warmup,
        ref_adata=ref_adata,
        label_key="cell_type",
        k=k,
        min_confidence=min_confidence,
        distance_threshold=distance_thresh,
        device="cpu",
        seed=seed,
    )

    cpu_rep_times: list[float] = []
    cpu_mapping_df: Optional[pd.DataFrame] = None
    cpu_meta: Optional[dict[str, Any]] = None

    for rep in range(n_reps):
        t0 = time.perf_counter()
        df_res, meta = map_knn_reference(
            query_adata=query_adata,
            ref_adata=ref_adata,
            label_key="cell_type",
            k=k,
            min_confidence=min_confidence,
            distance_threshold=distance_thresh,
            device="cpu",
            seed=seed,
        )
        elapsed = time.perf_counter() - t0
        cpu_rep_times.append(elapsed)
        if rep == 0:
            cpu_mapping_df = df_res
            cpu_meta = meta

    cpu_peak_rss = get_peak_rss_mb()
    cpu_metrics = {
        "device": "cpu",
        "repetitions": n_reps,
        "rep_times_sec": [round(t, 4) for t in cpu_rep_times],
        "median_wall_sec": round(float(np.median(cpu_rep_times)), 4),
        "mean_wall_sec": round(float(np.mean(cpu_rep_times)), 4),
        "peak_rss_mb": round(cpu_peak_rss, 2),
        "peak_rss_scope": "cumulative_process_ru_maxrss_includes_source_load",
        "vram_mb": None,
        "vram_reason": "Not applicable for CPU backend",
        "metadata": cpu_meta,
    }
    assert cpu_mapping_df is not None
    cpu_mapping_df.to_csv(cpu_dir / "mapping.csv")
    (cpu_dir / "metrics.json").write_text(json.dumps(cpu_metrics, indent=2), encoding="utf-8")

    # Step 4: Run GPU Benchmark
    logger.info("Starting GPU Benchmark (repetitions=%d)...", n_reps)
    gpu_available = True
    gpu_error = None
    gpu_mapping_df: Optional[pd.DataFrame] = None
    gpu_meta: Optional[dict[str, Any]] = None
    gpu_rep_times: list[float] = []
    gpu_metrics: dict[str, Any] = {}

    try:
        # Warm-up
        map_knn_reference(
            query_adata=query_warmup,
            ref_adata=ref_adata,
            label_key="cell_type",
            k=k,
            min_confidence=min_confidence,
            distance_threshold=distance_thresh,
            device="gpu",
            seed=seed,
        )

        for rep in range(n_reps):
            t0 = time.perf_counter()
            df_res, meta = map_knn_reference(
                query_adata=query_adata,
                ref_adata=ref_adata,
                label_key="cell_type",
                k=k,
                min_confidence=min_confidence,
                distance_threshold=distance_thresh,
                device="gpu",
                seed=seed,
            )
            elapsed = time.perf_counter() - t0
            gpu_rep_times.append(elapsed)
            if rep == 0:
                gpu_mapping_df = df_res
                gpu_meta = meta

        gpu_peak_rss = get_peak_rss_mb()
        gpu_metrics = {
            "device": "gpu",
            "repetitions": n_reps,
            "rep_times_sec": [round(t, 4) for t in gpu_rep_times],
            "median_wall_sec": round(float(np.median(gpu_rep_times)), 4),
            "mean_wall_sec": round(float(np.mean(gpu_rep_times)), 4),
            "peak_rss_mb": round(gpu_peak_rss, 2),
            "peak_rss_scope": "cumulative_process_ru_maxrss_includes_source_load",
            "vram_used_after_query_mb": gpu_meta.get("device_info", {}).get("vram_used_after_query_mb") if gpu_meta else None,
            "vram_delta_mb": gpu_meta.get("device_info", {}).get("vram_delta_mb") if gpu_meta else None,
            "vram_reason": None if gpu_meta else "GPU metadata unavailable",
            "backend_residency": gpu_meta.get("device_info", {}).get("backend_residency") if gpu_meta else None,
            "transfer_overhead_sec": (
                float(gpu_meta.get("timings", {}).get("h2d_sec", 0.0))
                + float(gpu_meta.get("timings", {}).get("d2h_sec", 0.0))
            ) if gpu_meta else None,
            "metadata": gpu_meta,
        }
        assert gpu_mapping_df is not None
        gpu_mapping_df.to_csv(gpu_dir / "mapping.csv")
        (gpu_dir / "metrics.json").write_text(json.dumps(gpu_metrics, indent=2), encoding="utf-8")

    except Exception as exc:
        gpu_available = False
        gpu_error = str(exc)
        logger.warning("GPU execution failed or unavailable: %s", exc)
        gpu_metrics = {
            "device": "gpu",
            "status": "failed_or_unavailable",
            "error": gpu_error,
        }
        (gpu_dir / "metrics.json").write_text(json.dumps(gpu_metrics, indent=2), encoding="utf-8")

    # Step 5: Parity and OOD Evaluation
    ground_truth = query_adata.obs["cell_type"].astype(str)
    parity_summary: dict[str, Any] = {}
    ood_summary: dict[str, Any] = {}

    if gpu_available and gpu_mapping_df is not None:
        logger.info("Computing CPU vs GPU analytical parity metrics...")
        parity_summary = compute_parity_metrics(
            df_a=cpu_mapping_df,
            df_b=gpu_mapping_df,
            ground_truth=ground_truth,
            ood_label=ood_label,
        )
        (val_dir / "cpu_gpu_parity.json").write_text(json.dumps(parity_summary, indent=2), encoding="utf-8")
        logger.info("Parity Verdict: %s", parity_summary.get("parity_verdict"))

    # Held-out OOD metrics on CPU mapping
    ood_mask = ground_truth == ood_label
    n_ood = int(ood_mask.sum())
    known_mask = ~ood_mask
    n_known = int(known_mask.sum())

    ood_rej_rate = float((cpu_mapping_df["reference_assignment_status"].loc[ood_mask] != ASSIGNMENT_STATUS_ACCEPTED).mean())
    known_rej_rate = float((cpu_mapping_df["reference_assignment_status"].loc[known_mask] != ASSIGNMENT_STATUS_ACCEPTED).mean())
    known_coverage = 1.0 - known_rej_rate

    from sklearn.metrics import f1_score

    # Accepted-only F1 measures label precision among cells that passed the
    # rejector. The all-known score also counts rejected known cells as
    # ``Unknown`` and is the stricter metric used by the promotion gate.
    acc_known_mask = known_mask & (cpu_mapping_df["reference_assignment_status"] == ASSIGNMENT_STATUS_ACCEPTED)
    if acc_known_mask.sum() > 0:
        accepted_known_macro_f1 = float(f1_score(
            ground_truth[acc_known_mask],
            cpu_mapping_df.loc[acc_known_mask, "reference_cell_type"],
            average="macro",
            zero_division=0,
        ))
    else:
        accepted_known_macro_f1 = 0.0
    all_known_macro_f1 = float(f1_score(
        ground_truth[known_mask],
        cpu_mapping_df.loc[known_mask, "reference_cell_type"],
        average="macro",
        zero_division=0,
    ))

    delta_rejection = ood_rej_rate - known_rej_rate

    # Check the 4 predeclared falsifiability criteria
    gate_known_cov = bool(known_coverage >= 0.80)
    gate_known_f1 = bool(all_known_macro_f1 >= 0.70)
    gate_ood_recall = bool(ood_rej_rate >= 0.80)
    gate_delta_rej = bool(delta_rejection >= 0.30)
    all_ood_gates_passed = gate_known_cov and gate_known_f1 and gate_ood_recall and gate_delta_rej

    ood_summary = {
        "held_out_label": ood_label,
        "n_held_out_cells": n_ood,
        "n_known_cells": n_known,
        "held_out_ood_rejection_recall": round(ood_rej_rate, 4),
        "known_cells_rejection_rate": round(known_rej_rate, 4),
        "known_cells_acceptance_coverage": round(known_coverage, 4),
        "delta_rejection_rate_ood_vs_known": round(delta_rejection, 4),
        "known_accepted_macro_f1": round(accepted_known_macro_f1, 4),
        "known_all_cells_macro_f1_rejected_as_unknown": round(all_known_macro_f1, 4),
        "gates": {
            "known_coverage_ge_80": gate_known_cov,
            "known_macro_f1_ge_70": gate_known_f1,
            "ood_rejection_recall_ge_80": gate_ood_recall,
            "delta_rejection_ge_30": gate_delta_rej,
        },
        "all_falsifiability_gates_passed": all_ood_gates_passed,
        "verdict": "PASS" if all_ood_gates_passed else "FAIL_NOT_PROMOTED",
    }
    (val_dir / "heldout_ood_metrics.json").write_text(json.dumps(ood_summary, indent=2), encoding="utf-8")

    overall_passed = bool(
        all_ood_gates_passed
        and gpu_available
        and parity_summary.get("parity_gates_passed", False)
    )

    # Generate Markdown Summary
    summary_md = f"""# Technical Validation & Benchmark Summary: Reference Atlas OOD

**Generated**: {datetime.now(timezone.utc).isoformat()}
**Source**: `{source_h5ad}`
**Run Directory**: `{run_dir}`

---

## 1. Held-Out OOD Falsifiability Evidence (`Microglia` Negative Control)

- **Query Cells**: {len(query_adata)} (Held-out Tissue `HFT3`)
- **Reference Cells**: {len(ref_adata)} (Tissues `HFT7, HFT5, HFT6`, Microglia purged)
- **Held-Out OOD Cells**: {n_ood} (`Microglia`)
- **Calibrated Cosine Distance Threshold**: `{distance_thresh:.6f}` (mode: `reference_quantile`, q={cal_quantile})

| Metric | Target Gate | Observed Value | Status |
|---|---|---|---|
| Known-label Acceptance Coverage | $\\ge 80.0\\%$ | **{known_coverage*100:.2f}%** | {'PASS' if gate_known_cov else 'FAIL'} |
| Known-label Macro-F1, all known cells (rejected = Unknown; proxy target) | $\\ge 0.700$ | **{all_known_macro_f1:.4f}** | {'PASS' if gate_known_f1 else 'FAIL'} |
| Accepted-known Macro-F1 (descriptive; proxy target) | descriptive | **{accepted_known_macro_f1:.4f}** | N/A |
| Held-out Microglia Rejection Recall | $\\ge 80.0\\%$ | **{ood_rej_rate*100:.2f}%** | {'PASS' if gate_ood_recall else 'FAIL'} |
| Rejection Separation ($\\Delta$ OOD - Known) | $\\ge 30.0\\%$ | **{delta_rejection*100:.2f}%** | {'PASS' if gate_delta_rej else 'FAIL'} |

**OOD Falsifiability Verdict**: `{ood_summary['verdict']}`

---

## 2. CPU vs GPU Analytical Parity Benchmark

| Parity Gate | Tolerance / Gate | Observed Value | Verdict |
|---|---|---|---|
| Candidate Label Agreement | $\\ge 99.50\\%$ | **{parity_summary.get('candidate_label_agreement', 0)*100:.2f}%** | {'PASS' if parity_summary.get('candidate_label_agreement', 0) >= 0.995 else 'FAIL'} |
| Assignment Status Agreement | $\\ge 99.50\\%$ | **{parity_summary.get('assignment_status_agreement', 0)*100:.2f}%** | {'PASS' if parity_summary.get('assignment_status_agreement', 0) >= 0.995 else 'FAIL'} |
| Rejection Rate Delta | $\\le 0.50\\%$ | **{parity_summary.get('delta_rejection_rate', 0)*100:.4f}%** | {'PASS' if parity_summary.get('delta_rejection_rate', 0) <= 0.005 else 'FAIL'} |
| Known Macro-F1 Delta | $\\le 0.010$ | **{parity_summary.get('delta_known_macro_f1', 0):.4f}** | {'PASS' if parity_summary.get('delta_known_macro_f1', 0) <= 0.010 else 'FAIL'} |
| Distance Agreement | `rtol=1e-4, atol=1e-4` | **{parity_summary.get('distance_allclose_1e4', False)}** (Max $\\Delta$: {parity_summary.get('distance_max_abs_diff', 0):.2e}) | {'PASS' if parity_summary.get('distance_allclose_1e4', False) else 'FAIL'} |

**Overall Parity Verdict**: `{parity_summary.get('parity_verdict', 'N/A')}`

---

## 3. Performance & Resource Comparison

| Backend | Median Wall Time (s) | Mean Wall Time (s) | Peak RSS (MB) | Accelerator / VRAM |
|---|---|---|---|---|
| **CPU (sklearn)** | {cpu_metrics['median_wall_sec']:.4f} | {cpu_metrics['mean_wall_sec']:.4f} | {cpu_metrics['peak_rss_mb']} | Host RAM |
| **GPU (cuML RAPIDS)** | {gpu_metrics.get('median_wall_sec', 'N/A')} | {gpu_metrics.get('mean_wall_sec', 'N/A')} | {gpu_metrics.get('peak_rss_mb', 'N/A')} | {gpu_meta.get('device_info', {}).get('cuda_device', 'NVIDIA GPU') if gpu_meta else 'N/A'} |

`peak_rss_mb` is cumulative process `ru_maxrss` and includes source-H5AD loading; it is retained as a reproducible upper bound, not a lane-isolated allocation. GPU VRAM after query: `{gpu_metrics.get('vram_used_after_query_mb', 'N/A')}` MB; measured H2D+D2H overhead: `{gpu_metrics.get('transfer_overhead_sec', 'N/A')}` s.

---

## 4. Promotion & Claim Honesty Boundary
1. Trevino labels are technical pipeline proxies; this validation confirms OOD rejection mathematics and CPU/GPU parity.
2. SCANVI is recorded as `not_run_missing_real_compatible_model_artifact`.
"""
    (val_dir / "validation_summary.md").write_text(summary_md, encoding="utf-8")
    logger.info("Benchmark complete. Validation summary written to %s", val_dir / "validation_summary.md")

    return {
        "run_id": run_dir.name,
        "technical_verdict": "PASS" if overall_passed else "FAIL_NOT_PROMOTED",
        "cal_receipt": cal_receipt,
        "cpu_metrics": cpu_metrics,
        "gpu_metrics": gpu_metrics,
        "parity_summary": parity_summary,
        "ood_summary": ood_summary,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Reference Atlas OOD Benchmark & Parity Evaluation")
    parser.add_argument("--project-root", required=True, help="Project root directory")
    parser.add_argument("--run-id", default=None, help="Optional run ID")
    parser.add_argument("--source-h5ad", default=str(DEFAULT_TREVINO_SOURCE), help="Path to Trevino source H5AD")
    parser.add_argument("--k", type=int, default=15, help="K nearest neighbors")
    parser.add_argument("--min-confidence", type=float, default=0.6, help="Minimum vote confidence")
    parser.add_argument("--cal-quantile", type=float, default=0.95, help="OOD calibration quantile")
    parser.add_argument("--n-reps", type=int, default=3, help="Benchmark repetitions")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()
    try:
        result = run_benchmark(
            project_root=Path(args.project_root),
            run_id=args.run_id,
            source_h5ad=Path(args.source_h5ad),
            k=args.k,
            min_confidence=args.min_confidence,
            cal_quantile=args.cal_quantile,
            n_reps=args.n_reps,
            seed=args.seed,
        )
        return 0 if result["technical_verdict"] == "PASS" else 1
    except Exception as exc:
        logger.error("Benchmark execution failed: %s", exc, exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
