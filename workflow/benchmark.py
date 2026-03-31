from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

from workflow.standard import StandardWorkflowConfig, run_standard_workflow


@dataclass(frozen=True)
class DatasetSpec:
    """Benchmark dataset specification.

    Parameters
    ----------
    name
        Dataset name.
    tenx_dir
        10x matrix directory path.
    baseline_h5ad
        Existing baseline processed file.
    """

    name: str
    tenx_dir: Path
    baseline_h5ad: Path


def cluster_ari(new_h5ad: Path, baseline_h5ad: Path) -> float:
    """Compute ARI between new and baseline Leiden labels."""

    new = ad.read_h5ad(new_h5ad)
    try:
        base = ad.read_h5ad(baseline_h5ad)
    except Exception:
        return float("nan")
    if "leiden" not in new.obs or "leiden" not in base.obs:
        raise ValueError("Both datasets must contain `leiden` in obs.")
    common = np.intersect1d(new.obs_names.values, base.obs_names.values)
    if common.size == 0:
        return float("nan")
    y1 = new[common].obs["leiden"].astype(str).values
    y2 = base[common].obs["leiden"].astype(str).values
    return float(adjusted_rand_score(y1, y2))


def velocity_direction_consistency(reference_h5ad: Path, candidate_h5ad: Path) -> float:
    """Estimate velocity direction consistency between two runs.

    Notes
    -----
    Consistency is computed as weighted sign agreement on overlapping
    non-zero edges in `velocity_graph` matrices.
    """

    ref = ad.read_h5ad(reference_h5ad)
    cand = ad.read_h5ad(candidate_h5ad)
    if "velocity_graph" not in ref.uns or "velocity_graph" not in cand.uns:
        return float("nan")
    g1 = ref.uns["velocity_graph"]
    g2 = cand.uns["velocity_graph"]
    if not sparse.issparse(g1):
        g1 = sparse.csr_matrix(g1)
    if not sparse.issparse(g2):
        g2 = sparse.csr_matrix(g2)
    n = min(g1.shape[0], g2.shape[0])
    g1 = g1[:n, :n].tocsr()
    g2 = g2[:n, :n].tocsr()
    overlap = g1.multiply(g2)
    denom = np.sum(np.abs(overlap.data))
    if denom == 0:
        return float("nan")
    agree = np.sum(np.abs(overlap.data[overlap.data > 0]))
    return float(agree / denom)


def run_benchmark(root: Path | None = None) -> Path:
    """Run end-to-end benchmark on three datasets."""

    root = root or Path("/home/zerlinshen/singlecell_factory")
    out_root = root / "output" / "workflow_benchmark"
    out_root.mkdir(parents=True, exist_ok=True)
    def latest_baseline(prefix: str) -> Path:
        cands = sorted((root / "output" / "workflow_baseline").glob(f"{prefix}_*/data/processed.h5ad"))
        return cands[-1] if cands else Path("")

    specs = [
        DatasetSpec(
            "pbmc_1k_v3",
            root / "data" / "raw" / "pbmc_1k_v3_count" / "outs" / "filtered_feature_bc_matrix",
            latest_baseline("baseline_pbmc"),
        ),
        DatasetSpec(
            "lung_carcinoma_3k",
            root / "data" / "raw" / "lung_carcinoma_3k_count" / "outs" / "filtered_feature_bc_matrix",
            latest_baseline("baseline_lung"),
        ),
        DatasetSpec(
            "pbmc_backup",
            root / "data" / "raw" / "pbmc_1k_v3_count_backup" / "outs" / "filtered_feature_bc_matrix",
            latest_baseline("baseline_pbmc_backup"),
        ),
    ]

    rows = []
    for ds in specs:
        t0 = time.perf_counter()
        out = run_standard_workflow(
            StandardWorkflowConfig(
                project=f"{ds.name}_workflow",
                tenx_dir=ds.tenx_dir,
                output_dir=out_root,
                mode="full",
                use_cache=False,
            )
        )
        elapsed = time.perf_counter() - t0
        ari = cluster_ari(out, ds.baseline_h5ad) if ds.baseline_h5ad.exists() else float("nan")
        rows.append(
            {
                "dataset": ds.name,
                "output_h5ad": str(out),
                "elapsed_seconds": elapsed,
                "cluster_ari": ari,
            }
        )

    velocity_score = velocity_direction_consistency(
        root / "rna_velocity_pseudotime_analysis" / "runtime" / "results" / "lung_with_loom_v24" / "integrated.h5ad",
        root / "rna_velocity_pseudotime_analysis" / "runtime" / "results" / "lung_with_loom_v25" / "integrated.h5ad",
    )
    report = {
        "datasets": rows,
        "velocity_direction_consistency": velocity_score,
        "targets": {
            "cluster_ari_min": 0.95,
            "velocity_direction_consistency_min": 0.9,
        },
    }
    out_json = out_root / "benchmark_report.json"
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    pd.DataFrame(rows).to_csv(out_root / "benchmark_summary.csv", index=False)
    return out_json


if __name__ == "__main__":
    print(run_benchmark())
