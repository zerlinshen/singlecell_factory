from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path
import json
import os
import resource
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from workflow.modular.modules.pseudo_velocity import compute_pseudo_velocity_vectorized
from workflow.modular.perf_baseline import (
    baseline_clustering,
    baseline_pseudo_velocity,
    baseline_qc,
)


@dataclass
class RunMetrics:
    mode: str
    dataset_used: str
    wall_seconds: float
    cpu_seconds: float
    cpu_utilization_pct_all_cores: float
    peak_rss_mb: float
    n_cells: int
    n_genes: int
    n_clusters: int


def _optimized_qc(adata):
    gene_upper = adata.var_names.astype(str).str.upper()
    adata.var["mt"] = gene_upper.str.startswith("MT-")
    adata.var["ribo"] = gene_upper.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)
    mask = (
        (adata.obs["n_genes_by_counts"] >= 200)
        & (adata.obs["n_genes_by_counts"] <= 7000)
        & (adata.obs["pct_counts_mt"] <= 20.0)
        & (adata.obs["pct_counts_ribo"] <= 50.0)
    )
    adata = adata[mask, :].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    return adata


def _optimized_clustering(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if sparse.issparse(adata.X):
        adata.X = adata.X.astype(np.float32)
    else:
        adata.X = np.asarray(adata.X, dtype=np.float32)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000)
    sc.tl.pca(adata, n_comps=40, svd_solver="arpack", use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40, use_rep="X_pca")
    sc.tl.umap(adata, random_state=0)
    sc.tl.leiden(adata, resolution=0.8, flavor="igraph", directed=False, random_state=0)
    return adata


def _with_pseudotime(adata):
    sc.tl.diffmap(adata)
    adata.uns["iroot"] = int(np.argmin(adata.obsm["X_umap"][:, 1]))
    sc.tl.dpt(adata)
    return adata


def _load_dataset(preferred_tenx_dir: Path) -> tuple[Path, str]:
    if preferred_tenx_dir.exists():
        return preferred_tenx_dir, "A549_5K_Multiplex"
    fallback = Path(
        "/home/zerlinshen/singlecell_factory/data/raw/lung_carcinoma_3k_count/outs/filtered_feature_bc_matrix"
    )
    if fallback.exists():
        return fallback, "lung_carcinoma_3k_fallback"
    raise FileNotFoundError("No available 10x filtered_feature_bc_matrix found.")


def run_once(mode: str, tenx_dir: Path) -> tuple[RunMetrics, pd.DataFrame]:
    t0 = time.perf_counter()
    c0 = time.process_time()

    adata = sc.read_10x_mtx(str(tenx_dir), var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()

    if mode == "baseline":
        adata = baseline_qc(adata, 200, 7000, 20.0, 50.0, 3)
        adata = baseline_clustering(adata, 1e4, 3000, 15, 40, 0.8, 0)
    else:
        adata = _optimized_qc(adata)
        adata = _optimized_clustering(adata)

    adata = _with_pseudotime(adata)
    if mode == "baseline":
        velocity = baseline_pseudo_velocity(adata, n_neighbors=15)
    else:
        velocity = compute_pseudo_velocity_vectorized(
            adata.obsm["X_umap"], adata.obs["dpt_pseudotime"].to_numpy(), n_neighbors=15
        )
    adata.obsm["X_pseudo_velocity"] = velocity
    speed = np.sqrt((velocity**2).sum(axis=1))
    adata.obs["pseudo_velocity_speed"] = speed

    wall = time.perf_counter() - t0
    cpu = time.process_time() - c0
    cpu_util = (cpu / wall) / max(os.cpu_count() or 1, 1) * 100
    rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    metrics = RunMetrics(
        mode=mode,
        dataset_used=str(tenx_dir),
        wall_seconds=wall,
        cpu_seconds=cpu,
        cpu_utilization_pct_all_cores=cpu_util,
        peak_rss_mb=rss_mb,
        n_cells=int(adata.n_obs),
        n_genes=int(adata.n_vars),
        n_clusters=int(adata.obs["leiden"].nunique()),
    )
    labels = pd.DataFrame({"barcode": adata.obs_names, "leiden": adata.obs["leiden"].astype(str).values})
    return metrics, labels


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--a549-tenx-dir",
        default="/home/zerlinshen/singlecell_factory/data/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count/outs/filtered_feature_bc_matrix",
    )
    parser.add_argument(
        "--output-dir",
        default="/home/zerlinshen/singlecell_factory/output/performance_benchmark",
    )
    parser.add_argument("--single-mode", choices=["baseline", "optimized"], default="")
    parser.add_argument("--metrics-path", default="")
    parser.add_argument("--labels-path", default="")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tenx_dir, dataset_name = _load_dataset(Path(args.a549_tenx_dir))

    if args.single_mode:
        metrics, labels = run_once(args.single_mode, tenx_dir)
        if args.metrics_path:
            Path(args.metrics_path).write_text(
                json.dumps(asdict(metrics), indent=2, ensure_ascii=False), encoding="utf-8"
            )
        if args.labels_path:
            labels.to_csv(args.labels_path, index=False)
        print(args.metrics_path or "single-mode done")
        return

    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        base_m = td_path / "base_metrics.json"
        base_l = td_path / "base_labels.csv"
        opt_m = td_path / "opt_metrics.json"
        opt_l = td_path / "opt_labels.csv"

        py = sys.executable
        base_cmd = [
            py,
            str(Path(__file__).resolve()),
            "--a549-tenx-dir",
            str(tenx_dir),
            "--output-dir",
            str(out_dir),
            "--single-mode",
            "baseline",
            "--metrics-path",
            str(base_m),
            "--labels-path",
            str(base_l),
        ]
        opt_cmd = [
            py,
            str(Path(__file__).resolve()),
            "--a549-tenx-dir",
            str(tenx_dir),
            "--output-dir",
            str(out_dir),
            "--single-mode",
            "optimized",
            "--metrics-path",
            str(opt_m),
            "--labels-path",
            str(opt_l),
        ]
        subprocess.run(base_cmd, check=True)
        subprocess.run(opt_cmd, check=True)

        base_metrics = RunMetrics(**json.loads(base_m.read_text(encoding="utf-8")))
        opt_metrics = RunMetrics(**json.loads(opt_m.read_text(encoding="utf-8")))
        base_labels = pd.read_csv(base_l)
        opt_labels = pd.read_csv(opt_l)

    common = base_labels.merge(opt_labels, on="barcode", suffixes=("_base", "_opt"))
    ari = adjusted_rand_score(common["leiden_base"], common["leiden_opt"]) if len(common) else float("nan")
    speedup = (base_metrics.wall_seconds - opt_metrics.wall_seconds) / base_metrics.wall_seconds * 100
    mem_save = (base_metrics.peak_rss_mb - opt_metrics.peak_rss_mb) / base_metrics.peak_rss_mb * 100

    compare = {
        "dataset_name": dataset_name,
        "baseline": asdict(base_metrics),
        "optimized": asdict(opt_metrics),
        "cluster_label_consistency_ari": ari,
        "speed_improvement_pct": speedup,
        "memory_reduction_pct": mem_save,
        "meets_speed_target_50pct": bool(speedup >= 50.0),
        "meets_memory_target_30pct": bool(mem_save >= 30.0),
    }

    (out_dir / "benchmark_compare.json").write_text(
        json.dumps(compare, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    pd.DataFrame([asdict(base_metrics), asdict(opt_metrics)]).to_csv(
        out_dir / "benchmark_runs.csv", index=False
    )
    print(out_dir / "benchmark_compare.json")


if __name__ == "__main__":
    main()
