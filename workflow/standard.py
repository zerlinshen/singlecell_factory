from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal
import argparse
import time
import json

import scanpy as sc


@dataclass(frozen=True)
class StandardWorkflowConfig:
    """Standard single-cell workflow configuration.

    Parameters
    ----------
    project
        Project name used for output folder naming.
    tenx_dir
        Path to Cell Ranger `outs/filtered_feature_bc_matrix/`.
    output_dir
        Output base directory.
    mode
        Controls output density and performance tradeoffs.
    n_jobs
        Number of CPU threads for Scanpy/Numba-backed operations.
    use_cache
        Reuse processed output when available.
    """

    project: str
    tenx_dir: Path
    output_dir: Path
    mode: Literal["fast", "full"] = "full"
    n_jobs: int = 8
    use_cache: bool = True


def run_standard_workflow(cfg: StandardWorkflowConfig) -> Path:
    """Run the standard Scanpy workflow.

    Notes
    -----
    This is a thin wrapper around the existing optimized pipeline script,
    intended for decoupled execution. It does not depend on the velocity workflow.
    """

    if not cfg.tenx_dir.exists():
        raise FileNotFoundError(f"tenx_dir not found: {cfg.tenx_dir}")

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    out_project = cfg.output_dir / cfg.project
    out_project.mkdir(parents=True, exist_ok=True)
    out = out_project / "processed.h5ad"
    stats_path = out_project / "runtime_stats.json"
    if cfg.use_cache and out.exists():
        return out

    if hasattr(sc, "settings"):
        sc.settings.n_jobs = max(1, cfg.n_jobs)
    t0 = time.perf_counter()
    adata = sc.read_10x_mtx(str(cfg.tenx_dir), var_names="gene_symbols", cache=False)
    if "gene_symbols" in adata.var.columns:
        gene_symbols = adata.var["gene_symbols"].astype(str).str.upper()
    else:
        gene_symbols = adata.var_names.astype(str).str.upper()
    adata.var["mt"] = gene_symbols.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata = adata[adata.obs["n_genes_by_counts"] > 200, :].copy()
    adata = adata[adata.obs["pct_counts_mt"] < 5, :].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, random_state=0)
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph", directed=False, random_state=0)
    adata.write(out)
    elapsed = time.perf_counter() - t0
    stats_path.write_text(
        json.dumps(
            {
                "project": cfg.project,
                "elapsed_seconds": elapsed,
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "n_jobs": cfg.n_jobs,
                "mode": cfg.mode,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    return out


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for standard workflow."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--tenx-dir", required=True)
    parser.add_argument("--output-dir", default="/home/zerlinshen/singlecell_factory/output/workflow_standard")
    parser.add_argument("--mode", choices=["fast", "full"], default="full")
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--no-cache", action="store_true")
    return parser.parse_args()


def main() -> None:
    """CLI entrypoint for standard workflow."""

    args = parse_args()
    cfg = StandardWorkflowConfig(
        project=args.project,
        tenx_dir=Path(args.tenx_dir),
        output_dir=Path(args.output_dir),
        mode=args.mode,
        n_jobs=args.n_jobs,
        use_cache=not args.no_cache,
    )
    out = run_standard_workflow(cfg)
    print(out)


if __name__ == "__main__":
    main()
