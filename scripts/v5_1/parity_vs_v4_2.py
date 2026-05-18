"""MV1a parity check: compare new canonical run vs v4.2 baseline.

Usage:
    python scripts/v5_1/parity_vs_v4_2.py \
        --new-adata <path/to/final_adata.h5ad> \
        --out-dir <path/to/parity/> \
        [--v42-adata <path>]  # optional; defaults to v4.2 canonical run

Acceptance criteria (AC-V51-MV1a-2):
    - n_leiden_clusters: |new - 53| / 53 <= 0.05  (range [50, 56])
    - cell_type_ari (RNA-only ARI vs seurat_clusters): |new - 0.1683| <= 0.05 absolute
    - cell_count: new == 8981 (exact)
    - peak_gene_overlap_top1000_strict: Jaccard(new top-1000, v4.2 top-1000) == 1.0
      (exact replay of the v4.2 peak-gene ranking; Trevino S2F strict overlap
      remains 0.236 in the v4.2/v5.0 ledgers and is a different metric).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata
import numpy as np
from sklearn.metrics import adjusted_rand_score


V42_BASELINES = {
    "n_leiden_clusters": 53,
    "cell_type_ari": 0.1683,
    "cell_count": 8981,
    "peak_gene_overlap_top1000_strict": 1.0,
}

V42_ADATA_PATH = (
    "/home/zerlinshen/projects/wave5-trevino/runs/"
    "20260516T0931Z-d192836f1bb0/python/final_adata.h5ad"
)


def _compute_ari(adata: anndata.AnnData) -> float:
    if "leiden" not in adata.obs or "seurat_clusters" not in adata.obs:
        return float("nan")
    return float(adjusted_rand_score(
        adata.obs["seurat_clusters"].astype(str),
        adata.obs["leiden"].astype(str),
    ))


def _compute_peak_gene_overlap(adata: anndata.AnnData, v42_adata: anndata.AnnData | None) -> float | None:
    """Jaccard overlap of top-1000 peak-gene pairs between new and v4.2."""
    if "peak_to_gene_top1000" not in adata.uns:
        return None
    new_top = adata.uns["peak_to_gene_top1000"]
    if v42_adata is None or "peak_to_gene_top1000" not in v42_adata.uns:
        return None
    v42_top = v42_adata.uns["peak_to_gene_top1000"]

    def pair_set(df):
        return set(zip(df["peak_idx"].astype(int), df["gene"].astype(str)))

    new_pairs = pair_set(new_top)
    v42_pairs = pair_set(v42_top)
    if not new_pairs or not v42_pairs:
        return None
    overlap = len(new_pairs & v42_pairs) / len(new_pairs | v42_pairs)
    return round(overlap, 4)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--new-adata", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--v42-adata", type=Path, default=Path(V42_ADATA_PATH))
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading new adata: {args.new_adata}", flush=True)
    new = anndata.read_h5ad(args.new_adata)

    v42 = None
    if args.v42_adata.exists():
        print(f"Loading v4.2 adata: {args.v42_adata}", flush=True)
        v42 = anndata.read_h5ad(args.v42_adata)
    else:
        print(f"WARNING: v4.2 adata not found at {args.v42_adata}; peak_gene overlap will use uns only", flush=True)

    # Compute metrics
    new_n_clusters = int(new.obs["leiden"].nunique()) if "leiden" in new.obs else None
    new_cell_count = int(new.n_obs)
    new_ari = _compute_ari(new)
    new_overlap = _compute_peak_gene_overlap(new, v42)

    baselines = V42_BASELINES.copy()

    results = {}
    all_pass = True

    # cell_count: exact equal
    cell_count_pass = new_cell_count == baselines["cell_count"]
    results["cell_count"] = {
        "v42": baselines["cell_count"],
        "v51": new_cell_count,
        "diff_abs": new_cell_count - baselines["cell_count"],
        "threshold": "exact_equal",
        "pass": cell_count_pass,
    }
    if not cell_count_pass:
        all_pass = False

    # n_leiden_clusters: within 5%
    if new_n_clusters is not None:
        rel_diff = abs(new_n_clusters - baselines["n_leiden_clusters"]) / baselines["n_leiden_clusters"]
        cluster_pass = rel_diff <= 0.05
        results["n_leiden_clusters"] = {
            "v42": baselines["n_leiden_clusters"],
            "v51": new_n_clusters,
            "diff_pct": round(rel_diff * 100, 2),
            "threshold_pct": 5.0,
            "pass": cluster_pass,
        }
        if not cluster_pass:
            all_pass = False
    else:
        results["n_leiden_clusters"] = {"error": "leiden not in obs", "pass": False}
        all_pass = False

    # cell_type_ari: within 0.05 absolute
    if not np.isnan(new_ari):
        ari_diff = abs(new_ari - baselines["cell_type_ari"])
        ari_pass = ari_diff <= 0.05
        results["cell_type_ari"] = {
            "v42": baselines["cell_type_ari"],
            "v51": round(new_ari, 4),
            "diff_abs": round(ari_diff, 4),
            "threshold_abs": 0.05,
            "pass": ari_pass,
        }
        if not ari_pass:
            all_pass = False
    else:
        results["cell_type_ari"] = {"error": "leiden or seurat_clusters missing", "pass": False}
        all_pass = False

    # peak_gene_overlap: exact equal (input-bounded)
    if new_overlap is not None:
        overlap_pass = new_overlap == baselines["peak_gene_overlap_top1000_strict"]
        results["peak_gene_overlap_top1000_strict"] = {
            "v42": baselines["peak_gene_overlap_top1000_strict"],
            "v51": new_overlap,
            "threshold": "exact_equal",
            "pass": overlap_pass,
            "note": "Jaccard overlap of top-1000 peak-gene pairs; input-bounded",
        }
        if not overlap_pass:
            all_pass = False
    else:
        results["peak_gene_overlap_top1000_strict"] = {
            "v51": None,
            "note": "peak_to_gene_top1000 not in uns or v4.2 adata unavailable",
            "pass": False,
        }
        all_pass = False

    report = {
        "parity_verdict": "PASS" if all_pass else "FAIL",
        "v42_run_id": "20260516T0931Z-d192836f1bb0",
        "v42_conda_env_hash": "288712b5192eaad5294a0508f085710b5cd88e146e0bbb10a46987fbafa661a3",
        "metrics": results,
    }

    out_path = args.out_dir / "parity_vs_v4_2.json"
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nParity report written to {out_path}", flush=True)
    print(f"Verdict: {report['parity_verdict']}", flush=True)
    for k, v in results.items():
        status = "PASS" if v.get("pass") else "FAIL"
        print(f"  {k}: {status}  {v}", flush=True)

    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
