"""Compute per-substate DE evidence for AC-2.

For each substate surfaced by context_aware_annotation, run Mann-Whitney U
on its curated defining markers (substate cells vs rest), BH-adjust FDR
across markers, and write the DE CSV that T8 (compare_marker_intelligence.py)
needs for AC-2 verification.

Output schema (per CSV at differential_expression/substates/<substate>.csv):
    gene,log2fc,pvalue,fdr,pct_in_substate,pct_in_parent
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
import scipy.stats as stats
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger("substate_de")


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvalues)
    if n == 0:
        return np.array([])
    order = np.argsort(pvalues)
    ranks = np.arange(1, n + 1)
    adj = pvalues[order] * n / ranks
    # Enforce monotonicity from the largest
    adj_rev = adj[::-1].copy()
    adj_rev = np.minimum.accumulate(adj_rev)
    adj = adj_rev[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--checkpoint",
        default="/home/zerlinshen/projects/nc-reproduction/runs/2026-05-14T1738Z-d192836/python/NC2024_NSCLC_CANDIDATE_20260515_013812/.checkpoints/after_annotation.h5ad",
    )
    parser.add_argument(
        "--run-dir",
        default="/home/zerlinshen/projects/nc-reproduction/runs/2026-05-14T1738Z-d192836/python/NC2024_NSCLC_CANDIDATE_20260515_013812",
    )
    parser.add_argument(
        "--curated-yaml",
        default="/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ref/markers/curated/lung_NSCLC.yaml",
    )
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    candidate_json_path = run_dir / "metrics" / "marker_intelligence_candidate.json"
    candidate = json.loads(candidate_json_path.read_text(encoding="utf-8"))

    # Read curated YAML to get defining markers per substate
    curated = yaml.safe_load(Path(args.curated_yaml).read_text(encoding="utf-8"))
    curated_markers = curated.get("markers", {})

    # Build cluster -> substate mapping from candidate JSON
    per_cluster = candidate.get("per_cluster", [])
    cluster_substate: dict[str, str | None] = {}
    cluster_celltype: dict[str, str] = {}
    for entry in per_cluster:
        cid = str(entry["cluster_id"])
        cluster_substate[cid] = entry.get("substate")
        cluster_celltype[cid] = entry["cell_type"]

    surfaced_substates = sorted({v for v in cluster_substate.values() if v})
    if not surfaced_substates:
        logger.warning("no substates surfaced in candidate JSON — AC-2 will FAIL")
        return 0

    logger.info("surfaced substates: %s", surfaced_substates)
    logger.info("loading %s ...", args.checkpoint)
    adata = ad.read_h5ad(args.checkpoint)
    logger.info("loaded: shape=%s", adata.shape)

    cluster_labels = adata.obs["leiden"].astype(str).values

    # Build cell-level substate label
    n_cells = adata.shape[0]
    cell_substate = np.array([cluster_substate.get(c) or "" for c in cluster_labels])

    de_dir = run_dir / "differential_expression" / "substates"
    de_dir.mkdir(parents=True, exist_ok=True)

    var_names = adata.var_names
    gene_to_idx = {g: i for i, g in enumerate(var_names)}
    X = adata.X
    is_sparse = sp.issparse(X)

    substate_records: list[dict] = []
    for substate in surfaced_substates:
        in_mask = cell_substate == substate
        out_mask = ~in_mask
        n_in = int(in_mask.sum())
        n_out = int(out_mask.sum())
        logger.info("  substate=%s: %d cells (vs %d others)", substate, n_in, n_out)

        # Find parent cell type (substate name format: "<parent> <suffix>")
        parent_ct = next(
            (entry["cell_type"] for entry in per_cluster if entry.get("substate") == substate),
            "Unknown",
        )

        # Defining markers from curated YAML
        defining = curated_markers.get(substate, {}).get("positive", []) if isinstance(curated_markers.get(substate, {}), dict) else []
        defining_in_dataset = [g for g in defining if g in gene_to_idx]
        if not defining_in_dataset:
            logger.warning("  no curated markers for %s found in dataset; skipping DE", substate)
            continue

        # Compute per-marker MWU
        rows = []
        pvals = []
        for gene in defining_in_dataset:
            j = gene_to_idx[gene]
            if is_sparse:
                col = np.asarray(X[:, j].todense()).flatten().astype(np.float32)
            else:
                col = np.asarray(X[:, j], dtype=np.float32).flatten()
            in_vals = col[in_mask]
            out_vals = col[out_mask]
            try:
                u, p = stats.mannwhitneyu(in_vals, out_vals, alternative="greater")
            except ValueError:
                # constant values etc.
                u, p = 0.0, 1.0
            mean_in = float(in_vals.mean())
            mean_out = float(out_vals.mean())
            # log2 fold-change (add pseudocount)
            log2fc = float(np.log2((mean_in + 1e-9) / (mean_out + 1e-9)))
            pct_in = float((in_vals > 0).sum()) / max(1, n_in)
            pct_out = float((out_vals > 0).sum()) / max(1, n_out)
            rows.append((gene, log2fc, p, mean_in, mean_out, pct_in, pct_out))
            pvals.append(p)

        # BH FDR adjust across this substate's markers
        fdrs = bh_fdr(np.array(pvals))
        min_fdr = float(fdrs.min()) if len(fdrs) else 1.0

        # Write CSV
        csv_path = de_dir / f"{substate}.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["gene", "log2fc", "pvalue", "fdr", "pct_in_substate", "pct_in_parent"])
            for (gene, log2fc, p, mean_in, mean_out, pct_in, pct_out), fdr in zip(rows, fdrs):
                w.writerow([gene, f"{log2fc:.6f}", f"{p:.6e}", f"{fdr:.6e}", f"{pct_in:.4f}", f"{pct_out:.4f}"])
        logger.info("  wrote %s (min_fdr=%.3e)", csv_path, min_fdr)

        substate_records.append({
            "name": substate,
            "parent_cell_type": parent_ct,
            "defining_markers": defining_in_dataset,
            "fdr": min_fdr,
        })

    # Update candidate JSON
    candidate["condition_specific_substates"] = substate_records
    candidate_json_path.write_text(json.dumps(candidate, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("updated candidate JSON with %d substates: %s", len(substate_records), candidate_json_path)

    # Summary
    print()
    print("=== AC-2 evidence summary ===")
    for rec in substate_records:
        passed = rec["fdr"] < args.fdr_threshold
        print(f"  {rec['name']:<28} parent={rec['parent_cell_type']:<8} min_fdr={rec['fdr']:.3e}  {'PASS' if passed else 'FAIL'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
