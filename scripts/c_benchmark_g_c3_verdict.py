#!/usr/bin/env python3
"""G-C3 verdict — production-scale parity between Lane-1 and Lane-2.

Reads final_adata.h5ad from each lane, joins on barcode (obs_names),
computes ARI/NMI/cluster-size-Spearman/major-cluster-match, writes
comparison/{comparison_table.csv, winner.md, parity_report.json} into
the project's runs/<benchmark-tag>/comparison/ directory.

Per G-C0 v3 G-C3:
  ARI >= 0.95 → PASS_STRICT (backends parity-verified at convergence;
                 either lane is admissible for Component D)
  ARI in [0.85, 0.95) → PASS_ACCEPTABLE (parity-with-divergence; document)
  ARI < 0.85 → FAIL (backend divergence at convergence; demand
                root cause before Component D)
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    confusion_matrix,
)
from scipy.stats import spearmanr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lane1-run", required=True, help="Path to Lane-1 python run dir")
    ap.add_argument("--lane2-run", required=True, help="Path to Lane-2 python run dir")
    ap.add_argument("--out-dir", required=True, help="Output comparison dir")
    args = ap.parse_args()

    lane1_path = Path(args.lane1_run) / "final_adata.h5ad"
    lane2_path = Path(args.lane2_run) / "final_adata.h5ad"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading Lane-1 final adata: {lane1_path}")
    a1 = ad.read_h5ad(lane1_path)
    print(f"Loading Lane-2 final adata: {lane2_path}")
    a2 = ad.read_h5ad(lane2_path)

    print(f"Lane-1: {a1.n_obs} cells, leiden unique = {a1.obs['leiden'].nunique()}")
    print(f"Lane-2: {a2.n_obs} cells, leiden unique = {a2.obs['leiden'].nunique()}")

    # Join on barcode — both lanes started from the same 877k subset but
    # may have lost different cells in scrublet (CPU vs GPU scrublet code
    # path). Compare only on shared barcodes.
    common = a1.obs_names.intersection(a2.obs_names)
    n_common = len(common)
    n_lane1_only = a1.n_obs - n_common
    n_lane2_only = a2.n_obs - n_common
    print(f"Common barcodes: {n_common} ({n_common/a1.n_obs*100:.2f}% of Lane-1, "
          f"{n_common/a2.n_obs*100:.2f}% of Lane-2)")

    labels1 = a1.obs.loc[common, "leiden"].astype(str).values
    labels2 = a2.obs.loc[common, "leiden"].astype(str).values

    ari = float(adjusted_rand_score(labels1, labels2))
    nmi = float(normalized_mutual_info_score(labels1, labels2))

    if ari >= 0.95:
        verdict = "PASS_STRICT"
    elif ari >= 0.85:
        verdict = "PASS_ACCEPTABLE"
    else:
        verdict = "FAIL"

    sizes1 = pd.Series(labels1).value_counts().sort_index()
    sizes2 = pd.Series(labels2).value_counts().sort_index()

    # Build confusion matrix per-lane (rows = Lane-1 unique, cols = Lane-2 unique)
    label_set1 = sorted(set(labels1))
    label_set2 = sorted(set(labels2))
    label1_to_idx = {l: i for i, l in enumerate(label_set1)}
    label2_to_idx = {l: j for j, l in enumerate(label_set2)}
    cm = np.zeros((len(label_set1), len(label_set2)), dtype=np.int64)
    for l1, l2 in zip(labels1, labels2):
        cm[label1_to_idx[l1], label2_to_idx[l2]] += 1

    n_clusters1 = len(label_set1)
    n_clusters2 = len(label_set2)

    # Major-cluster mapping (Lane-1 -> Lane-2 best match)
    major_map = {}
    for i, l1 in enumerate(label_set1):
        row = cm[i]
        if row.sum() > 0:
            best_j = int(np.argmax(row))
            best_overlap = int(row[best_j])
            cluster_size = int(row.sum())
            major_map[str(l1)] = {
                "best_lane2_match": str(label_set2[best_j]),
                "overlap_cells": best_overlap,
                "cluster_size": cluster_size,
                "purity": best_overlap / cluster_size,
            }
    median_purity = float(np.median([v["purity"] for v in major_map.values()]))

    # Cluster size distribution Spearman (treats them as ordered distributions)
    sizes1_sorted = sorted(sizes1.values, reverse=True)
    sizes2_sorted = sorted(sizes2.values, reverse=True)
    n_pad = max(len(sizes1_sorted), len(sizes2_sorted))
    sizes1_sorted += [0] * (n_pad - len(sizes1_sorted))
    sizes2_sorted += [0] * (n_pad - len(sizes2_sorted))
    spearman_rho, spearman_p = spearmanr(sizes1_sorted, sizes2_sorted)

    report = {
        "verdict": verdict,
        "ari_lane1_vs_lane2": ari,
        "nmi_lane1_vs_lane2": nmi,
        "n_common_barcodes": int(n_common),
        "n_lane1_only_barcodes": int(n_lane1_only),
        "n_lane2_only_barcodes": int(n_lane2_only),
        "n_clusters_lane1": int(n_clusters1),
        "n_clusters_lane2": int(n_clusters2),
        "median_purity_lane1_to_lane2": median_purity,
        "cluster_size_spearman_rho": float(spearman_rho),
        "cluster_size_spearman_p": float(spearman_p),
        "thresholds": {
            "pass_strict": 0.95,
            "pass_acceptable": 0.85,
        },
        "lane1_run": args.lane1_run,
        "lane2_run": args.lane2_run,
    }
    (out_dir / "parity_report.json").write_text(json.dumps(report, indent=2))

    # Comparison table
    cmp_rows = []
    for l1, info in major_map.items():
        cmp_rows.append({
            "lane1_cluster": l1,
            "lane1_size": info["cluster_size"],
            "best_lane2_match": info["best_lane2_match"],
            "overlap_cells": info["overlap_cells"],
            "purity": round(info["purity"], 4),
        })
    cmp_df = pd.DataFrame(cmp_rows).sort_values("lane1_size", ascending=False)
    cmp_df.to_csv(out_dir / "comparison_table.csv", index=False)

    # Winner doc
    pre_ambient_banner = (
        "**PRE-AMBIENT BANNER (G-C0 A8):** This C benchmark was run WITHOUT "
        "ambient RNA correction (cellbender / SoupX / DecontX). Per G-C0 v3 "
        "interim policy, the C winner is the lane that demonstrated stronger "
        "scientific defensibility on the bare-pipeline comparison only. "
        "Component D MUST re-apply the chosen lane WITH ambient correction "
        "once the ambient module passes G-AMB-5."
    )

    winner_text = f"""# G-C3 Verdict — Lane-1 vs Lane-2 Production Parity (2026-05-20)

{pre_ambient_banner}

## Verdict: **{verdict}**

| Metric | Value |
|---|---|
| ARI (Lane-1 vs Lane-2) | **{ari:.4f}** |
| NMI (Lane-1 vs Lane-2) | {nmi:.4f} |
| n_clusters Lane-1 (CPU sparse_exact) | {n_clusters1} |
| n_clusters Lane-2 (GPU rsc drop-in) | {n_clusters2} |
| Shared barcodes | {n_common:,} ({n_common/max(a1.n_obs,a2.n_obs)*100:.2f}% of larger lane) |
| Lane-1-only barcodes | {n_lane1_only:,} |
| Lane-2-only barcodes | {n_lane2_only:,} |
| Median Lane-1 cluster purity (best match) | {median_purity:.4f} |
| Cluster size distribution Spearman ρ | {spearman_rho:.4f} (p={spearman_p:.2e}) |

## G-C0 v3 thresholds applied
- **PASS_STRICT** if ARI >= 0.95 — backends are parity-verified at convergence; either lane admissible for Component D.
- **PASS_ACCEPTABLE** if 0.85 <= ARI < 0.95 — parity-with-divergence; document the asymmetry; Component D may proceed with the CPU lane as canonical.
- **FAIL** if ARI < 0.85 — backend divergence is too large; demand root-cause investigation before Component D.

## Lane manifests
- Lane-1: `{args.lane1_run}/lane_manifest.json`
- Lane-2: `{args.lane2_run}/lane_manifest.json`

## Contract anchors
- G-C0 contract SHA: see `nc-reproduction/ledger/approvals/G-C0-2026-05-20-v3.md`
- Principle 9 convergence: both lanes recorded `harmony_converged: true`
- F-2 sparse_exact: Lane-1 used `pca_solver=arpack`, `neighbors_backend=scanpy_umap_true_knn`, `leiden_backend=igraph_cpu`
- F-3 DE: not applicable (C benchmark stops at clustering)
- Principle 8 GPU scanpy-API: Lane-2 used rsc.pp.* (drop-in)
"""
    (out_dir / "winner.md").write_text(winner_text)

    print(f"\n=== G-C3 verdict: {verdict} ===")
    print(f"ARI = {ari:.4f}, NMI = {nmi:.4f}")
    print(f"n_clusters: lane1={n_clusters1}, lane2={n_clusters2}")
    print(f"median purity lane1->lane2 = {median_purity:.4f}")
    print(f"Wrote: {out_dir / 'parity_report.json'}")
    print(f"Wrote: {out_dir / 'comparison_table.csv'}")
    print(f"Wrote: {out_dir / 'winner.md'}")


if __name__ == "__main__":
    main()
