#!/usr/bin/env python3
"""Render Trevino Figure 3B public-matrix ATAC transferred pseudotime outputs.

This uses a bounded, auditable cross-modality transfer:
- reference: Figure 3A RNA velocity subset with velocity_pseudotime
- query: public ATAC gene-activity matrix from GSE162170
- features: fetal-brain marker genes from the project marker config
- method: cosine kNN in within-modality z-scored shared gene activity/expression

Boundary: this is public processed-matrix reproduction, not FASTQ/fragments-level
scATAC reprocessing and not the original author's exact integration model.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from scipy import sparse  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.neighbors import NearestNeighbors  # noqa: E402


TRAJECTORY_MARKERS = ["SOX2", "PAX6", "EOMES", "NEUROD1", "NEUROD6", "SLC17A7", "SATB2", "GAD2"]
EXPECTED_NEGATIVE = ["SOX2", "PAX6"]
EXPECTED_POSITIVE = ["NEUROD1", "NEUROD6", "SLC17A7", "SATB2"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "png": out_dir / f"{stem}.png",
        "svg": out_dir / f"{stem}.svg",
        "pdf": out_dir / f"{stem}.pdf",
    }
    fig.savefig(paths["png"], dpi=300, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {f"{k}_path": str(v) for k, v in paths.items()} | {
        f"{k}_sha256": sha256_file(v) for k, v in paths.items()
    }


def flatten_marker_config(path: Path) -> list[str]:
    payload = json.loads(path.read_text())
    genes: list[str] = []
    for values in payload.values():
        genes.extend(map(str, values))
    genes.extend(TRAJECTORY_MARKERS)
    # Stable order, deduplicated.
    return list(dict.fromkeys(genes))


def read_atac_gene_activity(path: Path, wanted_genes: Iterable[str]) -> tuple[list[str], list[str], np.ndarray, list[str]]:
    wanted = set(map(str, wanted_genes))
    found: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().strip().split()
        for line in fh:
            if not line:
                continue
            gene, rest = line.split(maxsplit=1)
            if gene not in wanted:
                continue
            values = np.fromstring(rest, sep=" ", dtype=np.float32)
            if values.shape[0] != len(header):
                raise ValueError(f"{gene} has {values.shape[0]} values, expected {len(header)}")
            found[gene] = values
            if len(found) == len(wanted):
                break
    genes = [g for g in wanted_genes if g in found]
    missing = [g for g in wanted_genes if g not in found]
    if not genes:
        raise ValueError("No requested marker genes found in ATAC gene activity matrix")
    matrix = np.column_stack([found[g] for g in genes]).astype(np.float32, copy=False)
    return header, genes, matrix, missing


def expression_matrix(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        raise ValueError("No selected transfer genes are present in RNA reference")
    matrix = adata[:, present].X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32), present


def zscore_columns(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float32)
    mean = np.nanmean(matrix, axis=0, keepdims=True)
    std = np.nanstd(matrix, axis=0, keepdims=True)
    std[std == 0] = 1.0
    out = (matrix - mean) / std
    out[~np.isfinite(out)] = 0.0
    return out.astype(np.float32, copy=False)


def choose_atac_cells(meta: pd.DataFrame, max_cells: int, seed: int) -> list[str]:
    if len(meta) <= max_cells:
        return list(map(str, meta.index))
    rng = np.random.default_rng(seed)
    chosen: list[str] = []
    group_col = "Age" if "Age" in meta.columns else "Iterative.LSI.Clusters"
    groups = meta.groupby(group_col, observed=True).indices
    per_group = max(100, max_cells // max(1, len(groups)))
    for _, idx in sorted(groups.items(), key=lambda kv: str(kv[0])):
        idx = np.asarray(idx)
        take = min(len(idx), per_group)
        chosen.extend(meta.index[rng.choice(idx, size=take, replace=False)].astype(str).tolist())
    if len(chosen) > max_cells:
        chosen = rng.choice(np.asarray(chosen), size=max_cells, replace=False).astype(str).tolist()
    elif len(chosen) < max_cells:
        remaining = np.setdiff1d(meta.index.astype(str), np.asarray(chosen), assume_unique=False)
        if len(remaining):
            chosen.extend(rng.choice(remaining, size=min(max_cells - len(chosen), len(remaining)), replace=False).astype(str).tolist())
    return chosen


def marker_direction_metrics(adata: ad.AnnData, genes: list[str]) -> tuple[dict[str, dict], list[dict]]:
    pt = np.asarray(adata.obs["transferred_pseudotime"], dtype=float)
    metrics: dict[str, dict] = {}
    checks: list[dict] = []
    for gene in genes:
        if gene not in adata.var_names:
            metrics[gene] = {"present": False}
            continue
        expr = adata[:, gene].X
        if sparse.issparse(expr):
            expr = expr.toarray()
        expr = np.asarray(expr).ravel().astype(float)
        valid = np.isfinite(pt) & np.isfinite(expr)
        rho, pval = spearmanr(pt[valid], expr[valid]) if valid.sum() >= 10 else (np.nan, np.nan)
        metrics[gene] = {
            "present": True,
            "valid_cells": int(valid.sum()),
            "spearman_rho": None if np.isnan(rho) else float(rho),
            "spearman_pvalue": None if np.isnan(pval) else float(pval),
        }
    for gene in EXPECTED_NEGATIVE:
        rho = metrics.get(gene, {}).get("spearman_rho")
        checks.append({"marker": gene, "expected": "negative", "rho": rho, "pass": rho is not None and rho < 0})
    for gene in EXPECTED_POSITIVE:
        rho = metrics.get(gene, {}).get("spearman_rho")
        checks.append({"marker": gene, "expected": "positive", "rho": rho, "pass": rho is not None and rho > 0})
    return metrics, checks


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rna-velocity-h5ad", required=True, type=Path)
    p.add_argument("--atac-gene-activity", required=True, type=Path)
    p.add_argument("--atac-metadata", required=True, type=Path)
    p.add_argument("--marker-json", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--max-atac-cells", type=int, default=12000)
    p.add_argument("--neighbors", type=int, default=15)
    p.add_argument("--seed", type=int, default=13)
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    marker_genes = flatten_marker_config(args.marker_json)
    rna = ad.read_h5ad(args.rna_velocity_h5ad)
    if "velocity_pseudotime" not in rna.obs.columns:
        raise KeyError("RNA reference missing obs['velocity_pseudotime']")
    atac_cells, atac_genes, atac_matrix, atac_missing = read_atac_gene_activity(args.atac_gene_activity, marker_genes)
    meta = pd.read_csv(args.atac_metadata, sep="\t", compression="gzip", dtype=str).set_index("Cell.ID")
    meta = meta.reindex(atac_cells)

    selected_cells = choose_atac_cells(meta, args.max_atac_cells, args.seed)
    selected_idx = meta.index.get_indexer(selected_cells)
    atac_matrix = atac_matrix[selected_idx, :]
    meta = meta.iloc[selected_idx].copy()

    rna_matrix, rna_genes = expression_matrix(rna, atac_genes)
    shared_genes = [g for g in atac_genes if g in rna_genes]
    rna_idx = [rna_genes.index(g) for g in shared_genes]
    atac_idx = [atac_genes.index(g) for g in shared_genes]
    if len(shared_genes) < 10:
        raise ValueError(f"Too few shared transfer genes: {len(shared_genes)}")
    rna_features = zscore_columns(rna_matrix[:, rna_idx])
    atac_features = zscore_columns(np.log1p(np.maximum(atac_matrix[:, atac_idx], 0.0)))
    pt = np.asarray(rna.obs["velocity_pseudotime"], dtype=float)
    finite = np.isfinite(pt)
    rna_features = rna_features[finite, :]
    pt = pt[finite]
    rna_cell_types = rna.obs.loc[finite, "cell_type"].astype(str).to_numpy() if "cell_type" in rna.obs else None

    nn = NearestNeighbors(n_neighbors=min(args.neighbors, len(pt)), metric="cosine")
    nn.fit(rna_features)
    distances, indices = nn.kneighbors(atac_features)
    weights = 1.0 / np.maximum(distances, 1e-6)
    transferred = (weights * pt[indices]).sum(axis=1) / weights.sum(axis=1)
    transferred_cell_type = []
    if rna_cell_types is not None:
        for neigh in indices:
            vals, counts = np.unique(rna_cell_types[neigh], return_counts=True)
            transferred_cell_type.append(str(vals[np.argmax(counts)]))

    atac = ad.AnnData(
        X=np.log1p(np.maximum(atac_matrix, 0.0)).astype(np.float32, copy=False),
        obs=meta,
        var=pd.DataFrame(index=pd.Index(atac_genes, name="gene")),
    )
    atac.obs["transferred_pseudotime"] = transferred
    atac.obs["mean_transfer_cosine_distance"] = distances.mean(axis=1)
    if transferred_cell_type:
        atac.obs["transferred_rna_cell_type"] = transferred_cell_type

    # Compute a reproducible public-matrix ATAC embedding from selected gene activities.
    sc.pp.scale(atac, max_value=10)
    sc.tl.pca(atac, n_comps=min(30, len(shared_genes) - 1), svd_solver="arpack", random_state=args.seed)
    sc.pp.neighbors(atac, n_neighbors=20, n_pcs=min(20, len(shared_genes) - 1), random_state=args.seed)
    sc.tl.umap(atac, random_state=args.seed)

    computed_h5ad = args.out_dir / "figure3b_atac_transferred_pseudotime_subset.h5ad"
    atac.write_h5ad(computed_h5ad, compression="gzip")

    figure_outputs: dict[str, dict[str, str]] = {}
    fig = sc.pl.embedding(
        atac,
        basis="umap",
        color="transferred_pseudotime",
        title="Trevino Fig. 3B reproduction — ATAC transferred pseudotime",
        show=False,
        return_fig=True,
        size=12,
        cmap="viridis",
    )
    figure_outputs["atac_transferred_pseudotime_umap"] = save_triple(
        fig, args.out_dir, "figure3b_atac_transferred_pseudotime_umap"
    )
    color = "Iterative.LSI.Clusters" if "Iterative.LSI.Clusters" in atac.obs.columns else None
    if color:
        fig = sc.pl.embedding(
            atac,
            basis="umap",
            color=color,
            title="Public ATAC gene-activity UMAP × iterative LSI clusters",
            show=False,
            return_fig=True,
            size=12,
        )
        figure_outputs["atac_cluster_umap"] = save_triple(fig, args.out_dir, "figure3b_atac_cluster_umap")
    if "Age" in atac.obs.columns:
        fig, ax = plt.subplots(figsize=(7, 4))
        order = sorted(atac.obs["Age"].dropna().unique())
        data = [atac.obs.loc[atac.obs["Age"] == age, "transferred_pseudotime"].astype(float).to_numpy() for age in order]
        ax.boxplot(data, labels=order, showfliers=False)
        ax.set_xlabel("ATAC sample age")
        ax.set_ylabel("transferred RNA velocity pseudotime")
        ax.set_title("Figure 3B check — transferred pseudotime by ATAC age")
        figure_outputs["age_pseudotime_boxplot"] = save_triple(fig, args.out_dir, "figure3b_age_pseudotime_boxplot")

    marker_metrics, direction_checks = marker_direction_metrics(atac, TRAJECTORY_MARKERS)
    validation_status = (
        "REPRODUCED_PUBLIC_MATRIX_SUBSET"
        if all(c["pass"] for c in direction_checks)
        else "GENERATED_PUBLIC_MATRIX_SUBSET_WEAK_MARKER_CONCORDANCE"
    )
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "truth_boundary": (
            "Public processed-matrix ATAC gene-activity transfer from 3A RNA velocity pseudotime; "
            "not FASTQ/fragments-level and not the original author integration model."
        ),
        "rna_velocity_h5ad": str(args.rna_velocity_h5ad),
        "atac_gene_activity": str(args.atac_gene_activity),
        "atac_metadata": str(args.atac_metadata),
        "marker_json": str(args.marker_json),
        "out_dir": str(args.out_dir),
        "computed_subset_h5ad": str(computed_h5ad),
        "parameters": {
            "max_atac_cells": args.max_atac_cells,
            "neighbors": args.neighbors,
            "seed": args.seed,
            "feature_scaling": "within-modality z-score on log1p ATAC gene activity and RNA reference expression",
            "transfer_metric": "cosine kNN, inverse-distance weighted pseudotime mean",
        },
        "status": validation_status,
        "atac_cells_total": len(atac_cells),
        "atac_cells_used": int(atac.n_obs),
        "rna_reference_cells_total": int(rna.n_obs),
        "rna_reference_cells_finite_pseudotime": int(finite.sum()),
        "requested_marker_genes": marker_genes,
        "atac_marker_genes_found": atac_genes,
        "atac_marker_genes_missing": atac_missing,
        "shared_transfer_genes": shared_genes,
        "n_shared_transfer_genes": len(shared_genes),
        "computed_shape": [int(atac.n_obs), int(atac.n_vars)],
        "transferred_pseudotime_summary": {
            "min": float(np.nanmin(transferred)),
            "max": float(np.nanmax(transferred)),
            "mean": float(np.nanmean(transferred)),
            "finite_cells": int(np.isfinite(transferred).sum()),
        },
        "mean_transfer_cosine_distance_summary": {
            "min": float(np.nanmin(distances.mean(axis=1))),
            "max": float(np.nanmax(distances.mean(axis=1))),
            "mean": float(np.nanmean(distances.mean(axis=1))),
        },
        "marker_pseudotime_metrics": marker_metrics,
        "marker_direction_checks": direction_checks,
        "all_marker_direction_checks_pass": all(c["pass"] for c in direction_checks),
        "figure_outputs": figure_outputs,
        "software": {"scanpy": sc.__version__},
    }
    (args.out_dir / "figure3b_atac_pseudotime_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
