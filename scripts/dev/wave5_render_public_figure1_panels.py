#!/usr/bin/env python3
"""Render Trevino Figure 1 public-matrix panels 1D/1F/1G/1H.

Panels are generated from the current public RNA pipeline output plus public
GEO ATAC gene-activity resources. This is an honest public-matrix
reconstruction, not FASTQ/fragments-level CellRanger/ArchR parity.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402


MARKER_GROUPS: dict[str, list[str]] = {
    "Radial glia": ["SOX2", "PAX6", "HES1", "VIM", "HOPX"],
    "IPC/GluN": ["EOMES", "NEUROG2", "NEUROD1", "NEUROD2", "NEUROD6", "SATB2", "SLC17A7"],
    "IN": ["DLX1", "DLX2", "GAD1", "GAD2"],
    "OPC/Oligo": ["OLIG1", "OLIG2", "SOX10", "PDGFRA"],
    "Astro/MG/Vascular": ["GFAP", "AQP4", "PTPRC", "CX3CR1", "PECAM1", "CLDN5", "PDGFRB", "RGS5"],
}
MARKERS = [g for genes in MARKER_GROUPS.values() for g in genes]
AGE_ORDER = ["pcw16", "pcw20", "pcw21", "pcw24"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    paths = {ext: out_dir / f"{stem}.{ext}" for ext in ("png", "svg", "pdf")}
    fig.savefig(paths["png"], dpi=300, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)
    return {f"{k}_path": str(v) for k, v in paths.items()} | {
        f"{k}_sha256": sha256_file(v) for k, v in paths.items()
    }


def ordered_categories(values: pd.Series, preferred: list[str] | None = None) -> list[str]:
    seen = values.astype(str).value_counts().index.tolist()
    if preferred:
        return [x for x in preferred if x in seen] + [x for x in seen if x not in preferred]
    return seen


def scatter_by_category(
    ax: plt.Axes,
    coords: np.ndarray,
    labels: pd.Series,
    title: str,
    preferred: list[str] | None = None,
    size: float = 2.0,
) -> None:
    cats = ordered_categories(labels, preferred)
    cmap = plt.get_cmap("tab20")
    for i, cat in enumerate(cats):
        mask = labels.astype(str).to_numpy() == cat
        ax.scatter(coords[mask, 0], coords[mask, 1], s=size, alpha=0.75, linewidths=0, color=cmap(i % 20), label=cat)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7, markerscale=3)


def matrix_to_dense(x) -> np.ndarray:
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def dotplot(
    mean_df: pd.DataFrame,
    pct_df: pd.DataFrame,
    title: str,
    out_dir: Path,
    stem: str,
    cmap: str = "viridis",
) -> dict[str, str]:
    clusters = list(mean_df.index)
    genes = list(mean_df.columns)
    fig_w = max(8, 0.28 * len(genes))
    fig_h = max(4, 0.35 * len(clusters))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    x, y, sizes, colors = [], [], [], []
    for iy, cl in enumerate(clusters):
        for ix, gene in enumerate(genes):
            x.append(ix)
            y.append(iy)
            sizes.append(15 + 180 * float(pct_df.loc[cl, gene]))
            colors.append(float(mean_df.loc[cl, gene]))
    sc = ax.scatter(x, y, s=sizes, c=colors, cmap=cmap, edgecolor="0.6", linewidth=0.2)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters, fontsize=7)
    ax.set_title(title)
    ax.set_xlabel("marker")
    ax.set_ylabel("cluster / annotation")
    fig.colorbar(sc, ax=ax, label="mean log activity / expression")
    return save_triple(fig, out_dir, stem)


def aggregate_by_group(matrix: np.ndarray, groups: pd.Series, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    cats = groups.astype(str).value_counts().index.tolist()
    mean_rows, pct_rows = [], []
    for cat in cats:
        mask = groups.astype(str).to_numpy() == cat
        sub = matrix[mask, :]
        mean_rows.append(np.nanmean(sub, axis=0))
        pct_rows.append(np.nanmean(sub > 0, axis=0))
    return pd.DataFrame(mean_rows, index=cats, columns=genes), pd.DataFrame(pct_rows, index=cats, columns=genes)


def read_atac_gene_activity_rows(path: Path, genes: list[str]) -> tuple[list[str], dict[str, np.ndarray]]:
    wanted = set(genes)
    rows: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().strip().split()
        for line in fh:
            if not line.strip():
                continue
            gene, rest = line.split(maxsplit=1)
            if gene not in wanted:
                continue
            rows[gene] = np.fromstring(rest, sep=" ", dtype=np.float32)
            if len(rows) == len(wanted):
                break
    return header, rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rna-h5ad", required=True, type=Path)
    p.add_argument("--atac-umap-h5ad", required=True, type=Path)
    p.add_argument("--atac-gene-activities", required=True, type=Path)
    p.add_argument("--atac-metadata", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rna = ad.read_h5ad(args.rna_h5ad, backed="r")
    rna_obs = rna.obs.copy()
    rna_umap = np.asarray(rna.obsm["X_umap"])
    rna_present = [g for g in MARKERS if g in set(map(str, rna.var_names))]
    rna_var_idx = rna.var_names.get_indexer(rna_present)
    rna_mat = matrix_to_dense(rna[:, rna_var_idx].X)
    rna.file.close()
    rna_mat = np.log1p(np.maximum(rna_mat.astype(np.float32, copy=False), 0))

    atac = ad.read_h5ad(args.atac_umap_h5ad, backed="r")
    atac_obs = atac.obs.copy()
    atac_umap = np.asarray(atac.obsm["X_umap"])
    atac.file.close()

    atac_meta = pd.read_csv(args.atac_metadata, sep="\t")
    atac_header, atac_rows = read_atac_gene_activity_rows(args.atac_gene_activities, MARKERS)
    atac_present = [g for g in MARKERS if g in atac_rows]
    header_index = pd.Index(atac_header)
    meta_idx = header_index.get_indexer(atac_meta["Cell.ID"].astype(str))
    keep = meta_idx >= 0
    atac_meta = atac_meta.loc[keep].copy()
    meta_idx = meta_idx[keep]
    atac_matrix = np.vstack([atac_rows[g][meta_idx] for g in atac_present]).T
    atac_matrix = np.log1p(np.maximum(atac_matrix.astype(np.float32, copy=False), 0))

    outputs: dict[str, dict[str, str]] = {}
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    scatter_by_category(axes[0], rna_umap, rna_obs["Age"], "1D RNA UMAP by gestational age", preferred=AGE_ORDER, size=1.0)
    scatter_by_category(axes[1], atac_umap, atac_obs["Age"], "1D ATAC UMAP by gestational age (12k subset)", preferred=AGE_ORDER, size=2.0)
    outputs["figure1d_rna_atac_umap_by_age"] = save_triple(fig, args.out_dir, "figure1d_rna_atac_umap_by_age")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    scatter_by_category(axes[0], rna_umap, rna_obs["cell_type"], "1F RNA UMAP by pipeline annotation", size=1.0)
    scatter_by_category(axes[1], atac_umap, atac_obs["Iterative.LSI.Clusters"], "1F ATAC UMAP by public LSI cluster (12k subset)", size=2.0)
    outputs["figure1f_rna_celltype_atac_cluster_umap"] = save_triple(fig, args.out_dir, "figure1f_rna_celltype_atac_cluster_umap")

    rna_mean, rna_pct = aggregate_by_group(rna_mat, rna_obs["cell_type"], rna_present)
    outputs["figure1g_rna_marker_dotplot"] = dotplot(
        rna_mean, rna_pct, "1G RNA marker expression by pipeline annotation", args.out_dir, "figure1g_rna_marker_dotplot"
    )
    rna_mean.to_csv(args.out_dir / "figure1g_rna_marker_mean.tsv", sep="\t")
    rna_pct.to_csv(args.out_dir / "figure1g_rna_marker_fraction.tsv", sep="\t")

    atac_mean, atac_pct = aggregate_by_group(atac_matrix, atac_meta["Iterative.LSI.Clusters"], atac_present)
    outputs["figure1h_atac_gene_activity_dotplot"] = dotplot(
        atac_mean,
        atac_pct,
        "1H ATAC gene activity markers by public LSI cluster",
        args.out_dir,
        "figure1h_atac_gene_activity_dotplot",
        cmap="magma",
    )
    atac_mean.to_csv(args.out_dir / "figure1h_atac_marker_mean.tsv", sep="\t")
    atac_pct.to_csv(args.out_dir / "figure1h_atac_marker_fraction.tsv", sep="\t")

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "GENERATED_PUBLIC_MATRIX_FIGURE1_1D_1F_1G_1H",
        "truth_boundary": (
            "Public processed RNA and ATAC gene-activity matrices were used. ATAC UMAP panels use the existing "
            "12k Figure 3B subset embedding; ATAC marker dotplot uses all public 31,304 ATAC metadata-aligned cells. "
            "This is not FASTQ/fragments-level or exact author embedding parity."
        ),
        "inputs": {
            "rna_h5ad": str(args.rna_h5ad),
            "atac_umap_h5ad": str(args.atac_umap_h5ad),
            "atac_gene_activities": str(args.atac_gene_activities),
            "atac_metadata": str(args.atac_metadata),
        },
        "out_dir": str(args.out_dir),
        "metrics": {
            "rna_cells": int(rna_obs.shape[0]),
            "rna_genes": int(len(rna_present)),
            "rna_age_groups": rna_obs["Age"].astype(str).value_counts().to_dict(),
            "rna_cell_type_groups": rna_obs["cell_type"].astype(str).value_counts().to_dict(),
            "atac_umap_cells": int(atac_obs.shape[0]),
            "atac_dotplot_cells": int(atac_meta.shape[0]),
            "atac_marker_genes": int(len(atac_present)),
            "atac_age_groups_subset": atac_obs["Age"].astype(str).value_counts().to_dict(),
            "atac_lsi_clusters_full": atac_meta["Iterative.LSI.Clusters"].astype(str).value_counts().to_dict(),
            "rna_missing_markers": [g for g in MARKERS if g not in rna_present],
            "atac_missing_markers": [g for g in MARKERS if g not in atac_present],
        },
        "outputs": outputs,
        "checks": {
            "rna_cells_ge_50000": int(rna_obs.shape[0]) >= 50000,
            "rna_umap_present_finite": bool(np.isfinite(rna_umap).all() and rna_umap.shape[1] == 2),
            "rna_age_groups_ge_4": int(rna_obs["Age"].nunique()) >= 4,
            "rna_markers_all_present": len(rna_present) == len(MARKERS),
            "atac_umap_cells_ge_10000": int(atac_obs.shape[0]) >= 10000,
            "atac_umap_present_finite": bool(np.isfinite(atac_umap).all() and atac_umap.shape[1] == 2),
            "atac_dotplot_cells_ge_30000": int(atac_meta.shape[0]) >= 30000,
            "atac_markers_all_present": len(atac_present) == len(MARKERS),
            "atac_clusters_ge_10": int(atac_meta["Iterative.LSI.Clusters"].nunique()) >= 10,
            "all_output_hashes_recorded": all("png_sha256" in v and "svg_sha256" in v and "pdf_sha256" in v for v in outputs.values()),
        },
    }
    summary_path = args.out_dir / "figure1_public_matrix_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
