#!/usr/bin/env python3
"""Render/audit public-resource Trevino Figure 2F/2G/2H/2I panels.

This script deliberately keeps the evidence boundary explicit:

* Figure 2F is a data/QC summary for the public author PCW21 multiome SCEs.
* Figure 2G uses the author-provided multiome RNA UMAP from
  ``Multiome_RNA_SCE.RDS`` and a bounded nearest-neighbor projection of
  multiome ATAC gene activity onto the previously generated public ATAC
  marker manifold. It does not claim exact author uwot/peak-LSI projection
  parity.
* Figure 2H is an honest resource-gap audit: the local/public inputs do not
  expose the exact multiome significant peak-gene linkage table needed for the
  paper Venn diagram.
* Figure 2I recomputes same-cell multiome gene-activity/RNA correlations for
  the public/count-aligned GPC genes and compares them with the public author
  singleome GA/RNA correlations used for Figure 2D/2E.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    for ext in ("png", "svg", "pdf"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight", dpi=300)
        outputs[ext] = str(path)
        outputs[f"{ext}_sha256"] = sha256_file(path)
    plt.close(fig)
    return outputs


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def matrix_to_dense(x: Any) -> np.ndarray:
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def read_peak_gene_link_counts(path: Path) -> dict[str, Any]:
    """Return compact counts for public peak_gene_links.tsv.gz."""
    rows = 0
    sig = 0
    sig_pairs: set[tuple[str, str]] = set()
    sig_peaks: set[str] = set()
    sig_genes: set[str] = set()
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows += 1
            peak, gene, significant = line.rstrip("\n").split("\t")[:3]
            if significant.upper() == "TRUE":
                sig += 1
                sig_pairs.add((peak, gene))
                sig_peaks.add(peak)
                sig_genes.add(gene)
    return {
        "header": header,
        "rows": rows,
        "significant_true_rows": sig,
        "significant_unique_pairs": len(sig_pairs),
        "significant_unique_peaks": len(sig_peaks),
        "significant_unique_genes": len(sig_genes),
    }


def write_gene_list_for_atac(single_atac_h5ad: Path, out_dir: Path) -> Path:
    adata = sc.read_h5ad(single_atac_h5ad, backed="r")
    genes = [str(g) for g in adata.var_names]
    out = out_dir / "single_atac_marker_genes.txt"
    out.write_text("\n".join(genes) + "\n", encoding="utf-8")
    return out


def run_r_export(
    *,
    out_dir: Path,
    rscript: str,
    multiome_rna_sce: Path,
    multiome_atac_sce: Path,
    multiome_atac_gene_activity: Path,
    cluster_names: Path,
    gpc_table: Path,
    singleome_ga_rna_table: Path,
    atac_marker_genes: Path,
) -> dict[str, Path]:
    r_code = r"""
args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[[1]]
multiome_rna_sce <- args[[2]]
multiome_atac_sce <- args[[3]]
multiome_atac_gene_activity <- args[[4]]
cluster_names <- args[[5]]
gpc_table <- args[[6]]
singleome_ga_rna_table <- args[[7]]
atac_marker_genes <- args[[8]]

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

rna <- readRDS(multiome_rna_sce)
atac <- readRDS(multiome_atac_sce)
ga <- readRDS(multiome_atac_gene_activity)

common_cells <- Reduce(intersect, list(colnames(rna), colnames(atac), colnames(ga)))
rna <- rna[, common_cells]
atac <- atac[, common_cells]
ga <- ga[, common_cells]

cd <- as.data.frame(colData(rna))
clusters <- read.table(cluster_names, sep="\t", header=TRUE, stringsAsFactors=FALSE)
rna_map <- clusters[clusters$Assay == "Multiome RNA", c("Cluster.ID", "Cluster.Name")]
names(rna_map) <- c("seurat_clusters", "rna_cluster_name")
cd <- merge(cd, rna_map, by="seurat_clusters", all.x=TRUE, sort=FALSE)
rownames(cd) <- cd$Cell.ID

um <- reducedDim(rna, "UMAP")
um_df <- data.frame(Cell.ID=rownames(um), UMAP_1=um[,1], UMAP_2=um[,2], stringsAsFactors=FALSE)
um_df <- merge(um_df, cd[, c("Cell.ID", "seurat_clusters", "rna_cluster_name", "Sample.ID",
                             "RNA.Counts", "RNA.Features", "CR_ATAC.TSS.enrichment.score",
                             "CR_ATAC.Median.high.quality.fragments.per.cell",
                             "CR_Feature.linkages.detected", "CR_Linked.genes",
                             "CR_Linked.peaks")],
               by="Cell.ID", all.x=TRUE, sort=FALSE)
write.table(um_df, file.path(out_dir, "multiome_rna_umap_metadata.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

qc_cols <- c("Cell.ID", "Sample.ID", "Sample.Age", "seurat_clusters", "rna_cluster_name",
             "RNA.Counts", "RNA.Features", "percentMT", "percentRibo",
             "CR_ATAC.TSS.enrichment.score",
             "CR_ATAC.Median.high.quality.fragments.per.cell",
             "CR_ATAC.Fraction.of.high.quality.fragments.overlapping.peaks",
             "CR_Feature.linkages.detected", "CR_Linked.genes", "CR_Linked.peaks",
             "DF_classification")
qc_cols <- qc_cols[qc_cols %in% colnames(cd)]
write.table(cd[, qc_cols], file.path(out_dir, "multiome_qc_metadata.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

write.table(as.data.frame(sort(table(cd$seurat_clusters), decreasing=TRUE)),
            file.path(out_dir, "multiome_cluster_counts.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Export multiome ATAC gene activity for the marker-gene reference used in the
# public singleome ATAC subset.
marker_genes <- unique(scan(atac_marker_genes, what=character(), quiet=TRUE))
marker_genes <- marker_genes[marker_genes %in% rownames(ga)]
ga_assay <- assay(ga, "gene.activity")
ga_marker <- t(as.matrix(ga_assay[marker_genes, common_cells, drop=FALSE]))
ga_marker_df <- data.frame(Cell.ID=rownames(ga_marker), ga_marker, check.names=FALSE)
ga_marker_con <- gzfile(file.path(out_dir, "multiome_atac_gene_activity_for_projection.tsv.gz"), "wt")
write.table(ga_marker_df, ga_marker_con, sep="\t", quote=FALSE, row.names=FALSE)
close(ga_marker_con)

# Figure 2I: same-cell multiome gene-activity/RNA correlations for public GPCs
# and all genes in the singleome GA/RNA table.
single_tbl <- read.table(singleome_ga_rna_table, sep="\t", header=TRUE, stringsAsFactors=FALSE)
gpc_tbl <- read.table(gpc_table, sep="\t", header=TRUE, stringsAsFactors=FALSE)
target_genes <- unique(c(single_tbl$gene.symbol, gpc_tbl$gene.symbol))
rna_gene_names <- as.character(rowData(rna)$gene_name)
names(rna_gene_names) <- seq_along(rna_gene_names)
rna_first <- !is.na(rna_gene_names) & !duplicated(rna_gene_names)
rna_index_by_symbol <- setNames(which(rna_first), rna_gene_names[rna_first])
ga_index_by_symbol <- setNames(seq_len(nrow(ga)), rownames(ga))
common_genes <- intersect(target_genes, intersect(names(rna_index_by_symbol), names(ga_index_by_symbol)))
rna_log <- assay(rna, "logcounts")

res <- lapply(common_genes, function(g) {
  ri <- rna_index_by_symbol[[g]]
  gi <- ga_index_by_symbol[[g]]
  rv <- as.numeric(rna_log[ri, common_cells])
  gv <- as.numeric(ga_assay[gi, common_cells])
  ok <- is.finite(rv) & is.finite(gv)
  if (sum(ok) < 10 || sd(rv[ok]) == 0 || sd(gv[ok]) == 0) {
    rho <- NA_real_
  } else {
    rho <- suppressWarnings(cor(rv[ok], gv[ok], method="spearman"))
  }
  data.frame(gene.symbol=g, multiome_spearman.cor=rho, multiome_n.cells=sum(ok),
             stringsAsFactors=FALSE)
})
res <- do.call(rbind, res)
merged <- merge(single_tbl, res, by="gene.symbol", all.x=FALSE, all.y=FALSE)
merged$is_count_aligned_gpc <- merged$gene.symbol %in% gpc_tbl$gene.symbol
write.table(merged, file.path(out_dir, "figure2i_singleome_multiome_correlations.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

summary <- data.frame(
  metric=c("common_cells", "rna_features", "atac_peaks", "ga_genes",
           "marker_genes_for_atac_projection", "figure2i_common_genes",
           "figure2i_gpc_genes"),
  value=c(length(common_cells), nrow(rna), nrow(atac), nrow(ga),
          length(marker_genes), nrow(merged), sum(merged$is_count_aligned_gpc))
)
write.table(summary, file.path(out_dir, "multiome_export_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
"""
    script = out_dir / "export_multiome_figure2fghi_inputs.R"
    script.write_text(r_code, encoding="utf-8")
    log = out_dir / "export_multiome_figure2fghi_inputs.log"
    cmd = [
        rscript,
        str(script),
        str(out_dir),
        str(multiome_rna_sce),
        str(multiome_atac_sce),
        str(multiome_atac_gene_activity),
        str(cluster_names),
        str(gpc_table),
        str(singleome_ga_rna_table),
        str(atac_marker_genes),
    ]
    with log.open("w", encoding="utf-8") as fh:
        subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, check=True)
    return {
        "r_export_script": script,
        "r_export_log": log,
        "rna_umap": out_dir / "multiome_rna_umap_metadata.tsv",
        "qc": out_dir / "multiome_qc_metadata.tsv",
        "cluster_counts": out_dir / "multiome_cluster_counts.tsv",
        "atac_ga_projection_matrix": out_dir
        / "multiome_atac_gene_activity_for_projection.tsv.gz",
        "figure2i_correlations": out_dir / "figure2i_singleome_multiome_correlations.tsv",
        "export_summary": out_dir / "multiome_export_summary.tsv",
    }


def render_figure2f(qc: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes = axes.ravel()
    axes[0].hist(qc["RNA.Counts"].dropna(), bins=60, color="#4C72B0", alpha=0.85)
    axes[0].set_title("RNA UMI counts per multiome cell")
    axes[0].set_xlabel("RNA.Counts")
    axes[0].set_ylabel("cells")

    axes[1].hist(qc["RNA.Features"].dropna(), bins=60, color="#55A868", alpha=0.85)
    axes[1].set_title("RNA detected genes per multiome cell")
    axes[1].set_xlabel("RNA.Features")

    axes[2].hist(
        qc["CR_ATAC.TSS.enrichment.score"].dropna(),
        bins=60,
        color="#C44E52",
        alpha=0.85,
    )
    axes[2].set_title("ATAC TSS enrichment")
    axes[2].set_xlabel("TSS enrichment")
    axes[2].set_ylabel("cells")

    counts = qc["rna_cluster_name"].fillna(qc["seurat_clusters"]).value_counts()
    counts.sort_values(ascending=True).plot.barh(ax=axes[3], color="#8172B3")
    axes[3].set_title("PCW21 multiome cluster composition")
    axes[3].set_xlabel("cells")
    fig.suptitle(
        "Figure 2F public multiome QC/design summary\n"
        "Public author SCEs: same-cell RNA + ATAC profiles after joint filtering",
        fontsize=13,
    )
    fig.tight_layout()
    return save_triple(fig, out_dir, "figure2f_multiome_qc_summary")


def _plot_projection_panel(
    ax: plt.Axes,
    reference_umap: np.ndarray,
    projected: pd.DataFrame,
    x_col: str,
    y_col: str,
    label_col: str,
    title: str,
) -> None:
    ax.scatter(
        reference_umap[:, 0],
        reference_umap[:, 1],
        s=2,
        c="#D0D0D0",
        alpha=0.25,
        linewidths=0,
        rasterized=True,
        label="singleome reference",
    )
    labels = projected[label_col].fillna(projected["seurat_clusters"]).astype(str)
    uniq = sorted(labels.unique())
    cmap = plt.get_cmap("tab20")
    for i, label in enumerate(uniq):
        mask = labels == label
        ax.scatter(
            projected.loc[mask, x_col],
            projected.loc[mask, y_col],
            s=8,
            alpha=0.8,
            linewidths=0,
            color=cmap(i % 20),
            label=label,
        )
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")


def render_figure2g(
    *,
    single_rna_h5ad: Path,
    single_atac_h5ad: Path,
    multiome_rna_umap: pd.DataFrame,
    multiome_atac_ga: pd.DataFrame,
    out_dir: Path,
) -> tuple[dict[str, str], dict[str, Any]]:
    single_rna = sc.read_h5ad(single_rna_h5ad, backed="r")
    single_atac = sc.read_h5ad(single_atac_h5ad)
    rna_ref_umap = np.asarray(single_rna.obsm["X_umap"])
    atac_ref_umap = np.asarray(single_atac.obsm["X_umap"])

    # Bounded ATAC projection onto the existing 70-gene public ATAC marker
    # manifold from Figure 3B. This is not the author peak-LSI/uwot transform.
    ref_genes = [g for g in single_atac.var_names if g in multiome_atac_ga.columns]
    ref_x = matrix_to_dense(single_atac[:, ref_genes].X).astype(float)
    multi_x = multiome_atac_ga[ref_genes].to_numpy(dtype=float)
    scaler = StandardScaler(with_mean=True, with_std=True)
    ref_scaled = scaler.fit_transform(ref_x)
    multi_scaled = scaler.transform(multi_x)
    n_comp = max(2, min(20, len(ref_genes) - 1, ref_scaled.shape[0] - 1))
    pca = PCA(n_components=n_comp, random_state=13)
    ref_pca = pca.fit_transform(ref_scaled)
    multi_pca = pca.transform(multi_scaled)
    nn = NearestNeighbors(n_neighbors=5, metric="euclidean")
    nn.fit(ref_pca)
    dist, idx = nn.kneighbors(multi_pca)
    projected_umap = atac_ref_umap[idx].mean(axis=1)

    atac_projected = pd.DataFrame(
        {
            "Cell.ID": multiome_atac_ga["Cell.ID"].to_numpy(),
            "ATAC_UMAP_1": projected_umap[:, 0],
            "ATAC_UMAP_2": projected_umap[:, 1],
            "mean_nn_distance": dist.mean(axis=1),
        }
    )
    atac_projected = atac_projected.merge(
        multiome_rna_umap[["Cell.ID", "seurat_clusters", "rna_cluster_name"]],
        on="Cell.ID",
        how="left",
    )
    atac_projected.to_csv(
        out_dir / "figure2g_multiome_atac_marker_nn_projection.tsv.gz",
        sep="\t",
        index=False,
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    _plot_projection_panel(
        axes[0],
        rna_ref_umap,
        multiome_rna_umap,
        "UMAP_1",
        "UMAP_2",
        "rna_cluster_name",
        "Fig.2G RNA: author public multiome UMAP over current RNA reference",
    )
    _plot_projection_panel(
        axes[1],
        atac_ref_umap,
        atac_projected,
        "ATAC_UMAP_1",
        "ATAC_UMAP_2",
        "rna_cluster_name",
        "Fig.2G ATAC: bounded marker-gene-activity NN projection",
    )
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=7,
        frameon=False,
    )
    fig.suptitle(
        "Figure 2G public multiome projection evidence\n"
        "Gray: current public singleome manifolds; color: PCW21 multiome cluster labels",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 0.86, 1))
    outputs = save_triple(fig, out_dir, "figure2g_multiome_projection_public")
    metrics = {
        "rna_multiome_cells": int(multiome_rna_umap.shape[0]),
        "singleome_rna_reference_cells": int(rna_ref_umap.shape[0]),
        "singleome_atac_reference_cells": int(atac_ref_umap.shape[0]),
        "atac_projection_common_marker_genes": int(len(ref_genes)),
        "atac_projection_mean_nn_distance_median": float(
            np.nanmedian(atac_projected["mean_nn_distance"])
        ),
    }
    return outputs, metrics


def render_figure2h_gap(
    *,
    peak_gene_link_counts: dict[str, Any],
    rds_significant_rows: int,
    out_dir: Path,
) -> dict[str, str]:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axis("off")
    text = textwrap.dedent(
        f"""
        Figure 2H linkage-overlap audit

        Paper target:
          • singleome significant links: 64,878
          • overlap observed in same-cell multiome: 40,181
          • additional multiome links: 23,849

        Local/public resources found:
          • PeakGeneLinks_Significant.RDS rows: {rds_significant_rows:,}
          • peak_gene_links.tsv.gz significant=TRUE rows: {peak_gene_link_counts['significant_true_rows']:,}
          • no separate exact multiome significant-link table with peak-gene pairs

        Honest status:
          METHOD/RESOURCE GAP — the Venn diagram is not reproduced here.
          We can audit available link resources, but cannot assign the exact
          singleome-vs-multiome overlap without the multiome significant link list
          or a full rerun of the paper pseudobulk linkage workflow.
        """
    ).strip()
    ax.text(0.02, 0.98, text, va="top", ha="left", fontsize=11, family="DejaVu Sans Mono")
    return save_triple(fig, out_dir, "figure2h_linkage_overlap_resource_gap")


def render_figure2i(corr: pd.DataFrame, out_dir: Path) -> tuple[dict[str, str], dict[str, Any]]:
    gpc = corr[corr["is_count_aligned_gpc"].astype(bool)].copy()
    gpc = gpc.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["spearman.cor", "multiome_spearman.cor"]
    )
    pear = pearsonr(gpc["spearman.cor"], gpc["multiome_spearman.cor"])
    spear = spearmanr(gpc["spearman.cor"], gpc["multiome_spearman.cor"])
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(
        gpc["spearman.cor"],
        gpc["multiome_spearman.cor"],
        s=np.clip(gpc["linked_cre_count"].fillna(1), 5, 60),
        alpha=0.75,
        color="#4C72B0",
        edgecolor="white",
        linewidth=0.3,
    )
    for _, row in gpc.sort_values("linked_cre_count", ascending=False).head(12).iterrows():
        ax.text(
            row["spearman.cor"],
            row["multiome_spearman.cor"],
            str(row["gene.symbol"]),
            fontsize=7,
            alpha=0.85,
        )
    lim_min = float(np.nanmin([gpc["spearman.cor"].min(), gpc["multiome_spearman.cor"].min()]))
    lim_max = float(np.nanmax([gpc["spearman.cor"].max(), gpc["multiome_spearman.cor"].max()]))
    ax.plot([lim_min, lim_max], [lim_min, lim_max], "--", color="black", linewidth=1, alpha=0.5)
    ax.set_xlabel("Singleome GA/RNA Spearman correlation (public author table)")
    ax.set_ylabel("Multiome same-cell GA/RNA Spearman correlation (recomputed)")
    ax.set_title(
        "Figure 2I predictive chromatin correspondence\n"
        f"GPCs n={len(gpc)}; Pearson r={pear.statistic:.2f}; Spearman ρ={spear.statistic:.2f}"
    )
    fig.tight_layout()
    outputs = save_triple(fig, out_dir, "figure2i_predictive_chromatin_correlation")
    metrics = {
        "gpc_genes_with_singleome_and_multiome_correlations": int(len(gpc)),
        "all_genes_with_singleome_and_multiome_correlations": int(
            corr.dropna(subset=["spearman.cor", "multiome_spearman.cor"]).shape[0]
        ),
        "pearson_r_gpc": float(pear.statistic),
        "pearson_p_gpc": float(pear.pvalue),
        "spearman_rho_gpc": float(spear.statistic),
        "spearman_p_gpc": float(spear.pvalue),
    }
    return outputs, metrics


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--single-rna-h5ad", type=Path, required=True)
    ap.add_argument("--single-atac-h5ad", type=Path, required=True)
    ap.add_argument("--multiome-rna-sce", type=Path, required=True)
    ap.add_argument("--multiome-atac-sce", type=Path, required=True)
    ap.add_argument("--multiome-atac-gene-activity", type=Path, required=True)
    ap.add_argument("--cluster-names", type=Path, required=True)
    ap.add_argument("--gpc-table", type=Path, required=True)
    ap.add_argument("--singleome-ga-rna-table", type=Path, required=True)
    ap.add_argument("--peak-gene-links-tsv", type=Path, required=True)
    ap.add_argument("--rds-significant-link-count", type=int, default=64030)
    ap.add_argument("--rscript", default="/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    configure_matplotlib()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    marker_genes = write_gene_list_for_atac(args.single_atac_h5ad, args.out_dir)
    exported = run_r_export(
        out_dir=args.out_dir,
        rscript=args.rscript,
        multiome_rna_sce=args.multiome_rna_sce,
        multiome_atac_sce=args.multiome_atac_sce,
        multiome_atac_gene_activity=args.multiome_atac_gene_activity,
        cluster_names=args.cluster_names,
        gpc_table=args.gpc_table,
        singleome_ga_rna_table=args.singleome_ga_rna_table,
        atac_marker_genes=marker_genes,
    )

    qc = pd.read_csv(exported["qc"], sep="\t")
    multiome_rna_umap = pd.read_csv(exported["rna_umap"], sep="\t")
    multiome_atac_ga = pd.read_csv(exported["atac_ga_projection_matrix"], sep="\t")
    corr = pd.read_csv(exported["figure2i_correlations"], sep="\t")
    export_summary = pd.read_csv(exported["export_summary"], sep="\t")
    export_metrics = dict(zip(export_summary["metric"], export_summary["value"]))

    figure_outputs: dict[str, dict[str, str]] = {}
    figure_outputs["2F"] = render_figure2f(qc, args.out_dir)
    fig2g_outputs, fig2g_metrics = render_figure2g(
        single_rna_h5ad=args.single_rna_h5ad,
        single_atac_h5ad=args.single_atac_h5ad,
        multiome_rna_umap=multiome_rna_umap,
        multiome_atac_ga=multiome_atac_ga,
        out_dir=args.out_dir,
    )
    figure_outputs["2G"] = fig2g_outputs
    peak_gene_counts = read_peak_gene_link_counts(args.peak_gene_links_tsv)
    figure_outputs["2H"] = render_figure2h_gap(
        peak_gene_link_counts=peak_gene_counts,
        rds_significant_rows=args.rds_significant_link_count,
        out_dir=args.out_dir,
    )
    fig2i_outputs, fig2i_metrics = render_figure2i(corr, args.out_dir)
    figure_outputs["2I"] = fig2i_outputs

    outputs = {
        **{k: str(v) for k, v in exported.items()},
        "figure2g_atac_projection": str(
            args.out_dir / "figure2g_multiome_atac_marker_nn_projection.tsv.gz"
        ),
    }
    hashes = {
        Path(p).name: sha256_file(Path(p))
        for p in outputs.values()
        if Path(p).exists() and Path(p).is_file()
    }
    checks = {
        "multiome_cells_eq_8981": int(export_metrics.get("common_cells", 0)) == 8981,
        "figure2f_qc_rows_eq_8981": int(qc.shape[0]) == 8981,
        "figure2g_rna_cells_eq_8981": int(fig2g_metrics["rna_multiome_cells"]) == 8981,
        "figure2g_atac_projection_common_genes_ge_20": int(
            fig2g_metrics["atac_projection_common_marker_genes"]
        )
        >= 20,
        "figure2h_gap_records_public_link_counts": peak_gene_counts["rows"] > 0,
        "figure2i_gpc_genes_ge_150": int(
            fig2i_metrics["gpc_genes_with_singleome_and_multiome_correlations"]
        )
        >= 150,
        "all_output_hashes_recorded": bool(hashes),
    }
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "panels": ["2F", "2G", "2H", "2I"],
        "status": {
            "2F": "REPRODUCED_PUBLIC_AUTHOR_MULTIOME_QC_SUMMARY",
            "2G": "GENERATED_PUBLIC_AUTHOR_RNA_AND_BOUNDED_ATAC_MULTIOME_PROJECTION",
            "2H": "METHOD_RESOURCE_GAP_MULTIOME_LINKAGE_VENN_NOT_REPRODUCED",
            "2I": "GENERATED_PUBLIC_MULTIOME_GPC_CORRELATION_CORRESPONDENCE",
        },
        "truth_boundary": {
            "2F": "Public author PCW21 multiome SCE QC/cluster metadata; schematic elements are represented as data/QC summary.",
            "2G": "RNA uses the public author Multiome_RNA_SCE UMAP; ATAC uses a bounded 70-gene marker NN projection onto the existing public ATAC marker manifold, not the exact author peak-LSI/uwot model.",
            "2H": "Exact singleome-vs-multiome linkage Venn cannot be reproduced because no separate public multiome significant peak-gene pair table is available locally.",
            "2I": "Same-cell multiome GA/RNA correlations are recomputed from public SCEs and compared to public singleome GA/RNA correlations for count-aligned GPCs; not a full rerun of pseudobulk linkage calling.",
        },
        "inputs": {
            "single_rna_h5ad": str(args.single_rna_h5ad),
            "single_atac_h5ad": str(args.single_atac_h5ad),
            "multiome_rna_sce": str(args.multiome_rna_sce),
            "multiome_atac_sce": str(args.multiome_atac_sce),
            "multiome_atac_gene_activity": str(args.multiome_atac_gene_activity),
            "cluster_names": str(args.cluster_names),
            "gpc_table": str(args.gpc_table),
            "singleome_ga_rna_table": str(args.singleome_ga_rna_table),
            "peak_gene_links_tsv": str(args.peak_gene_links_tsv),
        },
        "outputs": outputs,
        "figure_outputs": figure_outputs,
        "metrics": {
            "export": export_metrics,
            "figure2g": fig2g_metrics,
            "figure2h_peak_gene_links_tsv": peak_gene_counts,
            "figure2i": fig2i_metrics,
        },
        "checks": checks,
        "output_hashes": hashes,
    }
    summary_path = args.out_dir / "figure2fghi_public_multiome_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path), "checks": checks, "status": summary["status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
