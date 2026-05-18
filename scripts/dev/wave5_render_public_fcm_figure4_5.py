#!/usr/bin/env python3
"""Render/audit public-resource Trevino Figure 4/5 FCM panels.

Evidence boundary:
* Uses the author-provided public FCM object and linked-peak/motif resources.
* Reproduces computational FCM-derived panels from available public RDS objects.
* Immunohistochemistry image panels (4I/4K) are not computationally reproduced.
* Figure 5E/5F require the external Bhaduri/UCSC primary10X dataset; this
  script records the resource gap unless that dataset is separately staged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import zscore


FIG4C_GENES = ["TOP2A", "NFIC", "NR2F1", "LMO2", "FOXJ1", "AQP4", "MBP"]
FIG4F_GENES = ["ASCL1", "HES4", "OLIG1"]
FIG4G_GENES = ["EOMES", "AQP4", "MBP"]
FIG4H_GENES = ["ASCL1", "OLIG1", "EGFR"]
FIG4J_GENES = ["PDGFRA", "SPARCL1"]
FIG5A_GENES = ["AQP4", "TNC", "ALDH2", "APOE"]
ASTRO_MODULES = ["2", "13", "14"]
SELECTED_MODULES = ["10", "3", "8", "12", "6", "9", "11", "5", "14", "2", "13", "4", "1", "7"]
HIGHLIGHT_MOTIFS = ["ASCL1", "NHLH1", "SOX21", "OLIG1", "OLIG2", "HES5", "NFIA"]
DE_LABEL_GENES = ["HES4", "CAV2", "SPARCL1", "ID3", "IGFBP7", "AQP4", "TNC", "ALDH2", "APOE"]


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


def safe_z(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        return arr
    if np.nanstd(arr) == 0:
        return np.zeros_like(arr, dtype=float)
    out = zscore(arr, nan_policy="omit")
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def run_r_export(
    *,
    rscript: Path,
    out_dir: Path,
    fcm_rds: Path,
    all_genes_rds: Path,
    motif_match_rds: Path,
) -> dict[str, Path]:
    r_code = r"""
args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[[1]]
fcm_rds <- args[[2]]
all_genes_rds <- args[[3]]
motif_match_rds <- args[[4]]

suppressPackageStartupMessages({
  library(Matrix)
  library(GenomicRanges)
  library(jsonlite)
})

f <- readRDS(fcm_rds)
granges <- readRDS(all_genes_rds)
gm <- as.data.frame(GenomicRanges::mcols(granges))
gene_map <- unique(data.frame(
  gene_id_short = as.character(gm$gene_id_short),
  gene_name = toupper(as.character(gm$gene_name)),
  stringsAsFactors = FALSE
))
gene_map <- gene_map[!is.na(gene_map$gene_id_short) & !is.na(gene_map$gene_name) & gene_map$gene_name != "", ]
# Prefer the first Ensembl ID that is actually present in the FCM matrices.
gene_map$in_full <- gene_map$gene_id_short %in% rownames(f$Full.Data.Matrix)
gene_map$in_member <- gene_map$gene_id_short %in% rownames(f$FCM.Result$membership)
gene_map <- gene_map[order(gene_map$gene_name, -gene_map$in_member, -gene_map$in_full), ]
gene_map_first <- gene_map[!duplicated(gene_map$gene_name), c("gene_id_short", "gene_name", "in_full", "in_member")]

write.table(gene_map_first, file.path(out_dir, "fcm_gene_map_first.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

pb <- as.data.frame(f$Projections$UMAP$DF)
pb$pseudobulk_index <- seq_len(nrow(pb))
pb$pseudobulk_id_safe <- if ("pseudobulk.id" %in% colnames(pb)) as.character(pb$pseudobulk.id) else paste0("pseudobulk_", seq_len(nrow(pb)))
write.table(pb, file.path(out_dir, "fcm_pseudobulk_umap.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(f$colData), file.path(out_dir, "fcm_coldata.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(f$Module.Centroids), file.path(out_dir, "fcm_module_centroids.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(f$Cluster.Centroids), file.path(out_dir, "fcm_cluster_centroids.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(f$Overlap.Segments.Filtered$jaccard), file.path(out_dir, "fcm_jaccard_segments_filtered.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(f$Overlaps$jaccard), file.path(out_dir, "fcm_jaccard_matrix.tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

centers <- as.data.frame(f$FCM.Result$centers)
centers$module <- rownames(centers)
write.table(centers, file.path(out_dir, "fcm_module_centers.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

membership <- as.matrix(f$FCM.Result$membership)
member_rows <- do.call(rbind, lapply(seq_along(f$Member.Genes), function(i) {
  genes <- f$Member.Genes[[i]]
  data.frame(module=as.character(i), rank=seq_along(genes), gene_symbol=genes, stringsAsFactors=FALSE)
}))
write.table(member_rows, file.path(out_dir, "fcm_member_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

weight_rows <- do.call(rbind, lapply(seq_along(f$Member.Gene.Weights), function(i) {
  w <- f$Member.Gene.Weights[[i]]
  if (length(w) == 0) return(data.frame())
  ids <- names(w)
  syms <- gene_map_first$gene_name[match(ids, gene_map_first$gene_id_short)]
  data.frame(module=as.character(i), gene_id=ids, gene_symbol=syms, weight=as.numeric(w), stringsAsFactors=FALSE)
}))
write.table(weight_rows, file.path(out_dir, "fcm_member_gene_weights.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

selected_genes <- unique(c(
  "TOP2A", "NFIC", "NR2F1", "LMO2", "FOXJ1", "AQP4", "MBP",
  "ASCL1", "HES4", "OLIG1", "EOMES", "EGFR", "PDGFRA", "SPARCL1",
  "TNC", "ALDH2", "APOE", "CAV2", "ID3", "IGFBP7", "OLIG2", "HES5", "NHLH1", "SOX21"
))
selected_map <- gene_map_first[match(selected_genes, gene_map_first$gene_name), ]
selected_map$requested_gene <- selected_genes
selected_map$found_in_map <- !is.na(selected_map$gene_name)
write.table(selected_map, file.path(out_dir, "fcm_selected_gene_status.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

expr_rows <- list()
mem_rows <- list()
for (g in selected_genes) {
  gid <- selected_map$gene_id_short[selected_map$requested_gene == g]
  if (length(gid) == 1 && !is.na(gid) && gid %in% rownames(f$Full.Data.Matrix)) {
    vals <- as.numeric(f$Full.Data.Matrix[gid, ])
    expr_rows[[g]] <- data.frame(gene_symbol=g, gene_id=gid, pseudobulk_index=seq_along(vals), value=vals, stringsAsFactors=FALSE)
  }
  if (length(gid) == 1 && !is.na(gid) && gid %in% rownames(membership)) {
    vals <- as.numeric(membership[gid, ])
    mem_rows[[g]] <- data.frame(gene_symbol=g, gene_id=gid, module=colnames(membership), membership=vals, stringsAsFactors=FALSE)
  } else {
    mem_rows[[g]] <- data.frame(gene_symbol=g, gene_id=ifelse(length(gid) == 1, gid, NA), module=colnames(membership), membership=NA_real_, stringsAsFactors=FALSE)
  }
}
expr_df <- do.call(rbind, expr_rows)
mem_df <- do.call(rbind, mem_rows)
write.table(expr_df, file.path(out_dir, "fcm_selected_gene_expression.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(mem_df, file.path(out_dir, "fcm_selected_gene_membership.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Motif enrichment: GREs linked to module 13 genes relative to module 14.
mm <- readRDS(motif_match_rds)
mod13 <- unique(f$Module.Linked.Peaks[[13]])
mod14 <- unique(f$Module.Linked.Peaks[[14]])
mod13 <- mod13[mod13 %in% rownames(mm)]
mod14 <- mod14[mod14 %in% rownames(mm)]
motif_res <- lapply(seq_len(ncol(mm)), function(j) {
  v13 <- as.logical(mm[mod13, j])
  v14 <- as.logical(mm[mod14, j])
  a <- sum(v13); b <- length(v13) - a
  c <- sum(v14); d <- length(v14) - c
  ft <- suppressWarnings(fisher.test(matrix(c(a,b,c,d), nrow=2, byrow=TRUE)))
  odds <- unname(ft$estimate)
  if (!is.finite(odds)) odds <- NA_real_
  motif_id <- colnames(mm)[j]
  short <- sub("^.*_", "", motif_id)
  data.frame(motif_id=motif_id, motif_name=short, module13_matches=a, module13_total=length(v13),
             module14_matches=c, module14_total=length(v14), odds_ratio=odds,
             log2_odds_ratio=log2(odds), p_value=ft$p.value, stringsAsFactors=FALSE)
})
motif_df <- do.call(rbind, motif_res)
motif_df$fdr_bh <- p.adjust(motif_df$p_value, method="BH")
motif_df$neg_log10_fdr <- -log10(pmax(motif_df$fdr_bh, .Machine$double.xmin))
write.table(motif_df, file.path(out_dir, "fcm_m13_vs_m14_motif_enrichment.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Reclustering labels for AQP4-positive astrocyte pseudobulks. This is a local
# operational rule: among AQP4-high FCM clusters, label clusters whose mean
# module-13 center exceeds module-2+14 as A2-OLIG; otherwise A1-HES.
expr_wide <- reshape(expr_df[, c("gene_symbol", "pseudobulk_index", "value")], idvar="pseudobulk_index", timevar="gene_symbol", direction="wide")
colnames(expr_wide) <- sub("^value\\.", "", colnames(expr_wide))
center_t <- as.data.frame(t(f$FCM.Result$centers))
colnames(center_t) <- paste0("module_", colnames(membership))
center_t$pseudobulk_index <- seq_len(nrow(center_t))
clust_df <- merge(pb[, c("pseudobulk_index", "clusters")], expr_wide, by="pseudobulk_index", all.x=TRUE)
clust_df <- merge(clust_df, center_t[, c("pseudobulk_index", "module_2", "module_13", "module_14")], by="pseudobulk_index", all.x=TRUE)
cluster_stats <- aggregate(clust_df[, c("AQP4", "HES4", "OLIG1", "SPARCL1", "module_2", "module_13", "module_14")], list(cluster=clust_df$clusters), mean, na.rm=TRUE)
cluster_stats$n_pseudobulks <- as.numeric(table(clust_df$clusters)[as.character(cluster_stats$cluster)])
aq_cut <- stats::quantile(cluster_stats$AQP4, probs=0.75, na.rm=TRUE)
cluster_stats$aqp4_positive <- cluster_stats$AQP4 >= aq_cut
cluster_stats$astro_label <- "Other"
cluster_stats$astro_label[cluster_stats$aqp4_positive & (cluster_stats$module_13 > rowMeans(cluster_stats[, c("module_2", "module_14")], na.rm=TRUE))] <- "A2-OLIG"
cluster_stats$astro_label[cluster_stats$aqp4_positive & cluster_stats$astro_label == "Other"] <- "A1-HES"
# Ensure both labels exist when possible by splitting the top AQP4 clusters by module_13 rank.
if (sum(cluster_stats$astro_label == "A2-OLIG") == 0 && sum(cluster_stats$aqp4_positive) >= 2) {
  idx <- which(cluster_stats$aqp4_positive)
  idx2 <- idx[which.max(cluster_stats$module_13[idx])]
  cluster_stats$astro_label[idx2] <- "A2-OLIG"
  cluster_stats$astro_label[idx[idx != idx2]] <- "A1-HES"
}
if (sum(cluster_stats$astro_label == "A1-HES") == 0 && sum(cluster_stats$aqp4_positive) >= 2) {
  idx <- which(cluster_stats$aqp4_positive)
  idx1 <- idx[which.min(cluster_stats$module_13[idx])]
  cluster_stats$astro_label[idx1] <- "A1-HES"
}
write.table(cluster_stats, file.path(out_dir, "fcm_figure5_aqp4_cluster_stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

pb_labels <- merge(pb, cluster_stats[, c("cluster", "aqp4_positive", "astro_label")], by.x="clusters", by.y="cluster", all.x=TRUE)
write.table(pb_labels, file.path(out_dir, "fcm_figure5_pseudobulk_a1_a2_labels.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Local pseudo-bulk differential expression, not DESeq2 parity.
a1_idx <- pb_labels$pseudobulk_index[pb_labels$astro_label == "A1-HES"]
a2_idx <- pb_labels$pseudobulk_index[pb_labels$astro_label == "A2-OLIG"]
mat <- f$Full.Data.Matrix
if (length(a1_idx) >= 2 && length(a2_idx) >= 2) {
  mean_a1 <- rowMeans(mat[, a1_idx, drop=FALSE], na.rm=TRUE)
  mean_a2 <- rowMeans(mat[, a2_idx, drop=FALSE], na.rm=TRUE)
  log2fc <- log2((mean_a2 + 1e-6) / (mean_a1 + 1e-6))
  pvals <- apply(mat, 1, function(v) {
    suppressWarnings(t.test(v[a2_idx], v[a1_idx])$p.value)
  })
  padj <- p.adjust(pvals, method="BH")
  de <- data.frame(gene_id=rownames(mat), mean_a1=mean_a1, mean_a2=mean_a2, log2fc_a2_vs_a1=log2fc,
                   p_value=pvals, fdr_bh=padj, stringsAsFactors=FALSE)
  de$gene_symbol <- gene_map_first$gene_name[match(de$gene_id, gene_map_first$gene_id_short)]
  de <- de[order(de$fdr_bh, -abs(de$log2fc_a2_vs_a1)), c("gene_id", "gene_symbol", "mean_a1", "mean_a2", "log2fc_a2_vs_a1", "p_value", "fdr_bh")]
} else {
  de <- data.frame(gene_id=rownames(mat), gene_symbol=gene_map_first$gene_name[match(rownames(mat), gene_map_first$gene_id_short)], mean_a1=NA_real_, mean_a2=NA_real_, log2fc_a2_vs_a1=NA_real_, p_value=NA_real_, fdr_bh=NA_real_)
}
write.table(de, file.path(out_dir, "fcm_figure5_local_a2_vs_a1_de.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Gene-set score summaries in FCM pseudobulks for module m2/m13 and top DE genes,
# used only as an internal check when Bhaduri primary10X is absent.
module_sets <- list(m2=f$Member.Genes[[2]], m13=f$Member.Genes[[13]], m14=f$Member.Genes[[14]])
score_rows <- list()
for (set_name in names(module_sets)) {
  ids <- gene_map_first$gene_id_short[match(toupper(module_sets[[set_name]]), gene_map_first$gene_name)]
  ids <- ids[!is.na(ids) & ids %in% rownames(mat)]
  score <- if (length(ids) > 0) colMeans(mat[ids, , drop=FALSE], na.rm=TRUE) else rep(NA_real_, ncol(mat))
  score_rows[[set_name]] <- data.frame(set_name=set_name, pseudobulk_index=seq_along(score), mean_scaled_expression=score)
}
score_df <- do.call(rbind, score_rows)
write.table(score_df, file.path(out_dir, "fcm_local_module_gene_set_scores.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

summary <- list(
  fcm_cells_sc_rna = nrow(f$sc.colData),
  pseudobulks = nrow(pb),
  modules = nrow(f$Module.Centroids),
  variable_genes = nrow(f$Data.Matrix),
  full_genes = nrow(f$Full.Data.Matrix),
  member_genes_total = sum(vapply(f$Member.Genes, length, numeric(1))),
  jaccard_edges_gt_0_2 = nrow(f$Overlap.Segments.Filtered$jaccard),
  motif_matrix_peaks = nrow(mm),
  motif_matrix_motifs = ncol(mm),
  module13_linked_peaks = length(unique(f$Module.Linked.Peaks[[13]])),
  module14_linked_peaks = length(unique(f$Module.Linked.Peaks[[14]])),
  a1_hes_pseudobulks = sum(pb_labels$astro_label == "A1-HES", na.rm=TRUE),
  a2_olig_pseudobulks = sum(pb_labels$astro_label == "A2-OLIG", na.rm=TRUE),
  selected_genes_found_full = sum(selected_map$found_in_map & selected_map$in_full, na.rm=TRUE),
  selected_genes_found_membership = sum(selected_map$found_in_map & selected_map$in_member, na.rm=TRUE)
)
write_json(summary, file.path(out_dir, "fcm_export_summary.json"), pretty=TRUE, auto_unbox=TRUE)
"""
    r_file = out_dir / "export_fcm_figure4_5_tables.R"
    r_file.write_text(r_code, encoding="utf-8")
    cmd = [str(rscript), str(r_file), str(out_dir), str(fcm_rds), str(all_genes_rds), str(motif_match_rds)]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    (out_dir / "export_fcm_figure4_5_tables.stdout").write_text(proc.stdout, encoding="utf-8")
    (out_dir / "export_fcm_figure4_5_tables.stderr").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(
            "R FCM export failed with code "
            f"{proc.returncode}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return {
        "gene_map": out_dir / "fcm_gene_map_first.tsv",
        "pseudobulk_umap": out_dir / "fcm_pseudobulk_umap.tsv",
        "module_centroids": out_dir / "fcm_module_centroids.tsv",
        "jaccard_segments": out_dir / "fcm_jaccard_segments_filtered.tsv",
        "module_centers": out_dir / "fcm_module_centers.tsv",
        "selected_expression": out_dir / "fcm_selected_gene_expression.tsv",
        "selected_membership": out_dir / "fcm_selected_gene_membership.tsv",
        "selected_status": out_dir / "fcm_selected_gene_status.tsv",
        "member_genes": out_dir / "fcm_member_genes.tsv",
        "member_weights": out_dir / "fcm_member_gene_weights.tsv",
        "motif_enrichment": out_dir / "fcm_m13_vs_m14_motif_enrichment.tsv",
        "cluster_stats": out_dir / "fcm_figure5_aqp4_cluster_stats.tsv",
        "pb_labels": out_dir / "fcm_figure5_pseudobulk_a1_a2_labels.tsv",
        "de": out_dir / "fcm_figure5_local_a2_vs_a1_de.tsv",
        "module_scores": out_dir / "fcm_local_module_gene_set_scores.tsv",
        "summary": out_dir / "fcm_export_summary.json",
    }


def read_tables(paths: dict[str, Path]) -> dict[str, Any]:
    data: dict[str, Any] = {}
    for key, path in paths.items():
        if path.suffix == ".json":
            data[key] = json.loads(path.read_text())
        else:
            data[key] = pd.read_csv(path, sep="\t")
    return data


def pivot_expression(expr: pd.DataFrame) -> pd.DataFrame:
    wide = expr.pivot_table(index="gene_symbol", columns="pseudobulk_index", values="value", aggfunc="mean")
    return wide


def plot_figure4a(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    pb = data["pseudobulk_umap"]
    summary = data["summary"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    groups = sorted(pb["group"].astype(str).unique())
    cmap = plt.get_cmap("tab20")
    for i, grp in enumerate(groups):
        sub = pb[pb["group"].astype(str) == grp]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=7, alpha=0.7, color=cmap(i % 20), label=grp)
    ax.set_title("Fig. 4A public FCM pseudobulk embedding")
    ax.set_xlabel("FCM UMAP1")
    ax.set_ylabel("FCM UMAP2")
    ax.legend(markerscale=2, fontsize=6, ncol=2, frameon=False)
    axes[1].axis("off")
    txt = (
        "Author public FCM object\n\n"
        f"scRNA glial cells: {summary['fcm_cells_sc_rna']:,}\n"
        f"50-cell pseudobulks: {summary['pseudobulks']:,}\n"
        f"FCM modules: {summary['modules']}\n"
        f"variable genes: {summary['variable_genes']:,}\n"
        f"full matrix genes: {summary['full_genes']:,}\n\n"
        "This panel is a computational\n"
        "schematic/overview, not the\n"
        "paper's exact graphical schematic."
    )
    axes[1].text(0, 1, txt, va="top", ha="left", fontsize=10)
    fig.suptitle("Regulatory logic of glial specification: public FCM resource", y=1.02)
    return save_triple(fig, out_dir, "figure4A_fcm_public_overview")


def plot_heatmap(matrix: np.ndarray, row_labels: list[str], title: str, out_dir: Path, stem: str, *, cmap: str = "RdBu_r", vlim: float = 2.5) -> dict[str, str]:
    fig_h = max(3.2, 0.28 * len(row_labels) + 1.2)
    fig, ax = plt.subplots(figsize=(9, fig_h))
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=-vlim, vmax=vlim)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xticks([])
    ax.set_xlabel("pseudobulks ordered by glial cluster and pseudotime")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("row z-score")
    return save_triple(fig, out_dir, stem)


def plot_figure4b(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    pb = data["pseudobulk_umap"].copy()
    centers = data["module_centers"].copy()
    order = pb.sort_values(["group", "pseudotime", "pseudobulk_index"])["pseudobulk_index"].astype(int).to_numpy() - 1
    row_order = [m for m in SELECTED_MODULES if m in centers["module"].astype(str).values]
    centers["module"] = centers["module"].astype(str)
    cols = [c for c in centers.columns if c != "module"]
    mat = centers.set_index("module").loc[row_order, cols].to_numpy(dtype=float)[:, order]
    mat = np.vstack([safe_z(row) for row in mat])
    return plot_heatmap(mat, [f"m{m}" for m in row_order], "Fig. 4B module expression across public FCM pseudobulks", out_dir, "figure4B_module_expression_heatmap")


def plot_figure4c(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    pb = data["pseudobulk_umap"].copy()
    expr = pivot_expression(data["selected_expression"])
    order = pb.sort_values(["group", "pseudotime", "pseudobulk_index"])["pseudobulk_index"].astype(int).to_numpy()
    genes = [g for g in FIG4C_GENES if g in expr.index]
    mat = expr.loc[genes, order].to_numpy(dtype=float)
    mat = np.vstack([safe_z(row) for row in mat])
    return plot_heatmap(mat, genes, "Fig. 4C selected gene expression across public FCM pseudobulks", out_dir, "figure4C_selected_gene_heatmap", cmap="Reds", vlim=2.5)


def draw_expression_scatter(ax: plt.Axes, pb: pd.DataFrame, values: pd.Series, title: str, cmap: str = "viridis") -> None:
    vals = values.reindex(pb["pseudobulk_index"].astype(int)).to_numpy(dtype=float)
    vals = safe_z(vals)
    sc = ax.scatter(pb["UMAP1"], pb["UMAP2"], c=vals, s=8, cmap=cmap, alpha=0.85, linewidths=0)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)


def draw_membership(ax: plt.Axes, cent: pd.DataFrame, seg: pd.DataFrame, gene_mem: pd.DataFrame, title: str) -> None:
    if not seg.empty:
        for _, r in seg.iterrows():
            ax.plot([r["x1"], r["x2"]], [r["y1"], r["y2"]], color="0.75", linewidth=max(0.3, 2.5 * float(r["overlap"])), zorder=1)
    vals = gene_mem.set_index(gene_mem["module"].astype(str))["membership"]
    modules = cent["module.name"].astype(str)
    colors = vals.reindex(modules).to_numpy(dtype=float)
    if np.all(~np.isfinite(colors)):
        colors = np.zeros(len(modules))
    sc = ax.scatter(cent["UMAP1"], cent["UMAP2"], c=colors, s=260, cmap="magma", vmin=0, vmax=max(0.25, np.nanmax(colors)), edgecolor="black", linewidth=0.5, zorder=3)
    for _, r in cent.iterrows():
        ax.text(r["UMAP1"], r["UMAP2"], f"m{r['module.name']}", ha="center", va="center", fontsize=7, color="white", zorder=4)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)


def plot_figure4d(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    pb = data["pseudobulk_umap"]
    centers = data["module_centers"].copy()
    centers["module"] = centers["module"].astype(str)
    cols = [c for c in centers.columns if c != "module"]
    selected = ["10", "3", "8", "6", "12", "1", "13", "14"]
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.4), constrained_layout=True)
    for ax, module in zip(axes.flat, selected):
        row = centers.set_index("module").loc[module, cols]
        vals = pd.Series(row.to_numpy(dtype=float), index=np.arange(1, len(row) + 1))
        draw_expression_scatter(ax, pb, vals, f"module m{module}", cmap="plasma")
    fig.suptitle("Fig. 4D mean scaled module expression in FCM embedding", y=1.03)
    return save_triple(fig, out_dir, "figure4D_selected_module_umaps")


def plot_figure4e(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    pb = data["pseudobulk_umap"]
    cent = data["module_centroids"]
    seg = data["jaccard_segments"]
    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(pb["UMAP1"], pb["UMAP2"], c=pb["pseudotime"], cmap="viridis", s=9, alpha=0.55, linewidths=0)
    if not seg.empty:
        for _, r in seg.iterrows():
            ax.plot([r["x1"], r["x2"]], [r["y1"], r["y2"]], color="black", alpha=0.45, linewidth=max(0.4, 4 * float(r["overlap"])))
    ax.scatter(cent["UMAP1"], cent["UMAP2"], s=260, facecolor="white", edgecolor="black", linewidth=1.2, zorder=3)
    for _, r in cent.iterrows():
        ax.text(r["UMAP1"], r["UMAP2"], f"m{r['module.name']}", ha="center", va="center", fontsize=8, zorder=4)
    ax.set_title("Fig. 4E public FCM module centroids; Jaccard links >0.2")
    ax.set_xlabel("FCM UMAP1")
    ax.set_ylabel("FCM UMAP2")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("pseudotime")
    return save_triple(fig, out_dir, "figure4E_module_centroid_jaccard_network")


def plot_gene_membership_panel(data: dict[str, Any], out_dir: Path, genes: list[str], title: str, stem: str) -> dict[str, str]:
    pb = data["pseudobulk_umap"]
    cent = data["module_centroids"]
    seg = data["jaccard_segments"]
    expr = pivot_expression(data["selected_expression"])
    mem = data["selected_membership"]
    fig, axes = plt.subplots(len(genes), 2, figsize=(9.5, 3.1 * len(genes)), constrained_layout=True)
    if len(genes) == 1:
        axes = np.asarray([axes])
    for i, gene in enumerate(genes):
        if gene in expr.index:
            vals = expr.loc[gene]
        else:
            vals = pd.Series(np.zeros(len(pb)), index=pb["pseudobulk_index"].astype(int))
        draw_expression_scatter(axes[i, 0], pb, vals, f"{gene}: pseudobulk expression", cmap="Reds")
        gmem = mem[mem["gene_symbol"] == gene]
        draw_membership(axes[i, 1], cent, seg, gmem, f"{gene}: FCM module membership")
    fig.suptitle(title, y=1.02)
    return save_triple(fig, out_dir, stem)


def plot_figure5b(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    motif = data["motif_enrichment"].copy()
    motif = motif.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2_odds_ratio", "neg_log10_fdr"])
    motif["abs_x"] = motif["log2_odds_ratio"].abs()
    top = motif.sort_values(["neg_log10_fdr", "abs_x"], ascending=False).head(45)
    highlight = motif[motif["motif_name"].str.upper().isin(HIGHLIGHT_MOTIFS)]
    plot_df = pd.concat([top, highlight], ignore_index=True).drop_duplicates("motif_id")
    fig, ax = plt.subplots(figsize=(8, 5.2))
    ax.scatter(motif["log2_odds_ratio"], motif["neg_log10_fdr"], s=12, color="0.75", alpha=0.7, label="all motifs")
    ax.scatter(plot_df["log2_odds_ratio"], plot_df["neg_log10_fdr"], s=30, color="#2C7FB8", alpha=0.9, label="top/highlight")
    for _, r in plot_df.iterrows():
        if r["motif_name"].upper() in HIGHLIGHT_MOTIFS or r["neg_log10_fdr"] > plot_df["neg_log10_fdr"].quantile(0.82):
            ax.text(r["log2_odds_ratio"], r["neg_log10_fdr"], str(r["motif_name"]), fontsize=7)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("log2 odds ratio: motif in m13 GREs vs m14 GREs")
    ax.set_ylabel("-log10 BH FDR")
    ax.set_title("Fig. 5B motif enrichment in public FCM linked peaks (m13 vs m14)")
    ax.legend(frameon=False, fontsize=8)
    return save_triple(fig, out_dir, "figure5B_m13_vs_m14_motif_enrichment")


def plot_figure5c(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    labels = data["pb_labels"]
    colors = {"A1-HES": "#D95F02", "A2-OLIG": "#1B9E77", "Other": "0.78"}
    fig, ax = plt.subplots(figsize=(6.5, 5.7))
    for lab in ["Other", "A1-HES", "A2-OLIG"]:
        sub = labels[labels["astro_label"].fillna("Other") == lab]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=12 if lab != "Other" else 7, color=colors[lab], alpha=0.9 if lab != "Other" else 0.45, label=lab, linewidths=0)
    ax.set_title("Fig. 5C local AQP4-positive FCM reclustering labels")
    ax.set_xlabel("FCM UMAP1")
    ax.set_ylabel("FCM UMAP2")
    ax.legend(frameon=False)
    return save_triple(fig, out_dir, "figure5C_aqp4_positive_a1_a2_reclustering")


def plot_figure5d(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    de = data["de"].copy()
    de = de.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2fc_a2_vs_a1", "fdr_bh"])
    de["neg_log10_fdr"] = -np.log10(np.clip(de["fdr_bh"].astype(float), np.finfo(float).tiny, 1.0))
    sig = de["fdr_bh"] < 1e-6
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(de.loc[~sig, "log2fc_a2_vs_a1"], de.loc[~sig, "neg_log10_fdr"], s=6, color="0.75", alpha=0.55)
    ax.scatter(de.loc[sig, "log2fc_a2_vs_a1"], de.loc[sig, "neg_log10_fdr"], s=8, color="#2B8CBE", alpha=0.7)
    for gene in DE_LABEL_GENES:
        rows = de[de["gene_symbol"].astype(str).str.upper() == gene]
        if not rows.empty:
            r = rows.iloc[0]
            ax.text(r["log2fc_a2_vs_a1"], r["neg_log10_fdr"], gene, fontsize=8)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("log2 fold-change A2-OLIG vs A1-HES")
    ax.set_ylabel("-log10 BH FDR")
    ax.set_title("Fig. 5D local pseudo-bulk differential expression (not DESeq2 parity)")
    return save_triple(fig, out_dir, "figure5D_local_a2_vs_a1_de_volcano")


def plot_figure5ef_gap(data: dict[str, Any], out_dir: Path, bhaduri_dir: Path | None) -> dict[str, str]:
    summary = data["summary"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.axis("off")
    has_bhaduri = bool(bhaduri_dir and bhaduri_dir.exists() and any(bhaduri_dir.iterdir()))
    txt = (
        "Figure 5E/5F external Bhaduri primary10X gate\n\n"
        f"Bhaduri directory supplied: {str(bhaduri_dir) if bhaduri_dir else 'none'}\n"
        f"Local staged files detected: {has_bhaduri}\n\n"
        "The Trevino methods state that normalized data were downloaded\n"
        "from UCSC Cell Browser dataset organoidreportcard/primary10X.\n"
        "No local primary10X expression/metadata bundle is staged for this\n"
        "run yet, so this script does not claim Figure 5E/5F reproduction.\n\n"
        "Internal FCM checks completed instead:\n"
        f"A1-HES pseudobulks: {summary.get('a1_hes_pseudobulks', 'NA')}\n"
        f"A2-OLIG pseudobulks: {summary.get('a2_olig_pseudobulks', 'NA')}\n"
        "Next recovery path: download/stage UCSC Cell Browser primary10X\n"
        "and render Bhaduri astrocyte UMAP + module/DE gene-set scores."
    )
    ax.text(0.01, 0.98, txt, va="top", ha="left", fontsize=10)
    return save_triple(fig, out_dir, "figure5E_5F_bhaduri_primary10x_resource_gate")


def collect_outputs_hashes(out_dir: Path) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for path in sorted(out_dir.glob("figure*.png")):
        hashes[str(path)] = sha256_file(path)
    return hashes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--fcm-rds", default=Path("/home/zerlinshen/projects/wave5-trevino/inputs/brainchromatin_s3/rds/FCM_Object.RDS"), type=Path)
    ap.add_argument("--all-genes-rds", default=Path("/home/zerlinshen/projects/wave5-trevino/inputs/brainchromatin_s3/rds/AllGenes_GenomicRanges.RDS"), type=Path)
    ap.add_argument("--motif-match-rds", default=Path("/home/zerlinshen/projects/wave5-trevino/inputs/brainchromatin_s3/rds/scATAC_MotifMatchMatrix.RDS"), type=Path)
    ap.add_argument("--rscript", default=Path("/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript"), type=Path)
    ap.add_argument("--bhaduri-dir", default=None, type=Path)
    args = ap.parse_args()

    configure_matplotlib()
    run_dir = args.project_root / "runs" / args.run_id
    fig_dir = run_dir / "python" / "figures" / "public_matrix" / "figure_4_5_fcm_public_panels"
    evidence_dir = run_dir / "python" / "figure_reproduction_evidence"
    fig_dir.mkdir(parents=True, exist_ok=True)
    evidence_dir.mkdir(parents=True, exist_ok=True)

    required = [args.fcm_rds, args.all_genes_rds, args.motif_match_rds, args.rscript]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required resource(s): " + ", ".join(missing))

    paths = run_r_export(
        rscript=args.rscript,
        out_dir=fig_dir,
        fcm_rds=args.fcm_rds,
        all_genes_rds=args.all_genes_rds,
        motif_match_rds=args.motif_match_rds,
    )
    data = read_tables(paths)

    outputs: dict[str, dict[str, str]] = {}
    outputs["4A"] = plot_figure4a(data, fig_dir)
    outputs["4B"] = plot_figure4b(data, fig_dir)
    outputs["4C"] = plot_figure4c(data, fig_dir)
    outputs["4D"] = plot_figure4d(data, fig_dir)
    outputs["4E"] = plot_figure4e(data, fig_dir)
    outputs["4F"] = plot_gene_membership_panel(data, fig_dir, FIG4F_GENES, "Fig. 4F ASCL1/HES4/OLIG1 expression and FCM membership", "figure4F_ascl1_hes4_olig1_membership_expression")
    outputs["4G"] = plot_gene_membership_panel(data, fig_dir, FIG4G_GENES, "Fig. 4G EOMES/AQP4/MBP expression and FCM membership", "figure4G_eomes_aqp4_mbp_membership_expression")
    outputs["4H"] = plot_gene_membership_panel(data, fig_dir, FIG4H_GENES, "Fig. 4H ASCL1/OLIG1/EGFR expression and FCM membership", "figure4H_ascl1_olig1_egfr_membership_expression")
    outputs["4J"] = plot_gene_membership_panel(data, fig_dir, FIG4J_GENES, "Fig. 4J PDGFRA/SPARCL1 expression and FCM membership", "figure4J_pdgfra_sparcl1_membership_expression")
    outputs["5A"] = plot_gene_membership_panel(data, fig_dir, FIG5A_GENES, "Fig. 5A astrocyte-associated genes in public FCM resource", "figure5A_aqp4_tnc_aldh2_apoe_membership_expression")
    outputs["5B"] = plot_figure5b(data, fig_dir)
    outputs["5C"] = plot_figure5c(data, fig_dir)
    outputs["5D"] = plot_figure5d(data, fig_dir)
    outputs["5E_5F_gate"] = plot_figure5ef_gap(data, fig_dir, args.bhaduri_dir)

    selected_status = data["selected_status"]
    motif = data["motif_enrichment"]
    cluster_stats = data["cluster_stats"]
    de = data["de"]
    summary = data["summary"]

    def motif_row(name: str) -> dict[str, Any] | None:
        rows = motif[motif["motif_name"].astype(str).str.upper() == name]
        if rows.empty:
            return None
        r = rows.sort_values("fdr_bh").iloc[0]
        return {
            "motif_id": r["motif_id"],
            "log2_odds_ratio": float(r["log2_odds_ratio"]),
            "fdr_bh": float(r["fdr_bh"]),
            "neg_log10_fdr": float(r["neg_log10_fdr"]),
        }

    panel_status = {
        "4A": "GENERATED_AUTHOR_FCM_PUBLIC_PSEUDOBULK_OVERVIEW",
        "4B": "REPRODUCED_AUTHOR_FCM_MODULE_HEATMAP_FROM_CENTERS",
        "4C": "REPRODUCED_AUTHOR_FCM_SELECTED_GENE_HEATMAP",
        "4D": "GENERATED_AUTHOR_FCM_MODULE_SCORE_UMAPS",
        "4E": "REPRODUCED_AUTHOR_FCM_CENTROID_JACCARD_NETWORK",
        "4F": "REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION",
        "4G": "REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION",
        "4H": "REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION",
        "4I": "IMAGE_PANEL_NOT_COMPUTATIONAL_REPRODUCTION_LOCAL_PDF_ONLY",
        "4J": "REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION",
        "4K": "IMAGE_PANEL_NOT_COMPUTATIONAL_REPRODUCTION_LOCAL_PDF_ONLY",
        "5A": "REPRODUCED_AUTHOR_FCM_ASTRO_GENE_MEMBERSHIP_EXPRESSION",
        "5B": "GENERATED_AUTHOR_FCM_LINKED_PEAK_MOTIF_ENRICHMENT_M13_VS_M14",
        "5C": "GENERATED_AUTHOR_FCM_AQP4_POSITIVE_RECLUSTERING_LOCAL_RULE",
        "5D": "GENERATED_PUBLIC_FCM_PSEUDOBULK_DE_NOT_DESEQ2_PARITY",
        "5E": "METHOD_RESOURCE_GAP_BHADURI_PRIMARY10X_NOT_STAGED",
        "5F": "METHOD_RESOURCE_GAP_BHADURI_PRIMARY10X_NOT_STAGED",
    }

    evidence = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "project_root": str(args.project_root),
        "run_id": args.run_id,
        "output_dir": str(fig_dir),
        "resources": {
            "fcm_rds": str(args.fcm_rds),
            "fcm_rds_sha256": sha256_file(args.fcm_rds),
            "all_genes_rds": str(args.all_genes_rds),
            "all_genes_rds_sha256": sha256_file(args.all_genes_rds),
            "motif_match_rds": str(args.motif_match_rds),
            "motif_match_rds_sha256": sha256_file(args.motif_match_rds),
            "rscript": str(args.rscript),
            "bhaduri_dir": str(args.bhaduri_dir) if args.bhaduri_dir else None,
        },
        "summary": summary,
        "panel_status": panel_status,
        "checks": {
            "pseudobulks_eq_1267": summary.get("pseudobulks") == 1267,
            "modules_eq_14": summary.get("modules") == 14,
            "jaccard_edges_gt_0_2_eq_32": summary.get("jaccard_edges_gt_0_2") == 32,
            "motif_matrix_motifs_ge_400": summary.get("motif_matrix_motifs", 0) >= 400,
            "module13_and_14_have_linked_peaks": summary.get("module13_linked_peaks", 0) > 0 and summary.get("module14_linked_peaks", 0) > 0,
            "a1_and_a2_labels_nonempty": summary.get("a1_hes_pseudobulks", 0) > 0 and summary.get("a2_olig_pseudobulks", 0) > 0,
            "selected_fig4_5_genes_found_full_ge_18": summary.get("selected_genes_found_full", 0) >= 18,
            "all_output_hashes_recorded": True,
        },
        "selected_gene_status_table": str(paths["selected_status"]),
        "motif_highlights": {name: motif_row(name) for name in HIGHLIGHT_MOTIFS},
        "a1_a2_cluster_stats_table": str(paths["cluster_stats"]),
        "a1_a2_cluster_stats_preview": cluster_stats.to_dict(orient="records"),
        "local_de_table": str(paths["de"]),
        "local_de_top10": de.head(10).replace({np.nan: None}).to_dict(orient="records"),
        "outputs": outputs,
        "output_png_hashes": collect_outputs_hashes(fig_dir),
        "truth_boundaries": [
            "FCM-derived Figures 4A-H/J and 5A-D are rendered from author public FCM/linkage/motif RDS resources, not from FASTQ/fragments.",
            "Figure 4I/4K are immunohistochemistry image panels and are recorded as non-computational image panels rather than recomputed matrix outputs.",
            "Figure 5D uses a local Welch-test pseudo-bulk differential expression proxy; exact paper DESeq2/Table S5 parity is not claimed.",
            "Figure 5E/5F require UCSC Cell Browser organoidreportcard/primary10X Bhaduri data; no local staged bundle was detected by this script, so no reproduction is claimed for those panels yet.",
        ],
    }
    evidence_path = evidence_dir / "figure4_5_public_author_fcm_evidence.json"
    evidence_path.write_text(json.dumps(evidence, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps({"evidence": str(evidence_path), "output_dir": str(fig_dir), "panel_status": panel_status, "checks": evidence["checks"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
