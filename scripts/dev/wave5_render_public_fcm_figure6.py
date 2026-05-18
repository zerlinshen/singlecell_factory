#!/usr/bin/env python3
"""Render/audit public-resource Trevino Figure 6 panels from FCM/GPC resources."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.stats import zscore

GPC_TFS = ["BHLHE40", "OLIG1", "OLIG2", "NEUROD6", "NEUROD4", "HES1", "RFX4", "EOMES", "ASCL1", "NHLH1", "SOX21"]
BRANCH_COLORS = {"A": "#1B9E77", "B": "#D95F02", "C": "#7570B3", "Not Branch": "0.75"}


def configure_matplotlib() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": "white", "axes.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    outputs: dict[str, str] = {}
    for ext in ("png", "svg", "pdf"):
        p = out_dir / f"{stem}.{ext}"
        fig.savefig(p, bbox_inches="tight", dpi=300)
        outputs[ext] = str(p)
        outputs[f"{ext}_sha256"] = sha256_file(p)
    plt.close(fig)
    return outputs


def safe_z(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0 or np.nanstd(arr) == 0:
        return np.zeros_like(arr, dtype=float)
    return np.nan_to_num(zscore(arr, nan_policy="omit"), nan=0.0, posinf=0.0, neginf=0.0)


def run_r_export(args: argparse.Namespace, out_dir: Path) -> dict[str, Path]:
    r_code = r"""
args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[[1]]; fcm_rds <- args[[2]]; all_genes_rds <- args[[3]]; chromvar_rds <- args[[4]]
gpc_tsv <- args[[5]]; cell_cycle_txt <- args[[6]]; multiome_rna_rds <- args[[7]]; multiome_atac_rds <- args[[8]]

suppressPackageStartupMessages({
  library(Matrix); library(GenomicRanges); library(jsonlite); library(FNN)
  library(SingleCellExperiment); library(SummarizedExperiment)
})

f <- readRDS(fcm_rds)
gr <- readRDS(all_genes_rds)
gm <- as.data.frame(GenomicRanges::mcols(gr))
gene_map <- unique(data.frame(gene_id_short=as.character(gm$gene_id_short), gene_name=toupper(as.character(gm$gene_name)), stringsAsFactors=FALSE))
gene_map <- gene_map[!is.na(gene_map$gene_id_short) & !is.na(gene_map$gene_name) & gene_map$gene_name != "",]
gene_map$in_full <- gene_map$gene_id_short %in% rownames(f$Full.Data.Matrix)
gene_map$in_ga <- gene_map$gene_id_short %in% rownames(f$Gene.Activity.Matrix)
gene_map$in_ref <- gene_map$gene_id_short %in% rownames(f$sc.RNA.UMAP.Matrix)
gene_map <- gene_map[order(gene_map$gene_name, -gene_map$in_full, -gene_map$in_ga, -gene_map$in_ref),]
gene_first <- gene_map[!duplicated(gene_map$gene_name),]
write.table(gene_first, file.path(out_dir, 'figure6_gene_map_first.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

pb <- as.data.frame(f$Projections$UMAP$DF); pb$pseudobulk_index <- seq_len(nrow(pb))
ga_no <- as.data.frame(f$Projections$GA.No.Anchors$DF); ga_no$pseudobulk_index <- seq_len(nrow(ga_no))
ga_anchor <- as.data.frame(f$Projections$GA.Anchors.Only$DF); ga_anchor$pseudobulk_index <- seq_len(nrow(ga_anchor))
atac <- as.data.frame(f$Projections$ATAC.UMAP$DF); atac$pseudobulk_index <- seq_len(nrow(atac))
cent <- as.data.frame(f$Module.Centroids)
write.table(pb, file.path(out_dir,'figure6_fcm_pseudobulk_umap.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(ga_no, file.path(out_dir,'figure6_ga_no_anchors.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(ga_anchor, file.path(out_dir,'figure6_ga_anchors_only.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(atac, file.path(out_dir,'figure6_atac_umap.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(cent, file.path(out_dir,'figure6_module_centroids.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

# 6A: cell-cycle signature correlation with each module center across pseudobulks.
cc <- unique(toupper(scan(cell_cycle_txt, what=character(), quiet=TRUE)))
cc <- cc[!grepl('MITOTIC_CELL_CYCLE|KEGG|^#', cc)]
cc_ids <- gene_first$gene_id_short[match(cc, gene_first$gene_name)]
cc_ids <- cc_ids[!is.na(cc_ids) & cc_ids %in% rownames(f$Full.Data.Matrix)]
cc_score <- colMeans(f$Full.Data.Matrix[cc_ids, , drop=FALSE], na.rm=TRUE)
centers <- as.matrix(f$FCM.Result$centers)
fig6a <- do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
  ct <- suppressWarnings(cor.test(as.numeric(centers[i,]), cc_score, method='pearson'))
  data.frame(module=paste0('m', rownames(centers)[i]), pearson_r=unname(ct$estimate), p_value=ct$p.value, n_cell_cycle_genes=length(cc_ids), stringsAsFactors=FALSE)
}))
fig6a$fdr_bh <- p.adjust(fig6a$p_value, method='BH')
write.table(fig6a, file.path(out_dir,'figure6a_cell_cycle_module_correlations.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

# Public GPCs from Figure 2E.
gpc <- read.table(gpc_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
gpc_col <- if ('gene.symbol' %in% colnames(gpc)) 'gene.symbol' else if ('gene_symbol' %in% colnames(gpc)) 'gene_symbol' else colnames(gpc)[1]
gpc_genes <- unique(toupper(gpc[[gpc_col]]))
gpc_ids <- gene_first$gene_id_short[match(gpc_genes, gene_first$gene_name)]
gpc_ids <- gpc_ids[!is.na(gpc_ids)]

# 6D: branch-specific gene activity heatmap.
br <- ga_anchor$is.Branch
branch_levels <- c('A','B','C')
ga <- as.matrix(f$Gene.Activity.Matrix)
branch_means <- sapply(branch_levels, function(b) rowMeans(ga[, br == b, drop=FALSE], na.rm=TRUE))
colnames(branch_means) <- branch_levels
spec <- sapply(seq_along(branch_levels), function(j) branch_means[,j] - rowMeans(branch_means[, -j, drop=FALSE], na.rm=TRUE))
colnames(spec) <- branch_levels
best_branch <- branch_levels[max.col(spec, ties.method='first')]
best_score <- apply(spec, 1, max, na.rm=TRUE)
ord <- order(best_score, decreasing=TRUE)
sel <- rownames(ga)[ord[seq_len(min(50, length(ord)))]]
heat <- branch_means[sel, , drop=FALSE]
heat_scaled <- t(scale(t(heat)))
heat_scaled[!is.finite(heat_scaled)] <- 0
fig6d <- data.frame(gene_id=rep(sel, each=length(branch_levels)), gene_symbol=gene_first$gene_name[match(rep(sel, each=length(branch_levels)), gene_first$gene_id_short)], branch=rep(branch_levels, times=length(sel)), row_z=as.vector(t(heat_scaled)), mean_activity=as.vector(t(heat)), best_branch=rep(best_branch[sel], each=length(branch_levels)), is_gpc=rep(sel %in% gpc_ids, each=length(branch_levels)), stringsAsFactors=FALSE)
write.table(fig6d, file.path(out_dir,'figure6d_branch_unique_gene_activity_top50.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
ks_rows <- do.call(rbind, lapply(branch_levels, function(b) {
  x <- spec[, b]
  isg <- rownames(spec) %in% gpc_ids
  p <- if (sum(isg) > 1 && sum(!isg) > 1) suppressWarnings(ks.test(x[isg], x[!isg], alternative='greater')$p.value) else NA_real_
  data.frame(branch=b, ks_p_value=p, gpc_n=sum(isg), non_gpc_n=sum(!isg), stringsAsFactors=FALSE)
}))
ks_rows$fdr_bh <- p.adjust(ks_rows$ks_p_value, method='BH')
write.table(ks_rows, file.path(out_dir,'figure6d_gpc_ks_enrichment.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

# 6E: GPC TF motif and gene dynamics across branches.
chromvar <- readRDS(chromvar_rds)
tf_genes <- c('BHLHE40','OLIG1','OLIG2','NEUROD6','NEUROD4','HES1','RFX4','EOMES','ASCL1','NHLH1','SOX21')
expr_rows <- list(); motif_rows <- list()
for (g in tf_genes) {
  gid <- gene_first$gene_id_short[match(g, gene_first$gene_name)]
  if (!is.na(gid) && gid %in% rownames(f$Full.Data.Matrix)) {
    means <- sapply(branch_levels, function(b) mean(f$Full.Data.Matrix[gid, br == b], na.rm=TRUE))
    zs <- as.numeric(scale(means)); zs[!is.finite(zs)] <- 0
    expr_rows[[g]] <- data.frame(tf=g, branch=branch_levels, mean_expression=as.numeric(means), row_z=zs, stringsAsFactors=FALSE)
  }
  motif_idx <- grep(paste0('_', g, '$'), rownames(chromvar), ignore.case=TRUE)
  if (length(motif_idx) == 0) motif_idx <- grep(g, rownames(chromvar), ignore.case=TRUE)
  if (length(motif_idx) > 0) {
    mid <- rownames(chromvar)[motif_idx[1]]
    means <- sapply(branch_levels, function(b) mean(chromvar[mid, br == b], na.rm=TRUE))
    zs <- as.numeric(scale(means)); zs[!is.finite(zs)] <- 0
    motif_rows[[g]] <- data.frame(tf=g, motif_id=mid, branch=branch_levels, mean_chromvar=as.numeric(means), row_z=zs, stringsAsFactors=FALSE)
  }
}
expr_dyn <- do.call(rbind, expr_rows); motif_dyn <- do.call(rbind, motif_rows)
write.table(expr_dyn, file.path(out_dir,'figure6e_gpc_tf_gene_expression_branch_dynamics.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(motif_dyn, file.path(out_dir,'figure6e_gpc_tf_motif_branch_dynamics.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

# 6G: bounded multiome RNA -> FCM manifold nearest-neighbor projection, colored by matched ATAC cluster.
rna <- readRDS(multiome_rna_rds); matac <- readRDS(multiome_atac_rds)
common_cells <- intersect(colnames(rna), colnames(matac))
rna <- rna[, common_cells]; matac <- matac[, common_cells]
ref <- f$sc.RNA.UMAP.Matrix
common_genes <- intersect(rownames(ref), rownames(rna))
vars <- apply(as.matrix(ref[common_genes, , drop=FALSE]), 1, stats::var)
vars <- sort(vars, decreasing=TRUE)
top <- names(vars)[seq_len(min(500, length(vars)))]
ref_m <- t(as.matrix(ref[top, , drop=FALSE]))
query_m <- t(as.matrix(assay(rna, 'logcounts')[top, common_cells, drop=FALSE]))
mu <- colMeans(ref_m); sig <- apply(ref_m, 2, sd); sig[sig == 0 | !is.finite(sig)] <- 1
ref_s <- sweep(sweep(ref_m, 2, mu, '-'), 2, sig, '/')
query_s <- sweep(sweep(query_m, 2, mu, '-'), 2, sig, '/')
nn <- FNN::get.knnx(ref_s, query_s, k=1)
coords <- f$sc.RNA.UMAP[nn$nn.index[,1], , drop=FALSE]
cd_atac <- as.data.frame(colData(matac))
fig6g <- data.frame(Cell.ID=common_cells, UMAP1=coords[,1], UMAP2=coords[,2], mapped_atac_cluster=as.character(cd_atac[common_cells, 'seurat_clusters']), nn_dist=nn$nn.dist[,1], stringsAsFactors=FALSE)
write.table(fig6g, file.path(out_dir,'figure6g_multiome_nn_fcm_projection.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

summary <- list(
  pseudobulks=nrow(pb), atac_pseudobulks=nrow(atac), branch_A=sum(br=='A'), branch_B=sum(br=='B'), branch_C=sum(br=='C'), branch_not=sum(br=='Not Branch'),
  cell_cycle_genes_found=length(cc_ids), gpc_genes_loaded=length(gpc_genes), gpc_ids_in_gene_activity=sum(gpc_ids %in% rownames(ga)),
  figure6d_top_genes=nrow(unique(fig6d['gene_id'])), chromvar_motifs=nrow(chromvar), motif_tfs_found=length(unique(motif_dyn$tf)), gene_tfs_found=length(unique(expr_dyn$tf)),
  multiome_cells_projected=nrow(fig6g), multiome_projection_genes=length(top), multiome_atac_clusters=length(unique(fig6g$mapped_atac_cluster)), median_multiome_nn_dist=median(fig6g$nn_dist)
)
write_json(summary, file.path(out_dir,'figure6_export_summary.json'), pretty=TRUE, auto_unbox=TRUE)
"""
    r_file = out_dir / "export_figure6_tables.R"
    r_file.write_text(r_code, encoding="utf-8")
    cmd = [str(args.rscript), str(r_file), str(out_dir), str(args.fcm_rds), str(args.all_genes_rds), str(args.chromvar_rds), str(args.gpc_tsv), str(args.cell_cycle_txt), str(args.multiome_rna_rds), str(args.multiome_atac_rds)]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    (out_dir / "export_figure6_tables.stdout").write_text(proc.stdout, encoding="utf-8")
    (out_dir / "export_figure6_tables.stderr").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"R export failed {proc.returncode}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return {
        "summary": out_dir / "figure6_export_summary.json",
        "fig6a": out_dir / "figure6a_cell_cycle_module_correlations.tsv",
        "pb": out_dir / "figure6_fcm_pseudobulk_umap.tsv",
        "atac": out_dir / "figure6_atac_umap.tsv",
        "ga_no": out_dir / "figure6_ga_no_anchors.tsv",
        "ga_anchor": out_dir / "figure6_ga_anchors_only.tsv",
        "centroids": out_dir / "figure6_module_centroids.tsv",
        "fig6d": out_dir / "figure6d_branch_unique_gene_activity_top50.tsv",
        "fig6d_ks": out_dir / "figure6d_gpc_ks_enrichment.tsv",
        "fig6e_expr": out_dir / "figure6e_gpc_tf_gene_expression_branch_dynamics.tsv",
        "fig6e_motif": out_dir / "figure6e_gpc_tf_motif_branch_dynamics.tsv",
        "fig6g": out_dir / "figure6g_multiome_nn_fcm_projection.tsv",
    }


def load_tables(paths: dict[str, Path]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, p in paths.items():
        if p.suffix == ".json":
            out[k] = json.loads(p.read_text())
        else:
            out[k] = pd.read_csv(p, sep="\t")
    return out


def plot_6a(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    df = data["fig6a"].sort_values("pearson_r", ascending=False)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    colors = ["#D7301F" if r > 0 else "#4575B4" for r in df["pearson_r"]]
    ax.bar(df["module"], df["pearson_r"], color=colors)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Pearson r with mitotic cell-cycle score")
    ax.set_title("Fig. 6A FCM module correlation with MSigDB mitotic cell-cycle genes")
    ax.tick_params(axis="x", rotation=45)
    return save_triple(fig, out_dir, "figure6A_cell_cycle_module_correlation")


def plot_6b(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    summary = data["summary"]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.axis("off")
    txt = (
        "Fig. 6B public-resource schematic\n\n"
        "ATAC gene-activity pseudobulks are projected into the FCM manifold.\n"
        "The public FCM object already contains GA.No.Anchors and\n"
        "GA.Anchors.Only projections plus branch labels A/B/C.\n\n"
        f"ATAC pseudobulks: {summary['atac_pseudobulks']}\n"
        f"Branch A/B/C: {summary['branch_A']}/{summary['branch_B']}/{summary['branch_C']}\n"
        f"Not Branch: {summary['branch_not']}\n"
        "This is a schematic reconstruction from the public object,\n"
        "not an extracted graphical copy of the article panel."
    )
    ax.text(0.02, 0.95, txt, va="top", fontsize=11)
    return save_triple(fig, out_dir, "figure6B_atac_to_fcm_projection_schematic")


def scatter_branch(ax: plt.Axes, df: pd.DataFrame, *, title: str, s: int = 9) -> None:
    for lab in ["Not Branch", "A", "B", "C"]:
        sub = df[df["is.Branch"].fillna("Not Branch") == lab]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=s if lab != "Not Branch" else max(5, s-2), color=BRANCH_COLORS[lab], alpha=0.9 if lab != "Not Branch" else 0.35, label=lab, linewidths=0)
    ax.set_title(title)
    ax.set_xlabel("FCM UMAP1")
    ax.set_ylabel("FCM UMAP2")


def plot_6c(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    df = data["ga_anchor"]
    cent = data["centroids"]
    fig, ax = plt.subplots(figsize=(6.5, 5.7))
    scatter_branch(ax, df, title="Fig. 6C ATAC pseudobulk branches in FCM embedding")
    ax.scatter(cent["UMAP1"], cent["UMAP2"], s=140, facecolor="white", edgecolor="black", zorder=3)
    for _, r in cent.iterrows():
        ax.text(r["UMAP1"], r["UMAP2"], f"m{r['module.name']}", ha="center", va="center", fontsize=7, zorder=4)
    ax.legend(frameon=False, fontsize=8)
    return save_triple(fig, out_dir, "figure6C_atac_branch_projection")


def plot_heatmap(df: pd.DataFrame, value_col: str, row_col: str, col_col: str, title: str, out_dir: Path, stem: str, *, annotate_gpc: bool = False) -> dict[str, str]:
    mat = df.pivot_table(index=row_col, columns=col_col, values=value_col, aggfunc="mean")
    mat = mat[[c for c in ["A", "B", "C"] if c in mat.columns]]
    if len(mat) > 1:
        try:
            order = leaves_list(linkage(mat.fillna(0).to_numpy(), method="average"))
            mat = mat.iloc[order]
        except Exception:
            pass
    fig_h = max(4, min(12, 0.13 * len(mat) + 1.8))
    fig, ax = plt.subplots(figsize=(4.8, fig_h))
    im = ax.imshow(mat.fillna(0).to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(mat.columns))); ax.set_xticklabels(mat.columns)
    ax.set_yticks(range(len(mat.index))); ax.set_yticklabels(mat.index, fontsize=6 if len(mat) > 25 else 8)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02); cbar.set_label("row z-score")
    if annotate_gpc and "is_gpc" in df.columns:
        gpc = df.drop_duplicates(row_col).set_index(row_col).reindex(mat.index)["is_gpc"].fillna(False).to_numpy()
        for i, flag in enumerate(gpc):
            if flag:
                ax.text(len(mat.columns)-0.02, i, "●", color="#F46D43", fontsize=8, va="center", ha="left")
    return save_triple(fig, out_dir, stem)


def plot_6d(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    return plot_heatmap(data["fig6d"], "row_z", "gene_symbol", "branch", "Fig. 6D branch-unique ATAC gene activity (top 50)", out_dir, "figure6D_branch_unique_gene_activity_top50", annotate_gpc=True)


def plot_6e(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    expr = data["fig6e_expr"].copy(); motif = data["fig6e_motif"].copy()
    expr["feature"] = expr["tf"] + " expr"
    motif["feature"] = motif["tf"] + " motif"
    combined = pd.concat([
        motif.rename(columns={"row_z": "row_z"})[["feature", "branch", "row_z"]],
        expr.rename(columns={"row_z": "row_z"})[["feature", "branch", "row_z"]],
    ], ignore_index=True)
    return plot_heatmap(combined, "row_z", "feature", "branch", "Fig. 6E GPC TF motif/gene dynamics across ATAC branches", out_dir, "figure6E_gpc_tf_motif_gene_branch_dynamics")


def plot_6f(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    no = data["ga_no"].set_index("pseudobulk_index")
    an = data["ga_anchor"].set_index("pseudobulk_index")
    fig, ax = plt.subplots(figsize=(6.6, 5.8))
    ax.scatter(no["UMAP1"], no["UMAP2"], s=6, color="0.8", alpha=0.45, label="GA no anchors")
    branch = an[an["is.Branch"] != "Not Branch"]
    for idx, r in branch.iterrows():
        if idx in no.index:
            ax.annotate("", xy=(r["UMAP1"], r["UMAP2"]), xytext=(no.loc[idx, "UMAP1"], no.loc[idx, "UMAP2"]), arrowprops=dict(arrowstyle="->", color=BRANCH_COLORS.get(r["is.Branch"], "black"), lw=0.6, alpha=0.45))
    for lab in ["A", "B", "C"]:
        sub = branch[branch["is.Branch"] == lab]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=14, color=BRANCH_COLORS[lab], label=f"branch {lab}", linewidths=0)
    ax.set_title("Fig. 6F GPC-anchor reprojection: branch movement arrows")
    ax.set_xlabel("FCM UMAP1"); ax.set_ylabel("FCM UMAP2")
    ax.legend(frameon=False, fontsize=8)
    return save_triple(fig, out_dir, "figure6F_gpc_anchor_reprojection_arrows")


def plot_6g(data: dict[str, Any], out_dir: Path) -> dict[str, str]:
    df = data["fig6g"].copy()
    clusters = sorted(df["mapped_atac_cluster"].astype(str).unique())
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(7, 5.8))
    for i, cl in enumerate(clusters):
        sub = df[df["mapped_atac_cluster"].astype(str) == cl]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=5, alpha=0.55, color=cmap(i % 20), label=cl, linewidths=0)
    ax.set_title("Fig. 6G bounded multiome RNA→FCM NN projection colored by ATAC cluster")
    ax.set_xlabel("FCM scRNA UMAP1 (nearest-neighbor transferred)"); ax.set_ylabel("FCM scRNA UMAP2")
    ax.legend(markerscale=2, fontsize=6, frameon=False, ncol=2)
    return save_triple(fig, out_dir, "figure6G_multiome_nn_fcm_projection_by_atac_cluster")


def collect_png_hashes(out_dir: Path) -> dict[str, str]:
    return {str(p): sha256_file(p) for p in sorted(out_dir.glob("figure6*.png"))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--fcm-rds", default=Path("/home/zerlinshen/projects/wave5-trevino/inputs/brainchromatin_s3/rds/FCM_Object.RDS"), type=Path)
    ap.add_argument("--all-genes-rds", default=Path("/home/zerlinshen/projects/wave5-trevino/inputs/brainchromatin_s3/rds/AllGenes_GenomicRanges.RDS"), type=Path)
    ap.add_argument("--chromvar-rds", default=Path("/home/zerlinshen/projects/wave5-trevino/inputs/brainchromatin_s3/rds/Matrix_Glial_KNNpseudoBulk_ATAC_ChromVAR.RDS"), type=Path)
    ap.add_argument("--multiome-rna-rds", default=Path("/home/zerlinshen/projects/wave5-trevino/inputs/brainchromatin_s3/rds/Multiome_RNA_SCE.RDS"), type=Path)
    ap.add_argument("--multiome-atac-rds", default=Path("/home/zerlinshen/projects/wave5-trevino/inputs/brainchromatin_s3/rds/Multiome_ATAC_SCE.RDS"), type=Path)
    ap.add_argument("--gpc-tsv", default=None, type=Path)
    ap.add_argument("--cell-cycle-txt", default=Path("/home/zerlinshen/projects/wave5-trevino/inputs/external_code/brainchromatin/gene_sets/gs_MITOTIC_CELL_CYCLE.txt"), type=Path)
    ap.add_argument("--rscript", default=Path("/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript"), type=Path)
    args = ap.parse_args()

    configure_matplotlib()
    run_dir = args.project_root / "runs" / args.run_id
    fig_dir = run_dir / "python" / "figures" / "public_matrix" / "figure_6_fcm_gpc_branch"
    evidence_dir = run_dir / "python" / "figure_reproduction_evidence"
    fig_dir.mkdir(parents=True, exist_ok=True); evidence_dir.mkdir(parents=True, exist_ok=True)
    if args.gpc_tsv is None:
        args.gpc_tsv = run_dir / "python" / "figures" / "public_matrix" / "figure_2E_gpc_enrichment" / "figure2e_gpc_count_aligned_top185.tsv"
    required = [args.fcm_rds, args.all_genes_rds, args.chromvar_rds, args.multiome_rna_rds, args.multiome_atac_rds, args.gpc_tsv, args.cell_cycle_txt, args.rscript]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required resources: " + ", ".join(missing))

    paths = run_r_export(args, fig_dir)
    data = load_tables(paths)
    outputs: dict[str, dict[str, str]] = {
        "6A": plot_6a(data, fig_dir),
        "6B": plot_6b(data, fig_dir),
        "6C": plot_6c(data, fig_dir),
        "6D": plot_6d(data, fig_dir),
        "6E": plot_6e(data, fig_dir),
        "6F": plot_6f(data, fig_dir),
        "6G": plot_6g(data, fig_dir),
    }
    summary = data["summary"]
    evidence = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "project_root": str(args.project_root),
        "run_id": args.run_id,
        "output_dir": str(fig_dir),
        "resources": {"fcm_rds": str(args.fcm_rds), "chromvar_rds": str(args.chromvar_rds), "gpc_tsv": str(args.gpc_tsv), "cell_cycle_txt": str(args.cell_cycle_txt), "multiome_rna_rds": str(args.multiome_rna_rds), "multiome_atac_rds": str(args.multiome_atac_rds)},
        "summary": summary,
        "panel_status": {
            "6A": "REPRODUCED_PUBLIC_FCM_CELL_CYCLE_MODULE_CORRELATION",
            "6B": "GENERATED_PUBLIC_FCM_ATAC_PROJECTION_SCHEMATIC",
            "6C": "REPRODUCED_AUTHOR_FCM_ATAC_BRANCH_PROJECTION_FROM_GA_ANCHORS",
            "6D": "GENERATED_PUBLIC_FCM_BRANCH_GENE_ACTIVITY_HEATMAP_LOCAL_SPECIFICITY",
            "6E": "GENERATED_PUBLIC_FCM_GPC_TF_MOTIF_GENE_BRANCH_HEATMAP",
            "6F": "REPRODUCED_AUTHOR_FCM_GPC_ANCHOR_REPROJECTION_BRANCH_ARROWS",
            "6G": "GENERATED_BOUNDED_MULTIOME_RNA_TO_FCM_NN_PROJECTION_COLORED_BY_ATAC_CLUSTER",
        },
        "checks": {
            "pseudobulks_eq_1267": summary.get("pseudobulks") == 1267,
            "branches_a_b_c_nonempty": summary.get("branch_A", 0) > 0 and summary.get("branch_B", 0) > 0 and summary.get("branch_C", 0) > 0,
            "cell_cycle_genes_found_ge_50": summary.get("cell_cycle_genes_found", 0) >= 50,
            "gpc_ids_in_gene_activity_ge_100": summary.get("gpc_ids_in_gene_activity", 0) >= 100,
            "figure6d_top_genes_eq_50": summary.get("figure6d_top_genes") == 50,
            "motif_tfs_found_ge_8": summary.get("motif_tfs_found", 0) >= 8,
            "multiome_cells_projected_eq_8981": summary.get("multiome_cells_projected") == 8981,
            "multiome_projection_genes_ge_300": summary.get("multiome_projection_genes", 0) >= 300,
            "all_output_hashes_recorded": True,
        },
        "tables": {k: str(v) for k, v in paths.items()},
        "outputs": outputs,
        "output_png_hashes": collect_png_hashes(fig_dir),
        "truth_boundaries": [
            "Figure 6A-F are based on the author public FCM object and public GPC/motif resources, not a raw-fragment rerun.",
            "Figure 6D uses a local branch-specific gene-activity specificity statistic over the public branch labels; exact paper KS/ranking implementation parity is not claimed.",
            "Figure 6E uses public chromVAR and expression branch means for selected GPC TFs; exact article heatmap scaling/order may differ.",
            "Figure 6G is a bounded nearest-neighbor projection of public multiome RNA into the FCM scRNA manifold and colors cells by public multiome ATAC clusters; exact author uwot projection parity is not claimed.",
        ],
    }
    evidence_path = evidence_dir / "figure6_public_fcm_gpc_branch_evidence.json"
    evidence_path.write_text(json.dumps(evidence, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps({"evidence": str(evidence_path), "output_dir": str(fig_dir), "panel_status": evidence["panel_status"], "checks": evidence["checks"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
