#!/usr/bin/env python3
"""Render Trevino Figure 2C/2D public resource panels.

Figure 2C in the paper is a paired CRE-accessibility / RNA-expression heatmap
for significant CRE-gene links.  The public author resources available in this
run include the significant link table and matched pseudobulk RNA / ATAC gene
activity matrices, but the large ATAC accessibility pseudobulk matrix available
locally has no peak names, so exact CRE-row reconstruction is not claimed here.

This script therefore produces:

* Figure 2C proxy evidence: paired ATAC gene-activity and RNA pseudobulk
  heatmaps for genes with significant linked CREs.
* Figure 2D evidence: author GA-RNA correlation scatter joined to significant
  linked-CRE counts per gene, with TFs labeled.

All outputs are explicitly marked as author-resource/public-matrix evidence, not
FASTQ/fragments-level or exact 2C CRE-pair parity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.cluster import KMeans  # noqa: E402


EXPECTED_FIG2C_LINKS = 64878


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


def export_author_resource_tables(
    *,
    significant_links_rds: Path,
    ga_rna_correlations_rds: Path,
    ga_pseudobulk_rds: Path,
    rna_pseudobulk_rds: Path,
    out_dir: Path,
    rscript_bin: str,
    max_heatmap_genes: int,
    force: bool,
) -> dict[str, str | int | bool]:
    out_dir.mkdir(parents=True, exist_ok=True)
    sentinel = out_dir / "figure2cd_export_complete.json"
    if sentinel.exists() and not force:
        return json.loads(sentinel.read_text()) | {"reused_existing_export": True}

    code = f"""
    links <- readRDS({json.dumps(str(significant_links_rds))})
    cor.raw <- readRDS({json.dumps(str(ga_rna_correlations_rds))})
    ga <- readRDS({json.dumps(str(ga_pseudobulk_rds))})
    rna <- readRDS({json.dumps(str(rna_pseudobulk_rds))})
    req.links <- c("peak.name", "geneSymbol")
    if (any(!req.links %in% colnames(links))) stop("links missing required columns")
    req.cor <- c("gene.symbol", "gene.id", "spearman.cor", "ranked.genes", "n.cells")
    if (any(!req.cor %in% colnames(cor.raw))) stop("GA/RNA correlations missing required columns")
    cor.raw$n.cells.numeric <- suppressWarnings(as.numeric(as.character(cor.raw$n.cells)))
    cor <- cor.raw[order(cor.raw$gene.symbol, -cor.raw$n.cells.numeric, -abs(cor.raw$spearman.cor)), ]
    cor <- cor[!duplicated(cor$gene.symbol), ]

    link_counts <- aggregate(
      links$peak.name,
      by = list(gene.symbol = links$geneSymbol),
      FUN = function(x) length(unique(x))
    )
    colnames(link_counts)[2] <- "linked_cre_count"
    merged <- merge(cor, link_counts, by = "gene.symbol", all.x = TRUE, sort = FALSE)
    merged$linked_cre_count[is.na(merged$linked_cre_count)] <- 0
    merged$has_significant_link <- merged$linked_cre_count > 0

    shared <- merged[
      merged$gene.symbol %in% rownames(ga) &
      merged$gene.id %in% rownames(rna) &
      merged$linked_cre_count > 0,
    ]
    shared <- shared[order(-shared$linked_cre_count, -abs(shared$spearman.cor), shared$gene.symbol), ]
    shared <- shared[!duplicated(shared$gene.symbol), ]
    selected <- head(shared, {int(max_heatmap_genes)})
    if (nrow(selected) < 10) stop("too few selected shared heatmap genes")

    ga.sel <- ga[selected$gene.symbol, , drop = FALSE]
    rna.sel <- rna[selected$gene.id, , drop = FALSE]
    rownames(rna.sel) <- selected$gene.symbol
    colnames(ga.sel) <- sprintf("PB_%04d", seq_len(ncol(ga.sel)))
    colnames(rna.sel) <- sprintf("PB_%04d", seq_len(ncol(rna.sel)))

    write.table(link_counts, {json.dumps(str(out_dir / "figure2d_linked_cre_counts.tsv"))}, sep="\\t", quote=FALSE, row.names=FALSE)
    write.table(merged, {json.dumps(str(out_dir / "figure2d_ga_rna_correlation_with_link_counts.tsv"))}, sep="\\t", quote=FALSE, row.names=FALSE)
    write.table(selected, {json.dumps(str(out_dir / "figure2c_selected_linked_genes.tsv"))}, sep="\\t", quote=FALSE, row.names=FALSE)
    write.table(ga.sel, gzfile({json.dumps(str(out_dir / "figure2c_selected_atac_gene_activity.tsv.gz"))}, "wt"), sep="\\t", quote=FALSE)
    write.table(rna.sel, gzfile({json.dumps(str(out_dir / "figure2c_selected_rna_expression.tsv.gz"))}, "wt"), sep="\\t", quote=FALSE)

    info <- list(
      significant_links_rows = nrow(links),
      significant_unique_peaks = length(unique(links$peak.name)),
      significant_unique_genes = length(unique(links$geneSymbol)),
      ga_rna_correlation_rows = nrow(cor.raw),
      ga_rna_correlation_unique_genes = nrow(cor),
      ga_rna_correlation_ncells_selected = max(cor$n.cells.numeric, na.rm = TRUE),
      genes_with_link_counts = nrow(link_counts),
      shared_heatmap_candidates = nrow(shared),
      selected_heatmap_genes = nrow(selected),
      pseudobulk_columns = ncol(ga.sel),
      ga_pseudobulk_rows = nrow(ga),
      rna_pseudobulk_rows = nrow(rna)
    )
    cat(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE))
    """
    proc = subprocess.run([rscript_bin, "-e", code], check=True, text=True, capture_output=True)
    info = json.loads(proc.stdout)
    export_info = {
        "rscript_bin": rscript_bin,
        "significant_links_rds": str(significant_links_rds),
        "ga_rna_correlations_rds": str(ga_rna_correlations_rds),
        "ga_pseudobulk_rds": str(ga_pseudobulk_rds),
        "rna_pseudobulk_rds": str(rna_pseudobulk_rds),
        "stderr_chars": len(proc.stderr),
        **{k: int(v) for k, v in info.items()},
    }
    sentinel.write_text(json.dumps(export_info, indent=2, ensure_ascii=False) + "\n")
    return export_info


def row_zscore(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float)
    mean = np.nanmean(arr, axis=1, keepdims=True)
    std = np.nanstd(arr, axis=1, keepdims=True)
    std[std == 0] = 1.0
    z = (arr - mean) / std
    z[~np.isfinite(z)] = 0.0
    return pd.DataFrame(z, index=df.index, columns=df.columns)


def read_tf_symbols(path: Path) -> set[str]:
    out: set[str] = set()
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(">") or line.startswith("GO_"):
            continue
        out.add(line)
    return out


def render_figure2c_heatmap(ga: pd.DataFrame, rna: pd.DataFrame, selected: pd.DataFrame, out_dir: Path) -> tuple[dict, dict]:
    ga_z = row_zscore(ga)
    rna_z = row_zscore(rna)
    selected = selected.drop_duplicates("gene.symbol", keep="first").copy()
    common = selected["gene.symbol"].astype(str).tolist()
    common = [g for g in common if g in ga_z.index and g in rna_z.index]
    ga_z = ga_z.loc[common]
    rna_z = rna_z.loc[common]
    selected = selected.set_index("gene.symbol").loc[common].reset_index()

    concat = np.hstack([ga_z.to_numpy(), rna_z.to_numpy()])
    n_clusters = min(20, max(2, len(common) // 10))
    km = KMeans(n_clusters=n_clusters, random_state=13, n_init=10)
    clusters = km.fit_predict(concat)
    order_df = pd.DataFrame(
        {
            "gene": common,
            "cluster": clusters,
            "linked_cre_count": selected["linked_cre_count"].to_numpy(),
            "abs_cor": np.abs(selected["spearman.cor"].astype(float).to_numpy()),
        }
    ).sort_values(["cluster", "linked_cre_count", "abs_cor", "gene"], ascending=[True, False, False, True])
    order = order_df["gene"].tolist()
    ga_z = ga_z.loc[order]
    rna_z = rna_z.loc[order]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16, max(6, min(18, 0.025 * len(order) + 4))),
        constrained_layout=True,
        sharey=True,
    )
    for ax, data, title in [
        (axes[0], ga_z, "2C proxy — ATAC gene activity pseudobulks"),
        (axes[1], rna_z, "2C proxy — RNA expression pseudobulks"),
    ]:
        im = ax.imshow(data.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-2.5, vmax=2.5)
        ax.set_title(title)
        ax.set_xlabel("pseudobulk samples (author glial KNN resource)")
        ax.set_xticks([])
    if len(order) <= 80:
        ticks = range(len(order))
    else:
        ticks = np.linspace(0, len(order) - 1, 40).astype(int)
    axes[0].set_yticks(ticks)
    axes[0].set_yticklabels([order[i] for i in ticks], fontsize=6)
    axes[0].set_ylabel("selected linked genes ordered by k-means cluster")
    fig.colorbar(im, ax=axes, shrink=0.8, label="row-scaled value")
    outputs = save_triple(fig, out_dir, "figure2c_gene_activity_rna_heatmap_proxy")

    cluster_table = order_df.set_index("gene").loc[order].reset_index()
    cluster_path = out_dir / "figure2c_selected_gene_kmeans_order.tsv"
    cluster_table.to_csv(cluster_path, sep="\t", index=False)
    metrics = {
        "selected_genes_plotted": int(len(order)),
        "pseudobulk_columns": int(ga_z.shape[1]),
        "kmeans_clusters": int(n_clusters),
        "order_table": str(cluster_path),
        "ga_shape": [int(ga_z.shape[0]), int(ga_z.shape[1])],
        "rna_shape": [int(rna_z.shape[0]), int(rna_z.shape[1])],
    }
    return outputs, metrics


def render_figure2d_scatter(corr: pd.DataFrame, tf_symbols: set[str], out_dir: Path) -> tuple[dict, dict]:
    df = corr.copy()
    df["linked_cre_count"] = pd.to_numeric(df["linked_cre_count"], errors="coerce").fillna(0)
    df["spearman.cor"] = pd.to_numeric(df["spearman.cor"], errors="coerce")
    df = df[np.isfinite(df["spearman.cor"])].copy()
    df["abs_cor_for_dedup"] = df["spearman.cor"].abs()
    df = (
        df.sort_values(["linked_cre_count", "abs_cor_for_dedup", "gene.symbol"], ascending=[False, False, True])
        .drop_duplicates("gene.symbol", keep="first")
        .copy()
    )
    df["is_tf"] = df["gene.symbol"].astype(str).isin(tf_symbols)
    df["log10_linked_cre_count_plus1"] = np.log10(df["linked_cre_count"] + 1)

    fig, ax = plt.subplots(figsize=(7.5, 5.5), constrained_layout=True)
    bg = df[~df["is_tf"]]
    tf = df[df["is_tf"]]
    ax.scatter(
        bg["spearman.cor"],
        bg["log10_linked_cre_count_plus1"],
        s=6,
        alpha=0.25,
        linewidths=0,
        color="0.55",
        label="non-TF genes",
    )
    ax.scatter(
        tf["spearman.cor"],
        tf["log10_linked_cre_count_plus1"],
        s=11,
        alpha=0.7,
        linewidths=0,
        color="#D62728",
        label="TF genes",
    )
    label_df = tf[tf["linked_cre_count"] > 0].copy()
    label_df["label_score"] = label_df["linked_cre_count"] * (0.25 + label_df["spearman.cor"].abs())
    label_df = label_df.sort_values("label_score", ascending=False).head(30)
    for _, row in label_df.iterrows():
        ax.text(
            row["spearman.cor"],
            row["log10_linked_cre_count_plus1"],
            str(row["gene.symbol"]),
            fontsize=6,
            alpha=0.85,
        )
    ax.set_xlabel("Spearman correlation: RNA expression vs ATAC gene activity")
    ax.set_ylabel("log10(significant linked CRE count + 1)")
    ax.set_title("Figure 2D — GA/RNA correlation and linked CRE count")
    ax.legend(frameon=False, loc="upper left")
    outputs = save_triple(fig, out_dir, "figure2d_ga_rna_correlation_linked_cre_scatter")
    metrics = {
        "genes_in_scatter": int(df.shape[0]),
        "tf_genes_in_scatter": int(tf.shape[0]),
        "tf_labels": label_df["gene.symbol"].astype(str).tolist(),
        "genes_with_significant_links": int((df["linked_cre_count"] > 0).sum()),
        "max_linked_cre_count": int(df["linked_cre_count"].max()),
        "spearman_cor_min": float(df["spearman.cor"].min()),
        "spearman_cor_max": float(df["spearman.cor"].max()),
    }
    return outputs, metrics


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--significant-links-rds", required=True, type=Path)
    p.add_argument("--ga-rna-correlations-rds", required=True, type=Path)
    p.add_argument("--ga-pseudobulk-rds", required=True, type=Path)
    p.add_argument("--rna-pseudobulk-rds", required=True, type=Path)
    p.add_argument("--tf-gene-set", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--rscript-bin", default="Rscript")
    p.add_argument("--max-heatmap-genes", type=int, default=400)
    p.add_argument("--force-export", action="store_true")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    export_info = export_author_resource_tables(
        significant_links_rds=args.significant_links_rds,
        ga_rna_correlations_rds=args.ga_rna_correlations_rds,
        ga_pseudobulk_rds=args.ga_pseudobulk_rds,
        rna_pseudobulk_rds=args.rna_pseudobulk_rds,
        out_dir=args.out_dir,
        rscript_bin=args.rscript_bin,
        max_heatmap_genes=args.max_heatmap_genes,
        force=args.force_export,
    )

    corr = pd.read_csv(args.out_dir / "figure2d_ga_rna_correlation_with_link_counts.tsv", sep="\t")
    selected = pd.read_csv(args.out_dir / "figure2c_selected_linked_genes.tsv", sep="\t")
    ga = pd.read_csv(args.out_dir / "figure2c_selected_atac_gene_activity.tsv.gz", sep="\t", index_col=0)
    rna = pd.read_csv(args.out_dir / "figure2c_selected_rna_expression.tsv.gz", sep="\t", index_col=0)
    tf_symbols = read_tf_symbols(args.tf_gene_set)

    figure_outputs: dict[str, dict] = {}
    figure2c_outputs, figure2c_metrics = render_figure2c_heatmap(ga, rna, selected, args.out_dir)
    figure_outputs["figure2c_gene_activity_rna_heatmap_proxy"] = figure2c_outputs
    figure2d_outputs, figure2d_metrics = render_figure2d_scatter(corr, tf_symbols, args.out_dir)
    figure_outputs["figure2d_ga_rna_correlation_linked_cre_scatter"] = figure2d_outputs

    significant_links_rows = int(export_info["significant_links_rows"])
    checks = {
        "figure2c_public_link_rows_ge_64000": significant_links_rows >= 64000,
        "figure2c_public_link_count_matches_paper_64878": significant_links_rows == EXPECTED_FIG2C_LINKS,
        "figure2c_heatmap_genes_ge_200": figure2c_metrics["selected_genes_plotted"] >= 200,
        "figure2c_ga_rna_shapes_match": figure2c_metrics["ga_shape"] == figure2c_metrics["rna_shape"],
        "figure2c_pseudobulk_columns_ge_1000": figure2c_metrics["pseudobulk_columns"] >= 1000,
        "figure2d_ga_rna_correlations_ge_19000": int(export_info["ga_rna_correlation_rows"]) >= 19000,
        "figure2d_unique_genes_ge_2700": int(export_info["ga_rna_correlation_unique_genes"]) >= 2700,
        "figure2d_genes_with_links_ge_1800": figure2d_metrics["genes_with_significant_links"] >= 1800,
        "figure2d_tf_labels_ge_20": len(figure2d_metrics["tf_labels"]) >= 20,
        "all_output_hashes_recorded": all(
            "png_sha256" in entry and "svg_sha256" in entry and "pdf_sha256" in entry
            for entry in figure_outputs.values()
        ),
    }
    figure2c_status = "GENERATED_AUTHOR_RESOURCE_GENE_ACTIVITY_RNA_HEATMAP_NOT_CRE_ACCESSIBILITY_PAIR_PARITY"
    figure2d_status = (
        "REPRODUCED_AUTHOR_RESOURCE_GA_RNA_CORRELATION_WITH_LINK_COUNT_DISCREPANCY"
        if checks["figure2d_ga_rna_correlations_ge_19000"] and checks["figure2d_genes_with_links_ge_1800"]
        else "GENERATED_AUTHOR_RESOURCE_GA_RNA_CORRELATION_WITH_VALIDATION_WARNINGS"
    )

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "panels": ["Figure 2C", "Figure 2D"],
        "status": {"2C": figure2c_status, "2D": figure2d_status},
        "truth_boundary": (
            "Figure 2D uses author GA_RNA_Correlations.RDS joined to public author PeakGeneLinks_Significant.RDS. "
            "Figure 2C is a proxy heatmap using author glial KNN pseudobulk ATAC gene activity and RNA matrices for "
            "genes with significant linked CREs. Exact Figure 2C CRE-accessibility pair rows are not claimed because "
            "the locally available accessibility pseudobulk RDS lacks peak names and therefore cannot be safely joined "
            "to peak.name without additional author row-order metadata."
        ),
        "inputs": {
            "significant_links_rds": str(args.significant_links_rds),
            "ga_rna_correlations_rds": str(args.ga_rna_correlations_rds),
            "ga_pseudobulk_rds": str(args.ga_pseudobulk_rds),
            "rna_pseudobulk_rds": str(args.rna_pseudobulk_rds),
            "tf_gene_set": str(args.tf_gene_set),
        },
        "outputs": {
            "out_dir": str(args.out_dir),
            "figure_outputs": figure_outputs,
            "exported_tables": {
                "linked_cre_counts": str(args.out_dir / "figure2d_linked_cre_counts.tsv"),
                "ga_rna_correlation_with_link_counts": str(
                    args.out_dir / "figure2d_ga_rna_correlation_with_link_counts.tsv"
                ),
                "selected_linked_genes": str(args.out_dir / "figure2c_selected_linked_genes.tsv"),
                "selected_atac_gene_activity": str(args.out_dir / "figure2c_selected_atac_gene_activity.tsv.gz"),
                "selected_rna_expression": str(args.out_dir / "figure2c_selected_rna_expression.tsv.gz"),
            },
        },
        "metrics": {
            **export_info,
            "figure2c": figure2c_metrics,
            "figure2d": figure2d_metrics,
            "paper_expected_figure2c_significant_linked_pairs": EXPECTED_FIG2C_LINKS,
            "significant_link_row_delta_vs_paper": int(significant_links_rows - EXPECTED_FIG2C_LINKS),
        },
        "checks": checks,
        "checks_all_generated_outputs_pass": all(v for k, v in checks.items() if k != "figure2c_public_link_count_matches_paper_64878"),
    }
    summary_path = args.out_dir / "figure2cd_public_resource_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
