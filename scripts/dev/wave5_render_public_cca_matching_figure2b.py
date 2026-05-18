#!/usr/bin/env python3
"""Render Trevino Figure 2B public CCA-matching panels.

Figure 2B is described as UMAPs of scRNA-seq and scATAC-seq cells colored by
cluster assignment of matched cells.  This script uses the author's public
``CCA_Matching.RDS`` resource for the cross-modality nearest-neighbor links and
the current public-matrix RNA/ATAC embeddings available in the Wave-5 run.

Boundary: this is a public processed-matrix reconstruction.  The matching table
is author-provided, but the RNA UMAP comes from the current pipeline output and
the ATAC UMAP is the bounded public gene-activity embedding generated in the
Figure 3B lane, not the exact author peak/fragments UMAP.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


RNA_TO_ATAC = "nearest_acc_for_expr"
ATAC_TO_RNA = "nearest_expr_for_acc"


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


def export_rds_to_tsv_gz(rds: Path, out_tsv_gz: Path, rscript_bin: str) -> dict[str, str | int]:
    """Export a data.frame RDS to a gzipped TSV using Rscript."""
    out_tsv_gz.parent.mkdir(parents=True, exist_ok=True)
    code = f"""
    x <- readRDS({json.dumps(str(rds))})
    stopifnot(is.data.frame(x))
    req <- c("cell_expr", "cell_acc", "dist", "sampleId", "mapping")
    missing <- setdiff(req, colnames(x))
    if (length(missing) > 0) stop(paste("missing columns:", paste(missing, collapse=",")))
    con <- gzfile({json.dumps(str(out_tsv_gz))}, "wt")
    write.table(x[, req], con, sep="\\t", quote=FALSE, row.names=FALSE)
    close(con)
    """
    proc = subprocess.run(
        [rscript_bin, "-e", code],
        check=True,
        text=True,
        capture_output=True,
    )
    return {
        "rscript_bin": rscript_bin,
        "rds": str(rds),
        "exported_tsv_gz": str(out_tsv_gz),
        "stderr_chars": len(proc.stderr),
        "stdout_chars": len(proc.stdout),
    }


def read_cca_matching(rds: Path, out_dir: Path, rscript_bin: str, force: bool = False) -> tuple[pd.DataFrame, dict]:
    exported = out_dir / "cca_matching_author_public.tsv.gz"
    export_info: dict
    if force or not exported.exists() or exported.stat().st_size == 0:
        export_info = export_rds_to_tsv_gz(rds, exported, rscript_bin)
    else:
        export_info = {
            "rscript_bin": rscript_bin,
            "rds": str(rds),
            "exported_tsv_gz": str(exported),
            "reused_existing_export": True,
        }
    cca = pd.read_csv(exported, sep="\t", compression="gzip")
    expected = {"cell_expr", "cell_acc", "dist", "sampleId", "mapping"}
    missing = expected.difference(cca.columns)
    if missing:
        raise ValueError(f"CCA export missing columns: {sorted(missing)}")
    return cca, export_info


def ordered_categories(values: pd.Series, preferred: list[str] | None = None) -> list[str]:
    counts = values.dropna().astype(str).value_counts()
    seen = counts.index.tolist()
    if preferred:
        return [v for v in preferred if v in seen] + [v for v in seen if v not in preferred]
    return seen


def scatter_by_category(
    ax: plt.Axes,
    coords: np.ndarray,
    labels: pd.Series,
    title: str,
    size: float,
    preferred: list[str] | None = None,
) -> None:
    label_values = labels.astype("object").where(labels.notna(), "unmatched").astype(str)
    cats = ordered_categories(label_values, preferred)
    cmap = plt.get_cmap("tab20")
    for i, cat in enumerate(cats):
        mask = label_values.to_numpy() == cat
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=size,
            alpha=0.78 if cat != "unmatched" else 0.25,
            linewidths=0,
            color="0.75" if cat == "unmatched" else cmap(i % 20),
            label=cat,
        )
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=6, markerscale=4)


def row_normalized_crosstab(rows: pd.Series, cols: pd.Series) -> pd.DataFrame:
    tab = pd.crosstab(rows.astype(str), cols.astype(str))
    denom = tab.sum(axis=1).replace(0, np.nan)
    return tab.div(denom, axis=0).fillna(0.0)


def plot_heatmap(ax: plt.Axes, matrix: pd.DataFrame, title: str) -> None:
    im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", vmin=0, vmax=max(0.2, float(matrix.to_numpy().max())))
    ax.set_title(title)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=6)
    ax.set_xlabel("matched-cell cluster")
    ax.set_ylabel("own cluster")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="row fraction")


def majority_map(source_cluster: pd.Series, target_label: pd.Series) -> dict[str, str]:
    df = pd.DataFrame({"cluster": source_cluster.astype(str), "label": target_label.astype(str)})
    out: dict[str, str] = {}
    for cluster, sub in df.groupby("cluster", observed=True):
        counts = sub["label"].value_counts()
        if not counts.empty:
            out[str(cluster)] = str(counts.index[0])
    return out


def finite_umap(coords: np.ndarray) -> bool:
    return bool(coords.ndim == 2 and coords.shape[1] == 2 and np.isfinite(coords).all())


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rna-h5ad", required=True, type=Path)
    p.add_argument("--rna-metadata", required=True, type=Path)
    p.add_argument("--atac-umap-h5ad", required=True, type=Path)
    p.add_argument("--atac-metadata", required=True, type=Path)
    p.add_argument("--cca-matching-rds", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--rscript-bin", default="Rscript")
    p.add_argument("--force-export", action="store_true")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    cca, export_info = read_cca_matching(args.cca_matching_rds, args.out_dir, args.rscript_bin, args.force_export)
    mapping_counts = cca["mapping"].astype(str).value_counts().to_dict()
    rna_to_atac = (
        cca.loc[cca["mapping"].astype(str) == RNA_TO_ATAC, ["cell_expr", "cell_acc", "dist", "sampleId"]]
        .drop_duplicates("cell_expr")
        .set_index("cell_expr")
    )
    atac_to_rna = (
        cca.loc[cca["mapping"].astype(str) == ATAC_TO_RNA, ["cell_expr", "cell_acc", "dist", "sampleId"]]
        .drop_duplicates("cell_acc")
        .set_index("cell_acc")
    )

    rna = ad.read_h5ad(args.rna_h5ad, backed="r")
    rna_obs = rna.obs.copy()
    rna_umap = np.asarray(rna.obsm["X_umap"])
    rna.file.close()

    atac = ad.read_h5ad(args.atac_umap_h5ad, backed="r")
    atac_obs = atac.obs.copy()
    atac_umap = np.asarray(atac.obsm["X_umap"])
    atac.file.close()

    rna_meta = pd.read_csv(args.rna_metadata, sep="\t", compression="gzip", dtype=str).set_index("Cell.ID")
    atac_meta = pd.read_csv(args.atac_metadata, sep="\t", compression="gzip", dtype=str).set_index("Cell.ID")

    # RNA cells colored by the cluster of the matched ATAC cell.
    rna_labels = rna_obs.copy()
    rna_labels["cell_expr"] = rna_labels.index.astype(str)
    rna_labels = rna_labels.join(rna_to_atac.rename(columns={"dist": "cca_dist_rna_to_atac"}), how="left")
    rna_labels["matched_atac_cluster"] = rna_labels["cell_acc"].map(atac_meta["Iterative.LSI.Clusters"])
    rna_labels["matched_atac_age"] = rna_labels["cell_acc"].map(atac_meta["Age"])
    rna_labels["matched_atac_sample"] = rna_labels["cell_acc"].map(atac_meta["Sample.ID"])
    rna_labels["own_rna_seurat_cluster"] = rna_labels.get("seurat_clusters", pd.Series(index=rna_labels.index, dtype=object)).astype(str)
    rna_labels["own_rna_cell_type"] = rna_labels.get("cell_type", pd.Series(index=rna_labels.index, dtype=object)).astype(str)

    # ATAC cells colored by the cluster of the matched RNA cell.
    atac_labels = atac_obs.copy()
    atac_labels["cell_acc"] = atac_labels.index.astype(str)
    atac_labels = atac_labels.join(atac_to_rna.rename(columns={"dist": "cca_dist_atac_to_rna"}), how="left")
    atac_labels["matched_rna_seurat_cluster"] = atac_labels["cell_expr"].map(rna_meta["seurat_clusters"])
    atac_labels["matched_rna_age"] = atac_labels["cell_expr"].map(rna_meta["Age"])
    atac_labels["matched_rna_sample"] = atac_labels["cell_expr"].map(rna_meta["Sample.ID"])
    if "cell_type" in rna_obs.columns:
        atac_labels["matched_pipeline_rna_cell_type"] = atac_labels["cell_expr"].map(rna_obs["cell_type"].astype(str))
    atac_labels["own_atac_lsi_cluster"] = atac_labels.get(
        "Iterative.LSI.Clusters", pd.Series(index=atac_labels.index, dtype=object)
    ).astype(str)

    # Cell-type proxy agreement: assign each ATAC cluster its majority matched RNA
    # cell type, then compare RNA cells' own pipeline type with the matched ATAC
    # cluster's majority type.  This is a QC proxy, not an author metric.
    if "matched_pipeline_rna_cell_type" in atac_labels:
        atac_cluster_to_rna_type = majority_map(
            atac_labels["own_atac_lsi_cluster"],
            atac_labels["matched_pipeline_rna_cell_type"].fillna("unmapped"),
        )
        rna_labels["matched_atac_cluster_majority_rna_type"] = rna_labels["matched_atac_cluster"].map(atac_cluster_to_rna_type)
        valid_type = (
            rna_labels["own_rna_cell_type"].notna()
            & rna_labels["matched_atac_cluster_majority_rna_type"].notna()
            & (rna_labels["own_rna_cell_type"] != "nan")
        )
        celltype_proxy_agreement = (
            float(
                (
                    rna_labels.loc[valid_type, "own_rna_cell_type"].astype(str)
                    == rna_labels.loc[valid_type, "matched_atac_cluster_majority_rna_type"].astype(str)
                ).mean()
            )
            if int(valid_type.sum()) > 0
            else None
        )
    else:
        atac_cluster_to_rna_type = {}
        valid_type = pd.Series(False, index=rna_labels.index)
        celltype_proxy_agreement = None

    figure_outputs: dict[str, dict[str, str]] = {}
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), constrained_layout=True)
    scatter_by_category(
        axes[0],
        rna_umap,
        rna_labels["matched_atac_cluster"],
        "Figure 2B RNA UMAP colored by matched ATAC cluster",
        size=0.8,
    )
    scatter_by_category(
        axes[1],
        atac_umap,
        atac_labels["matched_rna_seurat_cluster"],
        "Figure 2B ATAC UMAP colored by matched RNA cluster",
        size=2.2,
    )
    figure_outputs["figure2b_matched_cluster_umaps"] = save_triple(fig, args.out_dir, "figure2b_matched_cluster_umaps")

    rna_heat = row_normalized_crosstab(rna_labels["own_rna_seurat_cluster"], rna_labels["matched_atac_cluster"])
    atac_heat = row_normalized_crosstab(atac_labels["own_atac_lsi_cluster"], atac_labels["matched_rna_seurat_cluster"])
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    plot_heatmap(axes[0], rna_heat, "RNA own cluster → matched ATAC cluster")
    plot_heatmap(axes[1], atac_heat, "ATAC own cluster → matched RNA cluster")
    figure_outputs["figure2b_match_concordance_heatmaps"] = save_triple(
        fig, args.out_dir, "figure2b_match_concordance_heatmaps"
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    axes[0].hist(rna_labels["cca_dist_rna_to_atac"].dropna().astype(float), bins=60, color="#4C78A8")
    axes[0].set_title("RNA→ATAC CCA nearest-neighbor distance")
    axes[0].set_xlabel("CCA distance")
    axes[0].set_ylabel("cells")
    axes[1].hist(atac_labels["cca_dist_atac_to_rna"].dropna().astype(float), bins=60, color="#F58518")
    axes[1].set_title("ATAC→RNA CCA nearest-neighbor distance")
    axes[1].set_xlabel("CCA distance")
    axes[1].set_ylabel("cells")
    figure_outputs["figure2b_cca_distance_qc"] = save_triple(fig, args.out_dir, "figure2b_cca_distance_qc")

    rna_label_path = args.out_dir / "figure2b_rna_cells_with_matched_atac_labels.tsv.gz"
    atac_label_path = args.out_dir / "figure2b_atac_cells_with_matched_rna_labels.tsv.gz"
    with gzip.open(rna_label_path, "wt") as fh:
        rna_labels.to_csv(fh, sep="\t", index=True)
    with gzip.open(atac_label_path, "wt") as fh:
        atac_labels.to_csv(fh, sep="\t", index=True)

    rna_dist = rna_labels["cca_dist_rna_to_atac"].dropna().astype(float)
    atac_dist = atac_labels["cca_dist_atac_to_rna"].dropna().astype(float)
    metrics = {
        "cca_rows": int(cca.shape[0]),
        "cca_mapping_counts": {str(k): int(v) for k, v in mapping_counts.items()},
        "rna_cells_plotted": int(rna_labels.shape[0]),
        "rna_cells_with_author_match": int(rna_labels["cell_acc"].notna().sum()),
        "rna_cells_with_matched_atac_cluster": int(rna_labels["matched_atac_cluster"].notna().sum()),
        "rna_matched_atac_cluster_coverage": float(rna_labels["matched_atac_cluster"].notna().mean()),
        "rna_matched_atac_clusters": rna_labels["matched_atac_cluster"].dropna().astype(str).value_counts().to_dict(),
        "atac_cells_plotted": int(atac_labels.shape[0]),
        "atac_cells_with_author_match": int(atac_labels["cell_expr"].notna().sum()),
        "atac_cells_with_matched_rna_cluster": int(atac_labels["matched_rna_seurat_cluster"].notna().sum()),
        "atac_matched_rna_cluster_coverage": float(atac_labels["matched_rna_seurat_cluster"].notna().mean()),
        "atac_matched_rna_clusters": atac_labels["matched_rna_seurat_cluster"].dropna().astype(str).value_counts().to_dict(),
        "rna_to_atac_distance": {
            "min": float(rna_dist.min()) if not rna_dist.empty else None,
            "median": float(rna_dist.median()) if not rna_dist.empty else None,
            "mean": float(rna_dist.mean()) if not rna_dist.empty else None,
            "max": float(rna_dist.max()) if not rna_dist.empty else None,
        },
        "atac_to_rna_distance": {
            "min": float(atac_dist.min()) if not atac_dist.empty else None,
            "median": float(atac_dist.median()) if not atac_dist.empty else None,
            "mean": float(atac_dist.mean()) if not atac_dist.empty else None,
            "max": float(atac_dist.max()) if not atac_dist.empty else None,
        },
        "rna_atac_row_normalized_heatmap_shape": [int(rna_heat.shape[0]), int(rna_heat.shape[1])],
        "atac_rna_row_normalized_heatmap_shape": [int(atac_heat.shape[0]), int(atac_heat.shape[1])],
        "atac_cluster_to_majority_matched_pipeline_rna_cell_type": atac_cluster_to_rna_type,
        "rna_celltype_proxy_agreement_with_matched_atac_cluster_majority_type": celltype_proxy_agreement,
        "rna_celltype_proxy_agreement_cells": int(valid_type.sum()),
    }

    outputs = {
        "out_dir": str(args.out_dir),
        "cca_matching_export": str(args.out_dir / "cca_matching_author_public.tsv.gz"),
        "rna_cell_label_table": str(rna_label_path),
        "atac_cell_label_table": str(atac_label_path),
        "figure_outputs": figure_outputs,
    }
    checks = {
        "cca_has_expected_bidirectional_rows": mapping_counts.get(RNA_TO_ATAC, 0) >= 57000
        and mapping_counts.get(ATAC_TO_RNA, 0) >= 31000,
        "rna_umap_present_finite": finite_umap(rna_umap),
        "atac_umap_present_finite": finite_umap(atac_umap),
        "rna_cells_ge_50000": int(rna_labels.shape[0]) >= 50000,
        "atac_cells_ge_10000": int(atac_labels.shape[0]) >= 10000,
        "rna_match_cluster_coverage_ge_95pct": metrics["rna_matched_atac_cluster_coverage"] >= 0.95,
        "atac_match_cluster_coverage_ge_95pct": metrics["atac_matched_rna_cluster_coverage"] >= 0.95,
        "rna_matched_atac_clusters_ge_10": len(metrics["rna_matched_atac_clusters"]) >= 10,
        "atac_matched_rna_clusters_ge_10": len(metrics["atac_matched_rna_clusters"]) >= 10,
        "distance_values_finite": bool(np.isfinite(rna_dist).all() and np.isfinite(atac_dist).all()),
        "all_output_hashes_recorded": all(
            "png_sha256" in entry and "svg_sha256" in entry and "pdf_sha256" in entry
            for entry in figure_outputs.values()
        ),
    }
    status = "REPRODUCED_AUTHOR_CCA_PUBLIC_MATRIX_FIGURE2B" if all(checks.values()) else "GENERATED_WITH_VALIDATION_WARNINGS"

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "panel": "Figure 2B",
        "status": status,
        "truth_boundary": (
            "Author public CCA_Matching.RDS was used for RNA↔ATAC nearest-neighbor links. "
            "RNA UMAP is the current singlecell_factory public RNA pipeline embedding; ATAC UMAP is the bounded "
            "public gene-activity embedding produced for Figure 3B. This is not FASTQ/fragments-level reprocessing "
            "and not exact author peak-LSI UMAP parity."
        ),
        "inputs": {
            "rna_h5ad": str(args.rna_h5ad),
            "rna_metadata": str(args.rna_metadata),
            "atac_umap_h5ad": str(args.atac_umap_h5ad),
            "atac_metadata": str(args.atac_metadata),
            "cca_matching_rds": str(args.cca_matching_rds),
        },
        "export_info": export_info,
        "outputs": outputs,
        "metrics": metrics,
        "checks": checks,
    }
    summary_path = args.out_dir / "figure2b_public_cca_matching_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
