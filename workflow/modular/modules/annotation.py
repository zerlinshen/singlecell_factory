from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from ..context import PipelineContext


DEFAULT_MARKERS = {
    "Tumor epithelial": ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1"],
    "T cell": ["CD3D", "CD3E", "TRAC", "CD4", "CD8A", "IL7R"],
    "NK cell": ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    "B cell": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Myeloid/Macro": ["LYZ", "CD68", "CST3", "FCGR3A", "AIF1"],
    "Fibroblast": ["DCN", "LUM", "COL1A1", "COL3A1", "FAP"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "Plasma cell": ["IGHG1", "IGKC", "MZB1", "XBP1"],
    "Mast cell": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
    "Dendritic cell": ["CD1C", "CLEC9A", "FCER1A", "IRF8"],
}


class AnnotationModule:
    """Optional module: marker-scoring based cell-type annotation with confidence scores."""

    name = "annotation"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "leiden" not in adata.obs:
            raise ValueError("Annotation requires clustered AnnData.")

        marker_map = ctx.cfg.markers or DEFAULT_MARKERS
        available = {
            cell_type: [gene for gene in genes if gene in adata.raw.var_names]
            for cell_type, genes in marker_map.items()
        }
        available = {cell_type: genes for cell_type, genes in available.items() if genes}
        if not available:
            raise ValueError("No valid marker genes found in dataset.")

        for cell_type, genes in available.items():
            sc.tl.score_genes(adata, genes, score_name=f"score_{cell_type}", use_raw=True)

        score_cols = [f"score_{cell_type}" for cell_type in available]
        score_mat = adata.obs[score_cols].copy()
        score_mat.columns = list(available.keys())

        # Assign cell type with confidence scoring
        max_scores = score_mat.max(axis=1)
        vals = score_mat.values
        second_scores = pd.Series(
            np.partition(vals, -2, axis=1)[:, -2], index=score_mat.index,
        )
        confidence = max_scores - second_scores

        adata.obs["cell_type"] = score_mat.idxmax(axis=1).values
        adata.obs["annotation_confidence"] = confidence.values
        # Mark low-confidence assignments as Unknown
        low_conf = max_scores < 0.1
        adata.obs.loc[low_conf, "cell_type"] = "Unknown"

        ctx.metadata["annotation_unknown_pct"] = round(
            float(low_conf.sum()) / len(low_conf) * 100, 2
        )

        adata.obs[["leiden", "cell_type", "annotation_confidence"] + score_cols].to_csv(
            ctx.table_dir / "cell_type_annotation.csv"
        )

        cluster_summary = (
            adata.obs.groupby("leiden", observed=True)["cell_type"]
            .agg(lambda x: x.value_counts().idxmax())
        ).rename("majority_cell_type")
        cluster_summary.to_csv(ctx.table_dir / "cluster_majority_cell_type.csv")

        # --- Visualizations ---
        sc.pl.umap(adata, color=["cell_type"], show=False)
        plt.savefig(ctx.figure_dir / "umap_cell_type.png", dpi=160, bbox_inches="tight")
        plt.close()

        # Confidence UMAP
        sc.pl.umap(adata, color=["annotation_confidence"], show=False, cmap="viridis")
        plt.savefig(ctx.figure_dir / "umap_annotation_confidence.png", dpi=160, bbox_inches="tight")
        plt.close()

        # Stacked bar: cell type composition per cluster
        self._plot_composition(adata, ctx)

    @staticmethod
    def _plot_composition(adata, ctx: PipelineContext) -> None:
        ct_counts = adata.obs.groupby(["leiden", "cell_type"], observed=True).size().unstack(fill_value=0)
        ct_frac = ct_counts.div(ct_counts.sum(axis=1), axis=0)
        ct_frac.plot(kind="bar", stacked=True, figsize=(12, 5), colormap="tab20")
        plt.ylabel("Fraction")
        plt.xlabel("Leiden cluster")
        plt.title("Cell type composition per cluster")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "cell_type_composition.png", dpi=160, bbox_inches="tight")
        plt.close()
