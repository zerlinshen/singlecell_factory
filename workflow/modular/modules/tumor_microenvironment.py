from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from ..context import PipelineContext
from . import score_gene_sets

# Literature-derived TME scoring signatures.
TME_SIGNATURES = {
    # Cytolytic Activity Score — Rooney et al., Cell 2015
    "CYT": ["GZMA", "PRF1"],
    # T-cell Inflamed Signature — Ayers et al., JCI 2017 (18-gene)
    "TIS": [
        "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
        "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7",
        "PDCD1LG2", "PSMB10", "STAT1", "TIGIT", "TGFB1",
    ],
    # IFN-gamma signature — Ayers et al., JCI 2017
    "IFNg": [
        "IFNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11",
        "IDO1", "PRF1", "GZMA", "HLA-DRA",
    ],
    # Immunosuppression markers
    "immunosuppression": [
        "IDO1", "CD274", "PDCD1LG2", "HAVCR2", "LAG3",
        "TGFB1", "IL10", "VEGFA", "ARG1",
    ],
    # ESTIMATE-like immune signature — Yoshihara et al., Nat Commun 2013
    "immune_estimate": [
        "CD2", "CD3D", "CD3E", "CD3G", "CD8A", "CD19", "CD79A",
        "MS4A1", "LCK", "GZMB", "NKG7", "PRF1", "GNLY",
    ],
    # ESTIMATE-like stromal signature — Yoshihara et al., Nat Commun 2013
    "stromal_estimate": [
        "DCN", "LUM", "COL1A1", "COL3A1", "COL5A1", "FAP",
        "ACTA2", "VIM", "FN1", "SPARC", "POSTN", "THY1",
    ],
}

# Key immune checkpoint molecules for profiling
CHECKPOINT_GENES = {
    "PD-1": "PDCD1",
    "PD-L1": "CD274",
    "PD-L2": "PDCD1LG2",
    "CTLA-4": "CTLA4",
    "LAG-3": "LAG3",
    "TIM-3": "HAVCR2",
    "TIGIT": "TIGIT",
    "VISTA": "VSIR",
    "IDO1": "IDO1",
    "B7-H3": "CD276",
}


class TumorMicroenvironmentModule:
    """Optional module: tumor microenvironment (TME) scoring and checkpoint profiling.

    Computes established immunotherapy-relevant signatures:
    - CYT score (Rooney et al., Cell 2015) — cytolytic activity
    - TIS (Ayers et al., JCI 2017) — T-cell inflammation, predicts anti-PD-1 response
    - IFN-gamma signature — interferon gamma response
    - Immune/Stromal ESTIMATE scores (Yoshihara et al., Nat Commun 2013)
    - Checkpoint molecule expression profiling per cell type

    Requires cell type annotations from the annotation module.
    """

    name = "tumor_microenvironment"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "cell_type" not in adata.obs.columns:
            raise ValueError("TME scoring requires cell type annotations.")

        var_names = set(adata.var_names if adata.raw is None else adata.raw.var_names)

        # --- Score TME signatures ---
        scored_sigs = score_gene_sets(adata, TME_SIGNATURES, "tme")

        # CYT score: geometric mean of GZMA and PRF1 on raw expression
        # (per Rooney et al. definition, distinct from additive gene set score)
        self._compute_cyt_score(adata, var_names)

        if not scored_sigs:
            raise ValueError("Insufficient TME signature genes found in dataset.")

        # --- Checkpoint expression profiling ---
        checkpoint_data = self._profile_checkpoints(adata, var_names)

        # --- Save tables ---
        tme_cols = [c for c in adata.obs.columns if c.startswith("tme_")]
        if "cyt_score" in adata.obs:
            tme_cols = ["cyt_score"] + tme_cols
        adata.obs[tme_cols].to_csv(ctx.table_dir / "tme_scores_per_cell.csv")

        if "leiden" in adata.obs:
            cluster_scores = adata.obs[["leiden"] + tme_cols].groupby("leiden", observed=True).mean()
            cluster_scores.to_csv(ctx.table_dir / "tme_scores_per_cluster.csv")

        if checkpoint_data is not None:
            checkpoint_data.to_csv(ctx.table_dir / "checkpoint_expression.csv")

        # Metadata
        ctx.metadata["tme_signatures_scored"] = scored_sigs
        if "cyt_score" in adata.obs:
            ctx.metadata["tme_mean_cyt"] = round(float(adata.obs["cyt_score"].mean()), 4)

        # --- Visualizations ---
        self._plot_visualizations(adata, ctx, scored_sigs, checkpoint_data)

    @staticmethod
    def _compute_cyt_score(adata, var_names: set) -> None:
        """Compute CYT score as geometric mean of GZMA and PRF1."""
        if "GZMA" not in var_names or "PRF1" not in var_names:
            return
        expr = adata.raw.to_adata() if adata.raw else adata
        gzma_idx = list(expr.var_names).index("GZMA")
        prf1_idx = list(expr.var_names).index("PRF1")
        gzma = expr.X[:, gzma_idx]
        prf1 = expr.X[:, prf1_idx]
        if hasattr(gzma, "toarray"):
            gzma = gzma.toarray().flatten()
            prf1 = prf1.toarray().flatten()
        else:
            gzma = np.asarray(gzma).flatten()
            prf1 = np.asarray(prf1).flatten()
        adata.obs["cyt_score"] = np.sqrt(np.maximum(gzma, 0) * np.maximum(prf1, 0))

    @staticmethod
    def _profile_checkpoints(adata, var_names: set) -> pd.DataFrame | None:
        """Compute mean checkpoint expression per cell type."""
        expr = adata.raw.to_adata() if adata.raw else adata
        available = {label: gene for label, gene in CHECKPOINT_GENES.items() if gene in var_names}
        if not available:
            return None

        gene_to_idx = {g: i for i, g in enumerate(expr.var_names)}
        records = []
        for cell_type in adata.obs["cell_type"].unique():
            mask = (adata.obs["cell_type"] == cell_type).values
            if mask.sum() < 5:
                continue
            row = {"cell_type": cell_type, "n_cells": int(mask.sum())}
            for label, gene in available.items():
                idx = gene_to_idx[gene]
                x = expr.X[mask, idx]
                if hasattr(x, "toarray"):
                    x = x.toarray()
                row[label] = round(float(np.mean(x)), 4)
                row[f"{label}_pct"] = round(float(np.mean(np.asarray(x).flatten() > 0)) * 100, 2)
            records.append(row)
        return pd.DataFrame(records) if records else None

    def _plot_visualizations(
        self, adata, ctx: PipelineContext, scored_sigs: list, checkpoint_data,
    ) -> None:
        if "X_umap" not in adata.obsm:
            return

        # CYT score UMAP
        if "cyt_score" in adata.obs:
            sc.pl.umap(adata, color="cyt_score", show=False, cmap="YlOrRd")
            plt.savefig(ctx.figure_dir / "tme_cyt_umap.png", dpi=160, bbox_inches="tight")
            plt.close()

        # TIS score UMAP
        if "tme_TIS" in adata.obs:
            sc.pl.umap(adata, color="tme_TIS", show=False, cmap="YlOrRd")
            plt.savefig(ctx.figure_dir / "tme_tis_umap.png", dpi=160, bbox_inches="tight")
            plt.close()

        # TME signature heatmap per cluster
        if "leiden" in adata.obs and scored_sigs:
            tme_cols = [f"tme_{s}" for s in scored_sigs]
            cluster_means = adata.obs[["leiden"] + tme_cols].groupby("leiden", observed=True).mean()
            fig, ax = plt.subplots(figsize=(max(8, len(tme_cols) * 0.8), max(4, len(cluster_means) * 0.4)))
            im = ax.imshow(cluster_means.values, aspect="auto", cmap="RdBu_r")
            ax.set_xticks(range(cluster_means.shape[1]))
            ax.set_xticklabels([s.replace("tme_", "") for s in cluster_means.columns], rotation=45, fontsize=8)
            ax.set_yticks(range(cluster_means.shape[0]))
            ax.set_yticklabels(cluster_means.index)
            plt.colorbar(im, ax=ax, label="Mean score")
            ax.set_xlabel("TME Signature")
            ax.set_ylabel("Cluster")
            ax.set_title("TME signature scores per cluster")
            plt.tight_layout()
            plt.savefig(ctx.figure_dir / "tme_signature_heatmap.png", dpi=160, bbox_inches="tight")
            plt.close()

        # Checkpoint dotplot
        if checkpoint_data is not None and len(checkpoint_data) >= 2:
            self._plot_checkpoint_dotplot(checkpoint_data, ctx)

        # Immune vs stromal bar chart
        if "leiden" in adata.obs and "tme_immune_estimate" in adata.obs and "tme_stromal_estimate" in adata.obs:
            self._plot_immune_stromal_bar(adata, ctx)

    @staticmethod
    def _plot_checkpoint_dotplot(df: pd.DataFrame, ctx: PipelineContext) -> None:
        """Dot plot: checkpoint expression per cell type (size=pct, color=expr)."""
        pct_cols = [c for c in df.columns if c.endswith("_pct")]
        expr_cols = [c.replace("_pct", "") for c in pct_cols]
        expr_cols = [c for c in expr_cols if c in df.columns and c not in ("cell_type", "n_cells")]
        if not expr_cols:
            return

        cell_types = df["cell_type"].values
        fig, ax = plt.subplots(figsize=(max(8, len(expr_cols) * 1.0), max(4, len(cell_types) * 0.4)))
        df_indexed = df.set_index("cell_type")
        pct_matrix = df_indexed[[f"{g}_pct" for g in expr_cols]].values
        expr_matrix = df_indexed[expr_cols].values
        vmax = expr_matrix.max() or 1
        jj, ii = np.meshgrid(range(len(expr_cols)), range(len(cell_types)))
        ax.scatter(
            jj.ravel(), ii.ravel(),
            s=pct_matrix.ravel() * 3,
            c=expr_matrix.ravel(),
            cmap="Reds", vmin=0, vmax=vmax,
            edgecolors="grey", linewidth=0.5,
        )

        ax.set_xticks(range(len(expr_cols)))
        ax.set_xticklabels(expr_cols, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(cell_types)))
        ax.set_yticklabels(cell_types, fontsize=8)
        ax.set_title("Checkpoint molecule expression per cell type")
        ax.set_xlabel("Checkpoint")
        ax.set_ylabel("Cell type")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "checkpoint_dotplot.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_immune_stromal_bar(adata, ctx: PipelineContext) -> None:
        """Bar chart of immune and stromal ESTIMATE scores per cluster."""
        cluster_data = adata.obs[["leiden", "tme_immune_estimate", "tme_stromal_estimate"]].groupby("leiden", observed=True).mean()
        x = np.arange(len(cluster_data))
        width = 0.35
        fig, ax = plt.subplots(figsize=(max(8, len(x) * 0.6), 5))
        ax.bar(x - width / 2, cluster_data["tme_immune_estimate"], width, label="Immune", color="steelblue")
        ax.bar(x + width / 2, cluster_data["tme_stromal_estimate"], width, label="Stromal", color="coral")
        ax.set_xticks(x)
        ax.set_xticklabels(cluster_data.index)
        ax.set_xlabel("Leiden cluster")
        ax.set_ylabel("Mean ESTIMATE score")
        ax.set_title("Immune vs Stromal ESTIMATE scores per cluster")
        ax.legend()
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "tme_immune_stromal_bar.png", dpi=160, bbox_inches="tight")
        plt.close()
