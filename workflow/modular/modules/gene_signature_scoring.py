from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from ..context import PipelineContext
from . import score_gene_sets

# Built-in cancer-relevant gene signatures from high-impact publications.
BUILTIN_SIGNATURES = {
    # Proliferation — MKI67/TOP2A panel (standard in oncology)
    "proliferation": [
        "MKI67", "TOP2A", "PCNA", "MCM2", "MCM6", "CDK1", "CCNB1", "CCNB2",
    ],
    # Apoptosis resistance
    "apoptosis_resistance": [
        "BCL2", "BCL2L1", "MCL1", "BIRC5", "XIAP", "CFLAR",
    ],
    # Angiogenesis — VEGF/FLT pathway
    "angiogenesis": [
        "VEGFA", "VEGFB", "FLT1", "KDR", "PECAM1", "ANGPT2", "NRP1",
    ],
    # Invasion and metastasis — MMP/EMT markers
    "invasion_metastasis": [
        "MMP2", "MMP9", "MMP14", "SNAI1", "TWIST1", "VIM", "CDH2",
    ],
    # EMT mesenchymal signature — Tan et al., EMBO Mol Med 2014
    "EMT_mesenchymal": [
        "VIM", "CDH2", "FN1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "MMP2",
    ],
    # EMT epithelial signature — Tan et al., EMBO Mol Med 2014
    "EMT_epithelial": [
        "CDH1", "EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4", "OCLN",
    ],
    # Stemness — Malta et al., Cell 2018; Ben-Porath et al., Nat Genet 2008
    "stemness": [
        "POU5F1", "NANOG", "SOX2", "KLF4", "MYC", "LIN28A", "SALL4", "BMI1",
    ],
    # Hypoxia — Buffa et al., Br J Cancer 2010
    "hypoxia": [
        "VEGFA", "SLC2A1", "HK2", "LDHA", "PGK1", "CA9", "BNIP3", "ENO1",
    ],
    # DNA damage response
    "DDR": [
        "BRCA1", "BRCA2", "ATM", "ATR", "RAD51", "CHEK1", "CHEK2", "TP53",
    ],
    # Glycolysis — Warburg effect
    "glycolysis": [
        "HK2", "PFKP", "PKM", "LDHA", "ENO1", "GAPDH", "TPI1", "ALDOA",
    ],
}


class GeneSignatureScoringModule:
    """Optional module: score cells against gene signatures.

    Scores built-in cancer hallmark signatures and optional user-defined signatures
    loaded from a JSON file. Produces per-cell and per-cluster signature scores,
    heatmaps, UMAP overlays, and signature correlation analysis.

    Built-in signatures cover: proliferation, apoptosis resistance, angiogenesis,
    invasion, EMT (Tan et al. 2014), stemness (Malta et al. 2018), hypoxia
    (Buffa et al. 2010), DNA damage response, and glycolysis.
    """

    name = "gene_signature_scoring"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Gene signature scoring requires AnnData.")

        # Merge built-in and user signatures
        signatures = {}
        sig_cfg = ctx.cfg.gene_signature
        if sig_cfg.use_builtin:
            signatures.update(BUILTIN_SIGNATURES)
        if sig_cfg.signature_json and sig_cfg.signature_json.exists():
            user_sigs = json.loads(sig_cfg.signature_json.read_text(encoding="utf-8"))
            signatures.update({str(k): [str(g) for g in v] for k, v in user_sigs.items()})

        if not signatures:
            raise ValueError("No gene signatures available (built-in disabled and no JSON provided).")

        # Score each signature
        scored = score_gene_sets(adata, signatures, "sig")

        if not scored:
            raise ValueError("No gene signatures had sufficient genes in dataset.")

        # --- Save tables ---
        sig_cols = [f"sig_{s}" for s in scored]
        adata.obs[sig_cols].to_csv(ctx.table_dir / "gene_signature_scores.csv")

        if "leiden" in adata.obs:
            cluster_scores = adata.obs[["leiden"] + sig_cols].groupby("leiden", observed=True).mean()
            cluster_scores.to_csv(ctx.table_dir / "gene_signature_per_cluster.csv")

        ctx.metadata["gene_signatures_scored"] = scored

        # --- Visualizations ---
        self._plot_visualizations(adata, ctx, scored, sig_cols)

    def _plot_visualizations(self, adata, ctx, scored, sig_cols) -> None:
        # Signature heatmap per cluster
        if "leiden" in adata.obs and len(scored) >= 2:
            cluster_scores = adata.obs[["leiden"] + sig_cols].groupby("leiden", observed=True).mean()
            fig, ax = plt.subplots(
                figsize=(max(8, len(scored) * 0.7), max(4, len(cluster_scores) * 0.4))
            )
            im = ax.imshow(cluster_scores.values, aspect="auto", cmap="RdBu_r")
            ax.set_xticks(range(cluster_scores.shape[1]))
            ax.set_xticklabels(
                [s.replace("sig_", "") for s in cluster_scores.columns], rotation=45, ha="right", fontsize=8,
            )
            ax.set_yticks(range(cluster_scores.shape[0]))
            ax.set_yticklabels(cluster_scores.index)
            plt.colorbar(im, ax=ax, label="Mean score")
            ax.set_xlabel("Gene Signature")
            ax.set_ylabel("Cluster")
            ax.set_title("Gene signature scores per cluster")
            plt.tight_layout()
            plt.savefig(ctx.figure_dir / "signature_heatmap.png", dpi=160, bbox_inches="tight")
            plt.close()

        # UMAP of top 4 most variable signatures
        if "X_umap" in adata.obsm and sig_cols:
            var_scores = adata.obs[sig_cols].var().nlargest(min(4, len(sig_cols)))
            n = len(var_scores)
            fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
            if n == 1:
                axes = [axes]
            for ax, col in zip(axes, var_scores.index):
                sc.pl.umap(adata, color=col, ax=ax, show=False, cmap="RdBu_r")
                ax.set_title(col.replace("sig_", ""))
            plt.tight_layout()
            plt.savefig(ctx.figure_dir / "signature_umap.png", dpi=160, bbox_inches="tight")
            plt.close()

        # Signature correlation matrix
        if len(sig_cols) >= 3:
            corr = adata.obs[sig_cols].corr()
            fig, ax = plt.subplots(figsize=(max(6, len(sig_cols) * 0.5), max(5, len(sig_cols) * 0.5)))
            im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
            labels = [s.replace("sig_", "") for s in corr.columns]
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels, fontsize=7)
            plt.colorbar(im, ax=ax, label="Pearson r")
            ax.set_title("Gene signature correlation")
            plt.tight_layout()
            plt.savefig(ctx.figure_dir / "signature_correlation.png", dpi=160, bbox_inches="tight")
            plt.close()
