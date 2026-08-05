from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ._scanpy_compat import import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext
from . import score_gene_sets


__references__ = {
    "Tirosh_scoring_2016": {
        "title": "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq",
        "authors": "Tirosh et al.",
        "journal": "Science",
        "year": "2016",
        "doi": "10.1126/science.aad0501",
        "description": "Signature score methodology (mean expression minus matched control bin). sc.tl.score_genes implements this.",
    },
    "Mulder_mac_atlas_2021": {
        "title": "Cross-tissue single-cell landscape of human monocytes and macrophages in health and disease",
        "authors": "Mulder et al.",
        "journal": "Immunity",
        "year": "2021",
        "doi": "10.1016/j.immuni.2021.07.007",
        "description": (
            "Cross-tissue monocyte/macrophage programmes; tissue-resident / foetal-like "
            "macrophage panels (FOLR2/LYVE1/MRC1-axis) inform the foetal_like_mac builtin. "
            "Panel is literature-inspired, not an exact paper gene-list reimplementation."
        ),
    },
    "DeZuani_NSCLC_2024": {
        "title": "Single-cell and spatial transcriptomics analysis of non-small cell lung cancer",
        "authors": "De Zuani et al.",
        "journal": "Nature Communications",
        "year": "2024",
        "doi": "10.1038/s41467-024-48700-8",
        "description": (
            "STAB1+ tumour macrophages with oncofoetal reprogramming, iron export "
            "(SLC40A1/ferroportin), and cholesterol-export programmes (ABCA1/TREM2). "
            "Paper text names foetal-typical STAB1 signature genes: STAB1, FOLR2, "
            "SLC40A1, MERTK, GPR34, F13A1. Full 20-gene STAB1 signature is Fig 5I; "
            "de_zuani_stab1_signature_reconstructed is DEA-reconstructed from "
            "Supplementary Data 20/21 (AM∩AIM STAB1-up, padj≤0.05, |log2FC|≥1, top-20 "
            "by min|log2FC|), not a pixel OCR of Fig 5I."
        ),
    },
}


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
    # Foetal-like / tissue-resident macrophage (TRM) programme — Mulder-inspired
    # (Wave-6). Secondary to De Zuani paper-anchored panels below.
    "foetal_like_mac": [
        "FOLR2", "LYVE1", "MRC1", "SIGLEC1", "C1QA", "C1QB",
        "SELENOP", "RNASE1", "F13A1", "CD163", "MARCO",
    ],
    # Contrasting SPP1+ / inflammatory macrophage axis (companion panel)
    "spp1_mac": [
        "SPP1", "FABP5", "TREM2", "APOE", "CTSB", "CTSD", "LGALS3", "CHIT1",
    ],
    # De Zuani et al. 2024 Nat Commun — foetal-typical STAB1 signature genes
    # named in Results (oncofoetal reprogramming section). Prefer this panel for
    # NC2024-F3-02 directional claims over foetal_like_mac.
    "de_zuani_stab1_foetal": [
        "STAB1", "FOLR2", "SLC40A1", "MERTK", "GPR34", "F13A1",
    ],
    # DEA-reconstructed STAB1 signature (Supp Data 20 AM vs STAB1 + Data 21 AIM vs
    # STAB1): STAB1-up (negative LFC in other-vs-STAB1 tables), padj≤0.05,
    # |log2FC|≥1, intersection ranked by min|log2FC|, top 20. Approx. Fig 5I.
    "de_zuani_stab1_signature_reconstructed": [
        "SLC40A1", "SELENOP", "OLFML3", "FCGBP", "IGSF21", "IL2RA", "C3",
        "SIGLEC8", "F13A1", "ITGA9", "CTTNBP2", "STAB1", "ADGRG6", "ADAMDEC1",
        "CXCL12", "NCKAP5", "ENPP2", "SRGAP3", "AC079015.1", "IGF1",
    ],
    # Cholesterol export / lipid handling markers highlighted for tumour AMɸ/AIMɸ
    "de_zuani_cholesterol_export": [
        "ABCA1", "TREM2", "APOE", "APOC1", "SPP1", "FABP5",
    ],
}


class GeneSignatureScoringModule:
    """Optional module: score cells against gene signatures.

    Scores built-in cancer hallmark signatures and optional user-defined signatures
    loaded from a JSON file. Produces per-cell and per-cluster signature scores,
    heatmaps, UMAP overlays, and signature correlation analysis.

    Built-in signatures cover: proliferation, apoptosis resistance, angiogenesis,
    invasion, EMT (Tan et al. 2014), stemness (Malta et al. 2018), hypoxia
    (Buffa et al. 2010), DNA damage response, glycolysis, Mulder-inspired TRM
    (foetal_like_mac / spp1_mac), and De Zuani 2024 STAB1/oncofoetal panels
    (de_zuani_stab1_foetal, de_zuani_stab1_signature_reconstructed,
    de_zuani_cholesterol_export).
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
            plt.savefig(ctx.figure_dir / "signature_heatmap.png", bbox_inches="tight")
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
            plt.savefig(ctx.figure_dir / "signature_umap.png", bbox_inches="tight")
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
            plt.savefig(ctx.figure_dir / "signature_correlation.png", bbox_inches="tight")
            plt.close()
