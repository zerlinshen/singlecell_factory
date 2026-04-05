from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..context import PipelineContext


class GeneRegulatoryNetworkModule:
    """Optional module: transcription factor activity inference.

    Supports multiple backends:
    1. decoupler + DoRothEA — TF activity scoring per cell/cluster
    2. Fallback — manual TF-target overlap scoring

    Identifies master regulators driving each cell state/cluster.
    """

    name = "gene_regulatory_network"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("GRN inference requires AnnData.")

        try:
            self._run_decoupler(adata, ctx)
        except ImportError:
            self._run_fallback(adata, ctx)

    def _run_decoupler(self, adata, ctx: PipelineContext) -> None:
        """TF activity inference using decoupler + DoRothEA regulons."""
        import decoupler as dc

        # Get DoRothEA regulons (confidence A/B/C)
        net = dc.get_dorothea(organism="human", levels=["A", "B", "C"])

        dc.run_ulm(
            mat=adata,
            net=net,
            source="source",
            target="target",
            weight="weight",
            use_raw=False,
        )

        if "ulm_estimate" not in adata.obsm:
            raise ValueError("decoupler ULM did not produce results.")

        tf_acts = pd.DataFrame(
            adata.obsm["ulm_estimate"],
            index=adata.obs_names,
            columns=[c for c in adata.obsm["ulm_estimate"].columns]
            if hasattr(adata.obsm["ulm_estimate"], "columns")
            else [f"TF_{i}" for i in range(adata.obsm["ulm_estimate"].shape[1])],
        )

        # Per-cluster mean TF activity
        if "leiden" in adata.obs:
            tf_acts["leiden"] = adata.obs["leiden"].values
            cluster_tf = tf_acts.groupby("leiden", observed=True).mean()
            cluster_tf.to_csv(ctx.table_dir / "tf_activity_per_cluster.csv")

            # Top TFs per cluster
            top_tfs = {}
            for cluster in cluster_tf.index:
                top = cluster_tf.loc[cluster].nlargest(5)
                top_tfs[str(cluster)] = {tf: round(float(score), 3) for tf, score in top.items()}

            (ctx.table_dir / "tf_top_per_cluster.json").write_text(
                json.dumps(top_tfs, indent=2), encoding="utf-8",
            )
            ctx.metadata["grn_top_tfs"] = list(top_tfs.get("0", {}).keys())

            # Heatmap
            self._plot_tf_heatmap(cluster_tf, ctx)
        else:
            tf_acts.to_csv(ctx.table_dir / "tf_activity_per_cell.csv")

        # UMAP of top TFs
        if "X_umap" in adata.obsm:
            self._plot_tf_umap(adata, ctx)

    def _run_fallback(self, adata, ctx: PipelineContext) -> None:
        """Manual TF activity scoring using curated TF-target sets."""
        import scanpy as sc

        # Curated key TF-target sets for lung cancer / TME
        tf_targets = {
            "TP53": ["CDKN1A", "MDM2", "BAX", "FAS", "GADD45A", "BBC3", "PMAIP1"],
            "MYC": ["NCL", "NPM1", "LDHA", "ENO1", "CDK4", "PCNA", "RAN"],
            "STAT1": ["IRF1", "GBP1", "ISG15", "CXCL10", "TAP1", "B2M"],
            "NFKB1": ["NFKBIA", "TNF", "IL6", "CXCL8", "CCL2", "ICAM1"],
            "HIF1A": ["VEGFA", "SLC2A1", "LDHA", "PGK1", "CA9", "BNIP3"],
            "SOX2": ["NANOG", "POU5F1", "SALL4", "BMI1", "KLF4"],
            "FOXP3": ["IL2RA", "CTLA4", "TNFRSF18", "IKZF2", "LRRC32"],
            "TBX21": ["IFNG", "CXCR3", "IL12RB2", "STAT4", "CCL3"],
            "GATA3": ["IL4", "IL5", "IL13", "CCR4", "STAT6"],
            "TWIST1": ["VIM", "CDH2", "SNAI1", "MMP2", "FN1"],
            "SNAI1": ["VIM", "CDH2", "MMP9", "FN1", "TWIST1"],
            "NKX2-1": ["SFTPA1", "SFTPB", "SFTPC", "SFTPD", "SCGB1A1"],
        }

        var_names = set(adata.var_names)
        results = {}

        for tf, targets in tf_targets.items():
            valid_targets = [g for g in targets if g in var_names]
            if len(valid_targets) >= 2:
                sc.tl.score_genes(
                    adata, valid_targets, score_name=f"tf_{tf}", use_raw=False,
                )
                results[tf] = valid_targets

        if not results:
            raise ValueError("Insufficient TF target genes found in dataset.")

        score_cols = [f"tf_{tf}" for tf in results]
        scores = adata.obs[score_cols].copy()
        scores.columns = list(results.keys())

        if "leiden" in adata.obs:
            scores["leiden"] = adata.obs["leiden"].values
            cluster_tf = scores.groupby("leiden", observed=True).mean()
            cluster_tf.to_csv(ctx.table_dir / "tf_activity_per_cluster.csv")

            top_tfs = {}
            for cluster in cluster_tf.index:
                top = cluster_tf.loc[cluster].nlargest(5)
                top_tfs[str(cluster)] = {tf: round(float(s), 3) for tf, s in top.items()}

            (ctx.table_dir / "tf_top_per_cluster.json").write_text(
                json.dumps(top_tfs, indent=2), encoding="utf-8",
            )
            ctx.metadata["grn_top_tfs"] = list(top_tfs.get("0", {}).keys())

            self._plot_tf_heatmap(cluster_tf, ctx)

    @staticmethod
    def _plot_tf_heatmap(cluster_tf: pd.DataFrame, ctx: PipelineContext) -> None:
        """Heatmap of TF activities across clusters."""
        # Select most variable TFs
        tf_var = cluster_tf.var()
        top_tfs = tf_var.nlargest(min(30, len(tf_var))).index
        data = cluster_tf[top_tfs]

        fig, ax = plt.subplots(figsize=(max(8, len(top_tfs) * 0.5), max(4, len(data) * 0.4)))
        im = ax.imshow(data.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set_xticks(range(data.shape[1]))
        ax.set_xticklabels(data.columns, rotation=90, fontsize=8)
        ax.set_yticks(range(data.shape[0]))
        ax.set_yticklabels(data.index, fontsize=9)
        ax.set_xlabel("Transcription Factor")
        ax.set_ylabel("Cluster")
        plt.colorbar(im, ax=ax, label="Activity score")
        plt.title("TF activity per cluster")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "tf_activity_heatmap.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_tf_umap(adata, ctx: PipelineContext) -> None:
        """UMAP plots for top variable TFs."""
        import scanpy as sc

        tf_cols = [c for c in adata.obs.columns if c.startswith("tf_") or c.startswith("ulm_")]
        if not tf_cols:
            return

        # Pick top 4 most variable
        var_scores = adata.obs[tf_cols].var().nlargest(4)
        if var_scores.empty:
            return

        n = len(var_scores)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
        if n == 1:
            axes = [axes]
        for ax, col in zip(axes, var_scores.index):
            sc.pl.umap(adata, color=col, ax=ax, show=False, cmap="RdBu_r")
            ax.set_title(col.replace("tf_", "").replace("ulm_", ""))
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "tf_activity_umap.png", dpi=160, bbox_inches="tight")
        plt.close()
