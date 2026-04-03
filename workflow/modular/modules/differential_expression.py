from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from ..context import PipelineContext


class DifferentialExpressionModule:
    """Optional module: identify cluster marker genes with significance filtering and visualizations."""

    name = "differential_expression"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "leiden" not in adata.obs:
            raise ValueError("Differential expression requires clustered AnnData.")

        sc.tl.rank_genes_groups(
            adata, groupby="leiden", method="wilcoxon", use_raw=True, n_genes=500,
        )
        markers = sc.get.rank_genes_groups_df(adata, group=None)

        # Save unfiltered results
        markers.to_csv(ctx.table_dir / "marker_genes_all.csv", index=False)

        # Apply significance filtering
        sig_markers = markers[
            (markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"].abs() > 0.25)
        ].copy()
        sig_markers.to_csv(ctx.table_dir / "marker_genes.csv", index=False)

        top5 = (
            sig_markers.sort_values(["group", "scores"], ascending=[True, False])
            .groupby("group", observed=True)
            .head(5)
            .reset_index(drop=True)
        )
        top5.to_csv(ctx.table_dir / "marker_top5_by_cluster.csv", index=False)

        ctx.metadata["de_significant_genes"] = int(len(sig_markers))

        # --- Visualizations ---
        # Dot plot of top markers per cluster
        try:
            sc.pl.rank_genes_groups_dotplot(
                adata, n_genes=5, show=False, standard_scale="var",
            )
            plt.savefig(
                ctx.figure_dir / "de_dotplot_top5.png", dpi=160, bbox_inches="tight",
            )
            plt.close()
        except Exception:
            plt.close("all")

        # Heatmap of top markers
        try:
            sc.pl.rank_genes_groups_heatmap(
                adata, n_genes=5, show=False, show_gene_labels=True,
                use_raw=True, swap_axes=True, vmin=-3, vmax=3,
            )
            plt.savefig(
                ctx.figure_dir / "de_heatmap_top5.png", dpi=160, bbox_inches="tight",
            )
            plt.close()
        except Exception:
            plt.close("all")

        # Volcano plot (all clusters combined)
        self._volcano_plot(sig_markers, ctx)

    @staticmethod
    def _volcano_plot(markers: pd.DataFrame, ctx: PipelineContext) -> None:
        """Generate a combined volcano plot for all clusters."""
        if markers.empty:
            return
        fig, ax = plt.subplots(figsize=(10, 7))
        log_pval = -np.log10(markers["pvals_adj"].clip(lower=1e-300))
        lfc = markers["logfoldchanges"]

        scatter = ax.scatter(
            lfc, log_pval, c=lfc, cmap="RdBu_r", s=4, alpha=0.6,
            vmin=-3, vmax=3, edgecolors="none",
        )
        ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.8)
        ax.axvline(-0.25, color="grey", linestyle="--", linewidth=0.8)
        ax.axvline(0.25, color="grey", linestyle="--", linewidth=0.8)
        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-log10(adjusted p-value)")
        ax.set_title("Volcano Plot (all clusters)")
        plt.colorbar(scatter, ax=ax, label="Log2FC")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "de_volcano.png", dpi=160, bbox_inches="tight")
        plt.close()
