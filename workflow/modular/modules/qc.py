from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from ..context import PipelineContext


class QCModule:
    """Mandatory module: compute QC metrics, filter low-quality cells/genes, and generate QC plots."""

    name = "qc"
    required = True

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("QC requires loaded AnnData.")

        gene_upper = adata.var_names.astype(str).str.upper()
        adata.var["mt"] = gene_upper.str.startswith("MT-")
        adata.var["ribo"] = gene_upper.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = gene_upper.str.startswith(("HBA", "HBB"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True)

        # --- QC visualizations (before filtering) ---
        self._plot_qc_metrics(adata, ctx, suffix="pre_filter")

        qc = ctx.cfg.qc
        before = adata.n_obs
        # Use one combined mask to avoid multiple full AnnData copies.
        mask = (
            (adata.obs["n_genes_by_counts"] >= qc.min_genes)
            & (adata.obs["n_genes_by_counts"] <= qc.max_genes)
            & (adata.obs["total_counts"] >= qc.min_counts)
            & (adata.obs["total_counts"] <= qc.max_counts)
            & (adata.obs["pct_counts_mt"] <= qc.max_mito_pct)
            & (adata.obs["pct_counts_ribo"] <= qc.max_ribo_pct)
        )
        adata = adata[mask, :].copy()
        sc.pp.filter_genes(adata, min_cells=qc.min_cells)

        # --- QC visualizations (after filtering) ---
        self._plot_qc_metrics(adata, ctx, suffix="post_filter")

        ctx.adata = adata
        ctx.metadata["cells_after_qc"] = int(adata.n_obs)
        ctx.metadata["genes_after_qc"] = int(adata.n_vars)
        ctx.metadata["qc_removed_cells"] = int(before - adata.n_obs)

    @staticmethod
    def _plot_qc_metrics(adata, ctx: PipelineContext, suffix: str) -> None:
        """Generate violin + scatter QC plots."""
        # Violin plots
        fig, axes = plt.subplots(1, 4, figsize=(16, 4))
        for ax, key, title in zip(
            axes,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
            ["Genes per cell", "UMI counts per cell", "Mito %", "Ribo %"],
        ):
            sc.pl.violin(adata, key, ax=ax, show=False, stripplot=False)
            ax.set_title(title)
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / f"qc_violin_{suffix}.png", dpi=160, bbox_inches="tight")
        plt.close()

        # Scatter plots: genes vs counts, colored by mito%
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.scatter(
            adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt",
            ax=axes[0], show=False, title="Genes vs UMI (color=mito%)",
        )
        sc.pl.scatter(
            adata, x="total_counts", y="pct_counts_mt",
            ax=axes[1], show=False, title="Mito% vs UMI counts",
        )
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / f"qc_scatter_{suffix}.png", dpi=160, bbox_inches="tight")
        plt.close()
