from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from ._scanpy_compat import has_api, import_scanpy_or_stub

sc = import_scanpy_or_stub()

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
        use_scanpy_qc = has_api(sc, "pp.calculate_qc_metrics") and ctx.cfg.scale_mode != "massive"
        if use_scanpy_qc:
            try:
                sc.pp.calculate_qc_metrics(
                    adata,
                    qc_vars=["mt", "ribo", "hb"],
                    inplace=True,
                    log1p=False,
                    percent_top=None,
                )
            except Exception:
                self._calculate_qc_metrics_fallback(adata)
        else:
            self._calculate_qc_metrics_fallback(adata)

        # Drop verbose columns to save memory
        self._ensure_axis_tables_support_drop(adata)
        obs_drop = [c for c in adata.obs.columns if c.startswith("log1p") or c.startswith("total_counts_") or c.startswith("pct_counts_in_top_")]
        adata.obs.drop(columns=obs_drop, inplace=True, errors="ignore")
        var_drop = [c for c in adata.var.columns if c.startswith("log1p") or c.startswith("pct_dropout") or "total_counts" in c or c.startswith("mean_counts") or c.startswith("n_cells_by") or c in ["mt", "ribo", "hb"]]
        adata.var.drop(columns=var_drop, inplace=True, errors="ignore")

        # --- QC visualizations (before filtering) ---
        if ctx.cfg.scale_mode != "massive":
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
        if has_api(sc, "pp.filter_genes"):
            sc.pp.filter_genes(adata, min_cells=qc.min_cells)
        else:
            adata = self._filter_genes_fallback(adata, min_cells=qc.min_cells)

        # --- QC visualizations (after filtering) ---
        if ctx.cfg.scale_mode != "massive":
            self._plot_qc_metrics(adata, ctx, suffix="post_filter")

        ctx.adata = adata
        ctx.metadata["cells_after_qc"] = int(adata.n_obs)
        ctx.metadata["genes_after_qc"] = int(adata.n_vars)
        ctx.metadata["qc_removed_cells"] = int(before - adata.n_obs)

    @staticmethod
    def _plot_qc_metrics(adata, ctx: PipelineContext, suffix: str) -> None:
        """Generate violin + scatter QC plots."""
        if not has_api(sc, "pl.violin") or not has_api(sc, "pl.scatter"):
            return
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


    @staticmethod
    def _ensure_axis_tables_support_drop(adata) -> None:
        for axis in ("obs", "var"):
            table = getattr(adata, axis)
            if hasattr(table, "drop"):
                continue
            if hasattr(table, "to_memory"):
                table = table.to_memory()
            elif hasattr(table, "copy"):
                table = table.copy()
            setattr(adata, axis, table)

    @staticmethod
    def _materialize_matrix(x):
        if sparse.issparse(x):
            return x.tocsr()
        if hasattr(x, "to_memory"):
            try:
                x = x.to_memory()
            except Exception:
                pass
        if sparse.issparse(x):
            return x.tocsr()
        if hasattr(x, "compute"):
            try:
                x = x.compute()
            except Exception:
                pass
        return x

    @staticmethod
    def _calculate_qc_metrics_fallback(adata) -> None:
        x = QCModule._materialize_matrix(adata.X)
        if sparse.issparse(x):
            x = x.tocsr()
            total_counts = np.asarray(x.sum(axis=1)).ravel()
            n_genes = np.asarray((x > 0).sum(axis=1)).ravel()
            mt_counts = np.asarray(x[:, np.asarray(adata.var["mt"], dtype=bool)].sum(axis=1)).ravel()
            ribo_counts = np.asarray(x[:, np.asarray(adata.var["ribo"], dtype=bool)].sum(axis=1)).ravel()
            hb_counts = np.asarray(x[:, np.asarray(adata.var["hb"], dtype=bool)].sum(axis=1)).ravel()
        else:
            arr = np.asarray(x, dtype=float)
            total_counts = arr.sum(axis=1)
            n_genes = (arr > 0).sum(axis=1)
            mt_counts = arr[:, np.asarray(adata.var["mt"], dtype=bool)].sum(axis=1)
            ribo_counts = arr[:, np.asarray(adata.var["ribo"], dtype=bool)].sum(axis=1)
            hb_counts = arr[:, np.asarray(adata.var["hb"], dtype=bool)].sum(axis=1)
        denom = np.clip(total_counts, 1e-9, None)
        adata.obs["n_genes_by_counts"] = n_genes
        adata.obs["total_counts"] = total_counts
        adata.obs["pct_counts_mt"] = mt_counts / denom * 100.0
        adata.obs["pct_counts_ribo"] = ribo_counts / denom * 100.0
        adata.obs["pct_counts_hb"] = hb_counts / denom * 100.0

    @staticmethod
    def _filter_genes_fallback(adata, min_cells: int):
        x = QCModule._materialize_matrix(adata.X)
        if sparse.issparse(x):
            keep = np.asarray((x > 0).sum(axis=0)).ravel() >= int(min_cells)
        else:
            keep = (np.asarray(x) > 0).sum(axis=0) >= int(min_cells)
        if keep.ndim != 1:
            keep = np.asarray(keep).ravel()
        return adata[:, keep].copy()
