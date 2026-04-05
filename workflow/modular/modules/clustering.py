from __future__ import annotations

import logging
import os
import matplotlib
import numpy as np
from scipy import sparse

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from ..context import PipelineContext

logger = logging.getLogger(__name__)


def _gpu_available() -> bool:
    """Check if rapids-singlecell GPU backend is usable."""
    try:
        import rapids_singlecell  # noqa: F401
        import cupy  # noqa: F401
        cupy.cuda.runtime.getDeviceCount()
        return True
    except Exception:
        return False


class ClusteringModule:
    """Optional module: normalization, PCA/UMAP and Leiden clustering.

    Supports GPU acceleration via rapids-singlecell when available.
    Falls back to standard scanpy CPU backend automatically.
    """

    name = "clustering"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Clustering requires AnnData.")
        cfg = ctx.cfg.clustering
        sc.settings.n_jobs = max(1, os.cpu_count() or 1)

        use_gpu = _gpu_available()
        if use_gpu:
            logger.info("GPU detected — using rapids-singlecell for PCA/neighbors/UMAP")
            self._run_gpu(adata, cfg, ctx)
        else:
            self._run_cpu(adata, cfg, ctx)

        self._plot_umap_clusters(adata, ctx)
        ctx.adata = adata
        ctx.metadata["n_clusters"] = int(adata.obs["leiden"].nunique())
        ctx.metadata["clustering_backend"] = "gpu" if use_gpu else "cpu"

    def _run_cpu(self, adata, cfg, ctx) -> None:
        """Standard scanpy CPU clustering pipeline."""
        sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
        sc.pp.log1p(adata)
        if adata.raw is None:
            # Preserve normalized/log-transformed matrix for marker-level downstream tasks.
            adata.raw = adata
        is_sparse = sparse.issparse(adata.X)
        if is_sparse:
            adata.X = adata.X.astype(np.float32)
        else:
            adata.X = np.asarray(adata.X, dtype=np.float32)
        if cfg.scale_data:
            sc.pp.scale(adata, max_value=10, zero_center=not sparse.issparse(adata.X))
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=cfg.n_top_genes)
        sc.tl.pca(
            adata,
            n_comps=cfg.n_pcs,
            svd_solver="arpack" if is_sparse else "randomized",
            mask_var="highly_variable",
        )
        self._plot_pca_variance(adata, ctx, cfg.n_pcs)
        sc.pp.neighbors(
            adata,
            n_neighbors=cfg.n_neighbors,
            n_pcs=cfg.n_pcs,
            use_rep="X_pca",
            method="umap",
        )
        sc.tl.umap(adata, random_state=cfg.random_state)
        sc.tl.leiden(
            adata,
            resolution=cfg.leiden_resolution,
            flavor="igraph",
            directed=False,
            random_state=cfg.random_state,
        )

    def _run_gpu(self, adata, cfg, ctx) -> None:
        """GPU-accelerated clustering via rapids-singlecell."""
        import rapids_singlecell as rsc

        sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
        sc.pp.log1p(adata)
        if adata.raw is None:
            adata.raw = adata
        if sparse.issparse(adata.X):
            adata.X = adata.X.astype(np.float32)
        else:
            adata.X = np.asarray(adata.X, dtype=np.float32)
        if cfg.scale_data:
            sc.pp.scale(adata, max_value=10, zero_center=not sparse.issparse(adata.X))
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=cfg.n_top_genes)

        # GPU-accelerated PCA, neighbors, UMAP
        rsc.pp.pca(adata, n_comps=cfg.n_pcs, mask_var="highly_variable")
        self._plot_pca_variance(adata, ctx, cfg.n_pcs)
        rsc.pp.neighbors(adata, n_neighbors=cfg.n_neighbors, n_pcs=cfg.n_pcs, use_rep="X_pca")
        rsc.tl.umap(adata, random_state=cfg.random_state)
        rsc.tl.leiden(adata, resolution=cfg.leiden_resolution, random_state=cfg.random_state)

    @staticmethod
    def _plot_pca_variance(adata, ctx: PipelineContext, n_pcs: int) -> None:
        """Plot PCA variance explained (elbow plot + cumulative)."""
        var_ratio = adata.uns["pca"]["variance_ratio"]
        cumulative = np.cumsum(var_ratio)

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        axes[0].plot(range(1, len(var_ratio) + 1), var_ratio, "o-", markersize=3)
        axes[0].set_xlabel("Principal Component")
        axes[0].set_ylabel("Variance Ratio")
        axes[0].set_title("PCA Elbow Plot")

        axes[1].plot(range(1, len(cumulative) + 1), cumulative, "o-", markersize=3)
        axes[1].axhline(0.9, color="grey", linestyle="--", linewidth=0.8, label="90%")
        axes[1].set_xlabel("Principal Component")
        axes[1].set_ylabel("Cumulative Variance")
        axes[1].set_title("Cumulative Variance Explained")
        axes[1].legend()

        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pca_variance_explained.png", dpi=160, bbox_inches="tight")
        plt.close()

        # Record how many PCs needed for 90% variance
        pcs_for_90 = int(np.searchsorted(cumulative, 0.9) + 1)
        ctx.metadata["pca_cumulative_var_90pct"] = min(pcs_for_90, n_pcs)

    @staticmethod
    def _plot_umap_clusters(adata, ctx: PipelineContext) -> None:
        sc.pl.umap(adata, color=["leiden"], show=False, legend_loc="on data")
        plt.savefig(ctx.figure_dir / "umap_leiden.png", dpi=160, bbox_inches="tight")
        plt.close()
