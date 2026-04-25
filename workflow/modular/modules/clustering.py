from __future__ import annotations

import gc
import logging
import os
import matplotlib
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ._scanpy_compat import import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext
from ._gpu_utils import gpu_available

logger = logging.getLogger(__name__)


class ClusteringModule:
    @staticmethod
    def _materialize_matrix(x):
        if sparse.issparse(x):
            return x.tocsr().astype(np.float32)
        if hasattr(x, "to_memory"):
            try:
                x = x.to_memory()
            except Exception:
                pass
        if sparse.issparse(x):
            return x.tocsr().astype(np.float32)
        if hasattr(x, "compute"):
            try:
                x = x.compute()
            except Exception:
                pass
        if sparse.issparse(x):
            return x.tocsr().astype(np.float32)
        return np.asarray(x, dtype=np.float32)

    """Optional module: normalization, PCA/UMAP and Leiden clustering.

    Supports GPU acceleration via rapids-singlecell when available.
    Falls back to standard scanpy CPU backend automatically.
    """

    name = "clustering"
    provides_keys = {"obs": ["leiden"], "obsm": ["X_pca", "X_umap", "X_css"]}

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Clustering requires AnnData.")
        cfg = ctx.cfg.clustering
        sc.settings.n_jobs = max(1, os.cpu_count() or 1)

        if self._should_use_css(ctx, adata):
            logger.info("Massive scale-mode with sample labels detected — using CSS-style clustering representation")
            self._run_css(adata, cfg, ctx)
            use_gpu = False
        else:
            use_gpu = gpu_available(ctx.cfg.gpu_mode)
            if ctx.cfg.scale_mode == "massive" and ctx.cfg.gpu_mode == "auto":
                logger.info("Massive scale-mode: prefer CPU clustering path to avoid GPU copy overhead")
                ctx.metadata["clustering_gpu_disabled_reason"] = "massive_scale_mode"
                use_gpu = False
            if use_gpu:
                logger.info("GPU detected — using rapids-singlecell for PCA/neighbors/UMAP/Leiden")
                adata_gpu = adata.copy()
                try:
                    self._run_gpu(adata_gpu, cfg, ctx)
                    adata = adata_gpu
                except Exception as exc:
                    if ctx.cfg.gpu_mode == "force":
                        raise RuntimeError(f"GPU clustering failed in force mode: {exc}") from exc
                    logger.warning("GPU clustering failed (%s), trying hybrid GPU path", exc)
                    ctx.metadata["clustering_gpu_fallback_reason"] = str(exc)
                    use_gpu = False
                    if self._run_hybrid_gpu_graph(adata, cfg, ctx):
                        ctx.metadata["clustering_backend"] = "hybrid"
                    else:
                        self._run_cpu(adata, cfg, ctx)
                finally:
                    gc.collect()
            else:
                self._run_cpu(adata, cfg, ctx)

        self._plot_umap_clusters(adata, ctx)
        ctx.adata = adata
        ctx.metadata["n_clusters"] = int(adata.obs["leiden"].nunique())
        ctx.metadata["gpu_mode"] = ctx.cfg.gpu_mode
        if "clustering_backend" not in ctx.metadata:
            ctx.metadata["clustering_backend"] = "gpu" if use_gpu else "cpu"

    def _should_use_css(self, ctx, adata) -> bool:
        if os.environ.get("SC_CLUSTERING_ENGINE", "").lower() == "sparse_exact":
            ctx.metadata["css_status"] = "disabled_sparse_exact_engine"
            ctx.metadata["clustering_engine"] = "sparse_exact"
            return False
        if ctx.cfg.scale_mode != "massive":
            return False
        if "sample" not in adata.obs.columns:
            ctx.metadata["css_status"] = "unavailable_no_sample_labels"
            return False
        n_samples = int(pd.Series(adata.obs["sample"]).nunique())
        if n_samples < 2:
            ctx.metadata["css_status"] = "unavailable_single_sample"
            return False
        ctx.metadata["css_status"] = "enabled"
        ctx.metadata["css_n_samples"] = n_samples
        return True

    def _run_css(self, adata, cfg, ctx) -> None:
        sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
        sc.pp.log1p(adata)
        adata.X = self._materialize_matrix(adata.X)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=cfg.n_top_genes)
        hvg_mask = np.asarray(adata.var["highly_variable"]).astype(bool)
        X = adata[:, hvg_mask].X
        if sparse.issparse(X):
            X = X.tocsr().astype(np.float32)
            svd_dim = min(max(2, cfg.n_pcs), 64, max(2, X.shape[1] - 1))
            svd = TruncatedSVD(n_components=svd_dim, random_state=cfg.random_state)
            X_latent = svd.fit_transform(X).astype(np.float32, copy=False)
            ctx.metadata["css_latent_method"] = "truncated_svd"
            ctx.metadata["css_latent_dim"] = int(X_latent.shape[1])
        else:
            X_latent = np.asarray(X, dtype=np.float32)
            ctx.metadata["css_latent_method"] = "dense_hvg"
            ctx.metadata["css_latent_dim"] = int(X_latent.shape[1])
        del X

        sample_labels = pd.Series(adata.obs["sample"].astype(str), index=adata.obs_names)
        sample_to_clusters = defaultdict(list)
        centroid_blocks = []

        for sample in sample_labels.unique():
            idx = np.where(sample_labels.values == sample)[0]
            if len(idx) == 0:
                continue
            Xi_latent = X_latent[idx]
            n_clusters = max(2, min(12, len(idx) // 5000 if len(idx) >= 5000 else 2))
            if len(idx) < n_clusters:
                n_clusters = max(1, len(idx))
            km = MiniBatchKMeans(n_clusters=n_clusters, random_state=cfg.random_state, batch_size=min(10000, max(1024, len(idx))))
            labels = km.fit_predict(Xi_latent)
            centroids = km.cluster_centers_.astype(np.float32, copy=False)
            for c in range(centroids.shape[0]):
                sample_to_clusters[sample].append(len(centroid_blocks) + c)
            centroid_blocks.append(centroids)

        centroids = np.vstack(centroid_blocks).astype(np.float32, copy=False)
        centroid_norm = np.linalg.norm(centroids, axis=1, keepdims=True) + 1e-8
        centroids = centroids / centroid_norm

        chunk_size = 10000 if adata.n_obs > 500000 else 20000
        css = np.empty((adata.n_obs, centroids.shape[0]), dtype=np.float32)
        for start in range(0, adata.n_obs, chunk_size):
            stop = min(start + chunk_size, adata.n_obs)
            chunk = np.asarray(X_latent[start:stop], dtype=np.float32)
            chunk_norm = np.linalg.norm(chunk, axis=1, keepdims=True) + 1e-8
            chunk = chunk / chunk_norm
            sims = (chunk @ centroids.T).astype(np.float32, copy=False)
            for sample, cols in sample_to_clusters.items():
                block = sims[:, cols]
                if block.shape[1] == 1:
                    block[:] = 0.0
                else:
                    mu = block.mean(axis=1, keepdims=True)
                    sd = block.std(axis=1, keepdims=True) + 1e-6
                    block[:] = (block - mu) / sd
            css[start:stop, :] = sims
            del chunk, sims
            gc.collect()

        adata.obsm["X_css"] = css
        css_n_pcs = min(cfg.n_pcs, X_latent.shape[1], 50)
        adata.obsm["X_pca"] = X_latent[:, :css_n_pcs].copy() if X_latent.shape[1] >= css_n_pcs else X_latent.copy()
        sc.pp.neighbors(adata, n_neighbors=cfg.n_neighbors, use_rep="X_css", method="umap")
        sc.tl.umap(adata, random_state=cfg.random_state)
        sc.tl.leiden(adata, resolution=cfg.leiden_resolution, flavor="igraph", directed=False, random_state=cfg.random_state)
        ctx.metadata["css_n_reference_clusters"] = int(css.shape[1])

    def _run_hybrid_gpu_graph(self, adata, cfg, ctx) -> bool:
        """Use CPU PCA with GPU neighbors/UMAP/Leiden when GPU PCA is unstable."""
        try:
            import rapids_singlecell as rsc
        except Exception as exc:
            ctx.metadata["clustering_hybrid_skip_reason"] = f"rapids_singlecell_import_failed: {exc}"
            return False

        try:
            sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
            sc.pp.log1p(adata)
            if adata.raw is None and ctx.cfg.scale_mode != "massive":
                adata.raw = adata
            is_sparse = sparse.issparse(adata.X)
            adata.X = self._materialize_matrix(adata.X)
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
            rsc.pp.neighbors(adata, n_neighbors=cfg.n_neighbors, n_pcs=cfg.n_pcs, use_rep="X_pca")
            rsc.tl.umap(adata, random_state=cfg.random_state)
            rsc.tl.leiden(adata, resolution=cfg.leiden_resolution, random_state=cfg.random_state)
            ctx.metadata["clustering_hybrid_reason"] = "gpu_pca_unstable_cpu_pca_gpu_graph"
            return True
        except Exception as exc:
            ctx.metadata["clustering_hybrid_error"] = str(exc)
            return False

    def _run_cpu(self, adata, cfg, ctx) -> None:
        """Standard scanpy CPU clustering pipeline."""
        sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
        sc.pp.log1p(adata)
        if adata.raw is None and ctx.cfg.scale_mode != "massive":
            # Preserve normalized/log-transformed matrix for marker-level downstream tasks.
            adata.raw = adata
        is_sparse = sparse.issparse(adata.X)
        adata.X = self._materialize_matrix(adata.X)
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
        if adata.raw is None and ctx.cfg.scale_mode != "massive":
            adata.raw = adata
        adata.X = self._materialize_matrix(adata.X)
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
