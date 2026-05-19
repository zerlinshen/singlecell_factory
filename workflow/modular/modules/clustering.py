from __future__ import annotations

import gc
import logging
import os
import matplotlib
import numpy as np
import pandas as pd
from scipy import sparse

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ._scanpy_compat import import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext
from ._gpu_utils import gpu_available
from .._neighbors_cache import get_or_compute_neighbors


__references__ = {
    "scanpy": {
        "title": "SCANPY: large-scale single-cell gene expression data analysis",
        "authors": "Wolf, Angerer, Theis",
        "journal": "Genome Biology",
        "year": "2018",
        "doi": "10.1186/s13059-017-1382-0",
        "description": "scanpy pp.pca / pp.neighbors / tl.umap / tl.leiden stack.",
    },
    "McInnes_UMAP_2018": {
        "title": "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction",
        "authors": "McInnes, Healy, Melville",
        "journal": "arXiv preprint",
        "year": "2018",
        "doi": "arXiv:1802.03426",
        "description": "KNN-graph + UMAP layout used by scanpy.pp.neighbors / sc.tl.umap.",
    },
    "Traag_Leiden_2019": {
        "title": "From Louvain to Leiden: guaranteeing well-connected communities",
        "authors": "Traag, Waltman, van Eck",
        "journal": "Scientific Reports",
        "year": "2019",
        "doi": "10.1038/s41598-019-41695-z",
        "description": "Leiden community detection used as the clustering algorithm.",
    },
    "rapids_singlecell": {
        "title": "rapids-singlecell \u2014 GPU-accelerated scanpy",
        "authors": "scverse contributors",
        "journal": "Software (scverse)",
        "year": "2024",
        "doi": "https://github.com/scverse/rapids_singlecell",
        "description": "GPU acceleration path. Falls back to CPU scanpy when unavailable.",
    },
}


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

    @staticmethod
    def _clone_for_gpu_lite(adata):
        """US-007 Hotspot 1 mitigation: GPU clone with smaller peak RSS than adata.copy().

        Deep-copies X (log1p/scale mutate in place) and var (HVG mutates columns);
        shares obs/obsm/varm/uns/layers by reference. Design + safety predicate
        in docs/HOTSPOT1_DIAGNOSIS.md "US-007 mitigation".
        """
        import anndata as ad
        X = adata.X
        if X is None:
            X_clone = None
        elif sparse.issparse(X):
            X_clone = X.copy()
        else:
            X_clone = np.array(X, copy=True)
        return ad.AnnData(
            X=X_clone,
            obs=adata.obs,
            var=adata.var.copy(),
            obsm=dict(adata.obsm) if adata.obsm is not None else None,
            varm=dict(adata.varm) if adata.varm is not None else None,
            uns=dict(adata.uns) if adata.uns is not None else None,
            layers=dict(adata.layers) if adata.layers is not None else None,
        )

    @staticmethod
    def _log_rss(ctx, tag: str) -> None:
        """Log peak resident set size (MB) at tag points for Hotspot 1 diagnosis.

        Audit doc docs/SCIENTIFIC_AUDIT_2026-05-15.md identifies the GPU clustering
        path as the Wave 1 memory hotspot. This helper captures `ru_maxrss` (peak
        RSS since process start, in KB on Linux) and writes both to logger and
        ctx.metadata["clustering_rss_trace"]. Purely instrumentation — no
        behavioural change. Used by tests/test_hotspot1_instrumentation.py.
        """
        try:
            import resource
            rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            rss_mb = rss_kb / 1024.0  # Linux ru_maxrss is in KB
            logger.info("[hotspot1-trace] %s peak_rss=%.1fMB", tag, rss_mb)
            trace = ctx.metadata.setdefault("clustering_rss_trace", [])
            trace.append({"tag": tag, "peak_rss_mb": round(rss_mb, 1)})
        except Exception as exc:  # pragma: no cover — instrumentation must never break clustering
            logger.debug("_log_rss(%s) failed: %s", tag, exc)

    """Optional module: normalization, PCA/UMAP and Leiden clustering.

    Supports GPU acceleration via rapids-singlecell when available.
    Falls back to standard scanpy CPU backend automatically.
    """

    name = "clustering"
    provides_keys = {"obs": ["leiden"], "obsm": ["X_pca", "X_umap"]}

    def run(self, ctx: PipelineContext) -> None:
        if os.environ.get("SC_MEM_GUARD", "").lower() == "on":
            from .._mem_guard import MemoryGuard, MemoryGuardError
            try:
                with MemoryGuard(ctx, self.name) as mg:
                    mg.check("entry")
                    return self._run_impl(ctx)
            except MemoryGuardError as exc:
                logger.warning("MemoryGuard: %s", exc)
                ctx.metadata.setdefault("mem_warnings", []).append(
                    {"module": self.name, "error": str(exc)}
                )
        return self._run_impl(ctx)

    def _run_impl(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Clustering requires AnnData.")
        cfg = ctx.cfg.clustering
        # Canonical seed: ctx.random_state overrides per-module ClusteringConfig.random_state
        cfg.random_state = ctx.random_state
        sc.settings.n_jobs = max(1, os.cpu_count() or 1)

        # F-1: CSS removed from production path. Invoke _should_use_css for its
        # side effects (raises on engine=="css", sets css_status/clustering_engine
        # metadata). Return value is always False; the prior dispatch branch is gone.
        self._should_use_css(ctx, adata)
        # F-2 (Plan ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md, Principle 7):
        # `--clustering-engine=sparse_exact` is now a dedicated dispatch arm with
        # strict input validation (raises if input is not sparse) and forced
        # `svd_solver="arpack"` (no silent fallback to randomized SVD).
        # sparse_exact is the CPU credibility-reference lane for Component C's
        # benchmark; it does NOT route to GPU even when GPU is available.
        _engine_pre_gpu = (
            os.environ.get("SC_CLUSTERING_ENGINE", "").strip().lower()
            or getattr(ctx.cfg, "clustering_engine", "auto")
        )
        if _engine_pre_gpu == "sparse_exact":
            self._run_sparse_exact(adata, cfg, ctx)
            use_gpu = False
            self._plot_umap_clusters(adata, ctx)
            ctx.adata = adata
            ctx.metadata["n_clusters"] = int(adata.obs["leiden"].nunique())
            ctx.metadata["gpu_mode"] = ctx.cfg.gpu_mode
            ctx.metadata["clustering_random_state"] = ctx.random_state
            ctx.metadata["clustering_backend"] = "cpu_sparse_exact"
            return
        use_gpu = gpu_available(ctx.cfg.gpu_mode)
        engine = _engine_pre_gpu
        _cp_policy = (
            os.environ.get("SC_CHECKPOINT_POLICY", "").strip().lower()
            or getattr(ctx.cfg, "checkpoint_policy", "full")
        )
        if engine not in {"gpu"} and ctx.cfg.gpu_mode == "auto" and _cp_policy != "full":
            logger.info("Non-full checkpoint policy: prefer CPU clustering to avoid GPU copy overhead")
            ctx.metadata["clustering_gpu_disabled_reason"] = "checkpoint_policy_non_full"
            use_gpu = False
        if use_gpu:
            logger.info("GPU detected — using rapids-singlecell for PCA/neighbors/UMAP/Leiden")
            # Wave 3 / US-W3-1 — M2 mechanism (selected over M1 in US-W3-0 spike
            # due to RTX 5090 D v2 + cuSOLVER incompatibility on rsc.pp.pca; see
            # docs/HOTSPOT1_DIAGNOSIS.md). Site-1 elimination: skip the host-side
            # _clone_for_gpu_lite and run _run_gpu directly on the original adata.
            # Mandate adata.raw preservation so the GPU-failure fallback can
            # restore counts via adata.raw.to_adata() — safe even though
            # log1p/scale mutate adata.X in place (scanpy/_simple.py:364,
            # scanpy/_scale.py:200), because raw holds an immutable reference
            # to the pre-mutation counts.
            from .._contract_violation import (
                ClusteringContractViolation,
                poison_adata,
                resolve_gpu_failure_policy,
            )
            gpu_failure_policy = resolve_gpu_failure_policy(ctx.cfg)
            ctx.metadata["clustering_clone_strategy"] = "m2_inplace_raw_mandate"
            ctx.metadata["gpu_failure_policy"] = gpu_failure_policy
            self._log_rss(ctx, "hotspot1:before_inplace_gpu")
            # M2 invariant — raw MUST be preserved for the GPU path to enable
            # restore-on-fallback. Override the checkpoint-policy gate.
            if adata.raw is None:
                adata.raw = adata.copy()
                ctx.metadata["m2_raw_preserved"] = True
            self._log_rss(ctx, "hotspot1:after_raw_preserve")
            try:
                self._run_gpu(adata, cfg, ctx)
                self._log_rss(ctx, "hotspot1:after_inplace_gpu_success")
                gc.collect()
            except Exception as exc:
                if ctx.cfg.gpu_mode == "force":
                    raise RuntimeError(f"GPU clustering failed in force mode: {exc}") from exc
                logger.warning("GPU clustering failed (%s); applying policy=%s", exc, gpu_failure_policy)
                ctx.metadata["clustering_gpu_fallback_reason"] = str(exc)
                if gpu_failure_policy == "raise":
                    poison_adata(adata, f"GPU clustering failure: {exc}")
                    raise ClusteringContractViolation(
                        f"GPU clustering failed and SC_GPU_FAILURE_POLICY=raise; "
                        f"adata poisoned. Original error: {exc}"
                    ) from exc
                if gpu_failure_policy == "restore-cpu":
                    if adata.raw is None:
                        poison_adata(adata, "GPU failed + raw missing under restore-cpu")
                        raise ClusteringContractViolation(
                            "GPU clustering failed under policy=restore-cpu but adata.raw "
                            "is None — cannot restore counts. adata poisoned."
                        ) from exc
                    logger.info("Restoring adata from preserved raw before CPU fallback")
                    adata = adata.raw.to_adata()
                    ctx.adata = adata
                    ctx.metadata["clustering_restored_from_raw"] = True
                    use_gpu = False
                    if self._run_hybrid_gpu_graph(adata, cfg, ctx):
                        ctx.metadata["clustering_backend"] = "hybrid"
                    else:
                        self._run_cpu(adata, cfg, ctx)
                elif gpu_failure_policy == "reload-checkpoint":
                    checkpoint_path = getattr(ctx, "last_checkpoint_path", None)
                    if checkpoint_path is None or not checkpoint_path.exists():
                        poison_adata(adata, "GPU failed + no checkpoint under reload-checkpoint")
                        raise ClusteringContractViolation(
                            "GPU clustering failed under policy=reload-checkpoint but no "
                            "pre-clustering checkpoint exists. Re-run with --checkpoint. "
                            "adata poisoned."
                        ) from exc
                    import anndata as ad
                    logger.info("Reloading adata from %s before CPU fallback", checkpoint_path)
                    adata = ad.read_h5ad(str(checkpoint_path))
                    ctx.adata = adata
                    ctx.metadata["clustering_reloaded_from_checkpoint"] = str(checkpoint_path)
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
        ctx.metadata["clustering_random_state"] = ctx.random_state
        if "clustering_backend" not in ctx.metadata:
            ctx.metadata["clustering_backend"] = "gpu" if use_gpu else "cpu"

    def _should_use_css(self, ctx, adata) -> bool:
        # F-1 (Plan ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md, Principle 2):
        # CSS clustering is removed from the production science path. This gate
        # never returns True. Explicit CSS requests raise loudly; the prior
        # scale_mode=massive auto-activation branch is gone (config.py F-4 remap
        # makes massive route to clustering_engine=auto, but defense-in-depth here).
        engine = (
            os.environ.get("SC_CLUSTERING_ENGINE", "").strip().lower()
            or getattr(ctx.cfg, "clustering_engine", "auto")
        )
        if engine == "css":
            raise RuntimeError(
                "CSS clustering is removed from the production science path "
                "(Plan F-1, Principle 2). To benchmark CSS for historical "
                "comparison, check out a pre-2026-05-19 commit. "
                "See ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md"
            )
        ctx.metadata["css_status"] = "removed_per_plan_F1"
        if engine == "sparse_exact":
            ctx.metadata["clustering_engine"] = "sparse_exact"
        return False

    def _run_css(self, adata, cfg, ctx) -> None:
        # F-1 tombstone (Plan ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md, Principle 2):
        # CSS is removed from the production science path. The legacy implementation
        # body was physically deleted on 2026-05-20. This stub remains as a
        # documented gravestone — if any callsite still routes here (dispatcher
        # regression, copy-pasted old patch, etc.) it raises loudly instead of
        # silently no-op'ing or AttributeError'ing.
        raise RuntimeError(
            "_run_css is unreachable: CSS clustering is removed from the "
            "production science path (Plan F-1, Principle 2). "
            "See ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md"
        )

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
            if adata.raw is None and self._should_preserve_raw(ctx):
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
            get_or_compute_neighbors(
                adata,
                backend_id="rapids",
                method="umap",
                metric="euclidean",
                n_pcs=cfg.n_pcs,
                n_neighbors=cfg.n_neighbors,
                use_rep="X_pca",
                knn=True,
                random_state=cfg.random_state,
                compute_fn=lambda: rsc.pp.neighbors(adata, n_neighbors=cfg.n_neighbors, n_pcs=cfg.n_pcs, use_rep="X_pca"),
            )
            rsc.tl.umap(adata, random_state=cfg.random_state)
            rsc.tl.leiden(adata, resolution=cfg.leiden_resolution, random_state=cfg.random_state)
            ctx.metadata["clustering_hybrid_reason"] = "gpu_pca_unstable_cpu_pca_gpu_graph"
            return True
        except Exception as exc:
            ctx.metadata["clustering_hybrid_error"] = str(exc)
            return False

    @staticmethod
    def _should_preserve_raw(ctx) -> bool:
        """Return True when it is safe/desired to preserve adata.raw.

        Evaluated to preserve adata.raw when checkpoint policy permits.
        """
        policy = (
            os.environ.get("SC_CHECKPOINT_POLICY", "").strip().lower()
            or getattr(ctx.cfg, "checkpoint_policy", "full")
        )
        return policy == "full"

    def _run_sparse_exact(self, adata, cfg, ctx) -> None:
        """F-2 dedicated sparse-exact CPU clustering lane (Principle 7).

        Contract:
        - Input MUST be sparse on entry. If not, raise RuntimeError loudly —
          there is no silent dense fallback. The lane's name promises exactness
          on sparse data; we honor it at the gate rather than degrade to
          ``svd_solver="randomized"`` and pretend it is still "exact".
        - PCA solver is forced to ``arpack`` (scipy ARPACK iterative
          Lanczos-bidiagonalization eigensolver; the scanpy/scientific-Python
          canonical sparse exact SVD method, peer-reviewed since Lehoucq et al.
          1998). No randomized SVD on the principal stage.
        - Neighbors use scanpy ``sc.pp.neighbors`` (true KNN via pynndescent /
          UMAP), Leiden uses ``sc.tl.leiden(flavor="igraph", directed=False)``
          on the true graph — Wolf 2018 + Traag 2019 reference path.
        - Full provenance is recorded in ``ctx.metadata`` so the
          ``lane_manifest.json`` contract (Plan G-C0) can attest the lane
          identity without inspecting executor code.
        """
        is_sparse_in = sparse.issparse(adata.X)
        if not is_sparse_in:
            raise RuntimeError(
                "_run_sparse_exact requires a sparse input matrix; received "
                f"{type(adata.X).__name__}. Either pass a sparse AnnData or "
                "select --clustering-engine=auto to use the standard dispatch. "
                "Silent densification + randomized SVD is forbidden under "
                "Plan F-2 / Principle 7 (no no-op functionality). "
                "See ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md"
            )

        sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
        sc.pp.log1p(adata)
        if adata.raw is None and self._should_preserve_raw(ctx):
            adata.raw = adata
        # F-2: keep sparse — _materialize_matrix is sparse-safe for sparse inputs.
        adata.X = self._materialize_matrix(adata.X)
        if not sparse.issparse(adata.X):
            raise RuntimeError(
                "_run_sparse_exact lost sparsity after _materialize_matrix; "
                "input was sparse on entry but is dense after materialization. "
                "This indicates a contract regression in _materialize_matrix."
            )
        if cfg.scale_data:
            # zero_center=False to preserve sparsity (scanpy.pp.scale contract).
            sc.pp.scale(adata, max_value=10, zero_center=False)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=cfg.n_top_genes)
        sc.tl.pca(
            adata,
            n_comps=cfg.n_pcs,
            svd_solver="arpack",  # F-2: forced; no randomized fallback.
            mask_var="highly_variable",
            random_state=cfg.random_state,
        )
        self._plot_pca_variance(adata, ctx, cfg.n_pcs)
        get_or_compute_neighbors(
            adata,
            backend_id="scanpy",
            method="umap",
            metric="euclidean",
            n_pcs=cfg.n_pcs,
            n_neighbors=cfg.n_neighbors,
            use_rep="X_pca",
            knn=True,
            random_state=cfg.random_state,
            compute_fn=lambda: sc.pp.neighbors(
                adata,
                n_neighbors=cfg.n_neighbors,
                n_pcs=cfg.n_pcs,
                use_rep="X_pca",
                method="umap",
            ),
        )
        sc.tl.umap(adata, random_state=cfg.random_state)
        sc.tl.leiden(
            adata,
            resolution=cfg.leiden_resolution,
            flavor="igraph",
            directed=False,
            random_state=cfg.random_state,
        )

        # F-2 provenance — lane_manifest.json (G-C0) reads these.
        ctx.metadata["clustering_engine"] = "sparse_exact"
        ctx.metadata["clustering_lane"] = "scanpy_cpu_sparse_exact"
        ctx.metadata["pca_solver"] = "arpack"
        ctx.metadata["pca_input_was_sparse"] = True
        ctx.metadata["neighbors_backend"] = "scanpy_umap_true_knn"
        ctx.metadata["leiden_backend"] = "igraph_cpu"
        ctx.metadata["leiden_flavor"] = "igraph"
        ctx.metadata["leiden_directed"] = False
        ctx.metadata["sparse_exact_principle_8_compatible"] = True

    def _run_cpu(self, adata, cfg, ctx) -> None:
        """Standard scanpy CPU clustering pipeline."""
        sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
        sc.pp.log1p(adata)
        if adata.raw is None and self._should_preserve_raw(ctx):
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
        get_or_compute_neighbors(
            adata,
            backend_id="scanpy",
            method="umap",
            metric="euclidean",
            n_pcs=cfg.n_pcs,
            n_neighbors=cfg.n_neighbors,
            use_rep="X_pca",
            knn=True,
            random_state=cfg.random_state,
            compute_fn=lambda: sc.pp.neighbors(
                adata,
                n_neighbors=cfg.n_neighbors,
                n_pcs=cfg.n_pcs,
                use_rep="X_pca",
                method="umap",
            ),
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
        # M2 mandates adata.raw upstream in _run_impl, so the legacy
        # preserve-raw gate that lived here is unreachable on the M2 path.
        # Tags retained for downstream provenance comparison.
        self._log_rss(ctx, "hotspot1:_run_gpu:before_preserve_raw")
        self._log_rss(ctx, "hotspot1:_run_gpu:after_preserve_raw")
        adata.X = self._materialize_matrix(adata.X)
        self._log_rss(ctx, "hotspot1:_run_gpu:after_materialize_X")
        if cfg.scale_data:
            sc.pp.scale(adata, max_value=10, zero_center=not sparse.issparse(adata.X))
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=cfg.n_top_genes)

        # GPU-accelerated PCA, neighbors, UMAP
        rsc.pp.pca(adata, n_comps=cfg.n_pcs, mask_var="highly_variable")
        self._plot_pca_variance(adata, ctx, cfg.n_pcs)
        get_or_compute_neighbors(
            adata,
            backend_id="rapids",
            method="umap",
            metric="euclidean",
            n_pcs=cfg.n_pcs,
            n_neighbors=cfg.n_neighbors,
            use_rep="X_pca",
            knn=True,
            random_state=cfg.random_state,
            compute_fn=lambda: rsc.pp.neighbors(adata, n_neighbors=cfg.n_neighbors, n_pcs=cfg.n_pcs, use_rep="X_pca"),
        )
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
