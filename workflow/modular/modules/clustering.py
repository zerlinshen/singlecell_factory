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
from ._gpu_utils import bind_cuda_context, gpu_available
from .._neighbors_cache import get_or_compute_neighbors


# Transient stash of pre-mutation counts for the GPU-failure fallback. Kept out
# of adata.raw so that adata.raw means exactly one thing in every lane:
# normalized + log1p expression, the scanpy convention that trajectory.py and
# tumor_microenvironment.py rely on. Deleted once clustering finishes.
_GPU_RESTORE_LAYER = "_gpu_restore_counts"


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
    def _run_hvg(adata, cfg, ctx) -> None:
        """Run sc.pp.highly_variable_genes honoring cfg.hvg_flavor.

        Default flavor is "seurat" (GT-validated on the suite's hgmm species-mixing
        benchmark: seurat ARI 0.293 >= seurat_v3 at res<=0.5; real-run 2026-06-27).
        "seurat_v3" (Stuart 2019; Heumos 2023) selects HVGs on RAW COUNTS via a
        variance-stabilizing transform, so it MUST run on a raw-count source. By HVG
        time adata.X is already normalized+log1p'd, so seurat_v3 uses layer="counts"
        when present. If seurat_v3 is requested but no raw-count layer is available,
        fall back to "seurat" on the log data and record a loud hvg_flavor_fallback
        status — never silently feed log data to seurat_v3 (Principle 9: no silent
        compromise). batch_key is passed when a batch column with >1 level exists so
        HVG selection is batch-aware.
        """
        flavor = str(getattr(cfg, "hvg_flavor", "seurat"))
        kwargs = {"n_top_genes": cfg.n_top_genes}

        # Batch-aware HVG: only when the run has a batch column with >1 level.
        batch_key = getattr(getattr(ctx.cfg, "batch", None), "batch_key", None)
        if (
            batch_key
            and batch_key in adata.obs.columns
            and adata.obs[batch_key].nunique() > 1
        ):
            kwargs["batch_key"] = batch_key
            ctx.metadata["hvg_batch_key"] = batch_key

        if flavor == "seurat_v3":
            raw_layer = "counts" if "counts" in getattr(adata, "layers", {}) else None
            if raw_layer is None:
                # No raw-count source -> do NOT run seurat_v3 on log data.
                ctx.metadata["hvg_flavor_requested"] = "seurat_v3"
                ctx.metadata["hvg_flavor_actually_used"] = "seurat"
                ctx.metadata["hvg_flavor_fallback"] = (
                    "seurat_v3 requested but no raw-count layer ('counts') available "
                    "at HVG time; fell back to 'seurat' on log data."
                )
                ctx.status("clustering", "warning", ctx.metadata["hvg_flavor_fallback"])
                sc.pp.highly_variable_genes(adata, flavor="seurat", **kwargs)
                return
            sc.pp.highly_variable_genes(adata, flavor="seurat_v3", layer=raw_layer, **kwargs)
            ctx.metadata["hvg_flavor_actually_used"] = "seurat_v3"
            ctx.metadata["hvg_seurat_v3_layer"] = raw_layer
            return

        sc.pp.highly_variable_genes(adata, flavor=flavor, **kwargs)
        ctx.metadata["hvg_flavor_actually_used"] = flavor

    @staticmethod
    def _rapids_leiden_n_iterations(cfg) -> int:
        """Resolve n_iterations for the rapids/cuGraph Leiden lane.

        CPU scanpy/igraph interprets n_iterations=-1 as "run to convergence"
        (Traag 2019), but cuGraph's Leiden requires a positive max_iter (passing
        -1 raises "can't convert negative value to size_t"). Map -1 to cuGraph's
        own default (100), which is effectively convergence for typical graphs, so
        the GPU lane still tracks the convergence intent without crashing.
        """
        n = int(getattr(cfg, "leiden_n_iterations", -1))
        return 100 if n < 0 else n

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
            self._write_resolution_sweep_audit(adata, ctx, cfg)
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
            # M2 invariant — the pre-mutation counts MUST be preserved so the
            # GPU-failure fallback can restore them.
            #
            # This used to overload `adata.raw` for that stash, which made
            # `adata.raw` MEAN DIFFERENT THINGS depending on which lane ran: raw
            # counts after a GPU run, but normalized+log1p data after any CPU lane
            # (which assigns adata.raw straight after log1p). Downstream consumers
            # read adata.raw as *expression* — trajectory.py correlates it against
            # pseudotime, and tumor_microenvironment.py computes the CYT score as
            # the geometric mean of GZMA/PRF1 (Rooney et al. 2015, *Cell*
            # 160:48-61, defined on normalized expression). So the same input
            # produced different biology depending on the host's GPU.
            #
            # The stash is now a dedicated counts layer, and adata.raw is assigned
            # inside _run_gpu at the same point the CPU lanes assign it.
            if _GPU_RESTORE_LAYER not in adata.layers:
                adata.layers[_GPU_RESTORE_LAYER] = adata.X.copy()
                ctx.metadata["m2_restore_counts_preserved"] = True
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
                    if _GPU_RESTORE_LAYER not in adata.layers:
                        poison_adata(adata, "GPU failed + counts stash missing under restore-cpu")
                        raise ClusteringContractViolation(
                            f"GPU clustering failed under policy=restore-cpu but layer "
                            f"{_GPU_RESTORE_LAYER!r} is missing — cannot restore counts. "
                            "adata poisoned."
                        ) from exc
                    logger.info("Restoring counts from preserved stash before CPU fallback")
                    # _run_gpu mutates X in place (normalize/log1p/scale), so the
                    # CPU lane must start from the pre-mutation counts or it would
                    # normalize already-normalized data.
                    adata.X = adata.layers[_GPU_RESTORE_LAYER].copy()
                    adata.raw = None
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

        self._write_resolution_sweep_audit(adata, ctx, cfg)
        self._plot_umap_clusters(adata, ctx)
        ctx.adata = adata
        ctx.metadata["n_clusters"] = int(adata.obs["leiden"].nunique())
        ctx.metadata["gpu_mode"] = ctx.cfg.gpu_mode
        ctx.metadata["clustering_random_state"] = ctx.random_state
        if "clustering_backend" not in ctx.metadata:
            ctx.metadata["clustering_backend"] = "gpu" if use_gpu else "cpu"
        # Drop the transient counts stash so it never reaches final_adata.h5ad,
        # and state what adata.raw holds so downstream consumers (and readers of
        # the manifest) never have to infer it from the backend that ran.
        if _GPU_RESTORE_LAYER in adata.layers:
            del adata.layers[_GPU_RESTORE_LAYER]
        ctx.metadata["raw_semantics"] = (
            "normalized_log1p_all_genes" if adata.raw is not None else "raw_absent"
        )
        self._record_batch_confounding_risk(adata, ctx)

    # obs columns that, with more than one level, indicate the run spans
    # multiple technical units and is therefore exposed to batch effects.
    _BATCH_CANDIDATE_COLUMNS = (
        "batch", "sample", "sample_id", "donor_id", "donor",
        "dataset", "study", "patient", "patient_id", "platform", "assay",
    )

    def _record_batch_confounding_risk(self, adata, ctx) -> None:
        """Flag clusters produced from multi-batch input with no integration.

        Without integration, clusters on multi-study/multi-platform input track the
        batch rather than the biology. On the LuCA LUSC cohort (8 studies, 5
        platforms, 87 samples) the unintegrated default scored iLISI 0.013 / kBET
        0.115 and split 24 published cell types into 47 clusters; enabling Harmony
        on the same input moved that to iLISI 0.051 / kBET 0.552 and exactly 24
        clusters, with ARI against the published labels rising 0.214 -> 0.455
        (governance/realrun_gt_concordance_lusc_2026-07-28.md).

        This does not change the analysis — it records, in the manifest, that a
        known confounder was present and unaddressed, so a batch-driven clustering
        can never be reported as an integrated one.
        """
        if "batch_correction" in (getattr(ctx.cfg, "optional_modules", None) or ()):
            return
        candidates = {
            col: int(adata.obs[col].nunique(dropna=True))
            for col in self._BATCH_CANDIDATE_COLUMNS
            if col in adata.obs.columns and adata.obs[col].nunique(dropna=True) > 1
        }
        if not candidates:
            return
        ctx.metadata["batch_confounding_risk"] = {
            "status": "multi_batch_input_without_integration",
            "candidate_batch_columns": candidates,
            "integration_ran": False,
        }
        logger.warning(
            "BATCH_CONFOUNDING_RISK: clustering ran with no batch correction, but the "
            "input spans multiple technical units (%s). Clusters may separate by batch "
            "rather than by biology; run with --optional-modules ...,batch_correction "
            "--batch-key <col> to integrate. Recorded in the manifest as "
            "batch_confounding_risk.",
            ", ".join(f"{k}={v}" for k, v in sorted(candidates.items())),
        )

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
            bind_cuda_context()  # ensure cuBLAS/cuSOLVER pre-warm before GPU ops (P1 fix; no-op if already warmed in probe)
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
            # HVG before scaling — see _run_cpu for the rationale.
            self._run_hvg(adata, cfg, ctx)
            if cfg.scale_data:
                sc.pp.scale(adata, max_value=10, zero_center=not sparse.issparse(adata.X))
            sc.tl.pca(
                adata,
                n_comps=cfg.n_pcs,
                svd_solver="arpack" if is_sparse else "randomized",
                mask_var="highly_variable",
                random_state=cfg.random_state,
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
            rsc.tl.leiden(
                adata,
                resolution=cfg.leiden_resolution,
                random_state=cfg.random_state,
                n_iterations=self._rapids_leiden_n_iterations(cfg),  # -1 -> cuGraph default 100; CPU/GPU parity (Traag 2019)
            )
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
        # HVG before scaling — see _run_cpu for the rationale.
        self._run_hvg(adata, cfg, ctx)
        if cfg.scale_data:
            # zero_center=False to preserve sparsity (scanpy.pp.scale contract).
            sc.pp.scale(adata, max_value=10, zero_center=False)
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
            n_iterations=int(getattr(cfg, "leiden_n_iterations", -1)),  # -1 = converge (Traag 2019)
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
        # HVG must be selected on normalized+log1p data, BEFORE scaling. The
        # "seurat" flavor bins genes by mean expression and ranks them by
        # normalized dispersion; sc.pp.scale sets every gene to mean 0 / variance
        # 1 (and clips at max_value), destroying exactly the mean-variance
        # relationship the selection depends on and making the ranking
        # degenerate. Canonical order is normalize -> log1p -> HVG -> scale
        # (Luecken & Theis 2019, *Mol Syst Biol* 15:e8746; Heumos et al. 2023,
        # *Nat Rev Genet* 24:550-572), and _run_hvg's own contract states that
        # adata.X is normalized+log1p at call time.
        self._run_hvg(adata, cfg, ctx)
        if cfg.scale_data:
            sc.pp.scale(adata, max_value=10, zero_center=not sparse.issparse(adata.X))
        # Both solvers are seeded: "randomized" draws a random projection and
        # ARPACK a random start vector, so without an explicit random_state PCA
        # silently used scanpy's own default (0) instead of the run's
        # --random-state, while the manifest recorded clustering_random_state as
        # if the seed had been honoured. Only _run_sparse_exact passed it.
        sc.tl.pca(
            adata,
            n_comps=cfg.n_pcs,
            svd_solver="arpack" if is_sparse else "randomized",
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
            n_iterations=int(getattr(cfg, "leiden_n_iterations", -1)),  # -1 = converge (Traag 2019)
        )

    def _run_gpu(self, adata, cfg, ctx) -> None:
        """GPU-accelerated clustering via rapids-singlecell."""
        import rapids_singlecell as rsc
        bind_cuda_context()  # rebind CUDA context poisoned by rapids import (P1 fix)

        sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
        sc.pp.log1p(adata)
        self._log_rss(ctx, "hotspot1:_run_gpu:before_preserve_raw")
        # Assign adata.raw at the SAME point the CPU lanes do, so adata.raw is
        # normalized+log1p expression regardless of which lane ran. The
        # pre-mutation counts the GPU fallback needs live in
        # _GPU_RESTORE_LAYER, not here.
        if adata.raw is None and self._should_preserve_raw(ctx):
            adata.raw = adata
        self._log_rss(ctx, "hotspot1:_run_gpu:after_preserve_raw")
        adata.X = self._materialize_matrix(adata.X)
        self._log_rss(ctx, "hotspot1:_run_gpu:after_materialize_X")
        # HVG before scaling — see _run_cpu for the rationale.
        self._run_hvg(adata, cfg, ctx)
        if cfg.scale_data:
            sc.pp.scale(adata, max_value=10, zero_center=not sparse.issparse(adata.X))

        # GPU-accelerated PCA, neighbors, UMAP
        rsc.pp.pca(
            adata, n_comps=cfg.n_pcs, mask_var="highly_variable",
            random_state=cfg.random_state,
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
        rsc.tl.leiden(
            adata,
            resolution=cfg.leiden_resolution,
            random_state=cfg.random_state,
            n_iterations=self._rapids_leiden_n_iterations(cfg),  # -1 -> cuGraph default 100; CPU/GPU parity (Traag 2019)
        )

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
        plt.savefig(ctx.figure_dir / "pca_variance_explained.png", bbox_inches="tight")
        plt.close()

        # Record how many PCs needed for 90% variance
        pcs_for_90 = int(np.searchsorted(cumulative, 0.9) + 1)
        ctx.metadata["pca_cumulative_var_90pct"] = min(pcs_for_90, n_pcs)

    @staticmethod
    def _write_resolution_sweep_audit(adata, ctx: PipelineContext, cfg) -> None:
        """Write optional Leiden resolution diagnostics without changing final labels."""
        sweep = tuple(getattr(cfg, "leiden_resolution_sweep", ()) or ())
        if not sweep:
            ctx.metadata["leiden_resolution_sweep_status"] = "skipped_not_requested"
            return
        if "neighbors" not in adata.uns:
            ctx.metadata["leiden_resolution_sweep_status"] = "skipped_missing_neighbors"
            return

        original = adata.obs["leiden"].copy() if "leiden" in adata.obs else None
        rows = []
        for resolution in sweep:
            key = f"leiden_res_{str(resolution).replace('.', 'p')}"
            try:
                sc.tl.leiden(
                    adata,
                    resolution=float(resolution),
                    key_added=key,
                    flavor="igraph",
                    directed=False,
                    random_state=cfg.random_state,
                )
            except TypeError:
                sc.tl.leiden(
                    adata,
                    resolution=float(resolution),
                    flavor="igraph",
                    directed=False,
                    random_state=cfg.random_state,
                )
                adata.obs[key] = adata.obs["leiden"].copy()
                if original is not None:
                    adata.obs["leiden"] = original.copy()
            labels = adata.obs[key].astype(str)
            sizes = labels.value_counts()
            small_cluster_threshold = max(10, int(round(adata.n_obs * 0.005)))
            small_cells = int(sizes[sizes < small_cluster_threshold].sum())
            rows.append(
                {
                    "resolution": float(resolution),
                    "n_clusters": int(sizes.shape[0]),
                    "min_cluster_size": int(sizes.min()),
                    "median_cluster_size": float(sizes.median()),
                    "max_cluster_size": int(sizes.max()),
                    "small_cluster_threshold_cells": int(small_cluster_threshold),
                    "small_cluster_cell_pct": round(
                        float(small_cells) / max(int(adata.n_obs), 1) * 100.0,
                        3,
                    ),
                }
            )

        if original is not None:
            adata.obs["leiden"] = original
        pd.DataFrame(rows).to_csv(ctx.table_dir / "leiden_resolution_sweep.csv", index=False)
        ctx.metadata["leiden_resolution_sweep_status"] = "completed"
        ctx.metadata["leiden_resolution_sweep_values"] = [float(v) for v in sweep]

    @staticmethod
    def _plot_umap_clusters(adata, ctx: PipelineContext) -> None:
        """Draw the cluster UMAP at journal geometry with a distinguishable palette.

        rcParams alone cannot fix this panel: ``sc.pl.umap`` supplies its own
        categorical palette, its own figure size and its own title, so a themed
        run still produced scanpy's default colours at an arbitrary aspect ratio
        titled with the obs key. At 40+ clusters that palette repeats visually
        similar hues, and on-data labels collide, so clusters cannot be told
        apart — the same class of defect as the plotting layer's silent colour
        wrap. Size to the double-column width, request distinct colours, and
        state n on the panel.
        """
        from .._figure_theme import journal_figure_size, qualitative_colors

        n_clusters = int(adata.obs["leiden"].nunique())
        figsize = journal_figure_size("double", height_mm=110.0)
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        sc.pl.umap(
            adata,
            color=["leiden"],
            palette=qualitative_colors("leiden_qualitative", n_clusters),
            show=False,
            legend_loc="on data",
            legend_fontsize=5,
            legend_fontoutline=1,   # keep labels legible over dense point clouds
            frameon=False,          # UMAP axes carry no units; a box implies scale
            size=3,
            ax=ax,
        )
        ax.set_title(
            f"Leiden clusters (n = {adata.n_obs:,} cells, {n_clusters} clusters, "
            f"resolution {ctx.cfg.clustering.leiden_resolution:g})"
        )
        fig.savefig(ctx.figure_dir / "umap_leiden.png")
        plt.close(fig)
