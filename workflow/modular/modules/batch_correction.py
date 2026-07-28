from __future__ import annotations

import inspect
import json
import logging
import os
import warnings
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ._scanpy_compat import import_scanpy_or_stub
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

sc = import_scanpy_or_stub()

from ..context import PipelineContext
from ._gpu_utils import bind_cuda_context, gpu_available
from .._neighbors_cache import get_or_compute_neighbors


__references__ = {
    "Korsunsky_Harmony_2019": {
        "title": "Fast, sensitive and accurate integration of single-cell data with Harmony",
        "authors": "Korsunsky et al.",
        "journal": "Nature Methods",
        "year": "2019",
        "doi": "10.1038/s41592-019-0619-0",
        "description": (
            "Reference algorithm. Two backends: CPU harmonypy via "
            "scanpy_external.pp.harmony_integrate, GPU rapids-singlecell via "
            "rsc.pp.harmony_integrate. Initial parity verified on a 50k/4-batch "
            "synthetic fixture (US-B3, ARI >= 0.95). Production-scale parity "
            "verification is part of the C benchmark (G-C0 v3)."
        ),
    },
    "Polanski_BBKNN_2020": {
        "title": "BBKNN: fast batch alignment of single cell transcriptomes",
        "authors": "Polanski et al.",
        "journal": "Bioinformatics",
        "year": "2020",
        "doi": "10.1093/bioinformatics/btz625",
        "description": "BBKNN backend.",
    },
    "Johnson_ComBat_2007": {
        "title": "Adjusting batch effects in microarray expression data using empirical Bayes methods",
        "authors": "Johnson, Li, Rabinovic",
        "journal": "Biostatistics",
        "year": "2007",
        "doi": "10.1093/biostatistics/kxj037",
        "description": "ComBat backend (sc.pp.combat).",
    },
    "Hie_Scanorama_2019": {
        "title": "Efficient integration of heterogeneous single-cell transcriptomes using Scanorama",
        "authors": "Hie, Bryson, Berger",
        "journal": "Nature Biotechnology",
        "year": "2019",
        "doi": "10.1038/s41587-019-0113-3",
        "description": "Scanorama backend.",
    },
    "Lopez_scVI_2018": {
        "title": "Deep generative modeling for single-cell transcriptomics",
        "authors": "Lopez et al.",
        "journal": "Nature Methods",
        "year": "2018",
        "doi": "10.1038/s41592-018-0229-2",
        "description": "scVI backend via scvi-tools.",
    },
    "Haghverdi_MNN_2018": {
        "title": "Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors",
        "authors": "Haghverdi et al.",
        "journal": "Nature Biotechnology",
        "year": "2018",
        "doi": "10.1038/nbt.4091",
        "description": "MNN / fastMNN backends.",
    },
}


logger = logging.getLogger(__name__)


def _to_cpu_value(value):
    if sparse.issparse(value):
        data = getattr(value, "data", None)
        if hasattr(data, "get"):
            value.data = data.get()
        return value
    if hasattr(value, "get"):
        return value.get()
    return value


def _ensure_cpu_batch_inputs(adata) -> None:
    adata.X = _to_cpu_value(adata.X)
    for layer_key in list(adata.layers.keys()):
        adata.layers[layer_key] = _to_cpu_value(adata.layers[layer_key])
    for rep_key in ["X_pca", "X_pca_harmony", "X_pca_combat", "X_scanorama", "X_scvi", "X_mnn", "X_fastmnn"]:
        if rep_key in adata.obsm:
            adata.obsm[rep_key] = _to_cpu_value(adata.obsm[rep_key])


class BatchCorrectionModule:
    """Optional module: batch effect correction for multi-sample integration.

    Supports multiple backends:
    1. Harmony (default) — fast, PCA-based correction
    2. BBKNN — batch-balanced k-nearest-neighbors
    3. Combat — linear model-based correction
    4. Scanorama — manifold alignment integration
    5. scVI — variational latent integration

    Should run after clustering (specifically after PCA), before downstream analysis.
    The batch key must be present in adata.obs.
    """

    name = "batch_correction"
    mutates_structure = True
    requires_keys = {"obs": ["leiden"]}
    provides_keys = {"obs": ["leiden"]}

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Batch correction requires AnnData.")

        cfg = ctx.cfg.batch
        # Sync canonical seed into clustering cfg so post-correction UMAP/Leiden are deterministic
        ctx.cfg.clustering.random_state = ctx.random_state
        batch_key = cfg.batch_key

        if batch_key not in adata.obs.columns:
            ctx.metadata["batch_correction_status"] = "skipped_missing_batch_key"
            ctx.status(
                self.name,
                "skipped",
                f"Batch key '{batch_key}' not found in adata.obs.",
            )
            return

        n_batches = adata.obs[batch_key].nunique()
        if n_batches < 2:
            ctx.metadata["batch_correction_status"] = "skipped_single_batch"
            ctx.status(self.name, "skipped", f"Only one batch found in '{batch_key}'.")
            return

        # The integration_select gate may set method="none" (no integration is
        # the discovery-safe choice when no candidate survives its gates). Honor
        # it as a clean skip rather than an unknown-method error.
        if cfg.method == "none":
            ctx.metadata["batch_correction_status"] = "skipped_method_none"
            ctx.status(
                self.name,
                "skipped",
                "cfg.batch.method='none' (no integration selected); "
                "leaving X_pca as the latent representation.",
            )
            return

        ctx.metadata["n_batches"] = n_batches
        ctx.metadata["batch_method"] = cfg.method
        ctx.metadata["gpu_mode"] = ctx.cfg.gpu_mode

        mixing_audit = {
            "method": cfg.method,
            "batch_key": batch_key,
            "interpretation": (
                "Higher normalized batch entropy and lower same-batch neighbor fraction after "
                "correction support better global batch mixing. Interpret this alongside "
                "cell-type markers; a batch remaining separated within the same cell type "
                "suggests residual batch effect, while different biological cell types should "
                "remain separable."
            ),
            "before": self._compute_batch_mixing_metrics(adata, batch_key, "X_pca"),
        }

        # Pre-correction UMAP (for comparison)
        if "X_umap" in adata.obsm:
            sc.pl.umap(adata, color=[batch_key], show=False)
            plt.savefig(
                ctx.figure_dir / "umap_batch_before.png", bbox_inches="tight",
            )
            plt.close()

        if cfg.method == "harmony":
            backend = self._run_harmony
        elif cfg.method == "bbknn":
            backend = self._run_bbknn
        elif cfg.method == "combat":
            backend = self._run_combat
        elif cfg.method == "scanorama":
            backend = self._run_scanorama
        elif cfg.method == "scvi":
            backend = self._run_scvi
        elif cfg.method == "mnn":
            backend = self._run_mnn
        elif cfg.method == "fastmnn":
            backend = self._run_fastmnn
        else:
            raise ValueError(f"Unknown batch correction method: {cfg.method}")
        _ensure_cpu_batch_inputs(adata)
        try:
            backend(adata, batch_key, ctx)
        except Exception as exc:
            if cfg.method not in {"scvi", "mnn", "fastmnn"}:
                raise
            # H-6 audit fix (2026-05-22): Harmony already fails loud on
            # backend failure. For scvi/mnn/fastmnn the user typed an
            # explicit --batch-method choice, so a silent skip masks the
            # fact that integration did NOT run. Require an explicit opt-in
            # env var (matching the SC_ALLOW_WELCH_FALLBACK / F-3 pattern)
            # before the silent skip is permitted.
            opt_in_env = "SC_ALLOW_BATCH_BACKEND_SKIP"
            if os.environ.get(opt_in_env, "").strip() != "1":
                raise RuntimeError(
                    f"batch_correction: explicitly-selected backend "
                    f"'{cfg.method}' failed ({exc!r}) and the silent skip "
                    f"is banned by default. Set {opt_in_env}=1 to opt in to "
                    "the skip (and accept that batch_correction did not run); "
                    "or fix the backend install."
                ) from exc
            msg = f"{cfg.method} backend unavailable or failed: {exc}"
            logger.error(
                "%s; %s=1 acknowledged — skipping batch correction.",
                msg, opt_in_env,
            )
            ctx.metadata["batch_correction_status"] = f"skipped_{cfg.method}_unavailable_or_failed_opt_in"
            ctx.metadata["batch_correction_method_actually_used"] = "none_opt_in_skip"
            ctx.metadata["batch_correction_skip_reason"] = str(exc)
            ctx.metadata["batch_correction_skip_opt_in_acknowledged"] = True
            ctx.status(self.name, "skipped", msg)
            return

        # Recompute UMAP + Leiden on corrected representation
        use_rep = "X_pca"
        if cfg.method == "harmony":
            use_rep = "X_pca_harmony"
        elif cfg.method == "combat":
            use_rep = "X_pca_combat"
        elif cfg.method == "scanorama":
            use_rep = "X_scanorama"
        elif cfg.method == "scvi":
            use_rep = "X_scvi"
        elif cfg.method == "mnn":
            use_rep = "X_mnn"
        elif cfg.method == "fastmnn":
            use_rep = "X_fastmnn"
        ctx.metadata["batch_correction_use_rep"] = use_rep

        use_gpu = gpu_available(ctx.cfg.gpu_mode)
        post_processing_ok = False
        if use_gpu:
            try:
                import rapids_singlecell as rsc
                bind_cuda_context()  # ensure cuBLAS/cuSOLVER pre-warm before GPU ops (P1 fix; no-op if already warmed in probe)
                logger.info("GPU batch post-processing: neighbors/UMAP/Leiden")
                if cfg.method != "bbknn":
                    # bbknn bypasses cache: incompatible neighbors API signature
                    get_or_compute_neighbors(
                        adata,
                        backend_id="rapids",
                        method="umap",
                        metric="euclidean",
                        n_pcs=None,
                        n_neighbors=ctx.cfg.clustering.n_neighbors,
                        use_rep=use_rep,
                        knn=True,
                        random_state=ctx.cfg.clustering.random_state,
                        compute_fn=lambda: rsc.pp.neighbors(adata, use_rep=use_rep),
                    )
                rsc.tl.umap(adata, random_state=ctx.cfg.clustering.random_state)
                rsc.tl.leiden(
                    adata,
                    resolution=ctx.cfg.clustering.leiden_resolution,
                    random_state=ctx.cfg.clustering.random_state,
                )
                post_processing_ok = True
            except Exception as exc:
                if ctx.cfg.gpu_mode == "force":
                    raise RuntimeError(f"GPU batch post-processing failed in force mode: {exc}") from exc
                logger.warning("GPU batch post-processing failed (%s), falling back to CPU", exc)
                ctx.metadata["batch_gpu_fallback_reason"] = str(exc)
                use_gpu = False

        if not use_gpu:
            try:
                if cfg.method != "bbknn":
                    # bbknn bypasses cache: incompatible neighbors API signature
                    get_or_compute_neighbors(
                        adata,
                        backend_id="scanpy",
                        method="umap",
                        metric="euclidean",
                        n_pcs=None,
                        n_neighbors=ctx.cfg.clustering.n_neighbors,
                        use_rep=use_rep,
                        knn=True,
                        random_state=ctx.cfg.clustering.random_state,
                        compute_fn=lambda: sc.pp.neighbors(adata, use_rep=use_rep),
                    )
                sc.tl.umap(adata, random_state=ctx.cfg.clustering.random_state)
                sc.tl.leiden(
                    adata,
                    resolution=ctx.cfg.clustering.leiden_resolution,
                    flavor="igraph",
                    directed=False,
                    random_state=ctx.cfg.clustering.random_state,
                )
                post_processing_ok = True
            except Exception as exc:
                # W14 fix: previously this swallowed the failure (warning + keep
                # PRE-correction UMAP/leiden) yet the module still reported
                # status="completed". That is a silent compromise: the corrected
                # representation (use_rep) was computed, but neighbors/UMAP/Leiden
                # on it failed, so the leiden labels + X_umap left behind are
                # STALE (computed on the pre-correction X_pca). Downstream
                # annotation/DE would then bind to clusters that do not match the
                # corrected embedding the run claims to have produced. Fail loud
                # by default; allow an explicit, audited opt-in via the SAME env
                # var the backend-failure path uses (SC_ALLOW_BATCH_BACKEND_SKIP).
                ctx.metadata["batch_post_error"] = str(exc)
                opt_in_env = "SC_ALLOW_BATCH_BACKEND_SKIP"
                if os.environ.get(opt_in_env, "").strip() != "1":
                    raise RuntimeError(
                        f"batch_correction: CPU post-processing "
                        f"(neighbors/UMAP/Leiden on use_rep='{use_rep}') failed "
                        f"({exc!r}). The corrected representation was computed but "
                        f"clustering on it did NOT run, so leiden/X_umap are STALE "
                        f"(pre-correction). Refusing to report 'completed' with "
                        f"stale clustering (Principle 9: no silent compromise). "
                        f"Set {opt_in_env}=1 to opt in to keeping the stale "
                        f"clustering (status will be recorded as "
                        f"'completed_with_stale_clustering_opt_in', NOT "
                        f"'completed'), or fix the post-processing failure."
                    ) from exc
                logger.error(
                    "CPU batch post-processing failed (%s); %s=1 acknowledged — "
                    "keeping STALE pre-correction UMAP/leiden. Clustering does NOT "
                    "reflect the corrected representation '%s'.",
                    exc, opt_in_env, use_rep,
                )
                ctx.metadata["batch_post_processing_stale_clustering"] = True
                ctx.metadata["batch_post_processing_skip_opt_in_acknowledged"] = True
        ctx.metadata["batch_post_backend"] = "gpu" if use_gpu else "cpu"
        # Non-destructive audit column: snapshot the post-correction leiden so
        # annotation can bind cell_type to the corrected clusters even after a
        # later module re-clusters. The `provides obs.leiden` overwrite above is
        # intentionally kept; leiden_corrected survives parallel merge-back as a
        # new obs column. See annotation leiden-race fix.
        if post_processing_ok and "leiden" in adata.obs:
            adata.obs["leiden_corrected"] = adata.obs["leiden"].copy()
            ctx.metadata["leiden_corrected_written"] = True
        ctx.metadata["n_clusters_after_batch"] = int(adata.obs["leiden"].nunique()) if "leiden" in adata.obs else 0
        # W14 fix: only report "completed" when the corrected representation was
        # actually re-clustered. A swallowed post-processing failure (reachable
        # only via the SC_ALLOW_BATCH_BACKEND_SKIP opt-in above; the default path
        # raises) leaves stale pre-correction clustering, so it gets a DISTINCT
        # status instead of the misleading "completed".
        if post_processing_ok:
            ctx.metadata["batch_correction_status"] = "completed"
        else:
            ctx.metadata["batch_correction_status"] = (
                "completed_with_stale_clustering_opt_in"
            )
        ctx.metadata["batch_correction_random_state"] = ctx.random_state
        if post_processing_ok:
            self._write_corrected_resolution_sweep(adata, ctx, use_rep)
        elif getattr(ctx.cfg.clustering, "leiden_resolution_sweep", ()):
            ctx.metadata["leiden_resolution_sweep_status"] = (
                "skipped_batch_post_processing_failed"
            )
        mixing_audit["after"] = self._compute_batch_mixing_metrics(adata, batch_key, use_rep)
        self._write_batch_mixing_audit(ctx, mixing_audit)

        # Post-correction UMAP
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        sc.pl.umap(adata, color=[batch_key], ax=axes[0], show=False)
        axes[0].set_title(f"After {cfg.method} — batch")
        sc.pl.umap(adata, color=["leiden"], ax=axes[1], show=False, legend_loc="on data")
        axes[1].set_title(f"After {cfg.method} — clusters")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "umap_batch_after.png", bbox_inches="tight")
        plt.close()

    @staticmethod
    def _compute_batch_mixing_metrics(
        adata,
        batch_key: str,
        use_rep: str,
        n_neighbors: int = 15,
    ) -> dict:
        if use_rep not in adata.obsm:
            return {"status": "skipped_missing_representation", "use_rep": use_rep}
        if batch_key not in adata.obs:
            return {"status": "skipped_missing_batch_key", "use_rep": use_rep}

        x = np.asarray(adata.obsm[use_rep])
        labels = adata.obs[batch_key].astype(str).to_numpy()
        finite = np.isfinite(x).all(axis=1)
        x = x[finite]
        labels = labels[finite]
        n_cells = int(labels.shape[0])
        batch_levels = np.array(sorted(set(labels.tolist())))
        n_batches = int(batch_levels.shape[0])
        if n_cells < 3 or n_batches < 2:
            return {
                "status": "skipped_insufficient_cells_or_batches",
                "use_rep": use_rep,
                "n_cells": n_cells,
                "n_batches": n_batches,
            }

        k = min(int(n_neighbors) + 1, n_cells)
        nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
        nn.fit(x)
        neighbor_idx = nn.kneighbors(x, return_distance=False)
        if neighbor_idx.shape[1] > 1 and np.array_equal(neighbor_idx[:, 0], np.arange(n_cells)):
            neighbor_idx = neighbor_idx[:, 1:]
        else:
            neighbor_idx = neighbor_idx[:, : min(int(n_neighbors), neighbor_idx.shape[1])]
        if neighbor_idx.shape[1] == 0:
            return {
                "status": "skipped_no_neighbors",
                "use_rep": use_rep,
                "n_cells": n_cells,
                "n_batches": n_batches,
            }

        neighbor_labels = labels[neighbor_idx]
        same_batch_fraction = float((neighbor_labels == labels[:, None]).mean())
        level_index = {label: idx for idx, label in enumerate(batch_levels)}
        entropies = []
        for row in neighbor_labels:
            counts = np.zeros(n_batches, dtype=float)
            for label in row:
                counts[level_index[label]] += 1.0
            probs = counts[counts > 0] / counts.sum()
            entropy = -float(np.sum(probs * np.log(probs))) / float(np.log(n_batches))
            entropies.append(entropy)
        return {
            "status": "ok",
            "use_rep": use_rep,
            "n_cells": n_cells,
            "n_batches": n_batches,
            "n_neighbors_effective": int(neighbor_idx.shape[1]),
            "same_batch_neighbor_fraction": round(same_batch_fraction, 6),
            "normalized_batch_entropy": round(float(np.mean(entropies)), 6),
        }

    @staticmethod
    def _write_batch_mixing_audit(ctx: PipelineContext, audit: dict) -> None:
        (ctx.table_dir / "batch_mixing_metrics.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        for phase in ("before", "after"):
            metrics = audit.get(phase, {})
            if metrics.get("status") != "ok":
                continue
            ctx.metadata[f"batch_mixing_{phase}_same_batch_neighbor_fraction"] = metrics[
                "same_batch_neighbor_fraction"
            ]
            ctx.metadata[f"batch_mixing_{phase}_normalized_batch_entropy"] = metrics[
                "normalized_batch_entropy"
            ]

    @staticmethod
    def _write_corrected_resolution_sweep(adata, ctx: PipelineContext, use_rep: str) -> None:
        from .clustering import ClusteringModule

        ClusteringModule._write_resolution_sweep_audit(adata, ctx, ctx.cfg.clustering)
        if ctx.metadata.get("leiden_resolution_sweep_status") == "completed":
            ctx.metadata["leiden_resolution_sweep_representation"] = use_rep

    @staticmethod
    def _resolve_harmony_backend(ctx) -> str:
        """Resolve the Harmony backend.

        Selection precedence (highest first):
        1. `SC_HARMONY_BACKEND` env var (`cpu` | `gpu` | `auto` | `direct`)
        2. `ctx.cfg.batch.harmony_backend` (`cpu` | `gpu` | `auto` | `direct`)
        3. Module default: `auto`

        `auto` resolves to `direct` on this env (DOCUMENTED ROUTING CHANGE,
        2026-05-27): both the `gpu` path (rapids CUBLAS_STATUS_NOT_INITIALIZED
        + CUDA-context corruption) and the `cpu` path (scanpy 1.12
        harmony_integrate mis-stores harmonypy 0.2.0's ``Z_corr`` with the wrong
        shape) are broken on the sc_gpu env / this hardware, whereas
        ``_run_harmony_direct`` (canonical Korsunsky-2019 via
        ``harmonypy.run_harmony`` + ``Z_corr.T``) WORKS. This is an EXPLICIT
        routing decision, NOT a silent try/except fallback around broken calls.
        Explicit `gpu` still raises if rsc is unavailable; explicit `cpu` still
        uses harmonypy via scanpy_external unchanged (both escape hatches
        preserved for when the upstream bugs are fixed). Operators can also
        force any backend via `SC_HARMONY_BACKEND` or `ctx.cfg.batch.harmony_backend`.

        `direct` bypasses scanpy/rapids' thin Harmony wrappers and calls
        ``harmonypy.run_harmony`` (the canonical Korsunsky 2019 reference)
        directly, transposing its (n_pcs, n_cells) ``Z_corr`` to
        (n_cells, n_pcs). The integration_select gate also selects this backend
        when it recommends Harmony.
        """
        choice = (
            os.environ.get("SC_HARMONY_BACKEND", "").strip().lower()
            or getattr(ctx.cfg.batch, "harmony_backend", "").strip().lower()
            or "auto"
        )
        if choice not in {"auto", "cpu", "gpu", "direct"}:
            raise ValueError(
                f"Invalid harmony_backend={choice!r}; expected one of "
                f"{{'auto', 'cpu', 'gpu', 'direct'}}."
            )
        if choice == "direct":
            return "direct"
        if choice == "auto":
            # DOCUMENTED ROUTING CHANGE (2026-05-27): auto -> direct because
            # both gpu (rapids CUBLAS) and cpu (scanpy 1.12 wrapper mis-stores
            # harmonypy Z_corr) are broken on this env; _run_harmony_direct
            # works. NOT a silent fallback — explicit gpu/cpu still honored.
            return "direct"
        if choice == "gpu":
            try:
                import rapids_singlecell  # noqa: F401
            except ImportError as exc:
                raise ImportError(
                    "harmony_backend=gpu requires rapids_singlecell. "
                    "Install rsc or switch to harmony_backend=cpu/auto."
                ) from exc
        return choice

    @staticmethod
    def _detect_non_convergence(caught_warnings) -> tuple[bool, list[str]]:
        """Scan captured warnings for Harmony non-convergence signals.

        Both rsc and harmonypy emit a Python warning when the iterative loop
        exits without convergence. Returns (converged, signals). `converged`
        is False iff at least one warning matches a non-convergence token.
        """
        tokens = ("did not converge", "maximum iter", "max iter", "reached max")
        signals = [
            str(w.message) for w in caught_warnings
            if any(t in str(w.message).lower() for t in tokens)
        ]
        return (len(signals) == 0, signals)

    @staticmethod
    def _run_harmony(adata, batch_key: str, ctx) -> None:
        """Run Harmony batch correction.

        Honors Principle 9 (No silent non-convergence) — if Harmony exits the
        iterative loop without converging at `harmony_max_iter`, the module
        raises `RuntimeError` unless the operator explicitly opts in via
        `SC_ALLOW_HARMONY_NON_CONVERGENCE=1`. Mirrors F-3's Welch opt-in.

        Backend selection: see `_resolve_harmony_backend`. The factory's
        default is `auto` (GPU if rapids_singlecell is importable, else CPU
        harmonypy). G-C0 v3 contract requires both lanes to converge —
        backend is allowed to differ across lanes provided both converge,
        because US-B3 + production C benchmark establish parity at
        convergence.
        """
        if "X_pca" not in adata.obsm:
            raise ValueError("Harmony requires PCA (run clustering first).")
        cfg_batch = ctx.cfg.batch
        backend = BatchCorrectionModule._resolve_harmony_backend(ctx)
        if backend == "direct":
            BatchCorrectionModule._run_harmony_direct(adata, batch_key, ctx)
            return
        ctx.metadata["harmony_backend"] = backend
        ctx.metadata["harmony_device"] = "gpu" if backend == "gpu" else "cpu"
        ctx.metadata["harmony_theta"] = cfg_batch.harmony_theta
        ctx.metadata["harmony_sigma"] = cfg_batch.harmony_sigma
        ctx.metadata["harmony_max_iter"] = cfg_batch.harmony_max_iter

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            if backend == "gpu":
                import rapids_singlecell as rsc
                bind_cuda_context()  # ensure cuBLAS/cuSOLVER pre-warm before GPU ops (P1 fix; no-op if already warmed in probe)
                try:
                    rsc.pp.harmony_integrate(
                        adata,
                        key=batch_key,
                        theta=cfg_batch.harmony_theta,
                        sigma=cfg_batch.harmony_sigma,
                        max_iter_harmony=cfg_batch.harmony_max_iter,
                    )
                except TypeError:
                    ctx.metadata["harmony_extra_args_unsupported"] = True
                    rsc.pp.harmony_integrate(adata, key=batch_key)
            else:  # backend == "cpu"
                # Use scanpy_external (harmonypy) — same param surface,
                # canonical CPU reference per Korsunsky 2019.
                import scanpy.external as sce
                try:
                    sce.pp.harmony_integrate(
                        adata,
                        key=batch_key,
                        theta=cfg_batch.harmony_theta,
                        sigma=cfg_batch.harmony_sigma,
                        max_iter_harmony=cfg_batch.harmony_max_iter,
                    )
                except TypeError:
                    ctx.metadata["harmony_extra_args_unsupported"] = True
                    sce.pp.harmony_integrate(adata, key=batch_key)

        converged, signals = BatchCorrectionModule._detect_non_convergence(caught)
        ctx.metadata["harmony_converged"] = converged
        if signals:
            ctx.metadata["harmony_non_convergence_signals"] = signals[:5]
        if not converged:
            opt_in = os.environ.get("SC_ALLOW_HARMONY_NON_CONVERGENCE", "").strip() == "1"
            if opt_in:
                logger.error(
                    "Harmony did not converge at max_iter=%d (backend=%s) — "
                    "SC_ALLOW_HARMONY_NON_CONVERGENCE=1 opt-in acknowledged; "
                    "downstream results may be iteration-trajectory-dependent. "
                    "Signals: %s",
                    cfg_batch.harmony_max_iter, backend, signals[:3],
                )
                ctx.metadata["harmony_non_convergence_opt_in_acknowledged"] = True
            else:
                raise RuntimeError(
                    f"Harmony did not converge at max_iter={cfg_batch.harmony_max_iter} "
                    f"(backend={backend}). Principle 9: no silent non-convergence. "
                    f"Raise --harmony-max-iter (recommended >= 50 for 75-batch tumor "
                    f"cohorts) or set SC_ALLOW_HARMONY_NON_CONVERGENCE=1 to override "
                    f"with audit-logged opt-in. Captured signals: {signals[:3]}"
                )

        # Ensure float32 for downstream consistency
        adata.obsm["X_pca_harmony"] = np.asarray(
            adata.obsm["X_pca_harmony"], dtype=np.float32
        )

    @staticmethod
    def _run_harmony_direct(adata, batch_key: str, ctx) -> None:
        """Harmony via ``harmonypy.run_harmony`` DIRECT — the proven-working path.

        DEVIATION (documented, authorized): both production scanpy/rapids Harmony
        backends are broken in the sc_gpu env on this hardware
        (``harmony_backend=gpu`` raises CUBLAS_STATUS_NOT_INITIALIZED;
        ``harmony_backend=cpu`` via scanpy.external.pp.harmony_integrate with
        harmonypy 0.2.0 stores ``X_pca_harmony`` with the wrong shape — (n_pcs,)
        instead of (n_cells, n_pcs)). ``harmonypy.run_harmony`` IS the canonical
        Korsunsky-2019 reference; this routes around scanpy 1.12's proven-buggy
        wrapper while running the SAME algorithm on the CPU.

        harmonypy returns ``Z_corr`` shape (n_pcs, n_cells); we transpose to
        (n_cells, n_pcs) before storing. Honors Principle 9 (no silent
        non-convergence) — non-convergence raises unless
        ``SC_ALLOW_HARMONY_NON_CONVERGENCE=1`` is set.
        """
        if "X_pca" not in adata.obsm:
            raise ValueError("Harmony requires PCA (run clustering first).")
        cfg_batch = ctx.cfg.batch
        ctx.metadata["harmony_backend"] = "direct"
        ctx.metadata["harmony_theta"] = cfg_batch.harmony_theta
        ctx.metadata["harmony_sigma"] = cfg_batch.harmony_sigma
        ctx.metadata["harmony_max_iter"] = cfg_batch.harmony_max_iter
        ctx.metadata["harmony_engine"] = "harmonypy.run_harmony direct (scanpy wrapper bypassed)"

        import harmonypy

        # Principle 9 (no silent non-convergence): harmonypy reports convergence
        # via the `logging` module, NOT `warnings` — it logs "Converged after N
        # iterations" on success and "Stopped before convergence" otherwise. The
        # warnings-based _detect_non_convergence alone would MISS a non-converged
        # direct run and silently report success. Capture the harmonypy logger to
        # detect convergence authoritatively (review MEDIUM-2).
        _hlog = logging.getLogger("harmonypy")

        class _ConvergenceCatcher(logging.Handler):
            def __init__(self) -> None:
                super().__init__()
                self.saw_record = False
                self.saw_converged = False
                self.signals: list[str] = []

            def emit(self, record) -> None:
                msg = record.getMessage()
                self.saw_record = True
                low = msg.lower()
                if "converged after" in low:
                    self.saw_converged = True
                if "stopped before convergence" in low or "did not converge" in low:
                    self.signals.append(msg)

        catcher = _ConvergenceCatcher()
        _prev_level = _hlog.level
        _hlog.addHandler(catcher)
        if _hlog.level == logging.NOTSET or _hlog.level > logging.INFO:
            _hlog.setLevel(logging.INFO)
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                pca = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
                harmony_kwargs = dict(
                    theta=cfg_batch.harmony_theta,
                    sigma=cfg_batch.harmony_sigma,
                    max_iter_harmony=cfg_batch.harmony_max_iter,
                    # Harmony initialises its soft-clustering with KMeans and is
                    # therefore stochastic (Korsunsky et al. 2019, *Nat Methods*
                    # 16:1289-1296). Scanorama, MNN and scVI already seed from
                    # ctx.random_state; Harmony — the default method — did not,
                    # so ctx.metadata["batch_correction_random_state"] recorded a
                    # seed the algorithm never received and reruns could drift.
                    random_state=ctx.random_state,
                )
                # harmonypy >= 0.2 grew a PyTorch backend whose `device=None`
                # default AUTO-SELECTS CUDA. Leaving it unset made this backend
                # log "Running Harmony (PyTorch on cuda)" and then die with
                # CUBLAS_STATUS_NOT_INITIALIZED on the very hardware whose broken
                # CUDA Harmony path is the documented reason `direct` exists —
                # while ctx.metadata still recorded harmony_device="cpu". Pin the
                # device explicitly so the recorded provenance is the truth.
                if "device" in inspect.signature(harmonypy.run_harmony).parameters:
                    harmony_kwargs["device"] = "cpu"
                    ctx.metadata["harmony_device"] = "cpu"
                    ctx.metadata["harmony_device_pinned"] = True
                else:
                    # Pre-0.2 harmonypy has no torch backend and is CPU-only.
                    ctx.metadata["harmony_device"] = "cpu"
                    ctx.metadata["harmony_device_pinned"] = False
                ho = harmonypy.run_harmony(pca, adata.obs, [batch_key], **harmony_kwargs)
        finally:
            _hlog.removeHandler(catcher)
            _hlog.setLevel(_prev_level)

        # Prefer the authoritative harmonypy log signal. If harmonypy emitted
        # records but never "Converged after", it did NOT converge. Only fall
        # back to the warnings-based check if harmonypy emitted no records at all
        # (anomalous logging) — recorded loud so the anomaly is visible.
        warn_converged, warn_signals = BatchCorrectionModule._detect_non_convergence(caught)
        if catcher.saw_record:
            converged = catcher.saw_converged
            signals = catcher.signals or warn_signals
        else:
            converged = warn_converged
            signals = warn_signals
            ctx.metadata["harmony_convergence_detection"] = (
                "harmonypy emitted no log records; fell back to warnings-based "
                "detection (anomalous — verify harmonypy logging is enabled)."
            )
        ctx.metadata["harmony_converged"] = converged
        if signals:
            ctx.metadata["harmony_non_convergence_signals"] = signals[:5]
        if not converged:
            opt_in = os.environ.get("SC_ALLOW_HARMONY_NON_CONVERGENCE", "").strip() == "1"
            if opt_in:
                logger.error(
                    "Harmony (direct) did not converge at max_iter=%d — "
                    "SC_ALLOW_HARMONY_NON_CONVERGENCE=1 opt-in acknowledged; "
                    "downstream results may be iteration-trajectory-dependent. "
                    "Signals: %s",
                    cfg_batch.harmony_max_iter, signals[:3],
                )
                ctx.metadata["harmony_non_convergence_opt_in_acknowledged"] = True
            else:
                raise RuntimeError(
                    f"Harmony (direct) did not converge at "
                    f"max_iter={cfg_batch.harmony_max_iter}. Principle 9: no "
                    f"silent non-convergence. Raise --harmony-max-iter "
                    f"(recommended >= 50 for 75-batch tumor cohorts) or set "
                    f"SC_ALLOW_HARMONY_NON_CONVERGENCE=1 to override with "
                    f"audit-logged opt-in. Captured signals: {signals[:3]}"
                )

        # harmonypy Z_corr is (n_pcs, n_cells); transpose to (n_cells, n_pcs).
        Z = np.asarray(ho.Z_corr, dtype=np.float32)
        if Z.shape[0] != adata.n_obs:
            Z = Z.T
        if Z.shape[0] != adata.n_obs:
            raise ValueError(
                f"harmonypy output shape {Z.shape} does not align to "
                f"n_obs={adata.n_obs} on either axis."
            )
        adata.obsm["X_pca_harmony"] = np.ascontiguousarray(Z, dtype=np.float32)

    @staticmethod
    def _run_bbknn(adata, batch_key: str, ctx) -> None:
        try:
            import bbknn
        except ImportError:
            raise ImportError(
                "bbknn is required for BBKNN batch correction. "
                "Install with: pip install bbknn"
            )
        if "X_pca" not in adata.obsm:
            raise ValueError("BBKNN requires PCA (run clustering first).")
        bbknn.bbknn(adata, batch_key=batch_key)

    @staticmethod
    def _run_combat(adata, batch_key: str, ctx) -> None:
        """ComBat correction (Johnson 2007) + corrected-embedding recompute.

        sc.pp.combat rewrites only adata.X. Previously the downstream
        neighbors/UMAP/Leiden and the 'after' mixing metric still ran on the
        stale X_pca, making ComBat a no-op for embedding/clustering. Recompute
        PCA on the corrected X into obsm['X_pca_combat'] so use_rep='X_pca_combat'
        flows through the post-processing recompute.
        """
        sc.pp.combat(adata, key=batch_key)
        n_pcs = int(getattr(ctx.cfg.clustering, "n_pcs", 50))
        n_pcs = max(2, min(n_pcs, adata.n_obs - 1, adata.n_vars - 1))
        # Recompute PCA on corrected X without clobbering the pre-combat
        # obsm['X_pca'] (kept for the 'before' baseline / audit).
        prev_pca = adata.obsm.get("X_pca")
        sc.pp.pca(adata, n_comps=n_pcs)
        adata.obsm["X_pca_combat"] = np.asarray(adata.obsm["X_pca"], dtype=np.float32)
        if prev_pca is not None:
            adata.obsm["X_pca"] = prev_pca
        ctx.metadata["combat_n_pcs"] = n_pcs

    @staticmethod
    def _run_scanorama(adata, batch_key: str, ctx) -> None:
        """Scanorama integration (Hie et al., Nature Biotechnology 2019)."""
        try:
            import scanorama
        except ImportError:
            raise ImportError(
                "scanorama is required for Scanorama batch correction. "
                "Install with: pip install scanorama"
            )

        if "X_pca" not in adata.obsm:
            raise ValueError("Scanorama requires PCA (run clustering first).")
        # Determinism: Scanorama uses randomized SVD internally (Hie 2019); seed
        # numpy from ctx.random_state for reproducible alignment.
        np.random.seed(ctx.random_state)
        ctx.metadata["scanorama_seed"] = int(ctx.random_state)
        batches = adata.obs[batch_key].unique().tolist()
        adatas = [adata[adata.obs[batch_key] == b].copy() for b in batches]
        scanorama.integrate_scanpy(adatas)
        # Reconstruct corrected embedding into original object
        adata.obsm["X_scanorama"] = np.zeros(
            (adata.n_obs, adatas[0].obsm["X_scanorama"].shape[1])
        )
        for b, sub in zip(batches, adatas):
            mask = adata.obs[batch_key] == b
            adata.obsm["X_scanorama"][mask] = sub.obsm["X_scanorama"]

    @staticmethod
    def _run_mnn(adata, batch_key: str, ctx) -> None:
        """MNN correction via scanpy.external.pp.mnn_correct (Haghverdi et al., 2018)."""
        if "X_pca" not in adata.obsm:
            raise ValueError("MNN requires PCA (run clustering first).")
        # Determinism: MNN/fastMNN PCA + matching are stochastic (Haghverdi 2018);
        # seed numpy from ctx.random_state for reproducible correction.
        np.random.seed(ctx.random_state)
        ctx.metadata["mnn_seed"] = int(ctx.random_state)
        try:
            _ = sc.external.pp.mnn_correct
        except (AttributeError, ImportError):
            raise ImportError(
                "mnn_correct is unavailable. Install mnnpy and ensure scanpy.external is available. "
                "Example: pip install mnnpy"
            )

        batches = adata.obs[batch_key].astype(str)
        batch_levels = list(dict.fromkeys(batches.tolist()))
        if len(batch_levels) < 2:
            raise ValueError("MNN requires at least two batches.")
        adatas = [adata[batches == b].copy() for b in batch_levels]

        try:
            corrected = sc.external.pp.mnn_correct(*adatas, batch_key=None, index_unique=None)[0]
        except TypeError:
            # Fallback for older mnnpy/scanpy signatures.
            corrected = sc.external.pp.mnn_correct(*adatas)[0]
        corrected = corrected[adata.obs_names, :].copy()
        sc.pp.pca(corrected, n_comps=ctx.cfg.clustering.n_pcs)
        adata.obsm["X_mnn"] = corrected.obsm["X_pca"].copy()

    @staticmethod
    def _run_fastmnn(adata, batch_key: str, ctx) -> None:
        """Fast-MNN style approximation using mnn_correct with reduced SVD dimension."""
        if "X_pca" not in adata.obsm:
            raise ValueError("fastMNN requires PCA (run clustering first).")
        try:
            _ = sc.external.pp.mnn_correct
        except (AttributeError, ImportError):
            raise ImportError(
                "mnn_correct is unavailable. Install mnnpy and ensure scanpy.external is available. "
                "Example: pip install mnnpy"
            )

        batches = adata.obs[batch_key].astype(str)
        batch_levels = list(dict.fromkeys(batches.tolist()))
        if len(batch_levels) < 2:
            raise ValueError("fastMNN requires at least two batches.")
        adatas = [adata[batches == b].copy() for b in batch_levels]

        svd_dim = max(10, int(ctx.cfg.clustering.n_pcs))
        try:
            corrected = sc.external.pp.mnn_correct(
                *adatas, batch_key=None, index_unique=None, svd_dim=svd_dim
            )[0]
        except TypeError:
            corrected = sc.external.pp.mnn_correct(*adatas)[0]
        corrected = corrected[adata.obs_names, :].copy()
        sc.pp.pca(corrected, n_comps=ctx.cfg.clustering.n_pcs)
        adata.obsm["X_fastmnn"] = corrected.obsm["X_pca"].copy()
        ctx.metadata["fastmnn_backend"] = "mnn_correct_svd_approx"
        ctx.metadata["fastmnn_svd_dim"] = svd_dim

    @staticmethod
    def _matrix_looks_integer_like(matrix, sample_size: int = 2048) -> bool:
        """Heuristic: check if matrix resembles raw integer counts."""
        if sparse.issparse(matrix):
            values = matrix.data
        else:
            values = np.asarray(matrix).ravel()
        if values.size == 0:
            return True
        if values.size > sample_size:
            idx = np.linspace(0, values.size - 1, sample_size, dtype=int)
            values = values[idx]
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return False
        if np.any(finite < 0):
            return False
        return bool(np.allclose(finite, np.round(finite), atol=1e-6))

    # W11 contract version. Bump when the counts-resolution contract changes so
    # prior runs are recognizable as having used a different (stale) contract.
    SCVI_COUNTS_CONTRACT = "scvi-counts-v3-ambient-corrected-rounded-int"

    @staticmethod
    def _ambient_correction_ran(adata) -> bool:
        """True iff ambient_correction (DecontX) actually rewrote ``adata.X``.

        Primary signal: ``adata.uns['ambient_correction']['decision'] ==
        'triggered'`` (set by ambient_correction.py only on the path that
        replaces ``.X`` with corrected counts). Secondary robust signal: the
        pre-correction snapshot ``layers['counts_raw_pre_decontx']`` is present
        (ambient_correction writes it on the SAME line it overwrites ``.X``).
        Either is sufficient — the snapshot survives even if ``uns`` is dropped
        by a downstream serialization round-trip.
        """
        uns = getattr(adata, "uns", {}) or {}
        amb = uns.get("ambient_correction")
        if isinstance(amb, dict) and amb.get("decision") == "triggered":
            return True
        layers = getattr(adata, "layers", {}) or {}
        return "counts_raw_pre_decontx" in layers

    @staticmethod
    def _run_scvi(adata, batch_key: str, ctx) -> None:
        """scVI latent integration (Lopez et al., Nature Methods 2018).

        W11 fix — counts-source consistency: clustering and DE operate on
        ``adata.X``, which ambient_correction OVERWRITES with the
        DecontX-corrected matrix when a trigger fires (preserving the
        pre-correction raw at ``layers['counts_raw_pre_decontx']``). It does NOT
        refresh ``layers['counts']``, so the previous "prefer layers['counts']"
        rule made scVI train on the STALE pre-ambient counts while
        clustering/DE used the corrected matrix — a silent cross-module
        inconsistency. The resolver below mirrors the clustering/DE path: when
        ambient ran, scVI consumes the ambient-corrected ``adata.X``; only when
        ambient did NOT run does it fall back to ``layers['counts']`` (raw UMI)
        or an integer-like ``adata.X``. No silent fallback: a non-count-like
        ``adata.X`` with no usable raw layer still RAISES (Principle 9)."""
        try:
            import scvi
        except ImportError:
            raise ImportError(
                "scvi-tools is required for scVI batch correction. "
                "Install with: pip install scvi-tools"
            )

        batch_cfg = ctx.cfg.batch
        ambient_ran = BatchCorrectionModule._ambient_correction_ran(adata)
        if ambient_ran:
            # Consume the SAME corrected signal clustering/DE consumed (the corrected .X),
            # but scVI's NB/ZINB likelihood requires INTEGER counts and DecontX corrected
            # counts are fractional (probabilistic redistribution). So round the corrected
            # matrix to integers into a dedicated layer for scVI: consistency (W11) + NB
            # count contract honored. We deliberately do NOT route to the stale
            # layers['counts'] (pre-ambient) here even when it exists.
            import numpy as _np
            from scipy import sparse as _sp
            _xc = adata.X
            if _sp.issparse(_xc):
                _xi = _xc.copy().astype(_np.float32)
                _xi.data = _np.rint(_xi.data)
            else:
                _xi = _np.rint(_np.asarray(_xc)).astype(_np.float32)
            adata.layers["counts_decontx_int"] = _xi
            layer = "counts_decontx_int"
            input_source = "adata.X_ambient_corrected_rounded_int"
        elif "counts" in adata.layers:
            layer = "counts"
            input_source = "layers.counts"
        else:
            layer = None
            if not BatchCorrectionModule._matrix_looks_integer_like(adata.X):
                raise ValueError(
                    "scVI requires raw count-like input. 'counts' layer is missing and adata.X "
                    "does not appear to be integer UMI counts."
                )
            input_source = "adata.X"

        # Real (non-mock) integrality guard on whatever scVI actually consumes — NB/ZINB
        # likelihood requires integer counts. Catches fractional input from any branch.
        _scvi_mat = adata.layers[layer] if layer is not None else adata.X
        if not BatchCorrectionModule._matrix_looks_integer_like(_scvi_mat):
            raise ValueError(
                f"scVI input ({input_source}) is not integer-like; NB likelihood requires "
                "integer counts. (DecontX corrected counts must be rounded before scVI.)")
        # Determinism: scVI is stochastic (Lopez 2018; scvi-tools reproducibility
        # docs). Seed from the canonical ctx.random_state before setup_anndata so
        # repeated runs with the same seed are reproducible. Guard on settings
        # presence so older scvi-tools layouts don't break.
        if hasattr(scvi, "settings"):
            scvi.settings.seed = int(ctx.random_state)
        ctx.metadata["scvi_seed"] = int(ctx.random_state)
        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key, layer=layer)
        model = scvi.model.SCVI(adata, n_latent=batch_cfg.scvi_n_latent)
        model.train(
            max_epochs=batch_cfg.scvi_max_epochs,
            early_stopping=batch_cfg.scvi_early_stopping,
        )
        adata.obsm["X_scvi"] = model.get_latent_representation()
        ctx.metadata["scvi_input_source"] = input_source
        # W11: record both the consumed-counts source and the contract version
        # so a post-run reviewer can tell which counts scVI actually trained on
        # and whether the run used the ambient-aware contract.
        ctx.metadata["scvi_consumed_counts_source"] = input_source
        ctx.metadata["scvi_counts_contract"] = BatchCorrectionModule.SCVI_COUNTS_CONTRACT
        ctx.metadata["scvi_ambient_correction_ran"] = bool(ambient_ran)
        ctx.metadata["scvi_corrected_counts_rounded_to_int"] = bool(ambient_ran)
        ctx.metadata["scvi_train_config"] = {
            "max_epochs": batch_cfg.scvi_max_epochs,
            "n_latent": batch_cfg.scvi_n_latent,
            "early_stopping": batch_cfg.scvi_early_stopping,
        }
