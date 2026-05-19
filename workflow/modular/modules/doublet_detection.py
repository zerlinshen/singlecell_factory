from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from ._scanpy_compat import import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext


__references__ = {
    "Wolock_Scrublet_2019": {
        "title": "Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data",
        "authors": "Wolock, Lopez, Klein",
        "journal": "Cell Systems",
        "year": "2019",
        "doi": "10.1016/j.cels.2018.11.005",
        "description": "Reference implementation; supports whole-dataset and per-sample (grouped) doublet rate calibration.",
    },
}


logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# GPU availability probe (evaluated once at import time)
# ---------------------------------------------------------------------------
try:
    import rapids_singlecell as _rsc  # noqa: F401
    _RSC_AVAILABLE = True
except ImportError:
    _RSC_AVAILABLE = False


class DoubletDetectionModule:
    """Mandatory module: detect and remove doublets using Scrublet.

    Should run after QC and before clustering. Doublets (two cells captured
    in one droplet) can create artificial intermediate clusters and corrupt
    downstream differential expression results.

    GPU path (preferred): uses ``rsc.pp.scrublet`` when rapids-singlecell is
    available.  Falls back to CPU ``scrublet.Scrublet`` only if rsc is not
    importable.
    """

    name = "doublet_detection"
    required = True

    @staticmethod
    def _n_prin_comps(n_obs: int, n_vars: int) -> int:
        return max(2, min(30, n_obs - 1, n_vars - 1))

    @staticmethod
    def _materialize_counts_matrix(x):
        """Return a CPU CSR matrix, materializing lazy/dask/GPU arrays."""
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
        if sparse.issparse(x):
            return x.tocsr()
        return np.asarray(x)

    @staticmethod
    def _resolve_doublet_strategy(ctx, adata) -> str:
        """Return the effective doublet strategy string.

        Priority: SC_DOUBLET_STRATEGY env var > cfg.doublet_strategy > "auto".
        "auto" → "grouped" when n_obs >= 100k and sample column with >= 2 labels exists,
                 otherwise "whole".
        """
        import os as _os
        import pandas as _pd
        flag = (
            _os.environ.get("SC_DOUBLET_STRATEGY", "").strip().lower()
            or getattr(ctx.cfg, "doublet_strategy", "auto")
        )
        if flag in {"grouped", "whole", "skip"}:
            return flag
        # "auto"
        has_sample = (
            "sample" in adata.obs.columns
            and _pd.Series(adata.obs["sample"]).nunique() >= 2
        )
        if adata.n_obs >= 100000 and has_sample:
            return "grouped"
        return "whole"

    # ------------------------------------------------------------------
    # GPU (rsc) paths
    # ------------------------------------------------------------------

    @staticmethod
    def _free_gpu_memory_pool() -> None:
        # HIGH fix (2026-05-19): release CuPy GPU memory after each RAPIDS
        # doublet call. Previously the GPU AnnData copy stayed resident until
        # Python GC eventually collected it, which accumulated across modules
        # (clustering, batch correction, DE) and pushed long pipelines toward
        # GPU OOM. CuPy's default mempool is process-wide so explicit free is
        # both safe and effective.
        try:
            import cupy as cp  # type: ignore
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()
        except Exception:
            pass

    def _run_rsc_scrublet_whole(self, adata, cfg, ctx) -> tuple[np.ndarray, np.ndarray, float | None]:
        """Run rsc.pp.scrublet on the full dataset."""
        import rapids_singlecell as rsc

        adata_gpu = adata.copy()
        try:
            rsc.get.anndata_to_GPU(adata_gpu)
            rsc.pp.scrublet(
                adata_gpu,
                expected_doublet_rate=cfg.expected_doublet_rate,
                n_prin_comps=self._n_prin_comps(adata.n_obs, adata.n_vars),
                random_state=ctx.random_state,
                verbose=False,
            )
            scores = self._col_to_numpy(adata_gpu.obs["doublet_score"]).astype(np.float32)
            predicted = self._col_to_numpy(adata_gpu.obs["predicted_doublet"]).astype(bool)
            threshold = adata_gpu.uns.get("scrublet", {}).get("threshold", None)
            ctx.metadata["doublet_method"] = "scrublet_gpu"
            return scores, predicted, threshold
        finally:
            del adata_gpu
            self._free_gpu_memory_pool()

    def _run_rsc_scrublet_grouped(self, adata, cfg, ctx) -> tuple[np.ndarray, np.ndarray, None]:
        """Run rsc.pp.scrublet with batch_key for per-sample doublet detection."""
        import rapids_singlecell as rsc

        sample_key = "sample" if "sample" in adata.obs.columns else None
        if sample_key is None:
            raise ValueError("grouped scrublet requested without sample labels")

        adata_gpu = adata.copy()
        try:
            rsc.get.anndata_to_GPU(adata_gpu)
            rsc.pp.scrublet(
                adata_gpu,
                batch_key=sample_key,
                expected_doublet_rate=cfg.expected_doublet_rate,
                n_prin_comps=self._n_prin_comps(adata.n_obs, adata.n_vars),
                random_state=ctx.random_state,
                verbose=False,
            )
            scores = self._col_to_numpy(adata_gpu.obs["doublet_score"]).astype(np.float32)
            predicted = self._col_to_numpy(adata_gpu.obs["predicted_doublet"]).astype(bool)
            ctx.metadata["doublet_grouped_key"] = sample_key
            ctx.metadata["doublet_method"] = "scrublet_gpu_grouped"
            return scores, predicted, None
        finally:
            del adata_gpu
            self._free_gpu_memory_pool()

    # ------------------------------------------------------------------
    # CPU fallback paths (used only when rsc is unavailable)
    # ------------------------------------------------------------------

    def _run_cpu_scrublet_grouped(self, adata, scrublet_cls, cfg, ctx) -> tuple[np.ndarray, np.ndarray, None]:
        sample_key = "sample" if "sample" in adata.obs.columns else None
        if sample_key is None:
            raise ValueError("grouped scrublet requested without sample labels")
        random_state = ctx.random_state
        labels = pd.Series(adata.obs[sample_key].astype(str), index=adata.obs_names)
        scores = np.zeros(adata.n_obs, dtype=np.float32)
        predicted = np.zeros(adata.n_obs, dtype=bool)
        thresholds = {}
        scrublet_params = {
            "min_counts": 2,
            "min_cells": 3,
            "min_gene_variability_pctl": 85,
            "n_prin_comps": self._n_prin_comps(adata.n_obs, adata.n_vars),
        }
        for group in labels.unique():
            idx = np.where(labels.values == group)[0]
            if len(idx) < 20:
                continue
            subX = self._materialize_counts_matrix(adata.X[idx])
            scrub = scrublet_cls(subX, expected_doublet_rate=cfg.expected_doublet_rate, random_state=random_state)
            try:
                s, p = scrub.scrub_doublets(**{**scrublet_params, "n_prin_comps": self._n_prin_comps(len(idx), adata.n_vars)})
                scores[idx] = s.astype(np.float32, copy=False)
                predicted[idx] = p.astype(bool, copy=False)
                thr = getattr(scrub, "threshold_", None)
                if thr is not None:
                    thresholds[str(group)] = float(thr)
            except Exception as exc:
                logger.warning("Grouped Scrublet failed for %s, falling back to singlets: %s", group, exc)
        ctx.metadata["doublet_grouped_key"] = sample_key
        ctx.metadata["doublet_grouped_thresholds"] = thresholds
        ctx.metadata["doublet_method"] = "scrublet_grouped"
        return scores, predicted, None

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _col_to_numpy(col) -> np.ndarray:
        """Convert pandas or cudf Series to a numpy array."""
        if hasattr(col, "to_numpy"):
            return col.to_numpy()
        return np.asarray(col)

    @staticmethod
    def _fallback_all_singlets(n_obs: int, reason: str) -> tuple[np.ndarray, np.ndarray]:
        scores = np.zeros(n_obs, dtype=np.float32)
        predicted = np.zeros(n_obs, dtype=bool)
        logger.info("Doublet detection fallback activated: %s", reason)
        return scores, predicted

    # ------------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------------

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Doublet detection requires loaded AnnData.")

        cfg = ctx.cfg.doublet
        random_state = ctx.random_state
        threshold = None

        # Scrublet can fail on tiny/degenerate datasets; keep pipeline usable by
        # falling back to "all singlets" rather than aborting mandatory stage.
        if adata.n_obs < 20 or adata.n_vars < 50:
            doublet_scores, predicted_doublets = self._fallback_all_singlets(
                adata.n_obs,
                reason=f"dataset too small (n_obs={adata.n_obs}, n_vars={adata.n_vars})",
            )
            ctx.metadata["doublet_method"] = "fallback_all_singlets"
        else:
            strategy = self._resolve_doublet_strategy(ctx, adata)
            use_grouped = strategy == "grouped"

            if _RSC_AVAILABLE:
                # --- GPU path ---
                try:
                    if use_grouped:
                        doublet_scores, predicted_doublets, threshold = self._run_rsc_scrublet_grouped(adata, cfg, ctx)
                    else:
                        doublet_scores, predicted_doublets, threshold = self._run_rsc_scrublet_whole(adata, cfg, ctx)
                    if threshold is not None:
                        logger.info("rsc.pp.scrublet auto-threshold: %.4f", threshold)
                except Exception as exc:
                    # HIGH-1 fix (2026-05-19): retry CPU scrublet before
                    # giving up. The previous behaviour silently dropped
                    # doublet detection entirely whenever the RAPIDS path
                    # raised (e.g. rsc.pp.scrublet sparse-dtype mismatch
                    # against int counts), keeping ~2% of cells that should
                    # have been removed. CPU scrublet is the canonical
                    # reference implementation (Wolock et al. 2019, Cell
                    # Systems) so retrying it is scientifically equivalent
                    # to the GPU path's intent.
                    logger.warning("rsc.pp.scrublet failed, retrying with CPU scrublet: %s", exc)
                    try:
                        import scrublet as scr  # type: ignore

                        if use_grouped:
                            doublet_scores, predicted_doublets, threshold = self._run_cpu_scrublet_grouped(
                                adata, scr.Scrublet, cfg, ctx
                            )
                        else:
                            scrub = scr.Scrublet(
                                self._materialize_counts_matrix(adata.X),
                                expected_doublet_rate=cfg.expected_doublet_rate,
                                random_state=random_state,
                            )
                            scrublet_params = {
                                "min_counts": 2,
                                "min_cells": 3,
                                "min_gene_variability_pctl": 85,
                                "n_prin_comps": self._n_prin_comps(adata.n_obs, adata.n_vars),
                            }
                            doublet_scores, predicted_doublets = scrub.scrub_doublets(**scrublet_params)
                            threshold = getattr(scrub, "threshold_", None)
                        ctx.metadata["doublet_method"] = "cpu_scrublet_fallback_from_rsc"
                        ctx.metadata["doublet_rsc_failure_reason"] = str(exc)
                    except Exception as cpu_exc:
                        logger.warning("CPU scrublet retry also failed, falling back to all singlets: %s", cpu_exc)
                        doublet_scores, predicted_doublets = self._fallback_all_singlets(
                            adata.n_obs, reason=f"rsc={exc}; cpu_scrublet={cpu_exc}",
                        )
                        ctx.metadata["doublet_method"] = "fallback_all_singlets"
                        ctx.metadata["doublet_rsc_failure_reason"] = str(exc)
                        ctx.metadata["doublet_cpu_failure_reason"] = str(cpu_exc)
            else:
                # --- CPU fallback path (rsc not installed) ---
                import scrublet as scr

                if use_grouped:
                    doublet_scores, predicted_doublets, threshold = self._run_cpu_scrublet_grouped(adata, scr.Scrublet, cfg, ctx)
                else:
                    scrub = scr.Scrublet(
                        self._materialize_counts_matrix(adata.X),
                        expected_doublet_rate=cfg.expected_doublet_rate,
                        random_state=random_state,
                    )
                    scrublet_params = {
                        "min_counts": 2,
                        "min_cells": 3,
                        "min_gene_variability_pctl": 85,
                        "n_prin_comps": self._n_prin_comps(adata.n_obs, adata.n_vars),
                    }
                    try:
                        doublet_scores, predicted_doublets = scrub.scrub_doublets(**scrublet_params)
                        threshold = getattr(scrub, "threshold_", None)
                        if threshold is not None:
                            logger.info("Scrublet auto-threshold: %.4f", threshold)
                        ctx.metadata["doublet_method"] = "scrublet"
                    except Exception as exc:
                        logger.warning("Scrublet failed, falling back to all singlets: %s", exc)
                        doublet_scores, predicted_doublets = self._fallback_all_singlets(
                            adata.n_obs, reason=str(exc)
                        )
                        ctx.metadata["doublet_method"] = "fallback_all_singlets"

        ctx.metadata["doublet_random_state"] = random_state

        adata.obs["doublet_score"] = doublet_scores
        adata.obs["predicted_doublet"] = predicted_doublets

        n_doublets = int(predicted_doublets.sum())
        ctx.metadata["doublets_detected"] = n_doublets
        ctx.metadata["doublet_rate_pct"] = round(n_doublets / adata.n_obs * 100, 2)

        # Visualize doublet score distribution
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(doublet_scores, bins=50, edgecolor="black", linewidth=0.5)
        if threshold is not None:
            ax.axvline(
                threshold, color="red", linestyle="--",
                label=f"Threshold ({threshold:.3f})",
            )
        ax.set_xlabel("Doublet score")
        ax.set_ylabel("Count")
        ax.set_title(f"Scrublet doublet scores (detected {n_doublets} doublets)")
        if threshold is not None:
            ax.legend()
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "doublet_scores.png", dpi=160, bbox_inches="tight")
        plt.close()

        if cfg.remove_doublets:
            before = adata.n_obs
            adata = adata[~adata.obs["predicted_doublet"]].copy()
            ctx.adata = adata
            ctx.metadata["cells_after_doublet_removal"] = int(adata.n_obs)
            ctx.metadata["doublets_removed"] = int(before - adata.n_obs)
