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

logger = logging.getLogger(__name__)


class DoubletDetectionModule:
    """Mandatory module: detect and remove doublets using Scrublet.

    Should run after QC and before clustering. Doublets (two cells captured
    in one droplet) can create artificial intermediate clusters and corrupt
    downstream differential expression results.
    """

    name = "doublet_detection"
    required = True

    @staticmethod
    def _scrublet_params(n_obs: int, n_vars: int) -> dict:
        return {
            "min_counts": 2,
            "min_cells": 3,
            "min_gene_variability_pctl": 85,
            "n_prin_comps": max(2, min(30, n_obs - 1, n_vars - 1)),
        }

    @staticmethod
    def _materialize_counts_matrix(x):
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

    def _run_grouped_scrublet(self, adata, scrublet_cls, cfg, ctx):
        sample_key = 'sample' if 'sample' in adata.obs.columns else None
        if sample_key is None:
            raise ValueError('grouped scrublet requested without sample labels')
        labels = pd.Series(adata.obs[sample_key].astype(str), index=adata.obs_names)
        scores = np.zeros(adata.n_obs, dtype=np.float32)
        predicted = np.zeros(adata.n_obs, dtype=bool)
        thresholds = {}
        for group in labels.unique():
            idx = np.where(labels.values == group)[0]
            if len(idx) < 20:
                continue
            subX = self._materialize_counts_matrix(adata.X[idx])
            scrub = scrublet_cls(subX, expected_doublet_rate=cfg.expected_doublet_rate)
            try:
                s, p = scrub.scrub_doublets(**self._scrublet_params(len(idx), adata.n_vars))
                scores[idx] = s.astype(np.float32, copy=False)
                predicted[idx] = p.astype(bool, copy=False)
                thr = getattr(scrub, 'threshold_', None)
                if thr is not None:
                    thresholds[str(group)] = float(thr)
            except Exception as exc:
                logger.warning('Grouped Scrublet failed for %s, falling back to singlets: %s', group, exc)
        ctx.metadata['doublet_grouped_key'] = sample_key
        ctx.metadata['doublet_grouped_thresholds'] = thresholds
        ctx.metadata['doublet_method'] = 'scrublet_grouped'
        return scores, predicted, None

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Doublet detection requires loaded AnnData.")

        import scrublet as scr

        cfg = ctx.cfg.doublet
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
            if use_grouped:
                doublet_scores, predicted_doublets, threshold = self._run_grouped_scrublet(adata, scr.Scrublet, cfg, ctx)
            else:
                scrub = scr.Scrublet(self._materialize_counts_matrix(adata.X), expected_doublet_rate=cfg.expected_doublet_rate)
                try:
                    doublet_scores, predicted_doublets = scrub.scrub_doublets(**self._scrublet_params(adata.n_obs, adata.n_vars))
                    threshold = getattr(scrub, "threshold_", None)
                    if threshold is not None:
                        logger.info("Scrublet auto-threshold: %.4f", threshold)
                    ctx.metadata["doublet_method"] = "scrublet"
                except Exception as exc:
                    logger.warning("Scrublet failed, falling back to all singlets: %s", exc)
                    doublet_scores, predicted_doublets = self._fallback_all_singlets(
                        adata.n_obs,
                        reason=str(exc),
                    )
                    ctx.metadata["doublet_method"] = "fallback_all_singlets"

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

    @staticmethod
    def _fallback_all_singlets(n_obs: int, reason: str) -> tuple[np.ndarray, np.ndarray]:
        scores = np.zeros(n_obs, dtype=np.float32)
        predicted = np.zeros(n_obs, dtype=bool)
        logger.info("Doublet detection fallback activated: %s", reason)
        return scores, predicted
