from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

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
            scrub = scr.Scrublet(adata.X, expected_doublet_rate=cfg.expected_doublet_rate)
            try:
                doublet_scores, predicted_doublets = scrub.scrub_doublets(
                    min_counts=2,
                    min_cells=3,
                    min_gene_variability_pctl=85,
                    n_prin_comps=max(2, min(30, adata.n_obs - 1, adata.n_vars - 1)),
                )
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
