from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from ..context import PipelineContext


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
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=cfg.expected_doublet_rate)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
        )

        adata.obs["doublet_score"] = doublet_scores
        adata.obs["predicted_doublet"] = predicted_doublets

        n_doublets = int(predicted_doublets.sum())
        ctx.metadata["doublets_detected"] = n_doublets
        ctx.metadata["doublet_rate_pct"] = round(n_doublets / adata.n_obs * 100, 2)

        # Visualize doublet score distribution
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(doublet_scores, bins=50, edgecolor="black", linewidth=0.5)
        ax.axvline(
            scrub.threshold_, color="red", linestyle="--",
            label=f"Threshold ({scrub.threshold_:.3f})",
        )
        ax.set_xlabel("Doublet score")
        ax.set_ylabel("Count")
        ax.set_title(f"Scrublet doublet scores (detected {n_doublets} doublets)")
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
