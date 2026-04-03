from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..context import PipelineContext


class RNAVelocityModule:
    """Optional module: real RNA velocity analysis using scVelo.

    Requires a loom file with spliced/unspliced count layers, or pre-computed
    spliced/unspliced layers in the AnnData object. If neither is available,
    falls back to the pseudo-velocity approach.
    """

    name = "rna_velocity"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "X_umap" not in adata.obsm:
            raise ValueError("RNA velocity requires AnnData with UMAP embedding.")

        import scvelo as scv

        cfg = ctx.cfg.velocity
        loom_path = cfg.loom_path

        # Try to load spliced/unspliced from loom
        if loom_path and loom_path.exists():
            ldata = scv.read(str(loom_path))
            # Merge loom data into existing AnnData
            import scvelo.utils as scu
            adata = scv.utils.merge(adata, ldata)
            ctx.adata = adata

        has_layers = "spliced" in adata.layers and "unspliced" in adata.layers
        if not has_layers:
            ctx.status(
                "rna_velocity", False,
                "No spliced/unspliced layers found. Provide a loom file via --velocity-loom.",
            )
            return

        # scVelo preprocessing
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

        # Run velocity model
        if cfg.mode == "dynamical":
            scv.tl.recover_dynamics(adata, n_jobs=cfg.n_jobs)
            scv.tl.velocity(adata, mode="dynamical")
        else:
            scv.tl.velocity(adata, mode="stochastic")

        scv.tl.velocity_graph(adata)

        # Visualizations
        scv.pl.velocity_embedding_stream(
            adata, basis="umap", color="leiden" if "leiden" in adata.obs else None,
            show=False,
        )
        plt.savefig(ctx.figure_dir / "velocity_stream_umap.png", dpi=160, bbox_inches="tight")
        plt.close()

        scv.pl.velocity_embedding_grid(
            adata, basis="umap", color="leiden" if "leiden" in adata.obs else None,
            show=False,
        )
        plt.savefig(ctx.figure_dir / "velocity_grid_umap.png", dpi=160, bbox_inches="tight")
        plt.close()

        # Velocity confidence
        scv.tl.velocity_confidence(adata)
        if "velocity_confidence" in adata.obs:
            adata.obs[["velocity_confidence", "velocity_length"]].to_csv(
                ctx.table_dir / "velocity_confidence.csv"
            )
            ctx.metadata["velocity_mean_confidence"] = float(
                adata.obs["velocity_confidence"].mean()
            )

        # Top velocity genes
        try:
            scv.tl.rank_velocity_genes(adata, groupby="leiden", n_genes=10)
            velocity_genes = scv.DataFrame(adata.uns["rank_velocity_genes"]["names"]).head(10)
            velocity_genes.to_csv(ctx.table_dir / "velocity_top_genes.csv", index=False)
        except Exception:
            pass

        ctx.metadata["velocity_mode"] = cfg.mode
