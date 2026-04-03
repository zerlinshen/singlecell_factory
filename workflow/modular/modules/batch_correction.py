from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from ..context import PipelineContext


class BatchCorrectionModule:
    """Optional module: batch effect correction for multi-sample integration.

    Supports multiple backends:
    1. Harmony (default) — fast, PCA-based correction
    2. BBKNN — batch-balanced k-nearest-neighbors
    3. Combat — linear model-based correction

    Should run after clustering (specifically after PCA), before downstream analysis.
    The batch key must be present in adata.obs.
    """

    name = "batch_correction"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Batch correction requires AnnData.")

        cfg = ctx.cfg.batch
        batch_key = cfg.batch_key

        if batch_key not in adata.obs.columns:
            raise ValueError(
                f"Batch key '{batch_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns[:20])}"
            )

        n_batches = adata.obs[batch_key].nunique()
        if n_batches < 2:
            ctx.status("batch_correction", True, f"Only 1 batch found in '{batch_key}', skipping")
            return

        ctx.metadata["n_batches"] = n_batches
        ctx.metadata["batch_method"] = cfg.method

        # Pre-correction UMAP (for comparison)
        if "X_umap" in adata.obsm:
            sc.pl.umap(adata, color=[batch_key], show=False)
            plt.savefig(
                ctx.figure_dir / "umap_batch_before.png", dpi=160, bbox_inches="tight",
            )
            plt.close()

        if cfg.method == "harmony":
            self._run_harmony(adata, batch_key, ctx)
        elif cfg.method == "bbknn":
            self._run_bbknn(adata, batch_key, ctx)
        elif cfg.method == "combat":
            self._run_combat(adata, batch_key, ctx)
        elif cfg.method == "scanorama":
            self._run_scanorama(adata, batch_key, ctx)
        else:
            raise ValueError(f"Unknown batch correction method: {cfg.method}")

        # Recompute UMAP + Leiden on corrected representation
        if cfg.method == "harmony":
            sc.pp.neighbors(adata, use_rep="X_pca_harmony")
        elif cfg.method == "scanorama":
            sc.pp.neighbors(adata, use_rep="X_scanorama")
        elif cfg.method != "bbknn":
            sc.pp.neighbors(adata, use_rep="X_pca")

        sc.tl.umap(adata, random_state=ctx.cfg.clustering.random_state)
        sc.tl.leiden(
            adata,
            resolution=ctx.cfg.clustering.leiden_resolution,
            flavor="igraph",
            directed=False,
            random_state=ctx.cfg.clustering.random_state,
        )
        ctx.metadata["n_clusters_after_batch"] = int(adata.obs["leiden"].nunique())

        # Post-correction UMAP
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        sc.pl.umap(adata, color=[batch_key], ax=axes[0], show=False)
        axes[0].set_title(f"After {cfg.method} — batch")
        sc.pl.umap(adata, color=["leiden"], ax=axes[1], show=False, legend_loc="on data")
        axes[1].set_title(f"After {cfg.method} — clusters")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "umap_batch_after.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _run_harmony(adata, batch_key: str, ctx) -> None:
        if "X_pca" not in adata.obsm:
            raise ValueError("Harmony requires PCA (run clustering first).")
        sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca")

    @staticmethod
    def _run_bbknn(adata, batch_key: str, ctx) -> None:
        import bbknn
        if "X_pca" not in adata.obsm:
            raise ValueError("BBKNN requires PCA (run clustering first).")
        bbknn.bbknn(adata, batch_key=batch_key)

    @staticmethod
    def _run_combat(adata, batch_key: str, ctx) -> None:
        sc.pp.combat(adata, key=batch_key)

    @staticmethod
    def _run_scanorama(adata, batch_key: str, ctx) -> None:
        """Scanorama integration (Hie et al., Nature Biotechnology 2019)."""
        import scanorama
        import numpy as np

        if "X_pca" not in adata.obsm:
            raise ValueError("Scanorama requires PCA (run clustering first).")
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
