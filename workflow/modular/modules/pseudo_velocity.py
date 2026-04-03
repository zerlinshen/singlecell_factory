from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from ..context import PipelineContext


class PseudoVelocityModule:
    """Optional module: pseudo-RNA velocity from local pseudotime gradients.

    Computes velocity vectors in UMAP space from k-NN pseudotime gradients,
    then produces arrow plots, stream plots, speed overlays, and per-cluster stats.
    """

    name = "pseudo_velocity"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "X_umap" not in adata.obsm:
            raise ValueError("Pseudo velocity requires UMAP embedding.")
        if "dpt_pseudotime" not in adata.obs:
            raise ValueError("Pseudo velocity requires pseudotime (run trajectory first).")

        umap = adata.obsm["X_umap"]
        pseudotime = adata.obs["dpt_pseudotime"].values.astype(np.float32)
        valid = np.isfinite(pseudotime)
        if valid.sum() < 50:
            raise ValueError("Too few cells with valid pseudotime.")

        velocity = compute_pseudo_velocity_vectorized(umap[valid], pseudotime[valid])
        speed = np.linalg.norm(velocity, axis=1)

        full_velocity = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
        full_speed = np.full(adata.n_obs, np.nan, dtype=np.float32)
        valid_idx = np.where(valid)[0]
        full_velocity[valid_idx] = velocity
        full_speed[valid_idx] = speed

        adata.obsm["X_pseudo_velocity"] = full_velocity
        adata.obs["pseudo_velocity_speed"] = full_speed

        # --- Tables ---
        pd.DataFrame({"barcode": adata.obs_names, "speed": full_speed}).to_csv(
            ctx.table_dir / "pseudo_velocity_speed.csv", index=False,
        )
        if "leiden" in adata.obs:
            stats = (
                adata.obs[["leiden", "pseudo_velocity_speed"]]
                .groupby("leiden", observed=True)["pseudo_velocity_speed"]
                .agg(["count", "mean", "median", "std"])
                .reset_index()
            )
            stats.to_csv(ctx.table_dir / "pseudo_velocity_per_cluster.csv", index=False)

        ctx.metadata["pseudo_velocity_mean_speed"] = round(float(np.nanmean(speed)), 4)

        # --- Visualizations ---
        umap_valid = umap[valid]
        self._plot_speed_umap(adata, ctx)
        self._plot_arrows(adata, umap_valid, velocity, ctx)
        self._plot_stream(adata, umap_valid, velocity, ctx)
        if "leiden" in adata.obs:
            self._plot_speed_boxplot(adata, ctx)

    @staticmethod
    def _plot_speed_umap(adata, ctx: PipelineContext) -> None:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        sc0 = axes[0].scatter(
            adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
            c=adata.obs["pseudo_velocity_speed"], cmap="YlOrRd", s=3, alpha=0.7,
        )
        plt.colorbar(sc0, ax=axes[0], label="Speed")
        axes[0].set_title("Pseudo-velocity Speed")
        axes[0].set_xlabel("UMAP1"); axes[0].set_ylabel("UMAP2")

        sc1 = axes[1].scatter(
            adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
            c=adata.obs["dpt_pseudotime"], cmap="viridis", s=3, alpha=0.7,
        )
        plt.colorbar(sc1, ax=axes[1], label="Pseudotime")
        axes[1].set_title("DPT Pseudotime")
        axes[1].set_xlabel("UMAP1"); axes[1].set_ylabel("UMAP2")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudo_velocity_speed_umap.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_arrows(adata, umap_valid, velocity, ctx: PipelineContext) -> None:
        n = len(umap_valid)
        step = max(1, n // 500)
        sel = slice(None, None, step)
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
                   c="lightgrey", s=1, alpha=0.3)
        ax.quiver(
            umap_valid[sel, 0], umap_valid[sel, 1],
            velocity[sel, 0], velocity[sel, 1],
            angles="xy", scale_units="xy", scale=3,
            color="black", alpha=0.6, width=0.002, headwidth=4,
        )
        ax.set_title("Pseudo-velocity Arrows")
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudo_velocity_arrows.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_stream(adata, umap_valid, velocity, ctx: PipelineContext) -> None:
        from scipy.spatial import cKDTree

        x_min, x_max = umap_valid[:, 0].min(), umap_valid[:, 0].max()
        y_min, y_max = umap_valid[:, 1].min(), umap_valid[:, 1].max()
        pad_x, pad_y = (x_max - x_min) * 0.05, (y_max - y_min) * 0.05
        ng = 40
        xi = np.linspace(x_min - pad_x, x_max + pad_x, ng)
        yi = np.linspace(y_min - pad_y, y_max + pad_y, ng)
        Xi, Yi = np.meshgrid(xi, yi)

        tree = cKDTree(umap_valid)
        grid_pts = np.column_stack([Xi.ravel(), Yi.ravel()])
        dist, idx = tree.query(grid_pts, k=10)
        dist = np.maximum(dist, 1e-6)
        w = 1.0 / dist
        w /= w.sum(axis=1, keepdims=True)
        Vx = (velocity[idx, 0] * w).sum(axis=1).reshape(ng, ng)
        Vy = (velocity[idx, 1] * w).sum(axis=1).reshape(ng, ng)
        min_d = dist[:, 0].reshape(ng, ng)
        thresh = np.percentile(dist[:, 0], 80)
        Vx[min_d > thresh] = np.nan
        Vy[min_d > thresh] = np.nan

        fig, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
                   c="lightgrey", s=1, alpha=0.3)
        try:
            speed_g = np.sqrt(np.nan_to_num(Vx)**2 + np.nan_to_num(Vy)**2)
            ax.streamplot(xi, yi, Vx, Vy, color=speed_g, cmap="coolwarm",
                          linewidth=1, density=1.5, arrowsize=1.2)
        except Exception:
            pass
        ax.set_title("Pseudo-velocity Stream Plot")
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudo_velocity_stream.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_speed_boxplot(adata, ctx: PipelineContext) -> None:
        fig, ax = plt.subplots(figsize=(max(8, adata.obs["leiden"].nunique() * 0.8), 5))
        clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x) if x.isdigit() else x)
        data = [adata.obs.loc[adata.obs["leiden"] == c, "pseudo_velocity_speed"].dropna().values
                for c in clusters]
        bp = ax.boxplot(data, labels=clusters, patch_artist=True, showfliers=False)
        for patch in bp["boxes"]:
            patch.set_facecolor("steelblue"); patch.set_alpha(0.6)
        ax.set_xlabel("Leiden Cluster"); ax.set_ylabel("Speed")
        ax.set_title("Pseudo-velocity Speed per Cluster")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudo_velocity_speed_boxplot.png", dpi=160, bbox_inches="tight")
        plt.close()


def compute_pseudo_velocity_vectorized(umap, pseudotime, n_neighbors: int = 15):
    """Compute pseudo velocity using vectorized neighbor operations."""
    knn = NearestNeighbors(n_neighbors=n_neighbors)
    knn.fit(umap)
    _, idx = knn.kneighbors(umap)
    nbr = idx[:, 1:]
    direction = umap[nbr] - umap[:, None, :]
    dt = pseudotime[nbr] - pseudotime[:, None]
    return (direction * dt[:, :, None]).mean(axis=1)
