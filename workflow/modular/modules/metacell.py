from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)


class MetacellModule:
    """Optional module: metacell aggregation to reduce noise and dataset size.

    This module is INFORMATIONAL -- it does NOT replace adata.X.
    The metacell AnnData is a separate output for optional downstream use.
    """

    name = "metacell"

    def run(self, ctx: PipelineContext) -> None:
        import anndata as ad
        from scipy import sparse

        adata = ctx.adata
        if adata is None:
            raise ValueError("Metacell module requires AnnData.")
        if "X_pca" not in adata.obsm:
            raise ValueError("Metacell module requires PCA embedding (run clustering first).")

        n_metacells = max(50, adata.n_obs // 75)
        logger.info(
            "Computing %d metacells for %d cells (ratio ~%d:1)",
            n_metacells, adata.n_obs, adata.n_obs // n_metacells,
        )

        labels = self._try_seacells(adata, n_metacells)
        if labels is None:
            labels = self._fallback_kmeans(adata, n_metacells)

        adata.obs["metacell"] = pd.Categorical(labels.astype(str))

        # Build assignment table
        assign_df = pd.DataFrame({
            "cell_barcode": adata.obs_names,
            "metacell_id": labels,
        })
        assign_df.to_csv(ctx.table_dir / "metacell_assignments.csv", index=False)

        # Aggregate into metacell AnnData
        mc_adata, summary_df = self._aggregate(adata, labels, n_metacells)
        summary_df.to_csv(ctx.table_dir / "metacell_summary.csv", index=False)
        mc_adata.write(ctx.table_dir / "metacells.h5ad")

        ctx.metadata["n_metacells"] = int(n_metacells)
        ctx.metadata["metacell_median_size"] = int(np.median(summary_df["n_cells"]))

        # Visualizations
        self._plot_umap(adata, ctx)
        self._plot_size_hist(summary_df, ctx)

    @staticmethod
    def _try_seacells(adata, n_metacells: int) -> np.ndarray | None:
        """Attempt metacell construction via SEACells."""
        try:
            from SEACells.core import SEACells as SEACellsModel

            model = SEACellsModel(
                adata, build_kernel_on="X_pca",
                n_SEACells=n_metacells, n_waypoint_eigs=10,
            )
            model.construct_kernel_matrix()
            model.initialize_archetypes()
            model.fit(min_iter=10, max_iter=50)
            labels = model.get_soft_assignments().idxmax(axis=1).values
            # Map SEACell identifiers to integer labels
            unique = np.unique(labels)
            mapping = {v: i for i, v in enumerate(unique)}
            return np.array([mapping[v] for v in labels])
        except Exception as exc:
            logger.info("SEACells unavailable or failed (%s), using MiniBatchKMeans fallback.", exc)
            return None

    @staticmethod
    def _fallback_kmeans(adata, n_metacells: int) -> np.ndarray:
        """MiniBatchKMeans metacell construction on PCA embedding."""
        from sklearn.cluster import MiniBatchKMeans

        km = MiniBatchKMeans(
            n_clusters=n_metacells, random_state=0,
            batch_size=min(1024, max(256, adata.n_obs // 10)),
        )
        labels = km.fit_predict(adata.obsm["X_pca"])
        return labels

    @staticmethod
    def _aggregate(adata, labels: np.ndarray, n_metacells: int):
        """Aggregate cells into metacell-level AnnData and summary table."""
        import anndata as ad
        from scipy import sparse

        X = adata.X
        records = []
        mc_X_rows = []

        for mc_id in range(n_metacells):
            mask = labels == mc_id
            n_cells = int(mask.sum())
            if n_cells == 0:
                continue

            # Sum raw counts for this metacell
            if sparse.issparse(X):
                row_sum = np.asarray(X[mask].sum(axis=0)).ravel()
            else:
                row_sum = X[mask].sum(axis=0).ravel()
            mc_X_rows.append(row_sum)

            sub = adata.obs.loc[mask]
            majority_ct = sub["cell_type"].mode().iloc[0] if "cell_type" in sub.columns else "NA"
            majority_leiden = sub["leiden"].mode().iloc[0] if "leiden" in sub.columns else "NA"

            rec = {
                "metacell_id": mc_id, "n_cells": n_cells,
                "majority_cell_type": majority_ct, "majority_leiden": majority_leiden,
            }
            records.append(rec)

        summary_df = pd.DataFrame(records)
        mc_X = np.vstack(mc_X_rows).astype(np.float32)
        mc_adata = ad.AnnData(
            X=mc_X,
            var=adata.var[[]].copy(),
            obs=summary_df.set_index(summary_df["metacell_id"].astype(str)),
        )
        mc_adata.obs_names = [f"MC_{i}" for i in summary_df["metacell_id"]]
        return mc_adata, summary_df

    @staticmethod
    def _plot_umap(adata, ctx: PipelineContext) -> None:
        """UMAP colored by metacell assignment (density-style for many groups)."""
        if "X_umap" not in adata.obsm:
            logger.warning("No UMAP embedding found; skipping metacell UMAP plot.")
            return
        n_mc = adata.obs["metacell"].nunique()
        fig, ax = plt.subplots(figsize=(8, 6))
        umap = adata.obsm["X_umap"]
        if n_mc > 30:
            # Too many metacells for a clean legend; show as scatter with continuous color
            mc_int = adata.obs["metacell"].cat.codes.values
            scatter = ax.scatter(
                umap[:, 0], umap[:, 1], c=mc_int,
                cmap="nipy_spectral", s=1, alpha=0.5, rasterized=True,
            )
            plt.colorbar(scatter, ax=ax, label="Metacell ID")
        else:
            for mc in adata.obs["metacell"].cat.categories:
                mask = adata.obs["metacell"] == mc
                ax.scatter(umap[mask, 0], umap[mask, 1], s=1, alpha=0.5, label=mc)
            ax.legend(fontsize=6, markerscale=3, bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title(f"Metacell assignments (n={n_mc})")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "metacell_umap.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_size_hist(summary_df: pd.DataFrame, ctx: PipelineContext) -> None:
        """Histogram of metacell sizes."""
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(summary_df["n_cells"], bins=min(50, len(summary_df) // 2 + 1), edgecolor="black", linewidth=0.5)
        ax.axvline(summary_df["n_cells"].median(), color="red", linestyle="--", label=f"median={summary_df['n_cells'].median():.0f}")
        ax.set_xlabel("Cells per metacell")
        ax.set_ylabel("Count")
        ax.set_title("Metacell size distribution")
        ax.legend()
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "metacell_size_hist.png", dpi=160, bbox_inches="tight")
        plt.close()
