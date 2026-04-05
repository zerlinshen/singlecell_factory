from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import diags, issparse

from ..context import PipelineContext

logger = logging.getLogger(__name__)


class CellFateModule:
    """Optional module: probabilistic cell fate mapping.

    Uses CellRank when available for kernel-based fate analysis.
    Falls back to a manual diffusion-based transition probability approach
    built on DPT pseudotime and the neighbor connectivity graph.
    """

    name = "cell_fate"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Cell fate analysis requires AnnData.")
        if "dpt_pseudotime" not in adata.obs:
            raise ValueError("Cell fate analysis requires dpt_pseudotime in adata.obs.")
        if "connectivities" not in adata.obsp:
            raise ValueError("Cell fate analysis requires a neighbors graph in adata.obsp.")

        try:
            terminal_df, fate_df = self._run_cellrank(adata, ctx)
        except ImportError:
            logger.info("CellRank not available, using manual fallback")
            terminal_df, fate_df = self._run_fallback(adata)
        except Exception as exc:
            logger.warning("CellRank failed (%s), using manual fallback", exc)
            terminal_df, fate_df = self._run_fallback(adata)

        # --- Save tables ---
        terminal_df.to_csv(ctx.table_dir / "terminal_states.csv", index=False)
        fate_df.to_csv(ctx.table_dir / "fate_probabilities.csv", index=False)

        n_terminal = len(terminal_df)
        ctx.metadata["n_terminal_states"] = n_terminal
        if n_terminal < 2:
            logger.warning("Fewer than 2 terminal states identified (%d)", n_terminal)

        # --- Figures ---
        terminal_names = [c for c in fate_df.columns if c != "cell"]
        if "X_umap" in adata.obsm and terminal_names:
            self._plot_fate_umaps(adata, fate_df, terminal_names, ctx)
        if "leiden" in adata.obs and terminal_names:
            self._plot_fate_heatmap(adata, fate_df, terminal_names, ctx)

    @staticmethod
    def _run_cellrank(adata, ctx: PipelineContext) -> tuple[pd.DataFrame, pd.DataFrame]:
        import cellrank as cr

        pk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
        pk.compute_transition_matrix()

        estimator = cr.estimators.GPCCA(pk)
        estimator.fit()
        estimator.predict_terminal_states()
        estimator.compute_fate_probabilities()

        # Terminal states table
        ts = adata.obs["terminal_states"].dropna()
        terminal_df = pd.DataFrame({
            "cell": ts.index,
            "terminal_state": ts.values,
        })

        # Fate probability matrix
        fate_probs = adata.obsm["lineages_fwd"]
        fate_df = pd.DataFrame(
            fate_probs.X if hasattr(fate_probs, "X") else np.asarray(fate_probs),
            columns=fate_probs.names if hasattr(fate_probs, "names") else [
                f"fate_{i}" for i in range(fate_probs.shape[1])
            ],
        )
        fate_df.insert(0, "cell", adata.obs_names)

        # UMAP plots via CellRank
        if "X_umap" in adata.obsm:
            try:
                estimator.plot_fate_probabilities(same_plot=False, basis="umap", save=False)
                plt.savefig(ctx.figure_dir / "fate_umap_cellrank.png", dpi=160, bbox_inches="tight")
                plt.close()
            except Exception:
                pass

        return terminal_df, fate_df

    @staticmethod
    def _run_fallback(adata) -> tuple[pd.DataFrame, pd.DataFrame]:
        pt = adata.obs["dpt_pseudotime"].values.astype(np.float64)
        clusters = adata.obs["leiden"] if "leiden" in adata.obs else None

        if clusters is None or clusters.nunique() < 2:
            logger.warning("Single cluster or no leiden labels; returning uniform fate.")
            fate_df = pd.DataFrame({"cell": adata.obs_names, "fate_0": 1.0})
            terminal_df = pd.DataFrame({"cluster": ["all"], "mean_pseudotime": [float(np.nanmean(pt))]})
            return terminal_df, fate_df

        # Row-normalize connectivities to get transition matrix
        conn = adata.obsp["connectivities"]
        if not issparse(conn):
            from scipy.sparse import csr_matrix
            conn = csr_matrix(conn)
        row_sums = np.array(conn.sum(axis=1)).flatten()
        T = diags(1.0 / np.maximum(row_sums, 1e-10)) @ conn

        # Identify terminal clusters: those with mean pseudotime in the top quartile
        cluster_labels = clusters.values
        unique_clusters = clusters.unique()
        cluster_mean_pt = {c: float(np.nanmean(pt[cluster_labels == c])) for c in unique_clusters}

        pt_threshold = np.percentile(list(cluster_mean_pt.values()), 75)
        terminal_clusters = [c for c, mpt in cluster_mean_pt.items() if mpt >= pt_threshold]

        if len(terminal_clusters) == 0:
            # Take the single highest-pseudotime cluster
            terminal_clusters = [max(cluster_mean_pt, key=cluster_mean_pt.get)]

        terminal_df = pd.DataFrame([
            {"cluster": c, "mean_pseudotime": round(cluster_mean_pt[c], 6),
             "n_cells": int((cluster_labels == c).sum())}
            for c in terminal_clusters
        ])

        # Compute fate probabilities: for each terminal cluster, sum transition
        # weights from each cell to cells belonging to that cluster.
        fate_matrix = np.zeros((adata.n_obs, len(terminal_clusters)), dtype=np.float64)
        for j, tc in enumerate(terminal_clusters):
            indicator = (cluster_labels == tc).astype(np.float64)
            fate_matrix[:, j] = T.dot(indicator)

        # Normalize rows so probabilities sum to 1
        row_totals = fate_matrix.sum(axis=1, keepdims=True)
        row_totals = np.maximum(row_totals, 1e-10)
        fate_matrix = fate_matrix / row_totals

        # Store in obsm
        fate_cols = [f"fate_{tc}" for tc in terminal_clusters]
        adata.obsm["fate_probabilities"] = pd.DataFrame(
            fate_matrix, index=adata.obs_names, columns=fate_cols,
        )

        fate_df = pd.DataFrame(fate_matrix, columns=fate_cols)
        fate_df.insert(0, "cell", adata.obs_names.values)

        return terminal_df, fate_df

    @staticmethod
    def _plot_fate_umaps(
        adata, fate_df: pd.DataFrame, terminal_names: list[str], ctx: PipelineContext,
    ) -> None:
        umap = adata.obsm["X_umap"]
        for name in terminal_names:
            fig, ax = plt.subplots(figsize=(6, 5))
            sc = ax.scatter(umap[:, 0], umap[:, 1], c=fate_df[name].values, cmap="viridis", s=1, rasterized=True)
            ax.set_title(f"Fate probability: {name}")
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            plt.colorbar(sc, ax=ax)
            plt.tight_layout()
            safe_name = name.replace("/", "_")
            plt.savefig(ctx.figure_dir / f"fate_umap_{safe_name}.png", dpi=160, bbox_inches="tight")
            plt.close()

    @staticmethod
    def _plot_fate_heatmap(
        adata, fate_df: pd.DataFrame, terminal_names: list[str], ctx: PipelineContext,
    ) -> None:
        clusters = adata.obs["leiden"].values
        unique_clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x) if x.isdigit() else x)

        mean_fate = np.zeros((len(unique_clusters), len(terminal_names)))
        for i, cl in enumerate(unique_clusters):
            mask = clusters == cl
            for j, tn in enumerate(terminal_names):
                mean_fate[i, j] = fate_df.loc[mask, tn].mean()

        fig, ax = plt.subplots(figsize=(max(6, len(terminal_names) * 1.2), max(4, len(unique_clusters) * 0.4)))
        im = ax.imshow(mean_fate, aspect="auto", cmap="YlOrRd")
        ax.set_xticks(range(len(terminal_names)))
        ax.set_xticklabels(terminal_names, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(unique_clusters)))
        ax.set_yticklabels(unique_clusters, fontsize=8)
        ax.set_xlabel("Terminal State")
        ax.set_ylabel("Leiden Cluster")
        ax.set_title("Mean Fate Probability per Cluster")
        plt.colorbar(im, ax=ax, label="Probability")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "fate_heatmap.png", dpi=160, bbox_inches="tight")
        plt.close()
