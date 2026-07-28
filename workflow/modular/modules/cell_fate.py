from __future__ import annotations

import logging
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import diags, issparse

from ..context import PipelineContext


_CELLRANK_FALLBACK_ENV = "SC_ALLOW_CELLRANK_FALLBACK"
_FALLBACK_ENGINE = "fallback_connectivity_diffusion"
_PRIMARY_ENGINE = "cellrank"


__references__ = {
    "Lange_CellRank_2022": {
        "title": "CellRank for directed single-cell fate mapping",
        "authors": "Lange et al.",
        "journal": "Nature Methods",
        "year": "2022",
        "doi": "10.1038/s41592-021-01346-6",
        "description": "Directed transition-probability based fate mapping \u2014 methodology this module's Pearson sampling approximates.",
    },
}


logger = logging.getLogger(__name__)


class CellFateModule:
    """Optional module: probabilistic cell fate mapping.

    Uses CellRank when available for kernel-based fate analysis.
    Falls back to a manual diffusion-based transition probability approach
    built on DPT pseudotime and the neighbor connectivity graph.
    """

    name = "cell_fate"
    requires_keys = {"obs": ["dpt_pseudotime"]}

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Cell fate analysis requires AnnData.")
        if "dpt_pseudotime" not in adata.obs:
            raise ValueError("Cell fate analysis requires dpt_pseudotime in adata.obs.")
        if "connectivities" not in adata.obsp:
            raise ValueError("Cell fate analysis requires a neighbors graph in adata.obsp.")

        engine_used = _PRIMARY_ENGINE
        fallback_reason: str | None = None
        cellrank_version: str | None = None

        try:
            terminal_df, fate_df = self._run_cellrank(adata, ctx)
            try:
                import cellrank as _cr

                cellrank_version = getattr(_cr, "__version__", None)
            except Exception:  # pragma: no cover - version probe only
                cellrank_version = None
        except ImportError as exc:
            # Canonical missing-dep path: always allowed because the cited
            # method (Lange 2022) requires the cellrank package. The fallback
            # is recorded as a separate engine so consumers cannot mistake it
            # for CellRank output.
            logger.info("CellRank not importable (%s); using connectivity-diffusion fallback", exc)
            engine_used = _FALLBACK_ENGINE
            fallback_reason = f"cellrank_import_error:{exc}"
            terminal_df, fate_df = self._run_fallback(adata)
        except Exception as exc:
            # Non-ImportError CellRank failure: this is a silent algorithmic
            # compromise risk (Principle 9). We refuse to silently swap the
            # cited Lange 2022 method for the project-local diffusion proxy
            # unless the operator explicitly opts in with SC_ALLOW_CELLRANK_FALLBACK=1.
            if os.environ.get(_CELLRANK_FALLBACK_ENV, "").strip() != "1":
                raise RuntimeError(
                    "cell_fate: CellRank computation failed and the silent "
                    "fallback to connectivity-diffusion is banned by default. "
                    f"Original CellRank error: {exc!r}. "
                    f"Set {_CELLRANK_FALLBACK_ENV}=1 to acknowledge that the "
                    "downstream fate_probabilities will NOT come from the cited "
                    "Lange 2022 method."
                ) from exc
            logger.error(
                "cell_fate: %s=1 OPT-IN ACKNOWLEDGED; replacing CellRank output "
                "with connectivity-diffusion fallback. Original error: %s",
                _CELLRANK_FALLBACK_ENV,
                exc,
            )
            engine_used = _FALLBACK_ENGINE
            fallback_reason = f"cellrank_runtime_error_opt_in:{exc}"
            ctx.metadata["cell_fate_fallback_opt_in_acknowledged"] = True
            terminal_df, fate_df = self._run_fallback(adata)

        # --- Provenance stamp (always written, before any side-effect tables) ---
        terminal_names = [c for c in fate_df.columns if c != "cell"]
        adata.uns["cell_fate"] = {
            "engine": engine_used,
            "primary_engine": _PRIMARY_ENGINE,
            "cellrank_version": cellrank_version,
            "fallback_reason": fallback_reason,
            "n_terminal_states": len(terminal_df),
            "terminal_state_names": list(terminal_names),
        }
        ctx.metadata["cell_fate_engine"] = engine_used
        ctx.metadata["cell_fate_method_actually_used"] = engine_used
        if fallback_reason is not None:
            ctx.metadata["cell_fate_fallback_reason"] = fallback_reason

        # --- Unified obsm shape across engines (numpy float32 + uns column names) ---
        if terminal_names:
            fate_matrix = fate_df[terminal_names].to_numpy(dtype=np.float32, copy=False)
            adata.obsm["fate_probabilities"] = fate_matrix
            adata.uns["cell_fate"]["fate_probability_columns"] = list(terminal_names)

        # --- Save tables ---
        terminal_df.to_csv(ctx.table_dir / "terminal_states.csv", index=False)
        fate_df.to_csv(ctx.table_dir / "fate_probabilities.csv", index=False)

        n_terminal = len(terminal_df)
        ctx.metadata["n_terminal_states"] = n_terminal
        if n_terminal < 2:
            logger.warning("Fewer than 2 terminal states identified (%d)", n_terminal)

        # --- Figures ---
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
                plt.savefig(ctx.figure_dir / "fate_umap_cellrank.png", bbox_inches="tight")
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

        fate_cols = [f"fate_{tc}" for tc in terminal_clusters]
        # Note: adata.obsm["fate_probabilities"] is written by run() in a
        # unified shape (np.ndarray + column names in adata.uns["cell_fate"])
        # to keep schema consistent across the cellrank and fallback engines.
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
            plt.savefig(ctx.figure_dir / f"fate_umap_{safe_name}.png", bbox_inches="tight")
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
        plt.savefig(ctx.figure_dir / "fate_heatmap.png", bbox_inches="tight")
        plt.close()
