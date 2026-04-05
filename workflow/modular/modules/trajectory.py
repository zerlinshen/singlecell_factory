from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from ..context import PipelineContext


class TrajectoryModule:
    """Optional module: comprehensive pseudotime trajectory inference.

    Performs:
    - PAGA (Partition-based graph abstraction) for trajectory topology
    - Diffusion pseudotime (DPT) for temporal ordering
    - Gene expression dynamics along pseudotime (top variable genes)
    - Per-cluster pseudotime statistics

    References:
    - Wolf et al., Genome Biology, 2019 (PAGA) DOI: 10.1186/s13059-019-1663-x
    - Haghverdi et al., Nature Methods, 2016 (DPT) DOI: 10.1038/nmeth.3971
    """

    name = "trajectory"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "X_umap" not in adata.obsm:
            raise ValueError("Trajectory requires UMAP embedding.")

        # --- PAGA trajectory graph ---
        if "leiden" in adata.obs:
            sc.tl.paga(adata, groups="leiden")

        # --- Diffusion map + DPT ---
        sc.tl.diffmap(adata)

        root_cluster = ctx.cfg.trajectory_root_cluster
        if root_cluster is not None and "leiden" in adata.obs:
            cluster_mask = adata.obs["leiden"] == str(root_cluster)
            if cluster_mask.any():
                dc1 = adata.obsm["X_diffmap"][:, 0]
                dc1_masked = np.where(cluster_mask, dc1, np.inf)
                adata.uns["iroot"] = int(np.argmin(np.abs(dc1_masked)))
            else:
                adata.uns["iroot"] = int(adata.obsm["X_diffmap"][:, 0].argmin())
        else:
            adata.uns["iroot"] = int(adata.obsm["X_diffmap"][:, 0].argmin())

        sc.tl.dpt(adata)

        # --- Tables ---
        adata.obs[["dpt_pseudotime"]].to_csv(ctx.table_dir / "dpt_pseudotime.csv")

        # Per-cluster pseudotime statistics
        if "leiden" in adata.obs:
            stats = (
                adata.obs[["leiden", "dpt_pseudotime"]]
                .groupby("leiden", observed=True)["dpt_pseudotime"]
                .agg(["count", "mean", "median", "std", "min", "max"])
                .reset_index()
            )
            stats.to_csv(ctx.table_dir / "pseudotime_per_cluster.csv", index=False)

        ctx.metadata["pseudotime_mean"] = round(float(adata.obs["dpt_pseudotime"].mean()), 4)

        # --- Visualizations ---
        self._plot_pseudotime_umap(adata, ctx)
        if "leiden" in adata.obs and "paga" in adata.uns:
            self._plot_paga(adata, ctx)
        self._plot_gene_trends(adata, ctx)
        self._plot_pseudotime_density(adata, ctx)

    @staticmethod
    def _plot_pseudotime_umap(adata, ctx: PipelineContext) -> None:
        """UMAP colored by pseudotime and clusters."""
        ncols = 3 if "leiden" in adata.obs else 2
        fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5))
        sc.pl.umap(adata, color="dpt_pseudotime", ax=axes[0], show=False, cmap="viridis")
        axes[0].set_title("DPT Pseudotime")

        # Diffusion components
        if "X_diffmap" in adata.obsm:
            sc.pl.embedding(
                adata, basis="diffmap", color="dpt_pseudotime",
                components="1,2", ax=axes[1], show=False, cmap="viridis",
            )
            axes[1].set_title("Diffusion Map (DC1 vs DC2)")

        if ncols == 3:
            sc.pl.umap(adata, color="leiden", ax=axes[2], show=False, legend_loc="on data")
            axes[2].set_title("Leiden Clusters")

        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudotime_dpt_umap.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_paga(adata, ctx: PipelineContext) -> None:
        """PAGA trajectory graph overlaid on UMAP."""
        try:
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            sc.pl.paga(
                adata, color="leiden", ax=axes[0], show=False,
                frameon=False, fontsize=9, node_size_scale=1.5,
            )
            axes[0].set_title("PAGA Graph (node = cluster)")
            sc.pl.umap(
                adata, color="dpt_pseudotime", ax=axes[1], show=False, cmap="viridis",
            )
            axes[1].set_title("DPT Pseudotime")
            plt.tight_layout()
            plt.savefig(ctx.figure_dir / "paga_trajectory.png", dpi=160, bbox_inches="tight")
            plt.close()
        except Exception:
            plt.close("all")

    @staticmethod
    def _plot_gene_trends(adata, ctx: PipelineContext) -> None:
        """Heatmap of top variable genes ordered by pseudotime."""
        if "dpt_pseudotime" not in adata.obs:
            return
        pt = adata.obs["dpt_pseudotime"]
        valid = pt[pt < pt.quantile(0.99)].index
        if len(valid) < 50:
            return

        # Find genes most correlated with pseudotime.
        # Limit candidate genes to keep memory/runtime bounded on large datasets.
        expr = adata.raw.to_adata() if adata.raw else adata
        expr_sub = expr[valid]
        pt_vals = pt.loc[valid].values.astype(np.float32)
        n_cells = float(len(pt_vals))
        if n_cells <= 1:
            return

        hv_mask = None
        if "highly_variable" in adata.var.columns:
            hv_mask = adata.var["highly_variable"].reindex(expr_sub.var_names, fill_value=False).to_numpy()
        candidate_idx = np.flatnonzero(hv_mask) if hv_mask is not None else np.arange(expr_sub.n_vars)
        if candidate_idx.size == 0:
            candidate_idx = np.arange(expr_sub.n_vars)

        max_candidates = min(3000, candidate_idx.size)
        X_all = expr_sub.X
        if candidate_idx.size > max_candidates:
            X_cand = X_all[:, candidate_idx]
            if sparse.issparse(X_cand):
                means = np.asarray(X_cand.mean(axis=0)).ravel()
                mean_sq = np.asarray(X_cand.power(2).mean(axis=0)).ravel()
                var = np.maximum(mean_sq - means**2, 0.0)
            else:
                X_cand_arr = np.asarray(X_cand, dtype=np.float32)
                var = X_cand_arr.var(axis=0)
            keep = np.argsort(var)[-max_candidates:]
            candidate_idx = candidate_idx[keep]

        X = X_all[:, candidate_idx]
        pt_mean = float(pt_vals.mean())
        pt_std = float(pt_vals.std())
        if pt_std == 0:
            return

        # Pearson correlation without densifying full cell x gene matrix.
        if sparse.issparse(X):
            X = X.tocsr()
            x_mean = np.asarray(X.mean(axis=0)).ravel()
            x_mean_sq = np.asarray(X.power(2).mean(axis=0)).ravel()
            x_std = np.sqrt(np.maximum(x_mean_sq - x_mean**2, 0.0))
            xy_mean = np.asarray(X.T.dot(pt_vals) / n_cells).ravel()
        else:
            X = np.asarray(X, dtype=np.float32)
            x_mean = X.mean(axis=0)
            x_std = X.std(axis=0)
            xy_mean = (X * pt_vals[:, None]).mean(axis=0)

        denom = x_std * pt_std
        denom[denom == 0] = np.inf
        corr = (xy_mean - x_mean * pt_mean) / denom
        corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)

        top_local_idx = np.argsort(np.abs(corr))[-30:][::-1]
        top_idx = candidate_idx[top_local_idx]
        top_genes = expr_sub.var_names[top_idx].tolist()

        # Save gene-pseudotime correlations
        corr_df = pd.DataFrame({
            "gene": expr_sub.var_names[top_idx],
            "pseudotime_correlation": corr[top_local_idx],
        })
        corr_df.to_csv(ctx.table_dir / "pseudotime_top_genes.csv", index=False)

        # Heatmap: cells ordered by pseudotime, top genes
        sort_idx = np.argsort(pt_vals)
        mat = expr_sub.X[:, top_idx]
        if sparse.issparse(mat):
            mat = mat.toarray()
        mat = np.asarray(mat, dtype=np.float32)
        mat = mat[sort_idx].T

        # Smooth for visualization (rolling mean, window=50)
        from scipy.ndimage import uniform_filter1d

        smoothed = uniform_filter1d(mat.astype(np.float64), size=min(50, mat.shape[1] // 5 or 1), axis=1)

        fig, ax = plt.subplots(figsize=(14, max(6, len(top_genes) * 0.25)))
        im = ax.imshow(smoothed, aspect="auto", cmap="RdBu_r", interpolation="none")
        ax.set_yticks(range(len(top_genes)))
        ax.set_yticklabels(top_genes, fontsize=7)
        ax.set_xlabel("Cells (ordered by pseudotime)")
        ax.set_title("Top 30 genes correlated with pseudotime")
        plt.colorbar(im, ax=ax, label="Smoothed expression", shrink=0.6)
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudotime_gene_heatmap.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_pseudotime_density(adata, ctx: PipelineContext) -> None:
        """Pseudotime distribution per cluster (violin/box)."""
        if "leiden" not in adata.obs or "dpt_pseudotime" not in adata.obs:
            return
        fig, ax = plt.subplots(figsize=(max(8, adata.obs["leiden"].nunique() * 0.8), 5))
        clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x) if x.isdigit() else x)
        data = [adata.obs.loc[adata.obs["leiden"] == c, "dpt_pseudotime"].dropna().values for c in clusters]
        parts = ax.violinplot(data, showmedians=True, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_alpha(0.6)
        ax.set_xticks(range(1, len(clusters) + 1))
        ax.set_xticklabels(clusters)
        ax.set_xlabel("Leiden Cluster")
        ax.set_ylabel("DPT Pseudotime")
        ax.set_title("Pseudotime Distribution per Cluster")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudotime_violin_per_cluster.png", dpi=160, bbox_inches="tight")
        plt.close()
