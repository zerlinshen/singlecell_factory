from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist

from ..context import PipelineContext


class EvolutionModule:
    """Optional module: tumor clonal evolution analysis.

    Integrates CNV profiles, pseudotime, and expression to reconstruct
    the clonal architecture and evolutionary trajectory of tumor cells.

    Produces:
    - CNV-based clonal clustering (hierarchical clustering on CNV profiles)
    - Clone frequency composition per Leiden cluster
    - Pseudotime-ordered clonal evolution timeline
    - Clone-specific gene expression signatures
    - Phylogenetic dendrogram of clones

    References:
    - Patel et al., Science, 2014. DOI: 10.1126/science.1254257
    - Tirosh et al., Science, 2016. DOI: 10.1126/science.aad0501
    - Gao et al., Nature Biotechnology, 2021 (CopyKAT). DOI: 10.1038/s41587-020-00795-2
    """

    name = "evolution"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Evolution analysis requires AnnData.")

        has_cnv = "X_cnv" in adata.obsm or "cnv_score" in adata.obs
        has_pt = "dpt_pseudotime" in adata.obs
        if not has_cnv:
            raise ValueError(
                "Evolution analysis requires CNV data. Run cnv_inference first."
            )

        # --- Step 1: Clonal clustering from CNV profiles ---
        clones, n_clones = self._cluster_clones(adata)
        adata.obs["clone"] = clones

        # --- Step 2: Tables ---
        # Clone composition
        clone_info = adata.obs[["clone"]].copy()
        if "leiden" in adata.obs:
            clone_info["leiden"] = adata.obs["leiden"]
        if "cell_type" in adata.obs:
            clone_info["cell_type"] = adata.obs["cell_type"]
        if has_pt:
            clone_info["dpt_pseudotime"] = adata.obs["dpt_pseudotime"]
        if "cnv_score" in adata.obs:
            clone_info["cnv_score"] = adata.obs["cnv_score"]
        clone_info.to_csv(ctx.table_dir / "evolution_clone_assignment.csv")

        # Per-clone statistics
        clone_stats = self._compute_clone_stats(adata, has_pt)
        clone_stats.to_csv(ctx.table_dir / "evolution_clone_stats.csv", index=False)

        # Clone-specific DE genes (top markers per clone)
        clone_markers = self._clone_markers(adata)
        if clone_markers is not None:
            clone_markers.to_csv(ctx.table_dir / "evolution_clone_markers.csv", index=False)

        ctx.metadata["n_clones"] = n_clones
        ctx.metadata["clone_sizes"] = adata.obs["clone"].value_counts().to_dict()

        # --- Step 3: Visualizations ---
        self._plot_clone_umap(adata, ctx)
        self._plot_clone_composition(adata, ctx)
        if has_pt:
            self._plot_evolution_timeline(adata, ctx)
        if "X_cnv" in adata.obsm:
            self._plot_phylo_dendrogram(adata, ctx)
        self._plot_cnv_by_clone(adata, ctx)

    @staticmethod
    def _cluster_clones(adata) -> tuple[np.ndarray, int]:
        """Cluster cells into clones based on CNV profiles."""
        if "X_cnv" in adata.obsm:
            cnv_mat = adata.obsm["X_cnv"]
            if hasattr(cnv_mat, "toarray"):
                cnv_mat = cnv_mat.toarray()
            cnv_mat = np.asarray(cnv_mat, dtype=np.float32)

            # Hierarchical clustering on CNV profiles
            if cnv_mat.shape[0] > 10000:
                # Subsample for distance computation, assign rest by nearest centroid
                rng = np.random.RandomState(42)
                sub_idx = rng.choice(cnv_mat.shape[0], 5000, replace=False)
                dist = pdist(cnv_mat[sub_idx], metric="correlation")
                Z = linkage(dist, method="ward")
                n_clones = min(6, max(2, int(np.sqrt(cnv_mat.shape[0] / 500))))
                sub_labels = fcluster(Z, t=n_clones, criterion="maxclust")

                # Assign remaining cells to nearest clone centroid
                centroids = np.array([
                    cnv_mat[sub_idx[sub_labels == c]].mean(axis=0)
                    for c in range(1, n_clones + 1)
                ])
                all_labels = np.zeros(cnv_mat.shape[0], dtype=int)
                all_labels[sub_idx] = sub_labels
                other_idx = np.setdiff1d(np.arange(cnv_mat.shape[0]), sub_idx)
                if len(other_idx) > 0:
                    dists = np.array([
                        np.linalg.norm(cnv_mat[other_idx] - c, axis=1)
                        for c in centroids
                    ])
                    all_labels[other_idx] = dists.argmin(axis=0) + 1
                labels = all_labels
            else:
                dist = pdist(cnv_mat, metric="correlation")
                Z = linkage(dist, method="ward")
                n_clones = min(6, max(2, int(np.sqrt(cnv_mat.shape[0] / 200))))
                labels = fcluster(Z, t=n_clones, criterion="maxclust")
        else:
            # Fallback: cluster by CNV score quantiles
            scores = adata.obs["cnv_score"].values
            n_clones = 3
            thresholds = np.percentile(scores, [33, 66])
            labels = np.digitize(scores, thresholds) + 1

        clone_names = np.array([f"Clone_{i}" for i in labels])
        return clone_names, int(np.max(labels))

    @staticmethod
    def _compute_clone_stats(adata, has_pt: bool) -> pd.DataFrame:
        records = []
        for clone in sorted(adata.obs["clone"].unique()):
            mask = adata.obs["clone"] == clone
            row = {
                "clone": clone,
                "n_cells": int(mask.sum()),
                "fraction": round(mask.sum() / adata.n_obs, 4),
            }
            if "cnv_score" in adata.obs:
                row["mean_cnv_score"] = round(float(adata.obs.loc[mask, "cnv_score"].mean()), 4)
            if has_pt:
                pt = adata.obs.loc[mask, "dpt_pseudotime"]
                pt = pt[np.isfinite(pt)]
                if len(pt) > 0:
                    row["mean_pseudotime"] = round(float(pt.mean()), 4)
                    row["median_pseudotime"] = round(float(pt.median()), 4)
            if "cell_type" in adata.obs:
                ct_counts = adata.obs.loc[mask, "cell_type"].value_counts()
                row["dominant_cell_type"] = ct_counts.index[0] if len(ct_counts) > 0 else "Unknown"
            records.append(row)
        return pd.DataFrame(records)

    @staticmethod
    def _clone_markers(adata) -> pd.DataFrame | None:
        """Find top DE genes per clone."""
        if adata.obs["clone"].nunique() < 2:
            return None
        try:
            sc.tl.rank_genes_groups(
                adata, groupby="clone", method="wilcoxon",
                use_raw=False, n_genes=50, key_added="rank_genes_clone",
            )
            markers = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_clone")
            return markers[markers["pvals_adj"] < 0.05].head(200)
        except Exception:
            return None

    @staticmethod
    def _plot_clone_umap(adata, ctx: PipelineContext) -> None:
        ncols = 2 + ("cnv_score" in adata.obs)
        fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5))
        sc.pl.umap(adata, color="clone", ax=axes[0], show=False)
        axes[0].set_title("Clone Assignment")
        if "leiden" in adata.obs:
            sc.pl.umap(adata, color="leiden", ax=axes[1], show=False, legend_loc="on data")
            axes[1].set_title("Leiden Clusters")
        if "cnv_score" in adata.obs and ncols > 2:
            sc.pl.umap(adata, color="cnv_score", ax=axes[ncols - 1], show=False, cmap="Reds")
            axes[ncols - 1].set_title("CNV Score")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "evolution_clone_umap.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_clone_composition(adata, ctx: PipelineContext) -> None:
        if "leiden" not in adata.obs:
            return
        ct = pd.crosstab(adata.obs["leiden"], adata.obs["clone"], normalize="index")
        ct.plot(kind="bar", stacked=True, figsize=(max(8, len(ct) * 0.8), 5), colormap="Set2")
        plt.ylabel("Fraction")
        plt.xlabel("Leiden Cluster")
        plt.title("Clone Composition per Cluster")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "evolution_clone_composition.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_evolution_timeline(adata, ctx: PipelineContext) -> None:
        """Pseudotime distribution per clone — reveals evolutionary ordering."""
        pt = adata.obs[["clone", "dpt_pseudotime"]].dropna()
        if len(pt) < 50:
            return

        clones = sorted(pt["clone"].unique())
        # Order clones by mean pseudotime
        clone_order = pt.groupby("clone")["dpt_pseudotime"].mean().sort_values().index.tolist()

        fig, axes = plt.subplots(1, 2, figsize=(16, 5))

        # Violin plot
        data = [pt.loc[pt["clone"] == c, "dpt_pseudotime"].values for c in clone_order]
        parts = axes[0].violinplot(data, showmedians=True, showextrema=False)
        for i, pc in enumerate(parts["bodies"]):
            pc.set_alpha(0.6)
        axes[0].set_xticks(range(1, len(clone_order) + 1))
        axes[0].set_xticklabels(clone_order, rotation=45, ha="right")
        axes[0].set_ylabel("DPT Pseudotime")
        axes[0].set_title("Pseudotime per Clone (evolutionary ordering)")

        # Density plot
        for c in clone_order:
            vals = pt.loc[pt["clone"] == c, "dpt_pseudotime"].values
            from scipy.stats import gaussian_kde
            try:
                kde = gaussian_kde(vals[np.isfinite(vals)])
                x = np.linspace(0, vals[np.isfinite(vals)].max(), 200)
                axes[1].plot(x, kde(x), label=c, linewidth=2)
            except Exception:
                pass
        axes[1].set_xlabel("DPT Pseudotime")
        axes[1].set_ylabel("Density")
        axes[1].set_title("Clone Density along Pseudotime")
        axes[1].legend(fontsize=8)

        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "evolution_timeline.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_phylo_dendrogram(adata, ctx: PipelineContext) -> None:
        """Hierarchical dendrogram of clone CNV profiles (phylogenetic tree proxy)."""
        cnv = adata.obsm["X_cnv"]
        if hasattr(cnv, "toarray"):
            cnv = cnv.toarray()
        cnv = np.asarray(cnv, dtype=np.float32)

        # Compute clone centroids
        clones = sorted(adata.obs["clone"].unique())
        centroids = np.array([
            cnv[adata.obs["clone"] == c].mean(axis=0) for c in clones
        ])

        if len(clones) < 2:
            return

        dist = pdist(centroids, metric="correlation")
        Z = linkage(dist, method="ward")

        fig, ax = plt.subplots(figsize=(max(6, len(clones) * 1.5), 5))
        dendrogram(Z, labels=clones, ax=ax, leaf_font_size=10)
        ax.set_title("Clone Phylogenetic Dendrogram (CNV-based)")
        ax.set_ylabel("Ward Distance")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "evolution_phylo_dendrogram.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_cnv_by_clone(adata, ctx: PipelineContext) -> None:
        """CNV score distribution per clone."""
        if "cnv_score" not in adata.obs:
            return
        fig, ax = plt.subplots(figsize=(8, 5))
        clones = sorted(adata.obs["clone"].unique())
        data = [adata.obs.loc[adata.obs["clone"] == c, "cnv_score"].values for c in clones]
        bp = ax.boxplot(data, labels=clones, patch_artist=True, showfliers=False)
        colors = plt.cm.Set2(np.linspace(0, 1, len(clones)))
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_xlabel("Clone")
        ax.set_ylabel("CNV Score")
        ax.set_title("CNV Score Distribution per Clone")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "evolution_cnv_by_clone.png", dpi=160, bbox_inches="tight")
        plt.close()
