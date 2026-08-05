from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ._scanpy_compat import import_scanpy_or_stub

sc = import_scanpy_or_stub()
from scipy import sparse

from .._gene_symbols import resolve_expression_axis
from ..context import PipelineContext

logger = logging.getLogger(__name__)



__references__ = {
    "Wolf_PAGA_2019": {
        "title": "PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells",
        "authors": "Wolf et al.",
        "journal": "Genome Biology",
        "year": "2019",
        "doi": "10.1186/s13059-019-1663-x",
        "description": "PAGA topology-preserving abstraction. sc.tl.paga used here.",
    },
    "Haghverdi_DPT_2016": {
        "title": "Diffusion pseudotime robustly reconstructs lineage branching",
        "authors": "Haghverdi, Buttner, Wolf et al.",
        "journal": "Nature Methods",
        "year": "2016",
        "doi": "10.1038/nmeth.3971",
        "description": "Diffusion pseudotime (sc.tl.dpt) used for ordering cells along trajectories.",
    },
    "Hao_WNN_2021": {
        "title": "Integrated analysis of multimodal single-cell data",
        "authors": "Hao et al.",
        "journal": "Cell",
        "year": "2021",
        "doi": "10.1016/j.cell.2021.04.048",
        "description": "WNN (Weighted Nearest Neighbor) joint embedding. X_wnn is preferred over X_pca for pseudotime when available, as it captures multi-modal cell state more faithfully.",
    },
    "Setty_Palantir_2019": {
        "title": "Characterization of cell fate probabilities in single-cell data with Palantir",
        "authors": "Setty et al.",
        "journal": "Nature Biotechnology",
        "year": "2019",
        "doi": "10.1038/s41587-019-0068-4",
        "description": "Pseudotime computed over a joint manifold embedding (X_wnn or X_pca). Supports the rationale for running DPT on the best available joint embedding rather than UMAP.",
    },
}



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
    requires_keys = {"obsm_any_of": ["X_wnn", "X_pca"]}
    provides_keys = {"obs": ["dpt_pseudotime"]}

    def run(self, ctx: PipelineContext) -> None:
        from .._contract_violation import ModuleContractError

        adata = ctx.adata
        if adata is None:
            raise ValueError("Trajectory requires AnnData.")

        if "X_wnn" in adata.obsm:
            embedding_key = "X_wnn"
        elif "X_pca" in adata.obsm:
            embedding_key = "X_pca"
        else:
            raise ModuleContractError(
                "trajectory: neither adata.obsm['X_wnn'] nor adata.obsm['X_pca'] is present. "
                "Run multimodal_integration (for X_wnn) or clustering (for X_pca) first."
            )

        ctx.metadata["trajectory_embedding_key"] = embedding_key

        # Neighbor graph must match the embedding used for PAGA/DPT.
        # Preferring X_wnn while reusing a PCA-neighbors graph (historical bug)
        # misattributes multiome geometry. Rebuild when WNN is selected.
        # PCA path keeps the existing clustering neighbors when present.
        neighbor_rep = embedding_key
        if embedding_key == "X_wnn":
            n_neighbors = int(getattr(ctx.cfg, "n_neighbors", 15) or 15)
            random_state = int(getattr(ctx, "random_state", 0) or 0)
            sc.pp.neighbors(
                adata,
                use_rep="X_wnn",
                n_neighbors=n_neighbors,
                random_state=random_state,
            )
            ctx.metadata["trajectory_neighbors_rebuilt_on"] = "X_wnn"
            neighbor_rep = "X_wnn"
        else:
            # Clustering normally left neighbors on X_pca; do not silently
            # rebuild unless the graph is missing (fail-soft rebuild on PCA).
            has_neighbors = (
                "neighbors" in getattr(adata, "uns", {})
                and "connectivities" in getattr(adata, "obsp", {})
            )
            if not has_neighbors:
                n_neighbors = int(getattr(ctx.cfg, "n_neighbors", 15) or 15)
                random_state = int(getattr(ctx, "random_state", 0) or 0)
                sc.pp.neighbors(
                    adata,
                    use_rep="X_pca",
                    n_neighbors=n_neighbors,
                    random_state=random_state,
                )
                ctx.metadata["trajectory_neighbors_rebuilt_on"] = "X_pca_missing_graph"
            else:
                ctx.metadata["trajectory_neighbors_rebuilt_on"] = "existing"
            neighbor_rep = "X_pca"
        ctx.metadata["trajectory_neighbor_rep"] = neighbor_rep

        # --- PAGA trajectory graph ---
        if "leiden" in adata.obs:
            sc.tl.paga(adata, groups="leiden")

        # --- Diffusion map + DPT ---
        sc.tl.diffmap(adata)

        root_cluster = ctx.cfg.trajectory_root_cluster
        root_justification = getattr(ctx.cfg, "trajectory_root_justification", None)
        root_justification = (
            str(root_justification).strip() if root_justification is not None else ""
        )
        root_cluster_found = False
        if root_cluster is not None and "leiden" in adata.obs:
            cluster_mask = adata.obs["leiden"] == str(root_cluster)
            if cluster_mask.any():
                dc1 = adata.obsm["X_diffmap"][:, 0]
                dc1_masked = np.where(cluster_mask, dc1, np.inf)
                adata.uns["iroot"] = int(np.argmin(np.abs(dc1_masked)))
                root_cluster_found = True
            else:
                adata.uns["iroot"] = int(adata.obsm["X_diffmap"][:, 0].argmin())
        else:
            adata.uns["iroot"] = int(adata.obsm["X_diffmap"][:, 0].argmin())

        sc.tl.dpt(adata)

        trajectory_claimable = bool(root_cluster_found and root_justification)
        if trajectory_claimable:
            claim_status = "claimable_biologically_justified_root"
            root_selection_method = "explicit_leiden_cluster"
        else:
            claim_status = "non_claimable_missing_biologically_justified_root"
            root_selection_method = "exploratory_diffmap_minimum_fallback"

        # M-3 audit fix (2026-05-22): persist trajectory provenance in
        # adata.uns so downstream consumers (cell_fate, pseudo_velocity,
        # bundle export) can tell which embedding was used to compute
        # dpt_pseudotime — PCA-only vs WNN joint embedding meaningfully
        # changes the interpretation of any fate or velocity claim.
        try:
            import scanpy as _sc

            scanpy_version = getattr(_sc, "__version__", None)
        except Exception:  # pragma: no cover - version probe only
            scanpy_version = None
        adata.uns["trajectory"] = {
            "engine": "scanpy",
            "method": "paga_dpt",
            "scanpy_version": scanpy_version,
            "embedding_key": embedding_key,
            "neighbor_rep": neighbor_rep,
            "neighbors_rebuilt_on": ctx.metadata.get("trajectory_neighbors_rebuilt_on"),
            "iroot": int(adata.uns.get("iroot", 0)),
            "root_cluster": (
                str(root_cluster) if root_cluster is not None else None
            ),
            "root_cluster_found": root_cluster_found,
            "root_justification": root_justification or None,
            "root_selection_method": root_selection_method,
            "claimable": trajectory_claimable,
            "claim_status": claim_status,
            "inference_class": (
                "biologically_rooted_dpt"
                if trajectory_claimable
                else "exploratory_unrooted_dpt"
            ),
            "n_obs": int(adata.n_obs),
        }
        ctx.metadata["trajectory_engine"] = "scanpy_paga_dpt"
        ctx.metadata["trajectory_method_actually_used"] = "scanpy_paga_dpt"
        ctx.metadata["trajectory_claimable"] = trajectory_claimable
        ctx.metadata["trajectory_claim_status"] = claim_status
        ctx.metadata["trajectory_root_selection_method"] = root_selection_method

        # --- Tables ---
        pseudotime_table = adata.obs[["dpt_pseudotime"]].copy()
        pseudotime_table["trajectory_claimable"] = trajectory_claimable
        pseudotime_table["trajectory_claim_status"] = claim_status
        pseudotime_table.to_csv(ctx.table_dir / "dpt_pseudotime.csv")

        # Per-cluster pseudotime statistics
        if "leiden" in adata.obs:
            stats = (
                adata.obs[["leiden", "dpt_pseudotime"]]
                .groupby("leiden", observed=True)["dpt_pseudotime"]
                .agg(["count", "mean", "median", "std", "min", "max"])
                .reset_index()
            )
            stats["trajectory_claimable"] = trajectory_claimable
            stats["trajectory_claim_status"] = claim_status
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
        claim_label = TrajectoryModule._claim_label(adata)
        ncols = 3 if "leiden" in adata.obs else 2
        fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5))
        sc.pl.umap(adata, color="dpt_pseudotime", ax=axes[0], show=False, cmap="viridis")
        axes[0].set_title(f"DPT Pseudotime{claim_label}")

        # Diffusion components
        if "X_diffmap" in adata.obsm:
            sc.pl.embedding(
                adata, basis="diffmap", color="dpt_pseudotime",
                components="1,2", ax=axes[1], show=False, cmap="viridis",
            )
            axes[1].set_title(f"Diffusion Map (DC1 vs DC2){claim_label}")

        if ncols == 3:
            sc.pl.umap(adata, color="leiden", ax=axes[2], show=False, legend_loc="on data")
            axes[2].set_title(f"Leiden Clusters{claim_label}")

        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudotime_dpt_umap.png", bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_paga(adata, ctx: PipelineContext) -> None:
        """PAGA trajectory graph overlaid on UMAP."""
        try:
            claim_label = TrajectoryModule._claim_label(adata)
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            sc.pl.paga(
                adata, color="leiden", ax=axes[0], show=False,
                frameon=False, fontsize=9, node_size_scale=1.5,
            )
            axes[0].set_title(f"PAGA Graph (node = cluster){claim_label}")
            sc.pl.umap(
                adata, color="dpt_pseudotime", ax=axes[1], show=False, cmap="viridis",
            )
            axes[1].set_title(f"DPT Pseudotime{claim_label}")
            plt.tight_layout()
            plt.savefig(ctx.figure_dir / "paga_trajectory.png", bbox_inches="tight")
            plt.close()
        except Exception:
            plt.close("all")

    @staticmethod
    def _hvg_mask_on_expression_axis(adata, expr_sub, ctx: PipelineContext):
        """Align ``adata.var["highly_variable"]`` onto the axis actually correlated.

        Real-data defect found 2026-08-02 on the 92,430-cell LUSC dataset. The gene
        trends are computed over ``adata.raw``, but the HVG flag was reindexed straight
        from ``adata.var`` — and those two axes are in different namespaces once the
        ingest boundary converts ``var_names`` Ensembl->symbol while leaving ``.raw``
        alone (see governance/raw_axis_namespace_divergence_2026-08-02.md).

        The reindex did not fail loudly, and it did not return empty either: the two
        axes overlap on the handful of features whose ``feature_name`` was itself an
        Ensembl ID (36 of 17,764 on the real file). So ``candidate_idx`` came back
        NON-empty with ~5 unnamed pseudogenes, the ``size == 0`` fallback never
        triggered, and "genes most correlated with pseudotime" was reported from that
        arbitrary handful instead of from the HVGs. A crash would have been kinder.

        Returns ``None`` when no trustworthy mask can be built, which makes the caller
        fall back to the full axis. That path is BOUNDED (candidates are variance-ranked
        down to ``max_candidates``) but it is a WEAKER selection criterion, not an
        equivalent one — on a raw counts matrix, variance ranking is dominated by the
        mean-variance relationship, so it selects the highest-EXPRESSED genes rather than
        the most variable. The refusal is therefore recorded, not silent.

        Why completeness rather than a sufficiency threshold
        ---------------------------------------------------
        The first version of this fix accepted the direct reindex whenever at least half
        the flagged HVGs survived it. That is backwards, and an independent review proved
        it: the accidental cross-namespace overlap consists of features whose
        ``feature_name`` was blank and therefore kept their Ensembl ID
        (``normalize_var_to_symbols`` does this deliberately). The more such features a
        file has, the LARGER that overlap — so above ~50% the sufficiency test accepted a
        mask containing ONLY the unnamed features, dropped every named HVG, and recorded
        nothing at all. Protection decreased as contamination increased, and the failure
        was quieter than the refusal it replaced. Measured on a constructed CELLxGENE-shaped
        file: at 52% unnamed it returned 52 of 100 HVGs, all of them the wrong ones, with
        an empty ``ctx.metadata``.

        So: compute BOTH alignments, take the more complete, and accept only a COMPLETE
        one. Anything short of every flagged HVG is a refusal with the shortfall recorded.
        Every reachable case in this pipeline aligns exactly — ``.raw`` is assigned as a
        copy or superset of ``var``, never a subset — so completeness is achievable rather
        than merely strict.
        """
        if "highly_variable" not in adata.var.columns:
            return None
        flags = adata.var["highly_variable"].astype(bool).to_numpy()
        if not flags.any():
            return None
        expr_names = np.asarray(list(expr_sub.var_names), dtype=object)

        def attempt(keys):
            """Join the flagged genes onto the expression axis by `keys`.

            Set-based on purpose. `Series.reindex` DUPLICATES a value across duplicated
            target labels, so one flagged gene facing a `.raw` axis that repeats its ID
            contributed several Trues and inflated the count to look complete while other
            flagged genes were never selected at all. Counting distinct KEYS on both
            sides removes that, and it also keeps the two counts in the same space --
            counting flagged rows in `var` (with duplicates) against selected labels in
            `raw` (deduplicated) made them incomparable, and produced permanent false
            refusals whenever two flagged genes shared an Ensembl ID.
            """
            keys = np.asarray(keys, dtype=object)
            flagged_keys = set(keys[flags].tolist())
            if not flagged_keys:
                return None, 0, 0
            mask = np.isin(expr_names, np.asarray(sorted(flagged_keys), dtype=object))
            present = flagged_keys & set(expr_names.tolist())
            return mask, len(present), len(flagged_keys)

        best = attempt(adata.var_names.astype(str))
        route = "direct"
        # Re-key through the Ensembl IDs, the namespace `.raw` is still in. Computed
        # unconditionally and compared on coverage, never gated behind the direct attempt
        # "looking good enough" -- that sufficiency test was the previous defect.
        if "ensembl_id" in adata.var.columns:
            # Only join through a key the ingest boundary vouched for. A shipped
            # `ensembl_id` that nothing verified is just a string column with the right
            # name -- and a wrong one is BIJECTIVE, so it satisfies the set-identity
            # check by construction and gets reported as a successful recovery. Measured
            # on the real file with a permuted shipped column: 2,000 HVGs accepted,
            # 290 correct, marked recovered_via_ensembl_id.
            from .._gene_symbols import TRUSTED_ENSEMBL_ID_SOURCES

            namespace = (getattr(ctx, "metadata", None) or {}).get("gene_namespace") or {}
            source = namespace.get("ensembl_id_source", "unknown")
            if source in TRUSTED_ENSEMBL_ID_SOURCES:
                alt = attempt(adata.var["ensembl_id"].astype(str))
                if alt[1] > best[1]:
                    best, route = alt, "recovered_via_ensembl_id"
            else:
                ctx.metadata["trajectory_hvg_ensembl_key_rejected"] = source
                logger.warning(
                    "TRAJECTORY: var['ensembl_id'] is present but its provenance is %r, "
                    "not one of %s. Refusing to re-key HVG flags through an unverified "
                    "join column -- a wrong one would be bijective and would pass every "
                    "completeness check while selecting the wrong genes.",
                    source, sorted(TRUSTED_ENSEMBL_ID_SOURCES),
                )

        mask, n_present, n_keys = best
        if mask is not None and n_present == n_keys:
            if route != "direct":
                ctx.metadata["trajectory_hvg_axis_alignment"] = route
                logger.info(
                    "TRAJECTORY: HVG flags re-keyed through var['ensembl_id'] to match "
                    "the expression axis (%d/%d distinct keys).", n_present, n_keys,
                )
            return mask

        ctx.metadata["trajectory_hvg_axis_alignment"] = "unavailable_axis_mismatch"
        ctx.metadata["trajectory_hvg_aligned_fraction"] = (
            round(n_present / n_keys, 4) if n_keys else 0.0)
        logger.warning(
            "TRAJECTORY: only %d of %d distinct highly-variable gene keys are present on "
            "the expression axis (best route: %s), so the HVG restriction is NOT applied "
            "and all %d genes are used as candidates, variance-ranked. On a counts matrix "
            "that favours the highest-expressed genes, not the most variable -- a weaker "
            "criterion, recorded rather than silently substituted. adata.var and the "
            "correlated matrix address genes in different namespaces.",
            n_present, n_keys, route, int(expr_sub.n_vars),
        )
        return None

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
        axis = resolve_expression_axis(adata, use_raw=True)
        expr = axis.expression
        expr_sub = expr[valid]
        pt_vals = pt.loc[valid].values.astype(np.float32)
        n_cells = float(len(pt_vals))
        if n_cells <= 1:
            return

        hv_mask = TrajectoryModule._hvg_mask_on_expression_axis(adata, expr_sub, ctx)
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

        # Save gene-pseudotime correlations.
        #
        # The two trajectory_claim_* columns describe the ROOT justification and say
        # nothing about which genes were searched. Without the candidate-set columns a
        # reader cannot tell "top 30 of the highly-variable genes" from "top 30 of the
        # 3,000 highest-variance columns of a counts matrix" -- and on raw counts,
        # variance ranking is dominated by the mean-variance relationship, so the second
        # is closer to "most highly expressed". Same claim-marking principle the rest of
        # the suite applies to engine substitution: the weaker selection is allowed,
        # unmarked is not.
        trajectory = adata.uns.get("trajectory", {})
        corr_df = pd.DataFrame({
            "gene": expr_sub.var_names[top_idx],
            "pseudotime_correlation": corr[top_local_idx],
            "trajectory_claimable": bool(trajectory.get("claimable", False)),
            "trajectory_claim_status": trajectory.get("claim_status", "non_claimable_unknown"),
            "candidate_selection": (
                "highly_variable_genes" if hv_mask is not None
                else "variance_ranked_full_axis"
            ),
            "n_candidates_searched": int(candidate_idx.size),
            "hvg_axis_alignment": ctx.metadata.get(
                "trajectory_hvg_axis_alignment", "direct"),
        })
        corr_df.to_csv(ctx.table_dir / "pseudotime_top_genes.csv", index=False)

        # Heatmap: cells ordered by pseudotime, top genes
        sort_idx = np.argsort(pt_vals)
        mat = expr_sub.X[:, top_idx]
        if sparse.issparse(mat):
            # densify-allowed: subset of top_k genes × pseudotime-sorted cells; bounded by n_top_genes (default ≤200)
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
        ax.set_title(
            f"Top 30 genes correlated with pseudotime{TrajectoryModule._claim_label(adata)}"
        )
        plt.colorbar(im, ax=ax, label="Smoothed expression", shrink=0.6)
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudotime_gene_heatmap.png", bbox_inches="tight")
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
        ax.set_title(
            f"Pseudotime Distribution per Cluster{TrajectoryModule._claim_label(adata)}"
        )
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudotime_violin_per_cluster.png", bbox_inches="tight")
        plt.close()

    @staticmethod
    def _claim_label(adata) -> str:
        """Return an explicit plot label when DPT cannot support biological claims."""
        trajectory = adata.uns.get("trajectory", {})
        if trajectory.get("claimable") is True:
            return ""
        return "\nNON-CLAIMABLE: missing biologically justified root"
