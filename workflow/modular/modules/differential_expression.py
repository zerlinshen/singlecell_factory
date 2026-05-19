from __future__ import annotations

import logging
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from ._scanpy_compat import has_api, import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext
from ._gpu_utils import gpu_available
from .._densify_policy import plan_densify, DensifyDecision


__references__ = {
    "scanpy": {
        "title": "SCANPY: large-scale single-cell gene expression data analysis",
        "authors": "Wolf, Angerer, Theis",
        "journal": "Genome Biology",
        "year": "2018",
        "doi": "10.1186/s13059-017-1382-0",
        "description": "scanpy.tl.rank_genes_groups implementation (Wilcoxon/t-test/MAST/ROC).",
    },
    "Soneson_DE_benchmark_2018": {
        "title": "Bias, robustness and scalability in single-cell differential expression analysis",
        "authors": "Soneson, Robinson",
        "journal": "Nature Methods",
        "year": "2018",
        "doi": "10.1038/nmeth.4612",
        "description": "Benchmark showing Wilcoxon competitive with bespoke single-cell DE methods.",
    },
    "Benjamini_Hochberg_1995": {
        "title": "Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing",
        "authors": "Benjamini, Hochberg",
        "journal": "Journal of the Royal Statistical Society B",
        "year": "1995",
        "doi": "10.1111/j.2517-6161.1995.tb02031.x",
        "description": "BH FDR correction applied to per-gene p-values and per-substate marker scoring.",
    },
}


logger = logging.getLogger(__name__)


class DifferentialExpressionModule:
    """Optional module: identify cluster marker genes with significance filtering and visualizations."""

    name = "differential_expression"
    requires_keys = {"obs": ["leiden"]}
    provides_keys = {"uns": ["rank_genes_groups"]}

    def run(self, ctx: PipelineContext) -> None:
        import os
        if os.environ.get("SC_MEM_GUARD", "").lower() == "on":
            from .._mem_guard import MemoryGuard, MemoryGuardError
            try:
                with MemoryGuard(ctx, self.name) as mg:
                    mg.check("entry")
                    return self._run_impl(ctx)
            except MemoryGuardError as exc:
                logger.warning("MemoryGuard: %s", exc)
                ctx.metadata.setdefault("mem_warnings", []).append(
                    {"module": self.name, "error": str(exc)}
                )
        return self._run_impl(ctx)

    def _run_impl(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "leiden" not in adata.obs:
            raise ValueError("Differential expression requires clustered AnnData.")

        method = ctx.cfg.de_method
        n_genes = max(1, int(ctx.cfg.de_n_genes))
        rank_kwargs = {"use_raw": False, "n_genes": n_genes, "pts": True}
        if method == "wilcoxon":
            rank_kwargs["tie_correct"] = True

        # HIGH-3 fix (2026-05-19): make the actual test used + correction method
        # explicit on both backends and record them in metadata so the contract
        # validator can detect a silent fallback to a different statistic.
        configured_corr_method = getattr(ctx.cfg, "de_correction", None) or "benjamini-hochberg"
        de_test_used = method  # the canonical Wilcoxon/t-test choice
        use_gpu = gpu_available(ctx.cfg.gpu_mode) and method in ("wilcoxon", "t-test", "t-test_overestim_var")
        if use_gpu:
            try:
                import rapids_singlecell as rsc
                logger.info("GPU DE: using rapids-singlecell rank_genes_groups (%s)", method)
                rsc.tl.rank_genes_groups(
                    adata,
                    groupby="leiden",
                    method=method,
                    corr_method=configured_corr_method,
                    **rank_kwargs,
                )
                ctx.metadata["de_backend"] = "gpu"
            except Exception as exc:
                if ctx.cfg.gpu_mode == "force":
                    raise RuntimeError(f"GPU DE failed in force mode: {exc}") from exc
                logger.warning("GPU DE failed (%s), falling back to CPU", exc)
                ctx.metadata["de_gpu_fallback_reason"] = str(exc)
                use_gpu = False

        if not use_gpu:
            markers = self._run_cpu_rank_genes_groups(
                adata=adata,
                method=method,
                rank_kwargs=rank_kwargs,
                n_genes=n_genes,
                corr_method=configured_corr_method,
            )
            ctx.metadata["de_backend"] = "cpu"
        elif has_api(sc, "get.rank_genes_groups_df"):
            markers = sc.get.rank_genes_groups_df(adata, group=None)
        else:
            # GPU wrote ranking state but scanpy accessor is unavailable in this environment.
            # The fallback uses Welch t-test + manual BH on its own — record the swap.
            logger.warning("scanpy get.rank_genes_groups_df unavailable; using Welch t-test fallback (HIGH-3)")
            de_test_used = "welch_t_test_fallback"
            markers = self._fallback_rank_genes_groups_df(
                adata=adata,
                groupby="leiden",
                n_genes=n_genes,
            )

        ctx.metadata["de_test_actually_used"] = de_test_used
        ctx.metadata["de_correction"] = configured_corr_method
        if "pvals_adj" not in markers.columns:
            markers["pvals_adj"] = 1.0
        if "logfoldchanges" not in markers.columns:
            markers["logfoldchanges"] = np.nan

        # Save unfiltered results
        markers.to_csv(ctx.table_dir / "marker_genes_all.csv", index=False)

        # Apply significance filtering
        sig_mask = markers["pvals_adj"] < ctx.cfg.de_pval_threshold
        if markers["logfoldchanges"].notna().any():
            sig_mask &= markers["logfoldchanges"].abs() > ctx.cfg.de_logfc_threshold
        sig_markers = markers[sig_mask].copy()
        sig_markers.to_csv(ctx.table_dir / "marker_genes.csv", index=False)

        top5 = (
            sig_markers.sort_values(["group", "pvals_adj", "scores"], ascending=[True, True, False])
            .groupby("group", observed=True)
            .head(5)
            .reset_index(drop=True)
        )
        top5.to_csv(ctx.table_dir / "marker_top5_by_cluster.csv", index=False)

        ctx.metadata["de_significant_genes"] = int(len(sig_markers))
        ctx.metadata["de_method"] = method
        ctx.metadata["gpu_mode"] = ctx.cfg.gpu_mode
        ctx.metadata["de_n_genes"] = int(n_genes)
        ctx.metadata["de_random_state"] = ctx.random_state

        # --- P1A.S3b: per-substate DE CSVs (AC-2 support) ---
        # Gated on column presence; AC-10 zero-regression when column absent.
        if "context_aware_substate" in adata.obs:
            self._write_substate_de(adata, markers, ctx)

        # --- Visualizations ---
        # Dot plot of top markers per cluster
        try:
            sc.pl.rank_genes_groups_dotplot(
                adata, n_genes=5, show=False, standard_scale="var",
            )
            plt.savefig(
                ctx.figure_dir / "de_dotplot_top5.png", dpi=160, bbox_inches="tight",
            )
            plt.close()
        except Exception:
            plt.close("all")

        # Heatmap of top markers
        try:
            sc.pl.rank_genes_groups_heatmap(
                adata, n_genes=5, show=False, show_gene_labels=True,
                use_raw=False, swap_axes=True, vmin=-3, vmax=3,
            )
            plt.savefig(
                ctx.figure_dir / "de_heatmap_top5.png", dpi=160, bbox_inches="tight",
            )
            plt.close()
        except Exception:
            plt.close("all")

        # Volcano plot (all clusters combined)
        self._volcano_plot(sig_markers, ctx)

    @staticmethod
    def _run_cpu_rank_genes_groups(
        adata,
        method: str,
        rank_kwargs: dict,
        n_genes: int,
        corr_method: str | None = None,
    ) -> pd.DataFrame:
        if has_api(sc, "tl.rank_genes_groups") and has_api(sc, "get.rank_genes_groups_df"):
            call_kwargs = dict(rank_kwargs)
            if corr_method is not None:
                import inspect as _inspect
                try:
                    sig = _inspect.signature(sc.tl.rank_genes_groups)
                    if "corr_method" in sig.parameters:
                        call_kwargs["corr_method"] = corr_method
                except Exception:
                    pass
            sc.tl.rank_genes_groups(
                adata,
                groupby="leiden",
                method=method,
                **call_kwargs,
            )
            return sc.get.rank_genes_groups_df(adata, group=None)
        logger.warning("scanpy DE APIs unavailable; using statistical fallback.")
        return DifferentialExpressionModule._fallback_rank_genes_groups_df(
            adata=adata,
            groupby="leiden",
            n_genes=n_genes,
        )

    @staticmethod
    def _fallback_rank_genes_groups_df(adata, groupby: str, n_genes: int) -> pd.DataFrame:
        if groupby not in adata.obs:
            adata.uns["rank_genes_groups"] = {"fallback": True, "groupby": groupby}
            return pd.DataFrame(columns=["group", "names", "scores", "pvals_adj", "logfoldchanges"])

        import os
        engine = os.environ.get("SC_DE_ENGINE", "dense").lower()
        X_raw = adata.X
        if engine == "sparse" and sparse.issparse(X_raw):
            return DifferentialExpressionModule._fallback_sparse_welch_df(
                adata=adata, groupby=groupby, n_genes=n_genes,
            )

        X = X_raw
        if sparse.issparse(X):
            decision = plan_densify(
                X.shape, float,
                reason="legacy dense Welch DE path; SC_DE_ENGINE != sparse",
            )
            if decision == DensifyDecision.ABORT:
                from .._mem_guard import MemoryGuardError
                raise MemoryGuardError(
                    f"differential_expression: matrix too large to densify "
                    f"({X.shape[0]} × {X.shape[1]}); set SC_DE_ENGINE=sparse to use sparse path"
                )
            if decision == DensifyDecision.CHUNK:
                logger.warning(
                    "differential_expression: DE densify in CHUNK range (%d × %d); proceeding",
                    X.shape[0], X.shape[1],
                )
            X = X.toarray()  # densify-allowed: guarded by plan_densify above
        X = np.asarray(X, dtype=float)
        if X.ndim != 2 or X.shape[1] == 0:
            adata.uns["rank_genes_groups"] = {"fallback": True, "groupby": groupby}
            return pd.DataFrame(columns=["group", "names", "scores", "pvals_adj", "logfoldchanges"])

        genes = np.asarray(adata.var_names.astype(str))
        groups = pd.Categorical(adata.obs[groupby].astype(str))
        rows: list[dict[str, object]] = []
        per_group_names: dict[str, list[str]] = {}

        for group in groups.categories:
            in_mask = np.asarray(groups == group)
            out_mask = ~in_mask
            if in_mask.sum() == 0 or out_mask.sum() == 0:
                continue
            in_x = X[in_mask]
            out_x = X[out_mask]
            mean_in = in_x.mean(axis=0) + 1e-9
            mean_out = out_x.mean(axis=0) + 1e-9
            logfc = np.log2(mean_in / mean_out)
            try:
                _, pvals = stats.ttest_ind(
                    in_x,
                    out_x,
                    axis=0,
                    equal_var=False,
                    nan_policy="omit",
                )
            except Exception:
                pvals = np.ones(X.shape[1], dtype=float)
            pvals = np.nan_to_num(np.asarray(pvals, dtype=float), nan=1.0, posinf=1.0, neginf=1.0)
            pvals_adj = DifferentialExpressionModule._benjamini_hochberg(pvals)
            order = np.lexsort((-np.abs(logfc), pvals_adj))
            keep = order[: max(1, min(n_genes, len(order)))]
            per_group_names[str(group)] = genes[keep].tolist()
            for idx in keep:
                rows.append(
                    {
                        "group": str(group),
                        "names": genes[idx],
                        "scores": float(logfc[idx]),
                        "pvals_adj": float(pvals_adj[idx]),
                        "logfoldchanges": float(logfc[idx]),
                    }
                )

        adata.uns["rank_genes_groups"] = {
            "fallback": True,
            "groupby": groupby,
            "names": per_group_names,
        }
        if not rows:
            return pd.DataFrame(columns=["group", "names", "scores", "pvals_adj", "logfoldchanges"])
        return pd.DataFrame(rows)

    @staticmethod
    def _fallback_sparse_welch_df(adata, groupby: str, n_genes: int) -> pd.DataFrame:
        from .._sparse_utils import sparse_welch_t

        X = adata.X
        if X.ndim != 2 or X.shape[1] == 0:
            adata.uns["rank_genes_groups"] = {"fallback": True, "groupby": groupby}
            return pd.DataFrame(columns=["group", "names", "scores", "pvals_adj", "logfoldchanges"])

        genes = np.asarray(adata.var_names.astype(str))
        groups = pd.Categorical(adata.obs[groupby].astype(str))
        rows: list[dict[str, object]] = []
        per_group_names: dict[str, list[str]] = {}

        for group in groups.categories:
            in_mask = np.asarray(groups == group)
            out_mask = ~in_mask
            if in_mask.sum() == 0 or out_mask.sum() == 0:
                continue
            try:
                _, pvals_raw, mean_in, mean_out = sparse_welch_t(X, in_mask, out_mask)
            except Exception as exc:
                logger.warning("sparse_welch_t failed for group %s: %s — skipping", group, exc)
                continue
            mean_in = np.asarray(mean_in, dtype=float) + 1e-9
            mean_out = np.asarray(mean_out, dtype=float) + 1e-9
            logfc = np.log2(mean_in / mean_out)
            pvals = np.nan_to_num(np.asarray(pvals_raw, dtype=float), nan=1.0, posinf=1.0, neginf=1.0)
            pvals_adj = DifferentialExpressionModule._benjamini_hochberg(pvals)
            order = np.lexsort((-np.abs(logfc), pvals_adj))
            keep = order[: max(1, min(n_genes, len(order)))]
            per_group_names[str(group)] = genes[keep].tolist()
            for idx in keep:
                rows.append(
                    {
                        "group": str(group),
                        "names": genes[idx],
                        "scores": float(logfc[idx]),
                        "pvals_adj": float(pvals_adj[idx]),
                        "logfoldchanges": float(logfc[idx]),
                    }
                )

        adata.uns["rank_genes_groups"] = {
            "fallback": True,
            "groupby": groupby,
            "engine": "sparse_welch",
            "names": per_group_names,
        }
        if not rows:
            return pd.DataFrame(columns=["group", "names", "scores", "pvals_adj", "logfoldchanges"])
        return pd.DataFrame(rows)

    @staticmethod
    def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
        p = np.clip(np.asarray(pvals, dtype=float), 0.0, 1.0)
        n = p.size
        if n == 0:
            return p
        order = np.argsort(p)
        ranked = p[order]
        adjusted = np.empty_like(ranked)
        prev = 1.0
        for i in range(n - 1, -1, -1):
            rank = i + 1
            val = ranked[i] * n / rank
            prev = min(prev, val)
            adjusted[i] = prev
        out = np.empty_like(adjusted)
        out[order] = np.clip(adjusted, 0.0, 1.0)
        return out

    @staticmethod
    def _volcano_plot(markers: pd.DataFrame, ctx: PipelineContext) -> None:
        """Generate a combined volcano plot for all clusters."""
        if markers.empty:
            return
        fig, ax = plt.subplots(figsize=(10, 7))
        log_pval = -np.log10(markers["pvals_adj"].clip(lower=1e-20))
        lfc = markers["logfoldchanges"]

        scatter = ax.scatter(
            lfc, log_pval, c=lfc, cmap="RdBu_r", s=4, alpha=0.6,
            vmin=-3, vmax=3, edgecolors="none",
        )
        ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.8)
        ax.axvline(-0.25, color="grey", linestyle="--", linewidth=0.8)
        ax.axvline(0.25, color="grey", linestyle="--", linewidth=0.8)
        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-log10(adjusted p-value)")
        ax.set_title("Volcano Plot (all clusters)")
        plt.colorbar(scatter, ax=ax, label="Log2FC")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "de_volcano.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _write_substate_de(adata, markers: pd.DataFrame, ctx: PipelineContext) -> None:
        """Write per-substate DE CSVs when context_aware_substate column is present (P1A.S3b).

        Output: runs/<run-id>/differential_expression/substates/<substate_name>.csv
        Columns: gene, log2fc, pvalue, fdr, pct_in_substate, pct_in_parent

        Gated strictly on column presence. When column is absent, this method is never
        called — preserving AC-10 zero-regression for default runs.
        """
        substate_col = adata.obs.get("context_aware_substate")
        if substate_col is None:
            return

        # Only process rows where substate is non-null
        valid_mask = substate_col.notna() & (substate_col.astype(str) != "None")
        if not valid_mask.any():
            return

        from scipy import sparse, stats as _stats

        substates = substate_col[valid_mask].unique()
        substates_dir = ctx.table_dir.parent / "differential_expression" / "substates"
        substates_dir.mkdir(parents=True, exist_ok=True)

        written: list[str] = []
        for substate in substates:
            substate_str = str(substate)
            # Find the parent cell type: everything before the last space-separated qualifier
            # or the full string if no qualifier. Use leiden-mapped cell type if available.
            substate_mask = (
                substate_col.notna() &
                (substate_col.astype(str) == substate_str)
            )
            # Parent = all cells that share the same leiden cluster (context_aware_celltype parent)
            # Fallback: all non-substate cells as parent reference
            parent_mask = ~substate_mask

            n_in = int(substate_mask.sum())
            n_parent = int(parent_mask.sum())
            if n_in < 3 or n_parent < 3:
                continue

            import scipy.sparse as sp
            X = adata.X
            if sp.issparse(X):
                sub_shape = (n_in, X.shape[1])
                par_shape = (n_parent, X.shape[1])
                sub_decision = plan_densify(
                    sub_shape,
                    np.float64,
                    reason="context-aware substate DE input group for Mann-Whitney test",
                )
                par_decision = plan_densify(
                    par_shape,
                    np.float64,
                    reason="context-aware substate DE parent group for Mann-Whitney test",
                )
                if DensifyDecision.ABORT in (sub_decision, par_decision):
                    from .._mem_guard import MemoryGuardError
                    raise MemoryGuardError(
                        "differential_expression substate DE: sparse subset too large "
                        f"to densify safely (substate={sub_shape}, parent={par_shape})"
                    )
                if DensifyDecision.CHUNK in (sub_decision, par_decision):
                    logger.warning(
                        "differential_expression substate DE densify in CHUNK range "
                        "(substate=%s, parent=%s); proceeding",
                        sub_shape,
                        par_shape,
                    )
                # densify-allowed: guarded substate subset for scipy Mann-Whitney fallback
                X_sub = np.asarray(X[substate_mask].todense(), dtype=np.float64)
                # densify-allowed: guarded parent subset for scipy Mann-Whitney fallback
                X_par = np.asarray(X[parent_mask].todense(), dtype=np.float64)
            else:
                X_sub = np.asarray(X[substate_mask], dtype=np.float64)
                X_par = np.asarray(X[parent_mask], dtype=np.float64)

            mean_sub = X_sub.mean(axis=0) + 1e-9
            mean_par = X_par.mean(axis=0) + 1e-9
            log2fc = np.log2(mean_sub / mean_par)
            pct_in = (X_sub > 0).mean(axis=0)
            pct_par = (X_par > 0).mean(axis=0)

            try:
                _, pvals = _stats.mannwhitneyu(X_sub, X_par, axis=0, alternative="two-sided")
                pvals = np.nan_to_num(np.asarray(pvals, dtype=float), nan=1.0)
            except Exception:
                pvals = np.ones(X_sub.shape[1], dtype=float)

            fdr = DifferentialExpressionModule._benjamini_hochberg(pvals)

            genes = np.asarray(adata.var_names.astype(str))
            result = pd.DataFrame({
                "gene": genes,
                "log2fc": log2fc.ravel(),
                "pvalue": pvals.ravel(),
                "fdr": fdr.ravel(),
                "pct_in_substate": pct_in.ravel(),
                "pct_in_parent": pct_par.ravel(),
            })
            result = result.sort_values("fdr").reset_index(drop=True)

            safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in substate_str)
            out_path = substates_dir / f"{safe_name}.csv"
            result.to_csv(out_path, index=False)
            written.append(substate_str)
            logger.info(
                "DE substates: wrote %d genes for substate '%s' → %s",
                len(result), substate_str, out_path,
            )

        ctx.metadata["de_substates_written"] = written
        ctx.metadata["de_substates_count"] = len(written)
