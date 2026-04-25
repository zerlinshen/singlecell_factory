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

        use_gpu = gpu_available(ctx.cfg.gpu_mode) and method in ("wilcoxon", "t-test", "t-test_overestim_var")
        if use_gpu:
            try:
                import rapids_singlecell as rsc
                logger.info("GPU DE: using rapids-singlecell rank_genes_groups (%s)", method)
                rsc.tl.rank_genes_groups(
                    adata,
                    groupby="leiden",
                    method=method,
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
            )
            ctx.metadata["de_backend"] = "cpu"
        elif has_api(sc, "get.rank_genes_groups_df"):
            markers = sc.get.rank_genes_groups_df(adata, group=None)
        else:
            # GPU wrote ranking state but scanpy accessor is unavailable in this environment.
            markers = self._fallback_rank_genes_groups_df(
                adata=adata,
                groupby="leiden",
                n_genes=n_genes,
            )
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
    ) -> pd.DataFrame:
        if has_api(sc, "tl.rank_genes_groups") and has_api(sc, "get.rank_genes_groups_df"):
            sc.tl.rank_genes_groups(
                adata,
                groupby="leiden",
                method=method,
                **rank_kwargs,
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
            X = X.toarray()
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
