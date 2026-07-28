from __future__ import annotations

from importlib import metadata as importlib_metadata
import logging
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from ._scanpy_compat import has_api, import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext
from ._gpu_utils import bind_cuda_context, gpu_available
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
        if os.environ.get("SC_MEM_GUARD", "").lower() == "on":
            from .._mem_guard import MemoryGuard, MemoryGuardError
            entered_impl = False
            try:
                with MemoryGuard(ctx, self.name) as mg:
                    mg.check("entry")
                    entered_impl = True
                    return self._run_impl(ctx)
            except MemoryGuardError as exc:
                logger.warning("MemoryGuard: %s", exc)
                ctx.metadata.setdefault("mem_warnings", []).append(
                    {"module": self.name, "error": str(exc)}
                )
                if entered_impl:
                    # MemoryGuardError was raised from inside _run_impl
                    # (substage abort): propagate as-is, no retry.
                    raise
                # else: raised by mg.check() (entry guard signal) ->
                # fall through to degraded _run_impl call below.
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
                bind_cuda_context()  # ensure cuBLAS/cuSOLVER pre-warm before GPU ops (P1 fix; no-op if already warmed in probe)
                rsc_version = self._rapids_singlecell_version()
                gpu_rank_genes = getattr(getattr(rsc, "tl", None), "rank_genes_groups", None)
                ctx.metadata["de_rapids_singlecell_version"] = rsc_version
                if gpu_rank_genes is None:
                    reason = (
                        f"rapids-singlecell {rsc_version} does not expose "
                        f"tl.rank_genes_groups for method {method}"
                    )
                    if ctx.cfg.gpu_mode == "force":
                        raise RuntimeError(reason)
                    logger.warning("GPU DE unavailable: %s; falling back to CPU", reason)
                    ctx.metadata["de_gpu_fallback_reason"] = reason
                    use_gpu = False
                else:
                    logger.info("GPU DE: using rapids-singlecell rank_genes_groups (%s)", method)
                    gpu_rank_genes(
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
            if self._should_use_sparse_cpu_de(ctx, adata):
                logger.info("CPU DE: using sparse Welch fallback for massive/sparse input")
                markers = self._fallback_sparse_welch_df(
                    adata=adata,
                    groupby="leiden",
                    n_genes=n_genes,
                    corr_method=configured_corr_method,
                )
                de_test_used = "sparse_welch_fallback"
                ctx.metadata["de_backend"] = "cpu_sparse"
                ctx.metadata["de_cpu_fallback_reason"] = "massive_or_requested_sparse_cpu_de"
            else:
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
            # F-3 (Plan ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md,
            # Principle 6): GPU wrote ranking state but scanpy accessor is
            # unavailable. The legacy fallback silently routed to Welch
            # t-test + manual BH — a hard scientific compromise. The route is
            # now banned unless SC_ALLOW_WELCH_FALLBACK=1 is set explicitly.
            if os.environ.get("SC_ALLOW_WELCH_FALLBACK", "").strip() != "1":
                raise RuntimeError(
                    "scanpy get.rank_genes_groups_df is unavailable in this "
                    "environment and the silent Welch fallback is banned for "
                    "final-claim runs (Plan F-3, Principle 6). Install / "
                    "upgrade scanpy to expose get.rank_genes_groups_df, or set "
                    "SC_ALLOW_WELCH_FALLBACK=1 to opt in to Welch explicitly "
                    "(loud warning + metadata flag). "
                    "See ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md"
                )
            logger.error(
                "Welch DE fallback OPT-IN via SC_ALLOW_WELCH_FALLBACK=1: "
                "scanpy get.rank_genes_groups_df is unavailable; running "
                "Welch t-test + manual BH. This is forbidden in final-claim "
                "runs (Plan F-3, Principle 6)."
            )
            de_test_used = "welch_t_test_fallback_OPT_IN"
            ctx.metadata["de_welch_opt_in_acknowledged"] = True
            markers = self._fallback_rank_genes_groups_df(
                adata=adata,
                groupby="leiden",
                n_genes=n_genes,
                corr_method=configured_corr_method,
            )

        ctx.metadata["de_test_actually_used"] = de_test_used
        ctx.metadata["de_correction"] = configured_corr_method
        ctx.metadata["de_correction_actually_used"] = self._normalise_correction_method(
            configured_corr_method
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
        # Minimum within-group detection filter (Seurat min.pct, Stuart 2019):
        # suppress genes detected in a negligible fraction of the group's cells
        # that can otherwise show extreme apparent LFC. pts=True populates
        # pct_nz_group (fraction 0-1); guard on column presence so fallback DE
        # paths (which lack it) are unaffected.
        de_min_pct = float(getattr(ctx.cfg, "de_min_pct", 0.10))
        if "pct_nz_group" in markers.columns and de_min_pct > 0:
            sig_mask &= markers["pct_nz_group"] >= de_min_pct
            ctx.metadata["de_min_pct_applied"] = de_min_pct
        else:
            ctx.metadata["de_min_pct_applied"] = None
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
        self._plot_marker_dotplot(adata, ctx)

        # Heatmap of top markers
        self._plot_marker_heatmap(adata, ctx)

        # Volcano plot (all clusters combined). The full tested set plus the
        # significance mask computed above are handed over so the panel can
        # contrast passing against non-passing tests and state both counts;
        # nothing is re-tested or re-thresholded here.
        self._volcano_plot(markers, ctx, sig_mask=sig_mask)

    @staticmethod
    def _rapids_singlecell_version() -> str:
        for dist_name in ("rapids-singlecell", "rapids_singlecell"):
            try:
                return importlib_metadata.version(dist_name)
            except importlib_metadata.PackageNotFoundError:
                continue
        return "unknown"

    @staticmethod
    def _should_use_sparse_cpu_de(ctx: PipelineContext, adata) -> bool:
        # F-3 (Plan ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md,
        # Principle 6): the sparse Welch t-test CPU fallback is a silent
        # algorithmic substitution (parametric Welch t-test for the
        # non-parametric Wilcoxon rank-sum gold standard per Soneson &
        # Robinson 2018). Both the prior auto-route on scale_mode=massive
        # AND the SC_DE_ENGINE=sparse direct request are now gated behind an
        # explicit opt-in env var SC_ALLOW_WELCH_FALLBACK=1. Default: return
        # False so the standard Wilcoxon CPU path runs.
        if not sparse.issparse(adata.X):
            return False
        engine = os.environ.get("SC_DE_ENGINE", "").strip().lower()
        if engine in {"scanpy", "dense", "cpu-scanpy"}:
            return False
        if engine == "sparse":
            if os.environ.get("SC_ALLOW_WELCH_FALLBACK", "").strip() != "1":
                raise RuntimeError(
                    "SC_DE_ENGINE=sparse requests the sparse Welch t-test "
                    "fallback, which is banned for final-claim runs "
                    "(Plan F-3, Principle 6). Set SC_ALLOW_WELCH_FALLBACK=1 "
                    "to opt in explicitly (loud warning + metadata flag), "
                    "or unset SC_DE_ENGINE to use the standard Wilcoxon CPU "
                    "path. See ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md"
                )
            logger.error(
                "Welch DE fallback OPT-IN via SC_ALLOW_WELCH_FALLBACK=1; "
                "this is forbidden in final-claim runs (Plan F-3, Principle 6). "
                "de_test_actually_used will record 'sparse_welch_fallback_OPT_IN'."
            )
            ctx.metadata["de_welch_opt_in_acknowledged"] = True
            return True
        # F-3: prior behavior auto-routed scale_mode=massive to Welch silently;
        # that route is removed. The standard Wilcoxon CPU path now runs even
        # on massive / sparse inputs.
        return False

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
            corr_method=corr_method,
        )

    @staticmethod
    def _fallback_rank_genes_groups_df(
        adata,
        groupby: str,
        n_genes: int,
        corr_method: str | None = None,
    ) -> pd.DataFrame:
        if groupby not in adata.obs:
            adata.uns["rank_genes_groups"] = {"fallback": True, "groupby": groupby}
            return pd.DataFrame(columns=["group", "names", "scores", "pvals_adj", "logfoldchanges"])

        import os
        engine = os.environ.get("SC_DE_ENGINE", "dense").lower()
        X_raw = adata.X
        if engine == "sparse" and sparse.issparse(X_raw):
            return DifferentialExpressionModule._fallback_sparse_welch_df(
                adata=adata,
                groupby=groupby,
                n_genes=n_genes,
                corr_method=corr_method,
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
            # Fallback-only path (scanpy unavailable / SC_ALLOW_WELCH_FALLBACK).
            # X is log1p-normalized at DE time, so the log fold change must be
            # taken on the LINEAR scale (expm1) to match scanpy's convention;
            # log2(mean_log_in / mean_log_out) is not a fold change. eps guards
            # against log2(0) for genes absent in one group.
            eps = 1e-9
            mean_in = in_x.mean(axis=0)
            mean_out = out_x.mean(axis=0)
            logfc = np.log2((np.expm1(mean_in) + eps) / (np.expm1(mean_out) + eps))
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
            pvals_adj = DifferentialExpressionModule._adjust_pvals(pvals, corr_method)
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
            "corr_method": DifferentialExpressionModule._normalise_correction_method(corr_method),
            "names": per_group_names,
        }
        if not rows:
            return pd.DataFrame(columns=["group", "names", "scores", "pvals_adj", "logfoldchanges"])
        return pd.DataFrame(rows)

    @staticmethod
    def _fallback_sparse_welch_df(
        adata,
        groupby: str,
        n_genes: int,
        corr_method: str | None = None,
    ) -> pd.DataFrame:
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
            # Fallback-only path (scanpy unavailable / SC_ALLOW_WELCH_FALLBACK).
            # X is log1p-normalized at DE time, so the log fold change must be
            # taken on the LINEAR scale (expm1) to match scanpy's convention;
            # log2 of a ratio of means of log values is not a fold change.
            eps = 1e-9
            mean_in = np.asarray(mean_in, dtype=float)
            mean_out = np.asarray(mean_out, dtype=float)
            logfc = np.log2((np.expm1(mean_in) + eps) / (np.expm1(mean_out) + eps))
            pvals = np.nan_to_num(np.asarray(pvals_raw, dtype=float), nan=1.0, posinf=1.0, neginf=1.0)
            pvals_adj = DifferentialExpressionModule._adjust_pvals(pvals, corr_method)
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
            "corr_method": DifferentialExpressionModule._normalise_correction_method(corr_method),
            "names": per_group_names,
        }
        if not rows:
            return pd.DataFrame(columns=["group", "names", "scores", "pvals_adj", "logfoldchanges"])
        return pd.DataFrame(rows)

    @staticmethod
    def _normalise_correction_method(corr_method: str | None) -> str:
        token = str(corr_method or "benjamini-hochberg").strip().lower()
        if token in {"bonferroni", "bonf"}:
            return "bonferroni"
        return "benjamini-hochberg"

    @staticmethod
    def _adjust_pvals(pvals: np.ndarray, corr_method: str | None) -> np.ndarray:
        method = DifferentialExpressionModule._normalise_correction_method(corr_method)
        p = np.clip(np.asarray(pvals, dtype=float), 0.0, 1.0)
        if method == "bonferroni":
            return np.clip(p * p.size, 0.0, 1.0)
        return DifferentialExpressionModule._benjamini_hochberg(p)

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

    # --- Figure tokens ---
    # The suite theme (plotting_factory/theme/nature_high_impact.yaml) builds
    # its continuous ramps out of ColorBrewer scales, so these matplotlib names
    # ARE the theme tokens rather than a second, drifting copy of the hexes:
    #   sequential_expression -> Blues  (#F7FBFF .. #6BAED6 .. #08306B)
    #   diverging_blue_red    -> RdBu_r (#2166AC .. #F7F7F7 .. #B2182B)
    # Both are colourblind-safe and neither is a rainbow ramp.
    _SEQUENTIAL_CMAP = "Blues"
    _DIVERGING_CMAP = "RdBu_r"
    # significance_highlight.{significant,nonsignificant} from the same theme
    # file; the significant hue is the one qc.py already uses for "this is the
    # cutoff decision", so the two panels speak with one vocabulary.
    _COLOUR_SIG = "#D55E00"
    _COLOUR_NS = "#9E9E9E"
    _COLOUR_GUIDE = "#666666"
    # Genes shown per cluster by the two top-5 panels (their filenames say so).
    _TOP_N_PER_CLUSTER = 5
    # Vertical pitch that keeps a 5 pt gene label (~1.8 mm) clear of its
    # neighbour, and the tallest figure a journal page can carry. Beyond the cap
    # the panel would be rescaled to fit anyway, shrinking the type with it.
    _GENE_ROW_PITCH_MM = 2.05
    _MAX_FIGURE_HEIGHT_MM = 247.0
    # The two-line title is placed above the canvas (see below), so the drawn
    # canvas has to give that band back or the saved image exceeds the cap.
    # Measured at 300 dpi, a 10 mm reserve still landed both panels at ~249 mm —
    # over the cap — because bbox_inches="tight" also captures the suptitle's
    # own padding. The typesetter then rescales by ~0.99, which drags the 5 pt
    # gene labels to ~4.95 pt, i.e. back under the legibility floor this sizing
    # logic exists to protect. 15 mm leaves both panels inside 247 mm.
    _TITLE_BAND_MM = 15.0

    _TEST_LABELS = {
        "wilcoxon": "Wilcoxon rank-sum",
        "t-test": "Welch t-test",
        "t-test_overestim_var": "t-test, overestimated variance",
        "logreg": "logistic regression",
        "sparse_welch_fallback": "Welch t-test (sparse fallback)",
        "sparse_welch_fallback_OPT_IN": "Welch t-test (sparse fallback)",
        "welch_t_test_fallback_OPT_IN": "Welch t-test (fallback)",
    }

    @classmethod
    def _test_label(cls, ctx: PipelineContext) -> str:
        """Human-readable name of the test that actually produced these numbers.

        Read from metadata rather than from ``cfg.de_method`` so a panel drawn
        after a fallback names the statistic that ran, not the one requested.
        """
        metadata = getattr(ctx, "metadata", {}) or {}
        token = str(
            metadata.get("de_test_actually_used")
            or getattr(ctx.cfg, "de_method", "")
            or "unknown test"
        )
        return cls._TEST_LABELS.get(token, token.replace("_", " "))

    @staticmethod
    def _correction_label(ctx: PipelineContext) -> str:
        metadata = getattr(ctx, "metadata", {}) or {}
        token = str(metadata.get("de_correction_actually_used") or "benjamini-hochberg")
        return "Bonferroni" if token == "bonferroni" else "Benjamini-Hochberg"

    @classmethod
    def _marker_row_count(cls, adata, n_genes: int) -> int:
        """Gene rows the top-N panels will draw (groups x N).

        The row count, not a fixed guess, is what decides whether the gene
        labels have room to sit next to each other, so it drives panel height.
        """
        try:
            names = adata.uns.get("rank_genes_groups", {}).get("names")
        except Exception:
            names = None
        groups = getattr(getattr(names, "dtype", None), "names", None)
        if not groups and isinstance(names, dict):
            groups = tuple(names)
        if not groups:
            groups = tuple(pd.Categorical(adata.obs["leiden"].astype(str)).categories)
        return max(1, len(groups) * n_genes)

    @staticmethod
    def _matrix_has_negative_values(X) -> bool:
        """Whether the plotted matrix is signed (z-scored) or non-negative.

        Only used to choose colour limits: a symmetric diverging scale is
        correct for signed input and a lie for log1p-normalised input, where the
        negative half of the bar names values the matrix cannot contain. Nothing
        is stored or reported from this read.
        """
        try:
            if sparse.issparse(X):
                data = X.data
            elif isinstance(X, np.ndarray):
                data = X
            else:
                # Backed / lazy matrix: sample the leading rows rather than
                # materialising the whole matrix just to pick a colour limit.
                data = np.asarray(X[: min(2000, int(X.shape[0]))])
            return bool(data.size) and bool(np.nanmin(data) < 0)
        except Exception:
            return False

    @classmethod
    def _plot_marker_dotplot(cls, adata, ctx: PipelineContext) -> None:
        """Top-5 marker dot plot at page geometry, with both channels named.

        rcParams cannot reach this panel: ``sc.pl.rank_genes_groups_dotplot``
        derives its own canvas from the category count, so 22 clusters x 5 genes
        came out 867 mm wide — 4.7x the double-column width. Whatever the
        typesetter does to fit that on a page takes the theme's 7 pt type down
        to roughly 1.5 pt with it. The colour bar also read "Mean expression in
        group" while ``standard_scale="var"`` was in force, so it named a
        quantity other than the one it showed.

        Genes move to the vertical axis because 110 gene labels cannot be
        resolved across a 183 mm width (1.7 mm each) while 22 cluster labels
        comfortably can.
        """
        from .._figure_theme import journal_figure_size

        try:
            n_rows = cls._marker_row_count(adata, cls._TOP_N_PER_CLUSTER)
            width_in, _ = journal_figure_size("double")
            height_mm = min(
                cls._MAX_FIGURE_HEIGHT_MM - cls._TITLE_BAND_MM,
                30.0 + cls._GENE_ROW_PITCH_MM * n_rows,
            )
            dot = sc.pl.rank_genes_groups_dotplot(
                adata,
                n_genes=cls._TOP_N_PER_CLUSTER,
                show=False,
                standard_scale="var",
                return_fig=True,
                figsize=(width_in, height_mm / 25.4),
            )
            dot.swap_axes()
            dot.style(
                cmap=cls._SEQUENTIAL_CMAP,
                dot_edge_lw=0.15,
                smallest_dot=1.0,
                largest_dot=30.0,
            )
            # Size is the reason this is a dot plot rather than a matrix plot,
            # so its legend is not optional; and the colour bar has to say that
            # standard_scale rescaled each gene, not report a bare mean.
            dot.legend(
                colorbar_title="Mean expression in cluster,\nmin-max scaled per gene",
                size_title="Cells in cluster with\nnon-zero detection (%)",
                width=1.25,
            )
            dot.make_figure()
            fig = dot.fig
            # scanpy lays the panel out inside matplotlib's default figure
            # margins, which on a page-tall canvas throws ~60 mm at blank
            # border. Reclaiming it is what gives each gene row the ~1.8 mm it
            # needs for a 5 pt label to clear its neighbour.
            fig.subplots_adjust(top=0.955, bottom=0.035)
            main_ax = (dot.get_axes() or {}).get("mainplot_ax")
            if main_ax is not None:
                main_ax.set_xlabel("Leiden cluster", fontsize=6)
                # Says what the row blocks and the brackets on the right mean;
                # without it the bracket numbers are unexplained.
                main_ax.set_ylabel(
                    "Marker genes, grouped by the cluster they rank for", fontsize=6
                )
                main_ax.tick_params(axis="y", labelsize=5)
            # DotPlot nests a grid inside its own grid, and the nested one
            # overruns the parent, so the dendrogram ends up above the canvas
            # top by an amount that depends on the figure's aspect. A title at
            # any fixed height therefore lands on the tree (or under the opaque
            # dendrogram axes, which is how it disappears). Place it above
            # whatever the tallest axes turned out to be; bbox_inches="tight"
            # grows the saved image to include it.
            top = max((a.get_position().y1 for a in fig.axes), default=1.0)
            fig.suptitle(
                f"Top {cls._TOP_N_PER_CLUSTER} marker genes per Leiden cluster "
                f"({cls._test_label(ctx)})\n"
                f"n = {adata.n_obs:,} cells, {n_rows // cls._TOP_N_PER_CLUSTER} clusters, "
                f"{adata.n_vars:,} genes tested per cluster",
                fontsize=7, y=max(1.0, top) + 0.008, va="bottom",
            )
            fig.savefig(ctx.figure_dir / "de_dotplot_top5.png", bbox_inches="tight")
            plt.close(fig)
        except Exception as exc:
            # Unchanged policy: a failed diagnostic never fails the module. Say
            # which one failed, so a missing panel is not silently invisible.
            logger.warning("DE dot plot skipped: %s", exc)
            plt.close("all")

    @classmethod
    def _plot_marker_heatmap(cls, adata, ctx: PipelineContext) -> None:
        """Top-5 marker heatmap on a colour scale the data can actually reach.

        The scale was ``vmin=-3, vmax=3`` on viridis. Those limits are written
        for z-scored input, but at DE time ``adata.X`` is log1p-normalised and
        strictly non-negative (0 .. ~7.9 on the LUSC cohort), so the lower half
        of the bar advertised negative expression that cannot occur while every
        real value above 3 saturated — and viridis puts a mid-tone at 0, which
        reads as "average" rather than "absent". The limits are now taken from
        the sign of the matrix: symmetric about zero on a diverging ramp only
        when the values really are signed, otherwise 0..3 on a sequential ramp
        whose lightest end is zero. The 3 is the module's existing clip
        magnitude, kept so only the impossible half of the range is dropped.
        """
        from .._figure_theme import journal_figure_size

        try:
            signed = cls._matrix_has_negative_values(adata.X)
            n_rows = cls._marker_row_count(adata, cls._TOP_N_PER_CLUSTER)
            width_in, _ = journal_figure_size("double")
            height_mm = min(
                cls._MAX_FIGURE_HEIGHT_MM - cls._TITLE_BAND_MM,
                34.0 + cls._GENE_ROW_PITCH_MM * n_rows,
            )
            axes = sc.pl.rank_genes_groups_heatmap(
                adata,
                n_genes=cls._TOP_N_PER_CLUSTER,
                show=False,
                show_gene_labels=True,
                use_raw=False,
                swap_axes=True,
                vmin=-3 if signed else 0,
                vmax=3,
                cmap=cls._DIVERGING_CMAP if signed else cls._SEQUENTIAL_CMAP,
                figsize=(width_in, height_mm / 25.4),
            )
            fig = plt.gcf()
            axes = axes if isinstance(axes, dict) else {}
            # scanpy lays the panel out inside matplotlib's default figure
            # margins, which on a page-tall canvas leaves ~60 mm of blank
            # border while the 110 gene labels overlap for want of 0.3 mm each.
            fig.subplots_adjust(top=0.955, bottom=0.035)
            heatmap_ax = axes.get("heatmap_ax")
            if heatmap_ax is not None:
                heatmap_ax.tick_params(axis="y", labelsize=5)
            groupby_ax = axes.get("groupby_ax")
            if groupby_ax is not None:
                # scanpy labels this axis with the obs key; "leiden" is not a
                # description of what the columns are. The cluster ticks are
                # spaced by cluster size, so the small clusters' labels collide
                # unless they are turned upright.
                groupby_ax.set_xlabel("Leiden cluster", fontsize=6)
                # Columns are cells, so a tick sits at each cluster's centre of
                # mass: the smallest clusters end up closer together than a
                # glyph is wide. Upright text on two alternating rows separates
                # them without dropping a cluster from the axis.
                from matplotlib.transforms import ScaledTranslation

                stagger = ScaledTranslation(0.0, -0.11, fig.dpi_scale_trans)
                for idx, label in enumerate(groupby_ax.get_xticklabels()):
                    label.set_rotation(90)
                    label.set_fontsize(5)
                    if idx % 2:
                        label.set_transform(label.get_transform() + stagger)
            # The colour bar is the one axis scanpy does not return; it is also
            # the only place the units of the matrix can be stated.
            returned = {id(a) for a in axes.values()}
            for cbar_ax in [a for a in fig.axes if id(a) not in returned]:
                # Outside the tick labels, not above the bar: the bar is 5 mm
                # wide and a title centred on it runs over the neighbouring
                # cluster colour strip.
                cbar_ax.yaxis.set_label_position("right")
                cbar_ax.set_ylabel(
                    "Expression (z-scored)" if signed
                    else "Expression (log1p-normalised)",
                    fontsize=5, labelpad=2,
                )
            fig.suptitle(
                f"Top {cls._TOP_N_PER_CLUSTER} marker genes per Leiden cluster, "
                f"single cells ({cls._test_label(ctx)})\n"
                f"n = {adata.n_obs:,} cells, {n_rows // cls._TOP_N_PER_CLUSTER} clusters, "
                f"{n_rows} gene rows (a gene recurs when it marks several clusters); "
                f"colour clipped at {'+/-3' if signed else '3'}",
                fontsize=7, y=1.0, va="bottom",
            )
            fig.savefig(ctx.figure_dir / "de_heatmap_top5.png", bbox_inches="tight")
            plt.close(fig)
        except Exception as exc:
            logger.warning("DE marker heatmap skipped: %s", exc)
            plt.close("all")

    @classmethod
    def _volcano_plot(
        cls,
        markers: pd.DataFrame,
        ctx: PipelineContext,
        sig_mask=None,
    ) -> None:
        """Volcano of every cluster-vs-rest test drawn in this run.

        What was wrong: a 10x7 in canvas that no journal column fits, whose
        colour channel re-encoded the x axis (log2 fold change mapped to a
        red/blue bar sitting next to an axis already labelled log2 fold change)
        while the one categorical distinction a volcano exists to show — passed
        vs did not pass — was never drawn. The guide lines were literal
        ``0.05`` / ``+/-0.25`` instead of the thresholds the run applied, so any
        run configured otherwise drew somebody else's cutoffs.

        ``sig_mask`` is the mask already computed in ``_run_impl``; no test,
        threshold or correction is recomputed here. When it is absent (a direct
        caller) every point is drawn neutrally and the pass count is left out
        rather than guessed.
        """
        if markers.empty:
            return
        from .._figure_theme import journal_figure_size

        p_cut = float(getattr(ctx.cfg, "de_pval_threshold", 0.05))
        lfc_cut = float(getattr(ctx.cfg, "de_logfc_threshold", 0.25))
        # The adjusted p is floored before the log so that exact zeros stay
        # finite. That floor is a drawing limit and becomes a visible ceiling of
        # stacked points, so the panel has to admit to it rather than let the
        # flat top read as data.
        p_floor = 1e-20
        lfc = np.asarray(markers["logfoldchanges"], dtype=float)
        log_pval = -np.log10(
            np.asarray(markers["pvals_adj"], dtype=float).clip(min=p_floor)
        )

        # Only pairs that can be positioned are drawn; a fallback DE frame with
        # no fold changes would otherwise contribute a legend entry and no dots.
        drawable = np.isfinite(lfc) & np.isfinite(log_pval)

        # Double column: this panel carries a three-line provenance caption that
        # is wider than any 89 mm axes, and at single-column width the caption set
        # the tight bbox and collapsed the axes to near-zero.
        fig, ax = plt.subplots(
            figsize=journal_figure_size("double", height_mm=88.0),
            constrained_layout=True,
        )
        # Colour carries the significance call and nothing else. Splitting it by
        # the sign of the fold change would repeat the x axis, which is the
        # defect this panel had; the pass/fail split is the one fact no axis
        # here shows, since it also depends on the detection-rate floor.
        if sig_mask is None:
            passed = np.zeros(len(markers), dtype=bool)
            series = [(drawable, cls._COLOUR_NS, None)]
        else:
            passed = np.asarray(sig_mask, dtype=bool)
            # Neutral cloud first so the passing tests sit on top of it.
            series = [
                (~passed & drawable, cls._COLOUR_NS, "not significant"),
                (passed & drawable, cls._COLOUR_SIG, "significant"),
            ]
        labelled = False
        for mask, colour, label in series:
            if not mask.any():
                continue  # a category with no points must not appear in the legend
            ax.scatter(
                lfc[mask], log_pval[mask], s=1.5, c=colour, alpha=0.65,
                edgecolors="none", rasterized=True, label=label,
            )
            labelled = labelled or label is not None

        # Guides are the thresholds this run applied, read from its config.
        ax.axhline(-np.log10(p_cut), color=cls._COLOUR_GUIDE,
                   linestyle="--", linewidth=0.5)
        # scanpy's ranking returns UP-regulated genes per cluster, so the frame
        # routinely contains no negative fold changes at all. Drawing the
        # symmetric -lfc_cut guide over a region that structurally cannot hold a
        # point invites the reading "not one cluster has a down-regulated
        # marker", which is an artefact of the ranking, not a result.
        one_sided = bool(drawable.any() and np.nanmin(lfc[drawable]) >= 0.0)
        bounds = (lfc_cut,) if one_sided else (-lfc_cut, lfc_cut)
        for bound in bounds:
            ax.axvline(bound, color=cls._COLOUR_GUIDE, linestyle="--", linewidth=0.5)
        if one_sided:
            ax.set_xlim(left=0.0)

        # `markers` is scanpy's RANKED table — the top de_n_genes rows per
        # cluster — not the test set. Calling its length "tests" overstated the
        # hit rate enormously (a 5,877/6,600 "significant" fraction reads as 89%
        # when the rate over all genes tested is ~1.5%), and a wrong n is the
        # single fastest way to lose a referee. Report what the rows actually
        # are, and derive the grouping from the frame instead of assuming it.
        group_col = next(
            (c for c in ("group", "cluster", "cell_type", "names_group") if c in markers.columns),
            None,
        )
        n_groups = int(markers[group_col].nunique()) if group_col else 0
        ranked_desc = f"{len(markers):,} ranked marker genes"
        if n_groups:
            ranked_desc += f" ({n_groups} clusters x top {len(markers) // n_groups})"

        ax.set_xlabel("log2 fold change (cluster vs all other cells)")
        ax.set_ylabel(f"-log10 adjusted P ({cls._correction_label(ctx)})")
        ax.set_title(
            f"Cluster marker genes ({cls._test_label(ctx)})\n"
            + ranked_desc
            + ("" if sig_mask is None else f"; {int(passed.sum()):,} pass the run's criteria"),
            fontsize=7,
        )
        # Legend and caption sit below the axes: the volcano's own empty corner
        # moves with the data, and anchoring them to the top collides with the
        # title. bbox_inches="tight" keeps both inside the saved canvas.
        if labelled:
            ax.legend(
                loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
                fontsize=5.5, handletextpad=0.2, columnspacing=1.0,
                borderpad=0.0, markerscale=4,
            )

        # The third significance criterion is a detection-rate floor, which has
        # no axis on a volcano; naming it keeps the pass count reconcilable with
        # the two cutoffs that are drawn.
        metadata = getattr(ctx, "metadata", {}) or {}
        min_pct = metadata.get("de_min_pct_applied")
        note = f"dashed guides: adjusted P < {p_cut:g} and |log2FC| > {lfc_cut:g}"
        if min_pct:
            # For most greys this floor, not either drawn guide, is the operative
            # criterion, so the visual grammar "grey = outside the dashed box"
            # would otherwise be wrong for nearly every grey point.
            note += (
                f"; also required, and not drawable here: detected in "
                f">= {float(min_pct):.0%} of the cluster"
            )
        if one_sided:
            note += (
                "\nranking returns up-regulated genes per cluster, so the table is "
                "one-sided; the absent left wing is a property of the ranking, not "
                "evidence that no gene is down-regulated"
            )
        n_floored = int((log_pval >= -np.log10(p_floor)).sum())
        if n_floored:
            # A dense line of points along the top edge is a drawing limit, not a
            # finding, and has to be named as one. It is a plotting clip, NOT
            # double-precision underflow: most clipped values are representable
            # and span many orders of magnitude below the floor.
            note += (
                f"\n{n_floored:,} of {len(markers):,} rows have adjusted P below the "
                f"{p_floor:g} drawing floor and are stacked on that ceiling"
            )
        import textwrap

        wrapped = "\n".join(
            "\n".join(textwrap.wrap(line, width=110)) for line in note.split("\n")
        )
        ax.text(
            0.5, -0.30, wrapped, transform=ax.transAxes, ha="center", va="top",
            fontsize=5, color=cls._COLOUR_GUIDE, linespacing=1.4,
        )

        fig.savefig(ctx.figure_dir / "de_volcano.png")
        plt.close(fig)

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
                    reason = (
                        "differential_expression substate DE: sparse subset too large "
                        f"to densify safely (substate={sub_shape}, parent={par_shape})"
                    )
                    logger.warning("%s — skipping substate '%s'", reason, substate_str)
                    ctx.metadata.setdefault("de_substate_skip_reasons", []).append(
                        {"substate": substate_str, "reason": reason}
                    )
                    continue
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

            fdr = DifferentialExpressionModule._adjust_pvals(
                pvals,
                getattr(ctx.cfg, "de_correction", None),
            )

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
