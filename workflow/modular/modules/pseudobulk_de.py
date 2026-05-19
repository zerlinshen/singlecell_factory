from __future__ import annotations

import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from ..context import PipelineContext


__references__ = {
    "Squair_pseudobulk_2021": {
        "title": "Confronting false discoveries in single-cell differential expression",
        "authors": "Squair et al.",
        "journal": "Nature Communications",
        "year": "2021",
        "doi": "10.1038/s41467-021-25960-2",
        "description": "Benchmark demonstrating pseudobulk DE outperforms per-cell DE \u2014 methodology applied here.",
    },
}


logger = logging.getLogger(__name__)
MIN_CELLS_PER_SAMPLE = 3
MIN_SAMPLES_PER_CONDITION = 2


class PseudobulkDEModule:
    """Pseudobulk DE with an explicit contrast contract.

    Supported modes:
    - confirmatory: explicit contrast(s) across a contrast column
    - exploratory: optional group-vs-rest summaries when explicitly enabled
    """

    name = "pseudobulk_de"

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
        if adata is None:
            raise ValueError("Pseudobulk DE requires AnnData.")

        cfg = getattr(ctx.cfg, "pseudobulk", None)
        min_cells = int(getattr(cfg, "min_cells_per_sample", MIN_CELLS_PER_SAMPLE))
        min_samples = int(getattr(cfg, "min_samples_per_condition", MIN_SAMPLES_PER_CONDITION))

        try:
            sample_col = self._resolve_sample_col(
                adata,
                getattr(cfg, "sample_col", None) or ctx.cfg.batch.batch_key,
            )
            group_col = self._resolve_group_col(
                adata,
                getattr(cfg, "group_col", "cell_type"),
            )
        except ValueError as exc:
            logger.warning("%s; skipping pseudobulk DE.", exc)
            ctx.metadata["pseudobulk_de_status"] = "skipped_missing_grouping_columns"
            ctx.metadata["pseudobulk_de_mode"] = "skipped"
            ctx.status(self.name, "skipped", str(exc))
            return

        if len(adata.obs[sample_col].unique()) < 2:
            msg = f"Only 1 sample in '{sample_col}'; pseudobulk DE needs >= 2. Skipping."
            logger.warning(msg)
            ctx.metadata["pseudobulk_de_status"] = "skipped_single_sample"
            ctx.metadata["pseudobulk_de_mode"] = "skipped"
            ctx.status(self.name, "skipped", msg)
            return

        if "counts" not in adata.layers:
            msg = "Missing adata.layers['counts']; pseudobulk DE requires raw UMI counts."
            logger.warning("%s; skipping pseudobulk DE.", msg)
            ctx.metadata["pseudobulk_de_status"] = "skipped_missing_counts_layer"
            ctx.metadata["pseudobulk_de_mode"] = "skipped"
            ctx.status(self.name, "skipped", msg)
            return

        explicit_contrasts, contrast_col, contract_message = self._resolve_confirmatory_contract(
            adata,
            ctx,
        )
        exploratory = bool(getattr(cfg, "exploratory_group_vs_rest", False))

        try:
            if explicit_contrasts:
                pb_counts, pb_meta = self._aggregate(
                    adata,
                    [sample_col, group_col, contrast_col],
                    min_cells_per_sample=min_cells,
                )
                results = self._run_confirmatory(
                    pb_counts,
                    pb_meta,
                    group_col,
                    contrast_col,
                    explicit_contrasts,
                    min_samples_per_condition=min_samples,
                )
                mode = "confirmatory"
                ctx.metadata["pseudobulk_de_contrast_contract"] = {
                    "mode": mode,
                    "sample_col": sample_col,
                    "group_col": group_col,
                    "contrast_col": contrast_col,
                    "contrasts": explicit_contrasts,
                }
            elif exploratory:
                pb_counts, pb_meta = self._aggregate(
                    adata,
                    [sample_col, group_col],
                    min_cells_per_sample=min_cells,
                )
                results = self._run_exploratory(
                    pb_counts,
                    pb_meta,
                    group_col,
                    min_samples_per_condition=min_samples,
                )
                mode = "exploratory"
                ctx.metadata["pseudobulk_de_contrast_contract"] = {
                    "mode": mode,
                    "sample_col": sample_col,
                    "group_col": group_col,
                    "contrast_policy": "group_vs_rest",
                }
            else:
                msg = contract_message or (
                    "No explicit pseudobulk contrast contract provided. Pass "
                    "--pseudobulk-contrast-col/--pseudobulk-contrast-a/--pseudobulk-contrast-b "
                    "for confirmatory pseudobulk or enable --pseudobulk-exploratory-group-vs-rest."
                )
                ctx.metadata["pseudobulk_de_status"] = "skipped_missing_contrast_contract"
                ctx.metadata["pseudobulk_de_mode"] = "skipped"
                ctx.status(self.name, "skipped", msg)
                return
        except ValueError as exc:
            logger.warning("%s; skipping pseudobulk DE.", exc)
            ctx.metadata["pseudobulk_de_status"] = "skipped_missing_counts_layer"
            ctx.metadata["pseudobulk_de_mode"] = "skipped"
            ctx.status(self.name, "skipped", str(exc))
            return

        if pb_counts.empty or pb_meta.empty:
            ctx.metadata["pseudobulk_de_status"] = "skipped_no_valid_pseudobulk_groups"
            ctx.metadata["pseudobulk_de_mode"] = "skipped"
            ctx.status(self.name, "skipped", "No pseudobulk groups passed minimum cell threshold.")
            return

        pb_counts.to_csv(ctx.table_dir / "pseudobulk_counts.csv", index=False)

        if results is None or results.empty:
            ctx.metadata["pseudobulk_de_status"] = "skipped_no_testable_contrasts"
            ctx.metadata["pseudobulk_de_mode"] = mode
            ctx.status(self.name, "skipped", "No pseudobulk contrasts met minimum sample requirements.")
            return

        results.to_csv(ctx.table_dir / "pseudobulk_de_results.csv", index=False)
        n_sig = int((results["padj"] < 0.05).sum()) if "padj" in results.columns else 0
        ctx.metadata["pseudobulk_de_significant_genes"] = n_sig
        ctx.metadata["pseudobulk_de_status"] = "completed"
        ctx.metadata["pseudobulk_de_mode"] = mode
        self._plot_volcano(results, ctx)
        self._plot_heatmap(results, pb_counts, pb_meta, sample_col, ctx)

    # -- contract helpers ------------------------------------------------
    @classmethod
    def _resolve_sample_col(cls, adata, preferred: str) -> str:
        return cls._resolve_col(
            adata,
            preferred,
            ("sample", "batch", "donor", "patient"),
            "sample/batch",
        )

    @classmethod
    def _resolve_group_col(cls, adata, preferred: str) -> str:
        if preferred:
            return cls._resolve_col(adata, preferred, ("cell_type", "leiden"), "grouping")
        return cls._resolve_col(adata, "cell_type", ("leiden",), "grouping")

    @classmethod
    def _resolve_optional_contrast_col(cls, adata, preferred: str | None) -> str | None:
        try:
            return cls._resolve_col(
                adata,
                preferred or "condition",
                ("condition", "disease", "group", "cohort"),
                "contrast",
            )
        except ValueError:
            return None

    def _resolve_confirmatory_contract(self, adata, ctx: PipelineContext):
        cfg = getattr(ctx.cfg, "pseudobulk", None)
        if cfg is None:
            return [], None, "No pseudobulk configuration found."

        contrast_specs: list[dict[str, str]] = []
        if getattr(cfg, "contrast_json", None):
            contrast_specs = self._load_contrast_json(Path(cfg.contrast_json))
        elif getattr(cfg, "contrast_a", None) and getattr(cfg, "contrast_b", None):
            contrast_specs = [
                {
                    "name": f"{cfg.contrast_a}_vs_{cfg.contrast_b}",
                    "contrast_a": str(cfg.contrast_a),
                    "contrast_b": str(cfg.contrast_b),
                }
            ]
        else:
            return [], None, "No explicit confirmatory pseudobulk contrast was provided."

        contrast_col = self._resolve_optional_contrast_col(adata, getattr(cfg, "contrast_col", None))
        if not contrast_col:
            return [], None, "Requested pseudobulk contrasts but no usable contrast column was found in adata.obs."

        for spec in contrast_specs:
            spec.setdefault("contrast_col", contrast_col)
        return contrast_specs, contrast_col, ""

    @staticmethod
    def _load_contrast_json(path: Path) -> list[dict[str, str]]:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, list):
            raise ValueError("pseudobulk contrast JSON must be a list of contrast specs.")
        specs: list[dict[str, str]] = []
        for item in payload:
            if not isinstance(item, dict):
                continue
            a = item.get("contrast_a")
            b = item.get("contrast_b")
            if not a or not b:
                continue
            specs.append(
                {
                    "name": str(item.get("name") or f"{a}_vs_{b}"),
                    "contrast_a": str(a),
                    "contrast_b": str(b),
                }
            )
        return specs

    # -- column helpers --------------------------------------------------
    @staticmethod
    def _resolve_col(adata, preferred: str, fallbacks: tuple, label: str) -> str:
        if preferred in adata.obs.columns:
            return preferred
        for c in fallbacks:
            if c in adata.obs.columns:
                logger.info("Using '%s' as %s column ('%s' not found).", c, label, preferred)
                return c
        raise ValueError(
            f"No {label} column found in adata.obs (tried '{preferred}' + {fallbacks})."
        )

    # -- aggregation -----------------------------------------------------
    @staticmethod
    def _aggregate(
        adata,
        group_cols: list[str] | str,
        group_col_legacy: str | None = None,
        min_cells_per_sample: int = MIN_CELLS_PER_SAMPLE,
    ):
        if isinstance(group_cols, str):
            cols = [group_cols]
            if group_col_legacy is not None:
                cols.append(group_col_legacy)
        else:
            cols = list(group_cols)
        raw = adata.layers["counts"]
        if hasattr(raw, "compute"):
            raw = raw.compute()
        elif hasattr(raw, "to_memory"):
            raw = raw.to_memory()
        if sparse.issparse(raw):
            raw = raw.tocsr()
        else:
            raw = sparse.csr_matrix(raw)
        if raw.shape != adata.shape:
            raise ValueError(
                f"counts layer shape {raw.shape} does not match adata shape {adata.shape}."
            )
        genes = adata.var_names.tolist()
        records, meta_rows = [], []
        for keys, idx in adata.obs.groupby(cols, observed=True).groups.items():
            if len(idx) < min_cells_per_sample:
                continue
            sub = raw[adata.obs.index.get_indexer(idx), :]
            sums = np.asarray(sub.sum(axis=0)).ravel()
            records.append(sums)
            if not isinstance(keys, tuple):
                keys = (keys,)
            row = {col: key for col, key in zip(cols, keys)}
            row["n_cells"] = len(idx)
            meta_rows.append(row)
        return pd.DataFrame(records, columns=genes), pd.DataFrame(meta_rows)

    def _run_confirmatory(
        self,
        pb_counts: pd.DataFrame,
        pb_meta: pd.DataFrame,
        group_col: str,
        contrast_col: str,
        contrasts: list[dict[str, str]],
        min_samples_per_condition: int,
    ) -> pd.DataFrame | None:
        results: list[pd.DataFrame] = []
        for group in pb_meta[group_col].unique():
            mask = pb_meta[group_col] == group
            group_meta = pb_meta.loc[mask].reset_index(drop=True)
            group_counts = pb_counts.loc[mask].reset_index(drop=True)
            for spec in contrasts:
                a = spec["contrast_a"]
                b = spec["contrast_b"]
                res = self._run_de(
                    group_counts,
                    group_meta,
                    contrast_col,
                    a,
                    b,
                    min_samples_per_condition=min_samples_per_condition,
                )
                if res is None:
                    continue
                res["group"] = f"{group}|{spec['name']}"
                res["contrast_col"] = contrast_col
                res["contrast_a"] = a
                res["contrast_b"] = b
                results.append(res)
        if not results:
            return None
        return pd.concat(results, ignore_index=True)

    def _run_exploratory(
        self,
        pb_counts: pd.DataFrame,
        pb_meta: pd.DataFrame,
        group_col: str,
        min_samples_per_condition: int,
    ) -> pd.DataFrame | None:
        groups = pb_meta[group_col].unique()
        if len(groups) < 2:
            return None

        all_results: list[pd.DataFrame] = []
        if len(groups) == 2:
            res = self._run_de(
                pb_counts,
                pb_meta,
                group_col,
                groups[0],
                groups[1],
                min_samples_per_condition=min_samples_per_condition,
            )
            if res is not None:
                res["group"] = f"{groups[0]}_vs_{groups[1]}"
                all_results.append(res)
        else:
            for grp in groups:
                m = pb_meta.copy()
                m["condition"] = np.where(m[group_col] == grp, grp, "rest")
                res = self._run_de(
                    pb_counts,
                    m,
                    "condition",
                    grp,
                    "rest",
                    min_samples_per_condition=min_samples_per_condition,
                )
                if res is not None:
                    res["group"] = f"{grp}_vs_rest"
                    all_results.append(res)
        if not all_results:
            return None
        return pd.concat(all_results, ignore_index=True)

    # -- DE dispatch -----------------------------------------------------
    def _run_de(
        self,
        counts,
        meta,
        cond_col,
        ca,
        cb,
        min_samples_per_condition: int = MIN_SAMPLES_PER_CONDITION,
    ):
        na, nb = (meta[cond_col] == ca).sum(), (meta[cond_col] == cb).sum()
        if na < min_samples_per_condition or nb < min_samples_per_condition:
            logger.info("Skipping %s vs %s: samples %d vs %d.", ca, cb, na, nb)
            return None
        try:
            return self._de_pydeseq2(counts, meta, cond_col, ca, cb)
        except ImportError:
            logger.info("pydeseq2 not available, falling back to Mann-Whitney U.")
        except Exception as exc:
            logger.warning("pydeseq2 failed (%s), falling back.", exc)
        try:
            from scipy.stats import mannwhitneyu

            return self._de_ranktest(
                counts,
                meta,
                cond_col,
                ca,
                cb,
                mannwhitneyu,
                two_sided_kw={"alternative": "two-sided"},
            )
        except Exception as exc:
            logger.warning("Mann-Whitney failed (%s), falling back to Wilcoxon.", exc)
        from scipy.stats import ranksums

        return self._de_ranktest(counts, meta, cond_col, ca, cb, ranksums)

    # -- backend: pydeseq2 ----------------------------------------------
    @staticmethod
    def _de_pydeseq2(counts, meta, cond_col, ca, cb):
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        ci = counts.round().astype(int)
        md = meta[[cond_col]].copy()
        md.index = ci.index = pd.RangeIndex(len(md))
        md[cond_col] = md[cond_col].astype(str)
        dds = DeseqDataSet(counts=ci, metadata=md, design=f"~{cond_col}")
        dds.deseq2()
        try:
            sr = DeseqStats(dds, contrast=[cond_col, str(ca), str(cb)])
        except TypeError:
            sr = DeseqStats(dds)
        sr.summary()
        r = sr.results_df.reset_index().rename(
            columns={"index": "gene", "log2FoldChange": "log2FC"}
        )
        return r[["gene", "log2FC", "pvalue", "padj"]].copy()

    # -- backend: generic rank test (MWU / ranksums) ---------------------
    @staticmethod
    def _de_ranktest(counts, meta, cond_col, ca, cb, test_fn, two_sided_kw=None):
        from .._mem_guard import MemoryGuard, MemoryAbortError
        kw = two_sided_kw or {}
        ia = meta[meta[cond_col] == ca].index
        ib = meta[meta[cond_col] == cb].index
        ma, mb = counts.loc[ia].values, counts.loc[ib].values

        # MEDIUM-4 fix (2026-05-19): the previous implementation computed
        # lfc as log2((mean(a)+1)/(mean(b)+1)) on raw pseudobulk counts,
        # which is biased against samples with deeper sequencing depth.
        # Normalize each pseudobulk to a per-sample CPM (counts per million)
        # before averaging so lfc is library-size-corrected. The rank test
        # itself is unaffected by per-sample monotone transforms, so the
        # p-value path stays equivalent.
        eps = 1.0
        def _cpm(matrix: np.ndarray) -> np.ndarray:
            totals = matrix.sum(axis=1, keepdims=True)
            totals = np.where(totals > 0, totals, 1.0)
            return matrix * 1e6 / totals

        ma_cpm = _cpm(ma.astype(np.float64))
        mb_cpm = _cpm(mb.astype(np.float64))

        genes, pvals, lfcs = [], [], []
        for j, g in enumerate(counts.columns):
            if MemoryGuard.abort_requested():
                raise MemoryAbortError("watchdog abort during pseudobulk DE gene loop")
            va, vb = ma[:, j], mb[:, j]
            va_cpm, vb_cpm = ma_cpm[:, j], mb_cpm[:, j]
            lfc = np.log2((va_cpm.mean() + eps) / (vb_cpm.mean() + eps))
            try:
                _, p = test_fn(va, vb, **kw)
            except ValueError:
                p = 1.0
            genes.append(g)
            pvals.append(p)
            lfcs.append(lfc)
        res = pd.DataFrame({"gene": genes, "log2FC": lfcs, "pvalue": pvals})
        res["padj"] = _bh_adjust(res["pvalue"].values)
        return res

    # -- plots -----------------------------------------------------------
    @staticmethod
    def _plot_volcano(results, ctx):
        if results.empty or "padj" not in results.columns:
            return
        fig, ax = plt.subplots(figsize=(10, 7))
        lp = -np.log10(results["padj"].clip(lower=1e-50))
        lfc = results["log2FC"]
        colors = np.where((results["padj"] < 0.05) & (lfc.abs() > 0.5), "tab:red", "grey")
        ax.scatter(lfc, lp, c=colors, s=6, alpha=0.6, edgecolors="none")
        ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.8)
        for v in (-0.5, 0.5):
            ax.axvline(v, color="grey", ls="--", lw=0.8)
        ax.set(
            xlabel="Log2 Fold Change",
            ylabel="-log10(adjusted p-value)",
            title="Pseudobulk DE -- Volcano Plot",
        )
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudobulk_volcano.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_heatmap(results, pb_counts, pb_meta, sample_col, ctx):
        if results.empty or "padj" not in results.columns:
            return
        top = (
            results[results["padj"] < 0.05]
            .sort_values("padj")
            .drop_duplicates("gene")
            .head(30)["gene"]
            .tolist()
        )
        avail = [g for g in top if g in pb_counts.columns]
        if len(avail) < 2:
            return
        mat = np.log2(pb_counts[avail] + 1)
        labels = pb_meta[sample_col].astype(str).values if sample_col in pb_meta.columns else np.arange(len(mat))
        fig, ax = plt.subplots(figsize=(max(8, len(avail) * 0.35), max(4, len(mat) * 0.4)))
        im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(avail)))
        ax.set_xticklabels(avail, rotation=90, fontsize=7)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set(ylabel="Sample", title="Top DE Genes -- Pseudobulk (log2 counts)")
        plt.colorbar(im, ax=ax, label="log2(count + 1)")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudobulk_heatmap.png", dpi=160, bbox_inches="tight")
        plt.close()


def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    if n == 0:
        return pvals.copy()
    order = np.argsort(pvals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n + 1)
    adj = pvals * n / ranks
    sorted_adj = adj[order]
    cummin = np.minimum.accumulate(sorted_adj[::-1])[::-1]
    adj[order] = cummin
    return np.clip(adj, 0.0, 1.0)
