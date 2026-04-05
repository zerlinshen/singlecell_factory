from __future__ import annotations

import logging

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from ..context import PipelineContext

logger = logging.getLogger(__name__)
MIN_CELLS_PER_SAMPLE = 3
MIN_SAMPLES_PER_CONDITION = 2


class PseudobulkDEModule:
    """Pseudobulk DE: aggregate counts by sample+group, then bulk-level DE test.

    Backends: 1) pydeseq2  2) Mann-Whitney U  3) Wilcoxon rank-sum.
    Does NOT modify adata.X or adata.obs (read-only / append-only metadata).
    """

    name = "pseudobulk_de"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Pseudobulk DE requires AnnData.")

        try:
            sample_col = self._resolve_col(
                adata,
                ctx.cfg.batch.batch_key,
                ("sample", "batch", "donor", "patient"),
                "sample/batch",
            )
            group_col = self._resolve_col(adata, "cell_type", ("leiden",), "grouping")
        except ValueError as exc:
            logger.warning("%s; skipping pseudobulk DE.", exc)
            ctx.metadata["pseudobulk_de_status"] = "skipped_missing_grouping_columns"
            return

        if len(adata.obs[sample_col].unique()) < 2:
            msg = f"Only 1 sample in '{sample_col}'; pseudobulk DE needs >= 2. Skipping."
            logger.warning(msg)
            ctx.metadata["pseudobulk_de_status"] = "skipped_single_sample"
            return

        pb_counts, pb_meta = self._aggregate(adata, sample_col, group_col)
        if pb_counts.empty or pb_meta.empty:
            ctx.metadata["pseudobulk_de_status"] = "skipped_no_valid_pseudobulk_groups"
            return
        pb_counts.to_csv(ctx.table_dir / "pseudobulk_counts.csv", index=False)

        groups = pb_meta[group_col].unique()
        all_results: list[pd.DataFrame] = []
        if len(groups) == 2:
            res = self._run_de(pb_counts, pb_meta, group_col, groups[0], groups[1])
            if res is not None:
                res["group"] = f"{groups[0]}_vs_{groups[1]}"
                all_results.append(res)
        else:
            for grp in groups:
                m = pb_meta.copy()
                m["condition"] = np.where(m[group_col] == grp, grp, "rest")
                res = self._run_de(pb_counts, m, "condition", grp, "rest")
                if res is not None:
                    res["group"] = f"{grp}_vs_rest"
                    all_results.append(res)

        if not all_results:
            ctx.metadata["pseudobulk_de_status"] = "no_testable_groups"
            return

        results = pd.concat(all_results, ignore_index=True)
        results.to_csv(ctx.table_dir / "pseudobulk_de_results.csv", index=False)
        n_sig = int((results["padj"] < 0.05).sum()) if "padj" in results.columns else 0
        ctx.metadata["pseudobulk_de_significant_genes"] = n_sig
        ctx.metadata["pseudobulk_de_status"] = "completed"
        self._plot_volcano(results, ctx)
        self._plot_heatmap(results, pb_counts, pb_meta, sample_col, ctx)

    # -- column helpers --------------------------------------------------
    @staticmethod
    def _resolve_col(adata, preferred: str, fallbacks: tuple, label: str) -> str:
        if preferred in adata.obs.columns:
            return preferred
        for c in fallbacks:
            if c in adata.obs.columns:
                logger.info("Using '%s' as %s column ('%s' not found).", c, label, preferred)
                return c
        raise ValueError(f"No {label} column found in adata.obs "
                         f"(tried '{preferred}' + {fallbacks}).")

    # -- aggregation -----------------------------------------------------
    @staticmethod
    def _aggregate(adata, sample_col: str, group_col: str):
        raw = adata.raw.X if adata.raw is not None else adata.X
        genes = (adata.raw.var_names if adata.raw is not None else adata.var_names).tolist()
        records, meta_rows = [], []
        for (smp, grp), idx in adata.obs.groupby(
                [sample_col, group_col], observed=True).groups.items():
            if len(idx) < MIN_CELLS_PER_SAMPLE:
                continue
            sub = raw[adata.obs.index.get_indexer(idx), :]
            sums = np.asarray(sub.sum(axis=0)).ravel()
            records.append(sums)
            meta_rows.append({sample_col: smp, group_col: grp, "n_cells": len(idx)})
        return pd.DataFrame(records, columns=genes), pd.DataFrame(meta_rows)

    # -- DE dispatch -----------------------------------------------------
    def _run_de(self, counts, meta, cond_col, ca, cb):
        na, nb = (meta[cond_col] == ca).sum(), (meta[cond_col] == cb).sum()
        if na < MIN_SAMPLES_PER_CONDITION or nb < MIN_SAMPLES_PER_CONDITION:
            logger.info("Skipping %s vs %s: samples %d vs %d.", ca, cb, na, nb)
            return None
        # Backend 1: pydeseq2
        try:
            return self._de_pydeseq2(counts, meta, cond_col)
        except ImportError:
            logger.info("pydeseq2 not available, falling back to Mann-Whitney U.")
        except Exception as exc:
            logger.warning("pydeseq2 failed (%s), falling back.", exc)
        # Backend 2: Mann-Whitney U
        try:
            from scipy.stats import mannwhitneyu
            return self._de_ranktest(counts, meta, cond_col, ca, cb, mannwhitneyu,
                                     two_sided_kw={"alternative": "two-sided"})
        except Exception as exc:
            logger.warning("Mann-Whitney failed (%s), falling back to Wilcoxon.", exc)
        # Backend 3: Wilcoxon rank-sum
        from scipy.stats import ranksums
        return self._de_ranktest(counts, meta, cond_col, ca, cb, ranksums)

    # -- backend: pydeseq2 ----------------------------------------------
    @staticmethod
    def _de_pydeseq2(counts, meta, cond_col):
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        ci = counts.round().astype(int)
        md = meta[[cond_col]].copy()
        md.index = ci.index = pd.RangeIndex(len(md))
        md[cond_col] = md[cond_col].astype(str)
        dds = DeseqDataSet(counts=ci, metadata=md, design=f"~{cond_col}")
        dds.deseq2()
        sr = DeseqStats(dds); sr.summary()
        r = sr.results_df.reset_index().rename(columns={
            "index": "gene", "log2FoldChange": "log2FC"})
        return r[["gene", "log2FC", "pvalue", "padj"]].copy()

    # -- backend: generic rank test (MWU / ranksums) ---------------------
    @staticmethod
    def _de_ranktest(counts, meta, cond_col, ca, cb, test_fn, two_sided_kw=None):
        kw = two_sided_kw or {}
        ia = meta[meta[cond_col] == ca].index
        ib = meta[meta[cond_col] == cb].index
        ma, mb = counts.loc[ia].values, counts.loc[ib].values
        genes, pvals, lfcs = [], [], []
        for j, g in enumerate(counts.columns):
            va, vb = ma[:, j], mb[:, j]
            lfc = np.log2((va.mean() + 1) / (vb.mean() + 1))
            try:
                _, p = test_fn(va, vb, **kw)
            except ValueError:
                p = 1.0
            genes.append(g); pvals.append(p); lfcs.append(lfc)
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
        ax.set(xlabel="Log2 Fold Change", ylabel="-log10(adjusted p-value)",
               title="Pseudobulk DE -- Volcano Plot")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudobulk_volcano.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_heatmap(results, pb_counts, pb_meta, sample_col, ctx):
        if results.empty or "padj" not in results.columns:
            return
        top = (results[results["padj"] < 0.05].sort_values("padj")
               .drop_duplicates("gene").head(30)["gene"].tolist())
        avail = [g for g in top if g in pb_counts.columns]
        if len(avail) < 2:
            return
        mat = np.log2(pb_counts[avail] + 1)
        labels = pb_meta[sample_col].astype(str).values
        fig, ax = plt.subplots(figsize=(max(8, len(avail) * 0.35), max(4, len(mat) * 0.4)))
        im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(avail))); ax.set_xticklabels(avail, rotation=90, fontsize=7)
        ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=8)
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
