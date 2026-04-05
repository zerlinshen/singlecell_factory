from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)


class CompositionModule:
    """Optional module: differential cell type composition analysis across samples."""

    name = "composition"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Composition analysis requires AnnData.")
        if "cell_type" not in adata.obs.columns:
            raise ValueError("Composition analysis requires 'cell_type' in adata.obs (run annotation first).")

        group_key = self._find_group_key(adata, ctx)
        logger.info("Composition analysis: cell_type grouped by '%s'", group_key)

        # Build count and proportion tables
        count_df = (
            adata.obs.groupby([group_key, "cell_type"], observed=True)
            .size()
            .unstack(fill_value=0)
        )
        prop_df = count_df.div(count_df.sum(axis=1), axis=0)

        count_df.to_csv(ctx.table_dir / "composition_counts.csv")
        prop_df.to_csv(ctx.table_dir / "composition_proportions.csv")
        ctx.metadata["composition_n_groups"] = int(count_df.shape[0])
        ctx.metadata["composition_n_cell_types"] = int(count_df.shape[1])

        # Try pertpy scCODA first, fall back to scipy tests
        test_results = self._try_pertpy(adata, group_key)
        if test_results is None:
            test_results = self._fallback_test(count_df, prop_df)

        test_results.to_csv(ctx.table_dir / "composition_test_results.csv", index=False)

        # Visualizations
        self._plot_barplot(prop_df, ctx)
        self._plot_boxplot(prop_df, ctx)

    @staticmethod
    def _find_group_key(adata, ctx: PipelineContext) -> str:
        """Detect the sample/condition column to group by."""
        batch_key = ctx.cfg.batch.batch_key
        if batch_key in adata.obs.columns and adata.obs[batch_key].nunique() > 1:
            return batch_key
        for candidate in ("sample", "batch", "donor", "patient", "condition"):
            if candidate in adata.obs.columns and adata.obs[candidate].nunique() > 1:
                return candidate
        # Single-sample fallback: use leiden clusters as groups
        if "leiden" in adata.obs.columns:
            logger.warning("No multi-sample grouping found; using leiden clusters as groups.")
            return "leiden"
        raise ValueError("No suitable grouping column found for composition analysis.")

    @staticmethod
    def _try_pertpy(adata, group_key: str) -> pd.DataFrame | None:
        """Attempt scCODA compositional analysis via pertpy."""
        try:
            import pertpy as pt

            sccoda = pt.tl.Sccoda()
            sccoda_data = sccoda.load(
                adata, type="cell_level", generate_sample_level=True,
                cell_type_identifier="cell_type", sample_identifier=group_key,
            )
            sccoda.prepare(sccoda_data, formula=f"C({group_key})", reference_cell_type="automatic")
            sccoda.run_nuts(sccoda_data, num_warmup=500, num_samples=1000)
            result = sccoda.credible_effects(sccoda_data)
            df = result.reset_index()
            df.columns = ["cell_type", "significant"]
            return df
        except Exception as exc:
            logger.info("pertpy scCODA unavailable or failed (%s), using fallback tests.", exc)
            return None

    @staticmethod
    def _fallback_test(count_df: pd.DataFrame, prop_df: pd.DataFrame) -> pd.DataFrame:
        """Statistical testing of composition differences using scipy."""
        from scipy import stats
        from statsmodels.stats.multitest import multipletests

        cell_types = prop_df.columns.tolist()
        n_groups = prop_df.shape[0]
        results = []

        if n_groups < 2:
            # Single group: permutation-based z-scores
            logger.warning("Only 1 group found; computing permutation-based enrichment z-scores.")
            total_cells = count_df.values.sum()
            observed_props = count_df.iloc[0] / count_df.iloc[0].sum()
            rng = np.random.default_rng(42)
            all_labels = []
            for ct, cnt in count_df.iloc[0].items():
                all_labels.extend([ct] * cnt)
            all_labels = np.array(all_labels)
            n_perms = 1000
            perm_props = np.zeros((n_perms, len(cell_types)))
            for i in range(n_perms):
                shuffled = rng.permutation(all_labels)
                half = shuffled[: len(shuffled) // 2]
                unique, counts = np.unique(half, return_counts=True)
                perm_map = dict(zip(unique, counts / len(half)))
                for j, ct in enumerate(cell_types):
                    perm_props[i, j] = perm_map.get(ct, 0.0)
            for j, ct in enumerate(cell_types):
                obs_p = float(observed_props.get(ct, 0.0))
                mean_p = perm_props[:, j].mean()
                std_p = perm_props[:, j].std()
                z = (obs_p - mean_p) / std_p if std_p > 0 else 0.0
                results.append({
                    "cell_type": ct, "proportion_mean": round(obs_p, 4),
                    "proportion_std": 0.0, "pvalue": np.nan,
                    "padj": np.nan, "z_score": round(z, 3), "significant": abs(z) > 1.96,
                })
            return pd.DataFrame(results)

        # Multi-group statistical tests
        pvals = []
        for ct in cell_types:
            values = prop_df[ct].values
            if n_groups == 2:
                _, p = stats.mannwhitneyu(
                    values[prop_df.index == prop_df.index[0]],
                    values[prop_df.index != prop_df.index[0]],
                    alternative="two-sided",
                )
            else:
                groups = [prop_df[ct][prop_df.index == g].values for g in prop_df.index.unique()]
                _, p = stats.kruskal(*groups)
            pvals.append(p)

        reject, padj, _, _ = multipletests(pvals, method="fdr_bh")
        for ct, pv, pa, sig in zip(cell_types, pvals, padj, reject):
            results.append({
                "cell_type": ct,
                "proportion_mean": round(float(prop_df[ct].mean()), 4),
                "proportion_std": round(float(prop_df[ct].std()), 4),
                "pvalue": round(pv, 6), "padj": round(pa, 6),
                "significant": bool(sig),
            })
        return pd.DataFrame(results)

    @staticmethod
    def _plot_barplot(prop_df: pd.DataFrame, ctx: PipelineContext) -> None:
        """Stacked bar chart of cell type proportions per group."""
        prop_df.plot(kind="bar", stacked=True, figsize=(max(8, len(prop_df) * 0.8), 5), colormap="tab20")
        plt.ylabel("Proportion")
        plt.xlabel("Sample / Group")
        plt.title("Cell type composition")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "composition_barplot.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_boxplot(prop_df: pd.DataFrame, ctx: PipelineContext) -> None:
        """Boxplot of proportions for each cell type across groups."""
        melted = prop_df.reset_index().melt(
            id_vars=prop_df.index.name or "index",
            var_name="cell_type", value_name="proportion",
        )
        fig, ax = plt.subplots(figsize=(max(8, len(prop_df.columns) * 0.7), 5))
        cell_types = prop_df.columns.tolist()
        positions = range(len(cell_types))
        for i, ct in enumerate(cell_types):
            vals = melted.loc[melted["cell_type"] == ct, "proportion"].values
            ax.boxplot(vals, positions=[i], widths=0.5, showfliers=True)
        ax.set_xticks(list(positions))
        ax.set_xticklabels(cell_types, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Proportion")
        ax.set_title("Cell type proportion distribution across groups")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "composition_boxplot.png", dpi=160, bbox_inches="tight")
        plt.close()
