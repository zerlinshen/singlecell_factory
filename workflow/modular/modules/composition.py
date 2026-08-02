from __future__ import annotations

import logging
import re
import textwrap
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..context import PipelineContext


__references__ = {
    "Buttner_scCODA_2021": {
        "title": "scCODA is a Bayesian model for compositional single-cell data analysis",
        "authors": "Buttner, Ostner, et al.",
        "journal": "Nature Communications",
        "year": "2021",
        "doi": "10.1038/s41467-021-27150-6",
        "description": "Bayesian compositional analysis; preferred path via pertpy. Chi-square fallback used when scCODA unavailable.",
    },
}


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CompositionDesign:
    """Resolved sample-level design used by composition inference."""

    sample_col: str
    condition_col: str | None
    contrast_a: str | None
    contrast_b: str | None
    covariates: tuple[str, ...]
    min_samples_per_condition: int
    formula: str | None
    sample_metadata: pd.DataFrame


class CompositionModule:
    """Optional module: differential cell type composition analysis across samples."""

    name = "composition"
    requires_keys = {"obs": ["cell_type"]}

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Composition analysis requires AnnData.")
        if "cell_type" not in adata.obs.columns:
            raise ValueError("Composition analysis requires 'cell_type' in adata.obs (run annotation first).")

        design = self._resolve_design(adata, ctx)
        analysis_adata = self._subset_to_contrast(adata, design)
        sample_key = design.sample_col
        logger.info(
            "Composition analysis: cell_type grouped by biological sample %r; "
            "condition=%r; formula=%r",
            sample_key,
            design.condition_col,
            design.formula,
        )

        # Build count and proportion tables
        count_df = (
            analysis_adata.obs.groupby([sample_key, "cell_type"], observed=True)
            .size()
            .unstack(fill_value=0)
        )
        prop_df = count_df.div(count_df.sum(axis=1), axis=0)

        count_df.to_csv(ctx.table_dir / "composition_counts.csv")
        prop_df.to_csv(ctx.table_dir / "composition_proportions.csv")
        ctx.metadata["composition_n_groups"] = int(count_df.shape[0])
        ctx.metadata["composition_n_cell_types"] = int(count_df.shape[1])

        ctx.metadata["composition_sample_col"] = sample_key
        ctx.metadata["composition_condition_col"] = design.condition_col
        ctx.metadata["composition_formula"] = design.formula
        ctx.metadata["composition_covariates"] = list(design.covariates)
        ctx.metadata["composition_contrast"] = (
            [design.contrast_a, design.contrast_b]
            if design.contrast_a is not None
            else None
        )

        # Try pertpy scCODA only for an explicit replicate-aware condition
        # contract. A sample column is an identifier, not an experimental
        # covariate; using C(sample) yields one coefficient per replicate and no
        # biological contrast while looking superficially confirmatory.
        #
        # CLAIMABILITY (added 2026-08-02). Cell-type proportions are COMPOSITIONAL: they
        # are constrained to sum to 1, so one cell type increasing mechanically decreases
        # the others. Testing each type independently with Mann-Whitney/chi-square
        # therefore produces spurious "significant" shifts in types that did not change --
        # which is exactly why scCODA exists (Buttner et al. 2021, Nat Commun 12:6876).
        # BH-FDR does not fix this: it corrects the multiplicity, not the invalid
        # independence assumption.
        #
        # Previously the fallback was taken at logger.info severity and recorded NOTHING,
        # so a run in the CPU env (where pertpy is absent) silently produced
        # compositionally-invalid statistics indistinguishable from the scCODA path.
        # pseudobulk_de.py already had the right pattern for this -- 86 claimability
        # markings vs 0 here -- so this module now follows it.
        if design.condition_col is None:
            test_results = pd.DataFrame(
                columns=["cell_type", "covariate", "significant"]
            )
            engine = "descriptive_only"
            inference_class = "sample_level_descriptive"
            inference_status = "descriptive_only_no_condition_contract"
            claimable = False
        else:
            test_results = self._try_pertpy(
                analysis_adata,
                design,
                random_state=int(getattr(ctx, "random_state", 42)),
            )

        if design.condition_col is not None and test_results is not None:
            engine = "sccoda_pertpy"
            inference_class = "compositional_bayesian"
            inference_status = "supported_confirmatory"
            claimable = True
        elif design.condition_col is not None:
            # Thread the canonical AC-10 seed (ctx.random_state) into the
            # permutation fallback so enrichment z-scores honor --random-state
            # (issue #20).
            random_state = int(getattr(ctx, "random_state", 42))
            condition_by_sample = (
                design.sample_metadata.set_index(sample_key)[design.condition_col]
                .astype(str)
                .reindex(count_df.index.astype(str))
            )
            condition_by_sample.index = count_df.index
            test_results = self._fallback_test(
                count_df,
                prop_df,
                random_state,
                condition_by_sample=condition_by_sample,
            )
            engine = "scipy_per_type_fallback"
            inference_class = "per_type_marginal_fallback"
            inference_status = "exploratory_nonclaimable_noncompositional_fallback"
            claimable = False
            logger.warning(
                "Composition: scCODA unusable (absent OR broken -- see the preceding "
                "warning for which) -- fell back to per-cell-type marginal "
                "tests. Proportions are compositional, so these p-values are NOT valid "
                "for a differential-abundance claim. Result marked non-claimable "
                "(inference_status=%s). Install pertpy/scCODA (declared in "
                "environment_gpu.yml) to obtain a claimable result.",
                inference_status,
            )

        ctx.metadata["composition_engine"] = engine
        ctx.metadata["composition_inference_class"] = inference_class
        ctx.metadata["composition_inference_status"] = inference_status
        ctx.metadata["composition_claimable"] = claimable

        test_results = test_results.copy()
        test_results["inference_class"] = inference_class
        test_results["inference_status"] = inference_status
        test_results["claimable"] = claimable
        test_results.to_csv(ctx.table_dir / "composition_test_results.csv", index=False)

        # Visualizations
        self._plot_barplot(prop_df, ctx)
        self._plot_boxplot(prop_df, ctx)

    @staticmethod
    def _resolve_sample_col(adata, ctx: PipelineContext) -> str:
        """Resolve the replicate identifier without treating a condition as one."""
        cfg = getattr(ctx.cfg, "composition", None)
        explicit = getattr(cfg, "sample_col", None)
        if explicit:
            if explicit not in adata.obs.columns:
                raise ValueError(
                    f"Composition sample column {explicit!r} is absent from adata.obs."
                )
            return str(explicit)

        batch_key = getattr(getattr(ctx.cfg, "batch", None), "batch_key", None)
        candidates = [batch_key, "sample", "donor", "patient", "batch"]
        seen: set[str] = set()
        for candidate in candidates:
            if not candidate or candidate in seen:
                continue
            seen.add(str(candidate))
            if candidate in adata.obs.columns and adata.obs[candidate].notna().any():
                return str(candidate)
        raise ValueError(
            "Composition analysis requires a biological sample column. Set "
            "--composition-sample-col; Leiden clusters and condition labels are not "
            "independent replicates."
        )

    @classmethod
    def _find_group_key(cls, adata, ctx: PipelineContext) -> str:
        """Backward-compatible name for the sample-column resolver."""
        return cls._resolve_sample_col(adata, ctx)

    @classmethod
    def _resolve_design(cls, adata, ctx: PipelineContext) -> CompositionDesign:
        cfg = getattr(ctx.cfg, "composition", None)
        sample_col = cls._resolve_sample_col(adata, ctx)
        condition_col = getattr(cfg, "condition_col", None)
        contrast_a = getattr(cfg, "contrast_a", None)
        contrast_b = getattr(cfg, "contrast_b", None)
        covariates = tuple(getattr(cfg, "covariates", ()) or ())
        min_replicates = int(getattr(cfg, "min_samples_per_condition", 2))

        if bool(contrast_a) != bool(contrast_b):
            raise ValueError(
                "Composition contrast requires both contrast_a and contrast_b."
            )
        if not condition_col:
            if contrast_a or contrast_b or covariates:
                raise ValueError(
                    "Composition contrasts/covariates require an explicit condition_col."
                )
            sample_metadata = (
                adata.obs[[sample_col]].drop_duplicates().astype(str).reset_index(drop=True)
            )
            return CompositionDesign(
                sample_col=sample_col,
                condition_col=None,
                contrast_a=None,
                contrast_b=None,
                covariates=(),
                min_samples_per_condition=min_replicates,
                formula=None,
                sample_metadata=sample_metadata,
            )

        condition_col = str(condition_col)
        if condition_col == sample_col:
            raise ValueError(
                "Composition sample_col and condition_col must be different; a sample "
                "identifier cannot be used as its own biological condition."
            )
        if sample_col in covariates:
            raise ValueError(
                "Composition sample identifier cannot be reintroduced as a model "
                "covariate."
            )
        if min_replicates < 2:
            raise ValueError(
                "Composition inference requires at least two biological samples per "
                "condition."
            )

        design_cols = [condition_col, *covariates]
        missing = [column for column in [sample_col, *design_cols]
                   if column not in adata.obs.columns]
        if missing:
            raise ValueError(
                f"Composition design columns are absent from adata.obs: {missing}"
            )
        if len(set(design_cols)) != len(design_cols):
            raise ValueError("Composition condition/covariate columns must be unique.")

        frame = adata.obs[[sample_col, *design_cols]].copy()
        if frame.isna().any(axis=None):
            raise ValueError("Composition design contains missing sample-level values.")
        for column in frame.columns:
            frame[column] = frame[column].astype(str)

        by_sample = frame.groupby(sample_col, observed=True, sort=False)
        conflicts = [
            column
            for column in design_cols
            if bool((by_sample[column].nunique(dropna=False) > 1).any())
        ]
        if conflicts:
            raise ValueError(
                "Each biological sample must map to one condition/covariate value; "
                f"conflicting columns: {conflicts}."
            )

        sample_metadata = frame.drop_duplicates(subset=[sample_col]).reset_index(drop=True)
        if contrast_a is not None:
            contrast_a, contrast_b = str(contrast_a), str(contrast_b)
            available = set(sample_metadata[condition_col])
            absent = [value for value in (contrast_a, contrast_b) if value not in available]
            if absent:
                raise ValueError(
                    f"Composition contrast labels are absent from {condition_col!r}: {absent}"
                )
            sample_metadata = sample_metadata[
                sample_metadata[condition_col].isin((contrast_a, contrast_b))
            ].reset_index(drop=True)

        condition_counts = sample_metadata[condition_col].value_counts()
        if len(condition_counts) < 2:
            raise ValueError(
                "Composition inference requires at least two modeled conditions."
            )
        underpowered = condition_counts[condition_counts < min_replicates]
        if not underpowered.empty:
            detail = ", ".join(
                f"{condition}={int(count)}"
                for condition, count in underpowered.items()
            )
            raise ValueError(
                "Composition inference has too few independent biological samples "
                f"(minimum {min_replicates} per condition; {detail})."
            )

        formula_columns = [condition_col, *covariates]
        invalid_formula_columns = [
            column for column in formula_columns
            if re.fullmatch(r"[A-Za-z_]\w*", str(column)) is None
        ]
        if invalid_formula_columns:
            raise ValueError(
                "Composition formula columns must be valid identifiers: "
                f"{invalid_formula_columns}"
            )
        formula = " + ".join([f"C({condition_col})", *covariates])
        return CompositionDesign(
            sample_col=sample_col,
            condition_col=condition_col,
            contrast_a=contrast_a,
            contrast_b=contrast_b,
            covariates=covariates,
            min_samples_per_condition=min_replicates,
            formula=formula,
            sample_metadata=sample_metadata,
        )

    @staticmethod
    def _subset_to_contrast(adata, design: CompositionDesign):
        if design.condition_col is None or design.contrast_a is None:
            return adata
        values = adata.obs[design.condition_col].astype(str)
        mask = values.isin((design.contrast_a, design.contrast_b)).to_numpy()
        return adata[mask].copy()

    @staticmethod
    def _try_pertpy(
        adata,
        design: CompositionDesign,
        *,
        random_state: int,
    ) -> pd.DataFrame | None:
        """Attempt scCODA compositional analysis via pertpy."""
        try:
            import pertpy as pt

            sccoda = pt.tl.Sccoda()
            sccoda_data = sccoda.load(
                adata, type="cell_level", generate_sample_level=True,
                cell_type_identifier="cell_type", sample_identifier=design.sample_col,
                covariate_obs=[design.condition_col, *design.covariates],
            )
            sccoda.prepare(
                sccoda_data,
                formula=design.formula,
                reference_cell_type="automatic",
            )
            sccoda.run_nuts(
                sccoda_data,
                num_warmup=500,
                num_samples=1000,
                rng_key=int(random_state),
            )
            result = sccoda.credible_effects(sccoda_data)

            # `credible_effects` returns a Series indexed by (Covariate, Cell Type), so
            # reset_index() yields THREE columns, not two. The previous code assigned
            # exactly two names and died with
            #   ValueError: Length mismatch: Expected axis has 3 elements, ...
            # That line had never actually executed: pertpy's import was broken by a
            # jax/numpyro incompatibility, so every run fell back long before reaching
            # here. Fixed 2026-08-02 after the import was repaired and a real-data run
            # (92,430 LUSC cells) surfaced it. Resolve columns by ROLE rather than by
            # position so a future pertpy layout change degrades to the honest
            # non-claimable fallback instead of mislabelling columns.
            df = result.reset_index()
            value_col = df.columns[-1]
            cell_type_col = next(
                (c for c in df.columns[:-1] if "cell" in str(c).lower()), None
            )
            if cell_type_col is None:
                raise ValueError(
                    "scCODA credible_effects has no identifiable cell-type column; "
                    f"got {list(df.columns)}"
                )
            covariate_cols = [c for c in df.columns[:-1] if c != cell_type_col]
            if not covariate_cols:
                raise ValueError(
                    "scCODA credible_effects has no identifiable model coefficient."
                )
            coefficient_values = df[covariate_cols[0]].astype(str)
            if not coefficient_values.str.contains(
                str(design.condition_col), regex=False
            ).any():
                raise ValueError(
                    "scCODA returned no coefficient for the explicit condition "
                    f"{design.condition_col!r}; got "
                    f"{sorted(coefficient_values.unique().tolist())}."
                )
            out = df[[cell_type_col, value_col]].copy()
            out.columns = ["cell_type", "significant"]
            # Keep the covariate so a multi-covariate design is not silently
            # collapsed into one row per cell type.
            out.insert(1, "covariate", coefficient_values)
            return out
        except Exception as exc:
            # Distinguish ABSENT from BROKEN. Real-data validation on 2026-08-02 showed
            # pertpy INSTALLED in sc_gpu yet failing on import, so a message saying
            # "unavailable" sent the reader looking for a missing package that was there.
            import importlib.util as _ilu
            try:
                _present = _ilu.find_spec("pertpy") is not None
            except (ImportError, ValueError):
                import sys as _sys
                _present = "pertpy" in _sys.modules
            logger.warning(
                "pertpy scCODA %s (%s: %s) -- falling back to per-cell-type marginal tests.",
                "is INSTALLED but FAILED" if _present else "is NOT INSTALLED",
                type(exc).__name__, exc,
            )
            return None

    @staticmethod
    def _fallback_test(
        count_df: pd.DataFrame,
        prop_df: pd.DataFrame,
        random_state: int = 42,
        *,
        condition_by_sample: pd.Series | None = None,
    ) -> pd.DataFrame:
        """Statistical testing of composition differences using scipy.

        ``random_state`` seeds the single-group permutation enrichment so the
        z-scores are reproducible under ``--random-state`` (issue #20).
        """
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
            rng = np.random.default_rng(random_state)
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

        # Multi-group statistical tests. On the production path each value is
        # one biological sample and the groups are experimental conditions.
        # ``condition_by_sample=None`` is retained only for backwards-compatible
        # direct unit calls; CompositionModule.run never treats sample IDs as
        # conditions.
        if condition_by_sample is not None:
            condition_by_sample = condition_by_sample.reindex(prop_df.index)
            if condition_by_sample.isna().any():
                missing = condition_by_sample[condition_by_sample.isna()].index.tolist()
                raise ValueError(
                    f"Composition condition metadata missing for samples: {missing}"
                )
            condition_levels = list(dict.fromkeys(condition_by_sample.astype(str)))
        else:
            condition_levels = [str(value) for value in prop_df.index.unique()]

        pvals = []
        for ct in cell_types:
            if condition_by_sample is not None:
                groups = [
                    prop_df.loc[
                        condition_by_sample.astype(str) == condition, ct
                    ].to_numpy()
                    for condition in condition_levels
                ]
            else:
                groups = [prop_df[ct][prop_df.index == g].values for g in prop_df.index.unique()]
            if len(groups) == 2:
                _, p = stats.mannwhitneyu(
                    groups[0], groups[1], alternative="two-sided"
                )
            else:
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
    def _present_cell_type_colours(prop_df: pd.DataFrame) -> tuple[list[str], dict[str, str]]:
        """Cell types that actually contribute data, and one stable colour each.

        Two things depend on this being shared by both panels. First, a cell type
        absent from every group draws nothing but is still listed by the legend,
        which hands the reader a colour key for data that is not on the figure.
        Second, the colours must be assigned from the same list in the same order
        in both figures, so a cell type keeps its colour between the stacked bars
        and the boxplot; ``qualitative_colors`` guarantees they are distinct
        rather than wrapping a fixed palette modulo its length.
        """
        from .._figure_theme import qualitative_colors

        present = [ct for ct in prop_df.columns if float(prop_df[ct].sum()) > 0]
        colours = qualitative_colors("cell_type_qualitative", len(present))
        return present, dict(zip(present, colours))

    @staticmethod
    def _compact_labels(values) -> tuple[list[str], bool]:
        """Group labels shortened only as far as they stay one-to-one.

        Multi-study sample IDs run to ~50 characters here, and rotated under 71
        bars they take more of the canvas than the bars do. Blind truncation is
        not an option — samples that differ only in a suffix would collapse onto
        one label — so the shortest middle-elided form that is still a unique
        relabelling wins, and the unabbreviated IDs stay in the CSV tables.
        """
        labels = [str(v) for v in values]
        distinct = len(set(labels))
        if max((len(label) for label in labels), default=0) <= 24:
            return labels, False
        for limit in (18, 24, 30, 36):
            elided = [
                label if len(label) <= limit else f"{label[:6]}…{label[7 - limit:]}"
                for label in labels
            ]
            if len(set(elided)) == distinct:
                return elided, True
        return labels, False

    # Characters of 7 pt text that fit one journal column (~1.23 mm per
    # character).
    _SUBTITLE_WRAP = {"single": 70, "double": 145}

    @classmethod
    def _wrap_subtitle(cls, text: str, column: str) -> str:
        """Wrap a caption to the column width, with the lines evenly filled.

        A caption wider than the figure is not clipped — savefig's tight bbox
        grows the canvas to hold it — so an unwrapped subtitle silently widens a
        figure past the column it was sized for. Wrapping greedily at the limit
        would then leave a one-word orphan line, so the text is divided evenly
        over the number of lines it needs.
        """
        limit = cls._SUBTITLE_WRAP[column]
        if len(text) <= limit:
            return text
        lines = -(-len(text) // limit)
        return textwrap.fill(text, width=-(-len(text) // lines) + 6)

    @staticmethod
    def _group_sizes(prop_df: pd.DataFrame, ctx: PipelineContext) -> pd.Series | None:
        """Cells per group for the n annotation, or ``None`` if unobtainable.

        A stacked proportion bar hides its own denominator: a bar built from 36
        cells is drawn exactly like one built from 6,009, so nothing on the panel
        tells the reader which composition estimates are essentially noise. The
        sizes are re-read from the same obs column the proportions were tabulated
        on (``prop_df.index.name`` is the grouping key) purely to annotate n —
        they feed no test, threshold, or table.
        """
        obs = getattr(getattr(ctx, "adata", None), "obs", None)
        key = prop_df.index.name
        if obs is None or not key or key not in getattr(obs, "columns", []):
            return None
        sizes = obs[key].value_counts().reindex(prop_df.index)
        if sizes.isna().any():
            return None
        return sizes.astype(int)

    @staticmethod
    def _plot_barplot(prop_df: pd.DataFrame, ctx: PipelineContext) -> None:
        """Stacked bar chart of cell type proportions per group.

        rcParams cannot repair this panel, because every defect was passed as an
        argument: the canvas grew with the cohort (0.8 in per group is a 57-inch
        figure at 71 samples, which the typesetter then scales to a column and in
        doing so scales 7 pt type down to about 1 pt), ``tab20`` is neither
        colourblind-safe nor large enough to stay distinct past 20 cell types,
        and the legend listed every column whether or not it drew a segment.
        The substantive omission was the denominator: proportions were shown with
        no n anywhere, so a 100 % bar from 36 cells read like a real result. The
        companion strip states cells per group on a log axis (36 to 6,009 in the
        LUSC cohort — a linear axis would flatten most of it to nothing).
        """
        from .._figure_theme import journal_figure_size

        present, colours = CompositionModule._present_cell_type_colours(prop_df)
        absent = int(prop_df.shape[1] - len(present))
        sizes = CompositionModule._group_sizes(prop_df, ctx)
        n_groups = int(prop_df.shape[0])
        group_key = prop_df.index.name or "group"
        x = np.arange(n_groups)

        # Rotated group labels, not the bars, drive the figure height, so the
        # height is budgeted explicitly: bars + size strip + chrome (two-line
        # title, legend, padding) + the longest label (~0.19 mm per character per
        # point, measured). Deriving it keeps the bars a usable ~45 mm tall for
        # any cohort instead of letting a 71-label axis squeeze them to a stripe.
        labels, abbreviated = CompositionModule._compact_labels(prop_df.index)
        label_size = 5.0 if n_groups <= 90 else 4.0
        longest = max((len(label) for label in labels), default=0)
        bars_mm, strip_ratio, chrome_mm = 45.0, 3.6, 26.0
        height_mm = (
            bars_mm + chrome_mm
            + (0.0 if sizes is None else bars_mm / strip_ratio)
            + min(48.0, longest * label_size * 0.19)
        )
        # A handful of bars stretched across 183 mm is dead space, not a figure.
        column = "single" if n_groups <= 8 else "double"

        nrows = 1 if sizes is None else 2
        fig, axes = plt.subplots(
            nrows, 1, sharex=True, constrained_layout=True,
            figsize=journal_figure_size(column, height_mm=height_mm),
            gridspec_kw=None if nrows == 1 else {"height_ratios": (1.0, strip_ratio)},
        )
        axes = np.atleast_1d(axes)
        ax_bar = axes[-1]

        if sizes is not None:
            ax_n = axes[0]
            ax_n.bar(x, sizes.to_numpy(), width=0.9, color="#666666", linewidth=0)
            ax_n.set_yscale("log")
            ax_n.minorticks_off()  # a log minor-tick ladder in a 12 mm strip is noise
            ax_n.set_ylabel("Cells")
            # 71 shared tick marks under a 12 mm strip read as a comb, not a scale.
            ax_n.tick_params(labelbottom=False, bottom=False)

        bottom = np.zeros(n_groups, dtype=float)
        for cell_type in present:
            values = prop_df[cell_type].to_numpy(dtype=float)
            ax_bar.bar(x, values, bottom=bottom, width=0.9, linewidth=0,
                       color=colours[cell_type], label=str(cell_type))
            bottom += values
        if not present:
            ax_bar.text(0.5, 0.5, "no cell type present in any group",
                        ha="center", va="center", transform=ax_bar.transAxes,
                        fontsize=6, color="#666666")
        ax_bar.set_xticks(x)
        ax_bar.set_xticklabels(labels, rotation=90, fontsize=label_size)
        ax_bar.set_xlim(-0.6, n_groups - 0.4)
        ax_bar.set_ylim(0.0, 1.0)
        ax_bar.set_ylabel("Proportion of cells")
        group_label = group_key.replace("_", " ").capitalize()
        ax_bar.set_xlabel(
            f"{group_label} (labels abbreviated; full IDs in composition_proportions.csv)"
            if abbreviated else group_label
        )

        subtitle = f"n = {n_groups} groups"
        if sizes is not None:
            subtitle = (
                f"n = {int(sizes.sum()):,} cells in {n_groups} groups "
                f"({int(sizes.min()):,}-{int(sizes.max()):,} cells per group, "
                f"median {round(float(sizes.median())):,}; top panel, log scale)"
            )
        if absent:
            subtitle += f"; {absent} cell type(s) absent from all groups, not shown"
        subtitle = CompositionModule._wrap_subtitle(subtitle, column)
        fig.suptitle(
            f"Cell type composition per {group_key.replace('_', ' ')}\n{subtitle}",
            fontsize=7,
        )
        fig.align_ylabels(axes)
        # Below the rotated labels: a legend inside the axes would cover bars, and
        # one to the right would take width from what is already a dense x axis.
        fig.legend(*ax_bar.get_legend_handles_labels(), loc="outside lower center",
                   ncol=min(3 if column == "single" else 4, max(1, len(present))),
                   fontsize=6, handlelength=1.2, columnspacing=1.2, borderaxespad=0.0)
        fig.savefig(ctx.figure_dir / "composition_barplot.png")
        plt.close(fig)

    @staticmethod
    def _plot_boxplot(prop_df: pd.DataFrame, ctx: PipelineContext) -> None:
        """Distribution of each cell type's proportion across groups.

        A box summarises a centre and a spread, so the panel has to say which
        ones: unlabelled boxes were previously drawn on an 8-inch canvas titled
        with no n, which leaves the reader unable to tell an IQR from a CI or to
        know how many groups each box rests on. Boxes are ordered by median and
        carry the same cell type colours as the stacked bars.
        """
        from .._figure_theme import journal_figure_size

        present, colours = CompositionModule._present_cell_type_colours(prop_df)
        absent = int(prop_df.shape[1] - len(present))
        n_groups = int(prop_df.shape[0])
        group_key = prop_df.index.name or "group"
        # Descending median: the informative comparison is between cell types,
        # and the tabulation order is arbitrary.
        order = sorted(present, key=lambda ct: float(prop_df[ct].median()), reverse=True)
        series = [prop_df[ct].to_numpy(dtype=float) for ct in order]
        x = np.arange(len(order))

        # A single column fits a handful of boxes without stretching them across
        # dead space; wider cohorts of cell types need the full width.
        column = "single" if len(order) <= 10 else "double"
        fig, ax = plt.subplots(
            figsize=journal_figure_size(column, height_mm=62.0), constrained_layout=True
        )

        if not order:
            ax.text(0.5, 0.5, "no cell type present in any group",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=6, color="#666666")
            ax.set_axis_off()
        elif n_groups < 2:
            # One group means there is no across-group distribution to summarise;
            # a box drawn from a single value would imply a spread that does not
            # exist, so plot the observed value and say so in the title (an
            # in-axes note would sit on top of the points at some proportions).
            ax.scatter(x, [values[0] for values in series], s=8,
                       color=[colours[ct] for ct in order], zorder=3)
        else:
            box = ax.boxplot(
                series, positions=x, widths=0.6, showfliers=True, patch_artist=True,
                medianprops={"color": "black", "linewidth": 0.8},
                whiskerprops={"color": "#333333", "linewidth": 0.5},
                capprops={"color": "#333333", "linewidth": 0.5},
                flierprops={"marker": "o", "markersize": 1.6, "markerfacecolor": "none",
                            "markeredgecolor": "#444444", "markeredgewidth": 0.3},
            )
            for patch, cell_type in zip(box["boxes"], order):
                patch.set_facecolor(colours[cell_type])
                patch.set_edgecolor("#333333")
                patch.set_linewidth(0.5)

        if order:
            ax.set_xticks(x)
            ax.set_xticklabels([str(ct) for ct in order], rotation=45, ha="right")
            ax.set_xlim(-0.7, len(order) - 0.3)
            ax.set_ylim(bottom=0.0)
            ax.set_ylabel(f"Proportion of cells per {group_key.replace('_', ' ')}")

        # The spread has to be named: an unlabelled box could be an IQR, an SD or
        # a CI, and the three support different claims.
        subtitle = (
            f"n = {n_groups} groups per box; box = median and IQR, "
            "whiskers 1.5x IQR, dots = outliers"
            if n_groups >= 2 else
            "single group: one observed value per cell type, no spread defined"
        )
        if not order:
            subtitle = "nothing to plot"
        if absent:
            subtitle += f"; {absent} absent cell type(s) not shown"
        subtitle = CompositionModule._wrap_subtitle(subtitle, column)
        # suptitle, not set_title: the caption is wider than the axes, and
        # centring it on the figure keeps it from hanging off the left edge.
        fig.suptitle(
            f"Cell type proportion by {group_key.replace('_', ' ')}\n{subtitle}",
            fontsize=7,
        )
        fig.savefig(ctx.figure_dir / "composition_boxplot.png")
        plt.close(fig)
