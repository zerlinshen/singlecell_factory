from __future__ import annotations

import json
import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from ._scanpy_compat import has_api, import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "scanpy": {
        "title": "SCANPY: large-scale single-cell gene expression data analysis",
        "authors": "Wolf, Angerer, Theis",
        "journal": "Genome Biology",
        "year": "2018",
        "doi": "10.1186/s13059-017-1382-0",
        "description": "scanpy.pp.calculate_qc_metrics / pp.filter_genes / pp.filter_cells semantics used here.",
    },
    "Luecken_Theis_2019": {
        "title": "Current best practices in single-cell RNA-seq analysis: a tutorial",
        "authors": "Luecken, Theis",
        "journal": "Molecular Systems Biology",
        "year": "2019",
        "doi": "10.15252/msb.20188746",
        "description": "Canonical QC threshold guidance (mito %, n_genes/cell, doublet detection).",
    },
}



class QCModule:
    """Mandatory module: compute QC metrics, filter low-quality cells/genes, and generate QC plots."""

    name = "qc"
    required = True
    # QC replaces ctx.adata with a cell/gene subset.  It must never execute on
    # a copy-on-branch appending lane because positional obsm data (including
    # ATAC peaks) have to be sliced on the canonical object.
    mutates_structure = True

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("QC requires loaded AnnData.")

        gene_upper = adata.var_names.astype(str).str.upper()
        adata.var["mt"] = gene_upper.str.startswith("MT-")
        adata.var["ribo"] = gene_upper.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = gene_upper.str.startswith(("HBA", "HBB"))
        self._assert_gene_classes_detectable(adata, ctx)
        use_scanpy_qc = has_api(sc, "pp.calculate_qc_metrics") and ctx.cfg.scale_mode != "massive"
        if use_scanpy_qc:
            try:
                sc.pp.calculate_qc_metrics(
                    adata,
                    qc_vars=["mt", "ribo", "hb"],
                    inplace=True,
                    log1p=False,
                    # percent_top must include 50 so pct_counts_in_top_50 is
                    # computed; ambient_correction._evaluate_triggers gates the T1
                    # DecontX trigger on it (Luecken & Theis 2019; scanpy percent_top QC).
                    percent_top=[20, 50],
                )
            except Exception:
                self._calculate_qc_metrics_fallback(adata)
        else:
            self._calculate_qc_metrics_fallback(adata)

        # scanpy names the percent_top column 'pct_counts_in_top_50_genes'; the
        # ambient_correction T1 trigger (and the whole repo) gates on the shorter
        # 'pct_counts_in_top_50'. Normalize so the trigger is live and the
        # keep-column logic below matches (the fallback already writes the short name).
        if (
            "pct_counts_in_top_50_genes" in adata.obs.columns
            and "pct_counts_in_top_50" not in adata.obs.columns
        ):
            adata.obs["pct_counts_in_top_50"] = adata.obs["pct_counts_in_top_50_genes"]

        # Drop verbose columns to save memory. KEEP pct_counts_in_top_50: the
        # ambient_correction T1 DecontX trigger gates on it (Luecken & Theis 2019).
        self._ensure_axis_tables_support_drop(adata)
        obs_drop = [
            c for c in adata.obs.columns
            if c.startswith("log1p")
            or c.startswith("total_counts_")
            or (c.startswith("pct_counts_in_top_") and c != "pct_counts_in_top_50")
        ]
        adata.obs.drop(columns=obs_drop, inplace=True, errors="ignore")
        var_drop = [c for c in adata.var.columns if c.startswith("log1p") or c.startswith("pct_dropout") or "total_counts" in c or c.startswith("mean_counts") or c.startswith("n_cells_by") or c in ["mt", "ribo", "hb"]]
        adata.var.drop(columns=var_drop, inplace=True, errors="ignore")

        # --- QC visualizations (before filtering) ---
        if ctx.cfg.scale_mode != "massive":
            self._plot_qc_metrics(adata, ctx, suffix="pre_filter")

        qc = ctx.cfg.qc
        before = adata.n_obs
        # Use one combined mask to avoid multiple full AnnData copies.
        criteria = {
            "n_genes_ge_min": adata.obs["n_genes_by_counts"] >= qc.min_genes,
            "n_genes_le_max": adata.obs["n_genes_by_counts"] <= qc.max_genes,
            "total_counts_ge_min": adata.obs["total_counts"] >= qc.min_counts,
            "total_counts_le_max": adata.obs["total_counts"] <= qc.max_counts,
            "pct_counts_mt_le_max": adata.obs["pct_counts_mt"] <= qc.max_mito_pct,
            "pct_counts_ribo_le_max": adata.obs["pct_counts_ribo"] <= qc.max_ribo_pct,
        }
        mask = np.logical_and.reduce([np.asarray(v, dtype=bool) for v in criteria.values()])
        threshold_audit = self._build_threshold_audit(adata, qc, criteria, mask)
        adata = adata[mask, :].copy()
        if has_api(sc, "pp.filter_genes"):
            sc.pp.filter_genes(adata, min_cells=qc.min_cells)
        else:
            adata = self._filter_genes_fallback(adata, min_cells=qc.min_cells)

        # --- QC visualizations (after filtering) ---
        if ctx.cfg.scale_mode != "massive":
            self._plot_qc_metrics(adata, ctx, suffix="post_filter")

        ctx.adata = adata
        ctx.metadata["cells_after_qc"] = int(adata.n_obs)
        ctx.metadata["genes_after_qc"] = int(adata.n_vars)
        ctx.metadata["qc_removed_cells"] = int(before - adata.n_obs)
        ctx.metadata["qc_retention_pct"] = round(float(adata.n_obs) / max(before, 1) * 100.0, 3)
        threshold_audit["cells_after_qc"] = int(adata.n_obs)
        threshold_audit["genes_after_qc"] = int(adata.n_vars)
        threshold_audit["cell_retention_pct"] = ctx.metadata["qc_retention_pct"]
        (ctx.table_dir / "qc_threshold_audit.json").write_text(
            json.dumps(threshold_audit, indent=2, sort_keys=True),
            encoding="utf-8",
        )

    @staticmethod
    def _metric_quantiles(values) -> dict[str, float]:
        arr = np.asarray(values, dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return {}
        qs = np.quantile(arr, [0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0])
        keys = ["min", "p01", "p05", "p25", "median", "p75", "p95", "p99", "max"]
        return {key: round(float(value), 6) for key, value in zip(keys, qs)}

    def _assert_gene_classes_detectable(self, adata, ctx: PipelineContext) -> None:
        """Refuse to report a QC threshold as applied when it cannot possibly fire.

        ``mt``/``ribo``/``hb`` are matched by gene-symbol prefix. If the gene axis
        is not in symbol space the flags are all-False, ``pct_counts_mt`` is 0.0
        for every cell, and ``pct_counts_mt <= max_mito_pct`` filters nothing
        while the manifest still lists the threshold. Mitochondrial fraction is a
        primary quality covariate (Luecken & Theis 2019, *Mol Syst Biol*
        15:e8746; Heumos et al. 2023, *Nat Rev Genet* 24:550-572), so a silently
        inert filter admits stressed and dying cells into every downstream
        result. Fail loudly instead.

        A dataset can legitimately contain no mitochondrial genes (pre-filtered
        atlas, targeted panel). That case is distinguished from a namespace
        failure by the ingest-recorded gene-namespace status and is a loud
        warning plus a machine-readable flag, not an error.
        """
        counts = {cls_: int(adata.var[cls_].sum()) for cls_ in ("mt", "ribo", "hb")}
        ctx.metadata["qc_gene_class_counts"] = counts

        namespace = ctx.metadata.get("gene_namespace") or {}
        namespace_status = namespace.get("status", "unknown")
        mito_filter_active = float(ctx.cfg.qc.max_mito_pct) < 100.0

        if counts["mt"] > 0:
            return

        ctx.metadata["qc_mito_filter_status"] = "inactive_no_mitochondrial_genes_detected"
        if namespace_status == "ensembl_unresolved" and mito_filter_active:
            raise ValueError(
                "QC contract violation: no mitochondrial genes are detectable, because "
                f"the gene axis is Ensembl-indexed and no symbol column was found "
                f"(gene_namespace={namespace!r}). pct_counts_mt would be 0.0 for every "
                f"cell, so --max-mito-pct={ctx.cfg.qc.max_mito_pct} would filter nothing "
                "while still being recorded as applied. Provide a symbol column in "
                "adata.var (e.g. 'feature_name' per the CZ CELLxGENE schema), or set "
                "--max-mito-pct 100 to declare explicitly that no mitochondrial filter "
                "is intended."
            )
        logger.warning(
            "QC: no mitochondrial genes detected in %d genes (gene_namespace=%s). "
            "pct_counts_mt is 0.0 for all cells and --max-mito-pct=%.1f will not "
            "remove any cell. Recorded as qc_mito_filter_status=%s.",
            adata.n_vars, namespace_status, ctx.cfg.qc.max_mito_pct,
            ctx.metadata["qc_mito_filter_status"],
        )

    @classmethod
    def _build_threshold_audit(cls, adata, qc, criteria: dict[str, object], mask) -> dict:
        failure_counts = {
            name: int((~np.asarray(value, dtype=bool)).sum())
            for name, value in criteria.items()
        }
        failure_counts["any_cell_filter"] = int((~np.asarray(mask, dtype=bool)).sum())
        metric_quantiles = {
            col: cls._metric_quantiles(adata.obs[col])
            for col in [
                "n_genes_by_counts",
                "total_counts",
                "pct_counts_mt",
                "pct_counts_ribo",
            ]
            if col in adata.obs
        }
        return {
            "threshold_mode": "fixed_configurable_audit_only",
            "adaptivity_status": (
                "not_adaptive; defaults are operator-configurable and must not be changed "
                "without a second dataset or a data-shape conditional"
            ),
            "thresholds": {
                "min_genes": int(qc.min_genes),
                "max_genes": int(qc.max_genes),
                "min_counts": int(qc.min_counts),
                "max_counts": int(qc.max_counts),
                "max_mito_pct": float(qc.max_mito_pct),
                "max_ribo_pct": float(qc.max_ribo_pct),
                "min_cells": int(qc.min_cells),
            },
            "cells_before_qc": int(adata.n_obs),
            "genes_before_qc": int(adata.n_vars),
            "candidate_cells_after_cell_filters": int(np.asarray(mask, dtype=bool).sum()),
            "failure_counts": failure_counts,
            "metric_quantiles": metric_quantiles,
        }

    @staticmethod
    def _plot_qc_metrics(adata, ctx: PipelineContext, suffix: str) -> None:
        """Generate violin + scatter QC plots."""
        if not has_api(sc, "pl.violin") or not has_api(sc, "pl.scatter"):
            return
        # Violin plots.
        #
        # The point of this panel is to justify the QC cutoffs, so it draws them.
        # Previously it was a 16x4 inch canvas (no journal column fits it) of
        # scanpy violins labelled "value", with counts on a linear axis where a
        # single 7e6-UMI outlier flattened the whole distribution to a line, and
        # with metrics that are identically zero (e.g. mitochondrial % on an atlas
        # whose curators removed MT genes) drawn as a fake flat "distribution" on
        # a +/-0.04 axis. A reader could not tell an inactive metric from a real
        # one, nor see which cells a threshold removes.
        from .._figure_theme import journal_figure_size

        qc_cfg = ctx.cfg.qc
        specs = [
            ("n_genes_by_counts", "Genes per cell", True,
             (float(qc_cfg.min_genes), float(qc_cfg.max_genes))),
            ("total_counts", "UMI counts per cell", True,
             (float(qc_cfg.min_counts), float(qc_cfg.max_counts))),
            ("pct_counts_mt", "Mitochondrial %", False, (None, float(qc_cfg.max_mito_pct))),
            ("pct_counts_ribo", "Ribosomal %", False, (None, float(qc_cfg.max_ribo_pct))),
        ]
        fig, axes = plt.subplots(
            1, len(specs), figsize=journal_figure_size("double", 58.0), constrained_layout=True
        )
        for idx, (ax, (key, title, log_scale, bounds)) in enumerate(zip(axes, specs)):
            ax.set_title(title)
            ax.set_xticks([])
            ax.text(-0.02, 1.10, chr(ord("a") + idx), transform=ax.transAxes,
                    ha="right", va="top", fontweight="bold", fontsize=8)
            if key not in adata.obs.columns:
                ax.text(0.5, 0.5, "not computed", ha="center", va="center",
                        transform=ax.transAxes, fontsize=6, color="#666666")
                ax.set_axis_off()
                continue
            values = np.asarray(adata.obs[key], dtype=float)
            values = values[np.isfinite(values)]
            # A metric with no spread carries no QC information; say so rather
            # than drawing a flat line on an autoscaled axis that implies one.
            if values.size == 0 or np.nanmax(values) == np.nanmin(values):
                constant = float(values[0]) if values.size else float("nan")
                ax.text(0.5, 0.5, f"no variation\n(constant {constant:g};\nthreshold inactive)",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=6, color="#666666")
                ax.set_axis_off()
                continue

            positive = values[values > 0]
            use_log = bool(log_scale and positive.size)
            plotted = np.log10(positive) if use_log else values
            ax.violinplot([plotted], showextrema=False, widths=0.7)
            ax.set_ylabel(f"log10({title.split(' per ')[0].lower()})" if use_log else "%")

            # Draw the cutoffs this run actually applies, and report how many
            # cells each one removes — that is the decision the panel supports.
            removed = 0
            for bound in bounds:
                if bound is None or not np.isfinite(bound):
                    continue
                y = np.log10(bound) if (use_log and bound > 0) else bound
                if use_log and bound <= 0:
                    continue
                ax.axhline(y, color="#D55E00", linewidth=0.6, linestyle="--")
            lo, hi = bounds
            if lo is not None:
                removed += int((values < lo).sum())
            if hi is not None:
                removed += int((values > hi).sum())
            ax.text(0.02, 0.02, f"{removed:,} outside", transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=5.5, color="#D55E00")

        fig.suptitle(
            f"QC metrics {suffix.replace('_', ' ')} — n = {adata.n_obs:,} cells; "
            "dashed lines = applied thresholds",
            fontsize=7,
        )
        fig.savefig(ctx.figure_dir / f"qc_violin_{suffix}.png")
        plt.close(fig)

        # Scatter panels: genes vs UMI counts, and mitochondrial fraction vs UMI.
        #
        # These were two scanpy scatters on a 12x5 inch canvas (no journal column
        # fits it) drawing one opaque point per cell on linear axes. At ~90k cells
        # that is a saturated blob: the overplotted core hides where the cells
        # actually lie, which is the only thing a QC scatter is for, and a handful
        # of 5e4-UMI cells compress the bulk of the data into the left tenth of the
        # axis. Panel b additionally drew a metric that is identically zero on
        # atlases whose curators removed MT genes as if it were a distribution.
        #
        # Draw the DENSITY instead (hexbin on log count axes) so the distribution
        # the panel claims to show is visible, and keep the axes box: unlike a UMAP
        # these axes carry real units.
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

        def _flag_uninformative(ax, message: str) -> None:
            """State that a panel carries no information rather than faking one."""
            ax.text(0.5, 0.5, message, ha="center", va="center",
                    transform=ax.transAxes, fontsize=6, color="#666666")
            ax.set_axis_off()

        def _finite_positive(*bounds) -> list[float]:
            return [b for b in bounds if b is not None and np.isfinite(b) and b > 0]

        # Decide up front whether the mitochondrial panel can say anything, because
        # that decides the layout: on an atlas whose curators removed MT genes
        # pct_counts_mt is constant 0, and giving a panel that carries one sentence
        # half of a 183 mm column is the same defect (dead canvas) as the 12x5 inch
        # figure this replaces. The slot is kept — the reader must still be told the
        # threshold is inert — but it is narrowed to what the sentence needs.
        mito_message = None
        if not {"pct_counts_mt", "total_counts"} <= set(adata.obs.columns):
            mito_message = "not computed"
        else:
            mito = np.asarray(adata.obs["pct_counts_mt"], dtype=float)
            finite_mito = mito[np.isfinite(mito)]
            # Same rule as the violins: a metric with no spread is reported as
            # inactive, never drawn as a density on an autoscaled axis.
            if finite_mito.size == 0 or finite_mito.max() == finite_mito.min():
                constant = float(finite_mito[0]) if finite_mito.size else float("nan")
                mito_message = f"no variation\n(constant {constant:g};\nthreshold inactive)"

        # Height follows the layout: two data panels sit side by side at roughly
        # square proportions, whereas one wide data panel needs a taller canvas or
        # the density is squashed into a letterbox that hides the low-count tail.
        fig, axes = plt.subplots(
            1, 2,
            figsize=journal_figure_size("double", 102.0 if mito_message else 78.0),
            constrained_layout=True,
            gridspec_kw={"width_ratios": [1.0, 0.24 if mito_message else 1.0]},
        )
        for idx, ax in enumerate(axes):
            ax.text(-0.02, 1.10, chr(ord("a") + idx), transform=ax.transAxes,
                    ha="right", va="top", fontweight="bold", fontsize=8)

        def _density_panel(ax, x, y, *, y_log: bool, title: str, xlabel: str,
                           ylabel: str, vlines=(), hlines=()) -> None:
            """Hexbin density of one QC pair, with the applied cutoffs drawn."""
            keep = np.isfinite(x) & np.isfinite(y) & (x > 0)
            if y_log:
                keep &= y > 0
            n_drawn = int(keep.sum())
            if n_drawn == 0:
                ax.set_title(title)
                _flag_uninformative(ax, "no cells with positive counts")
                return
            # Bin count follows the cell count: 46 bins across a ~90k-cell cloud
            # resolve the ridge, but the same grid over a few hundred cells is one
            # cell per hexagon — confetti, not a density.
            gridsize = int(np.clip(round(np.sqrt(n_drawn) / 3.0), 10, 46))
            hexes = ax.hexbin(
                x[keep], y[keep], gridsize=gridsize, mincnt=1, linewidths=0.0,
                xscale="log", yscale="log" if y_log else "linear",
                # Perceptually uniform and colourblind-safe (the same ramp as the
                # theme's continuous token), on a log colour scale because bin
                # occupancy spans several orders of magnitude — a linear colour
                # scale renders everything but the mode as a single shade.
                cmap="viridis", norm=LogNorm(),
            )
            bar = fig.colorbar(hexes, ax=ax, fraction=0.046, pad=0.02)
            bar.set_label("cells per hexagonal bin")
            bar.outline.set_visible(False)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            # Powers of ten alone put a single labelled tick on a 200-7000 gene
            # axis, so the reader cannot recover a value from the panel. Label the
            # 2x and 5x decade steps too, as plain numbers.
            for axis, axis_is_log in ((ax.xaxis, True), (ax.yaxis, y_log)):
                if not axis_is_log:
                    continue
                axis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 5.0)))
                axis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{value:g}"))
                axis.set_minor_formatter(NullFormatter())
            # n belongs on the panel, and it is the number of cells this panel can
            # actually draw: a log axis cannot show a zero, so cells dropped by it
            # are named rather than silently plotted away.
            hidden = int(x.size - n_drawn)
            ax.set_title(f"{title} (n = {n_drawn:,} cells)")
            if hidden:
                ax.text(0.98, 0.03,
                        f"{hidden:,} cell{'' if hidden == 1 else 's'} with zero "
                        "counts not shown",
                        transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=5.5, color="#666666",
                        # The corner can hold data; keep the note readable there
                        # rather than moving it somewhere it means less.
                        bbox={"facecolor": "white", "edgecolor": "none",
                              "alpha": 0.75, "pad": 1.0})
            # A cutoff outside the observed range must not rescale the panel, or
            # the cells themselves get squashed into a corner by an annotation.
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            for value in vlines:
                ax.axvline(value, color="#D55E00", linewidth=0.5, linestyle="--")
            for value in hlines:
                ax.axhline(value, color="#D55E00", linewidth=0.5, linestyle="--")
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        count_bounds = _finite_positive(float(qc_cfg.min_counts), float(qc_cfg.max_counts))

        if {"total_counts", "n_genes_by_counts"} <= set(adata.obs.columns):
            _density_panel(
                axes[0],
                np.asarray(adata.obs["total_counts"], dtype=float),
                np.asarray(adata.obs["n_genes_by_counts"], dtype=float),
                y_log=True,
                title="Genes vs UMI counts per cell",
                xlabel="UMI counts per cell (log scale)",
                ylabel="Genes detected per cell (log scale)",
                vlines=count_bounds,
                hlines=_finite_positive(float(qc_cfg.min_genes), float(qc_cfg.max_genes)),
            )
        else:
            axes[0].set_title("Genes vs UMI counts per cell")
            _flag_uninformative(axes[0], "not computed")

        if mito_message is not None:
            axes[1].set_title("Mitochondrial %")
            _flag_uninformative(axes[1], mito_message)
        else:
            _density_panel(
                axes[1],
                np.asarray(adata.obs["total_counts"], dtype=float),
                mito,
                y_log=False,
                title="Mitochondrial % vs UMI counts",
                xlabel="UMI counts per cell (log scale)",
                ylabel="Mitochondrial counts (%)",
                vlines=count_bounds,
                hlines=_finite_positive(float(qc_cfg.max_mito_pct)),
            )

        fig.suptitle(
            f"QC joint distributions {suffix.replace('_', ' ')} — "
            "dashed lines = applied thresholds",
            fontsize=7,
        )
        fig.savefig(ctx.figure_dir / f"qc_scatter_{suffix}.png")
        plt.close(fig)


    @staticmethod
    def _ensure_axis_tables_support_drop(adata) -> None:
        for axis in ("obs", "var"):
            table = getattr(adata, axis)
            if hasattr(table, "drop"):
                continue
            if hasattr(table, "to_memory"):
                table = table.to_memory()
            elif hasattr(table, "copy"):
                table = table.copy()
            setattr(adata, axis, table)

    @staticmethod
    def _materialize_matrix(x):
        if sparse.issparse(x):
            return x.tocsr()
        if hasattr(x, "to_memory"):
            try:
                x = x.to_memory()
            except Exception:
                pass
        if sparse.issparse(x):
            return x.tocsr()
        if hasattr(x, "compute"):
            try:
                x = x.compute()
            except Exception:
                pass
        return x

    @staticmethod
    def _calculate_qc_metrics_fallback(adata) -> None:
        x = QCModule._materialize_matrix(adata.X)
        if sparse.issparse(x):
            x = x.tocsr()
            total_counts = np.asarray(x.sum(axis=1)).ravel()
            n_genes = np.asarray((x > 0).sum(axis=1)).ravel()
            mt_counts = np.asarray(x[:, np.asarray(adata.var["mt"], dtype=bool)].sum(axis=1)).ravel()
            ribo_counts = np.asarray(x[:, np.asarray(adata.var["ribo"], dtype=bool)].sum(axis=1)).ravel()
            hb_counts = np.asarray(x[:, np.asarray(adata.var["hb"], dtype=bool)].sum(axis=1)).ravel()
        else:
            arr = np.asarray(x, dtype=float)
            total_counts = arr.sum(axis=1)
            n_genes = (arr > 0).sum(axis=1)
            mt_counts = arr[:, np.asarray(adata.var["mt"], dtype=bool)].sum(axis=1)
            ribo_counts = arr[:, np.asarray(adata.var["ribo"], dtype=bool)].sum(axis=1)
            hb_counts = arr[:, np.asarray(adata.var["hb"], dtype=bool)].sum(axis=1)
        denom = np.clip(total_counts, 1e-9, None)
        adata.obs["n_genes_by_counts"] = n_genes
        adata.obs["total_counts"] = total_counts
        adata.obs["pct_counts_mt"] = mt_counts / denom * 100.0
        adata.obs["pct_counts_ribo"] = ribo_counts / denom * 100.0
        adata.obs["pct_counts_hb"] = hb_counts / denom * 100.0
        # pct_counts_in_top_50: % of a cell's counts in its 50 highest-expressed
        # genes. Mirrors scanpy percent_top=[50]; the ambient_correction T1
        # DecontX trigger gates on this column (Luecken & Theis 2019).
        adata.obs["pct_counts_in_top_50"] = QCModule._pct_counts_in_top_n(x, total_counts, 50)

    @staticmethod
    def _pct_counts_in_top_n(x, total_counts, n_top: int) -> np.ndarray:
        """Per-cell percentage of counts contained in the top-n expressed genes."""
        denom = np.clip(np.asarray(total_counts, dtype=float), 1e-9, None)
        if sparse.issparse(x):
            x = x.tocsr()
            n_obs = x.shape[0]
            top_counts = np.zeros(n_obs, dtype=float)
            for i in range(n_obs):
                row = x.data[x.indptr[i]:x.indptr[i + 1]]
                if row.size == 0:
                    continue
                if row.size > n_top:
                    # Top n_top values (partition is O(nnz), not full sort).
                    top_counts[i] = np.partition(row, row.size - n_top)[row.size - n_top:].sum()
                else:
                    top_counts[i] = row.sum()
        else:
            arr = np.asarray(x, dtype=float)
            k = min(n_top, arr.shape[1])
            if k <= 0:
                top_counts = np.zeros(arr.shape[0], dtype=float)
            else:
                part = np.partition(arr, arr.shape[1] - k, axis=1)[:, arr.shape[1] - k:]
                top_counts = part.sum(axis=1)
        return top_counts / denom * 100.0

    @staticmethod
    def _filter_genes_fallback(adata, min_cells: int):
        x = QCModule._materialize_matrix(adata.X)
        if sparse.issparse(x):
            keep = np.asarray((x > 0).sum(axis=0)).ravel() >= int(min_cells)
        else:
            keep = (np.asarray(x) > 0).sum(axis=0) >= int(min_cells)
        if keep.ndim != 1:
            keep = np.asarray(keep).ravel()
        return adata[:, keep].copy()
