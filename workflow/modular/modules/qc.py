from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from ._scanpy_compat import has_api, import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext


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

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("QC requires loaded AnnData.")

        gene_upper = adata.var_names.astype(str).str.upper()
        adata.var["mt"] = gene_upper.str.startswith("MT-")
        adata.var["ribo"] = gene_upper.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = gene_upper.str.startswith(("HBA", "HBB"))
        use_scanpy_qc = has_api(sc, "pp.calculate_qc_metrics") and ctx.cfg.scale_mode != "massive"
        if use_scanpy_qc:
            try:
                sc.pp.calculate_qc_metrics(
                    adata,
                    qc_vars=["mt", "ribo", "hb"],
                    inplace=True,
                    log1p=False,
                    percent_top=None,
                )
            except Exception:
                self._calculate_qc_metrics_fallback(adata)
        else:
            self._calculate_qc_metrics_fallback(adata)

        # Drop verbose columns to save memory
        self._ensure_axis_tables_support_drop(adata)
        obs_drop = [c for c in adata.obs.columns if c.startswith("log1p") or c.startswith("total_counts_") or c.startswith("pct_counts_in_top_")]
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
        # Violin plots
        fig, axes = plt.subplots(1, 4, figsize=(16, 4))
        for ax, key, title in zip(
            axes,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
            ["Genes per cell", "UMI counts per cell", "Mito %", "Ribo %"],
        ):
            sc.pl.violin(adata, key, ax=ax, show=False, stripplot=False)
            ax.set_title(title)
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / f"qc_violin_{suffix}.png", dpi=160, bbox_inches="tight")
        plt.close()

        # Scatter plots: genes vs counts, colored by mito%
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.scatter(
            adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt",
            ax=axes[0], show=False, title="Genes vs UMI (color=mito%)",
        )
        sc.pl.scatter(
            adata, x="total_counts", y="pct_counts_mt",
            ax=axes[1], show=False, title="Mito% vs UMI counts",
        )
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / f"qc_scatter_{suffix}.png", dpi=160, bbox_inches="tight")
        plt.close()


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
