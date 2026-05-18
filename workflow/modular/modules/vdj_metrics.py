"""VDJ metrics (Wave 2B / P2.S16): Shannon, Gini, clonal expansion.

Consumes adata.obs["clonotype_id"] from vdj_ingest. Writes per-sample diversity
metrics and per-cell expansion class. Memory-safe (no AnnData duplication).

Outputs:
  adata.obs["clonal_expansion"]            string class: "singleton" | "small_2-5" | "medium_6-20" | "large_>20"
  adata.uns["vdj_diversity"]               per-sample DataFrame with Shannon + Gini + n_clonotypes + n_cells
  runs/<run-id>/vdj_metrics/vdj_diversity.csv
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "Shannon_1948": {
        "title": "A Mathematical Theory of Communication",
        "authors": "Shannon",
        "journal": "Bell System Technical Journal",
        "year": "1948",
        "doi": "10.1002/j.1538-7305.1948.tb01338.x",
        "description": "Shannon entropy used here for clonotype diversity per sample.",
    },
    "Gini_1912": {
        "title": "Variabilita e mutabilita",
        "authors": "Gini",
        "journal": "Studi Economico-Giuridici della Universita de Cagliari",
        "year": "1912",
        "doi": "https://en.wikipedia.org/wiki/Gini_coefficient",
        "description": "Gini coefficient used to quantify clonotype abundance inequality (1 - sum_i p_i (2 i - n - 1)/n).",
    },
    "Bagaev_VDJtools_2015": {
        "title": "VDJtools: Unifying Post-analysis of T Cell Receptor Repertoires",
        "authors": "Bagaev et al.",
        "journal": "PLoS Computational Biology",
        "year": "2015",
        "doi": "10.1371/journal.pcbi.1004503",
        "description": "Canonical repertoire diversity conventions adopted here (Shannon / Gini per sample; clonal expansion buckets).",
    },
}


_EXPANSION_BUCKETS = [
    (1, 1, "singleton"),
    (2, 5, "small_2-5"),
    (6, 20, "medium_6-20"),
    (21, float("inf"), "large_>20"),
]


def _expansion_class(n_cells: int) -> str:
    for lo, hi, label in _EXPANSION_BUCKETS:
        if lo <= n_cells <= hi:
            return label
    return "unknown"


def _shannon(counts: np.ndarray) -> float:
    """Shannon entropy H = -sum p_i log p_i (natural log)."""
    if counts.size == 0:
        return 0.0
    total = float(counts.sum())
    if total <= 0:
        return 0.0
    p = counts / total
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def _gini(counts: np.ndarray) -> float:
    """Gini coefficient over clonotype counts (0=uniform, 1=monopoly)."""
    if counts.size == 0:
        return 0.0
    counts = np.sort(counts.astype(np.float64))
    n = counts.size
    total = counts.sum()
    if total <= 0:
        return 0.0
    cumulative = (counts * np.arange(1, n + 1)).sum()
    return float((2.0 * cumulative) / (n * total) - (n + 1.0) / n)


class VDJMetricsModule:
    """Per-sample VDJ diversity + per-cell clonal expansion class."""

    name = "vdj_metrics"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {
        "obs": ["clonotype_id"],
    }
    provides_keys: dict[str, list[str]] = {
        "obs": ["clonal_expansion"],
        "uns": ["vdj_diversity"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata
        if "clonotype_id" not in adata.obs.columns:
            ctx.status(self.name, "skipped", "clonotype_id column absent — run vdj_ingest first")
            ctx.metadata["vdj_metrics_status"] = "skipped_no_clonotypes"
            return

        # Per-cell clonal expansion class (global, across all samples)
        clonotype_counts = adata.obs["clonotype_id"].value_counts()
        # Cells with empty clonotype_id (no VDJ) get "no_vdj"
        def _classify(cid: str) -> str:
            if cid == "" or pd.isna(cid):
                return "no_vdj"
            return _expansion_class(int(clonotype_counts.get(cid, 0)))
        adata.obs["clonal_expansion"] = adata.obs["clonotype_id"].map(_classify).astype(str)

        # Per-sample diversity. Use "sample" column if present, otherwise treat all cells as one sample.
        sample_col = "sample" if "sample" in adata.obs.columns else None
        rows: list[dict] = []
        if sample_col is None:
            # Aggregate one row
            cell_mask = adata.obs["clonotype_id"] != ""
            counts = adata.obs.loc[cell_mask, "clonotype_id"].value_counts().to_numpy()
            rows.append({
                "sample": "_ALL_",
                "n_cells_with_vdj": int(cell_mask.sum()),
                "n_clonotypes": int(counts.size),
                "shannon": _shannon(counts),
                "gini": _gini(counts),
            })
        else:
            for sample_id, group in adata.obs.groupby(sample_col, observed=True):
                mask = group["clonotype_id"] != ""
                counts = group.loc[mask, "clonotype_id"].value_counts().to_numpy()
                rows.append({
                    "sample": str(sample_id),
                    "n_cells_with_vdj": int(mask.sum()),
                    "n_clonotypes": int(counts.size),
                    "shannon": _shannon(counts),
                    "gini": _gini(counts),
                })

        diversity_df = pd.DataFrame(rows)
        adata.uns["vdj_diversity"] = diversity_df

        out_dir = ctx.run_dir / "vdj_metrics"
        out_dir.mkdir(parents=True, exist_ok=True)
        diversity_df.to_csv(out_dir / "vdj_diversity.csv", index=False)
        (out_dir / "vdj_diversity.summary.json").write_text(
            json.dumps({
                "n_samples": int(diversity_df.shape[0]),
                "shannon_mean": float(diversity_df["shannon"].mean()) if not diversity_df.empty else 0.0,
                "gini_mean": float(diversity_df["gini"].mean()) if not diversity_df.empty else 0.0,
                "expansion_counts": adata.obs["clonal_expansion"].value_counts().to_dict(),
            }, indent=2),
            encoding="utf-8",
        )

        ctx.metadata["vdj_metrics_status"] = "ok"
        ctx.metadata["vdj_n_samples_with_diversity"] = int(diversity_df.shape[0])
        logger.info(
            "%s: per-sample diversity computed for %d samples; expansion buckets: %s",
            self.name, ctx.metadata["vdj_n_samples_with_diversity"],
            adata.obs["clonal_expansion"].value_counts().to_dict(),
        )
