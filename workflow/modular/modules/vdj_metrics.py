"""VDJ metrics (Wave 2B / P2.S16): Shannon, Gini, clonal expansion.

Consumes adata.obs["clonotype_id"] from vdj_ingest. Writes per-sample diversity
metrics and per-cell expansion class. When ``obs.sample`` exists, clonal
expansion is keyed by ``(sample, clonotype_id)``; clonotype-only counting and
the ``_ALL_`` repertoire are used only when the sample column is absent.
Diversity remains per-sample with the established output columns. Memory-safe
(no AnnData duplication).

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
    """Compute per-sample diversity and sample-local clonal expansion classes.

    With ``obs.sample``, each expansion count is keyed by ``(sample,
    clonotype_id)`` so Cell Ranger-local clonotype IDs never aggregate across
    biological samples. Only an absent sample column permits clonotype-only
    counting in the ``_ALL_`` repertoire. Diversity retains the established
    per-sample ``sample``, ``n_cells_with_vdj``, ``n_clonotypes``, ``shannon``,
    and ``gini`` columns.
    """

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

        # Determine expansion unit and compute per-cell clonal expansion class
        sample_present = "sample" in adata.obs.columns
        vdj_mask = (adata.obs["clonotype_id"] != "") & (adata.obs["clonotype_id"].notna())

        if sample_present:
            expansion_unit = "obs.sample"
            counting_keys = "sample + clonotype_id"
            sample_local_clonotypes = True
            fallback_all_repertoire = False

            # Validate that every VDJ-bearing cell has a non-null, non-blank sample identifier
            if vdj_mask.any():
                vdj_samples = adata.obs.loc[vdj_mask, "sample"]
                if vdj_samples.isna().any() or (vdj_samples.astype(str).str.strip() == "").any():
                    raise ValueError("VDJ-bearing cells contain null, NaN, or blank sample identifiers.")

            vdj_obs = adata.obs.loc[vdj_mask, ["sample", "clonotype_id"]].copy()
            vdj_obs["sample_str"] = vdj_obs["sample"].astype(str)
            vdj_obs["clonotype_str"] = vdj_obs["clonotype_id"].astype(str)
            counts_series = vdj_obs.groupby(["sample_str", "clonotype_str"], observed=True).size()

            clonal_exp = pd.Series("no_vdj", index=adata.obs.index, dtype=object)
            if vdj_mask.any():
                keys = list(zip(vdj_obs["sample_str"], vdj_obs["clonotype_str"]))
                mapped_counts = [counts_series.get(k, 0) for k in keys]
                clonal_exp.loc[vdj_mask] = [_expansion_class(int(c)) for c in mapped_counts]
            adata.obs["clonal_expansion"] = clonal_exp.astype(str)

            # Per-sample diversity
            rows: list[dict[str, Any]] = []
            for sample_id, group in adata.obs.groupby("sample", observed=True):
                mask = (group["clonotype_id"] != "") & (group["clonotype_id"].notna())
                counts = group.loc[mask, "clonotype_id"].astype(str).value_counts().to_numpy()
                rows.append({
                    "sample": str(sample_id),
                    "n_cells_with_vdj": int(mask.sum()),
                    "n_clonotypes": int(counts.size),
                    "shannon": _shannon(counts),
                    "gini": _gini(counts),
                })
        else:
            expansion_unit = "_ALL_"
            counting_keys = "clonotype_id"
            sample_local_clonotypes = False
            fallback_all_repertoire = True

            counts_series = adata.obs.loc[vdj_mask, "clonotype_id"].astype(str).value_counts()
            clonal_exp = pd.Series("no_vdj", index=adata.obs.index, dtype=object)
            if vdj_mask.any():
                cids = adata.obs.loc[vdj_mask, "clonotype_id"].astype(str)
                clonal_exp.loc[vdj_mask] = [_expansion_class(int(counts_series.get(cid, 0))) for cid in cids]
            adata.obs["clonal_expansion"] = clonal_exp.astype(str)

            counts = counts_series.to_numpy()
            rows = [{
                "sample": "_ALL_",
                "n_cells_with_vdj": int(vdj_mask.sum()),
                "n_clonotypes": int(counts.size),
                "shannon": _shannon(counts),
                "gini": _gini(counts),
            }]

        diversity_df = pd.DataFrame(
            rows,
            columns=["sample", "n_cells_with_vdj", "n_clonotypes", "shannon", "gini"],
        )
        adata.uns["vdj_diversity"] = diversity_df

        vdj_provenance = {
            "expansion_unit": expansion_unit,
            "counting_keys": counting_keys,
            "sample_local_clonotypes": sample_local_clonotypes,
            "fallback_all_repertoire": fallback_all_repertoire,
            "sample_column_present": sample_present,
        }

        out_dir = ctx.run_dir / "vdj_metrics"
        out_dir.mkdir(parents=True, exist_ok=True)
        diversity_df.to_csv(out_dir / "vdj_diversity.csv", index=False)
        (out_dir / "vdj_diversity.summary.json").write_text(
            json.dumps({
                "n_samples": int(diversity_df.shape[0]),
                "shannon_mean": float(diversity_df["shannon"].mean()) if not diversity_df.empty else 0.0,
                "gini_mean": float(diversity_df["gini"].mean()) if not diversity_df.empty else 0.0,
                "expansion_counts": adata.obs["clonal_expansion"].value_counts().to_dict(),
                "provenance": vdj_provenance,
            }, indent=2),
            encoding="utf-8",
        )

        ctx.metadata["vdj_metrics_status"] = "ok"
        ctx.metadata["vdj_n_samples_with_diversity"] = int(diversity_df.shape[0])
        ctx.metadata["vdj_expansion_unit"] = expansion_unit
        ctx.metadata["vdj_counting_keys"] = counting_keys
        ctx.metadata["vdj_sample_local_clonotypes"] = sample_local_clonotypes
        ctx.metadata["vdj_fallback_all_repertoire"] = fallback_all_repertoire
        logger.info(
            "%s: per-sample diversity computed for %d samples (expansion_unit=%s); expansion buckets: %s",
            self.name, ctx.metadata["vdj_n_samples_with_diversity"], expansion_unit,
            adata.obs["clonal_expansion"].value_counts().to_dict(),
        )
