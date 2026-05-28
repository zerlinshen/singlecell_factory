"""Hi-C TAD boundary + A/B compartment scoring (Wave 2B / P2.S19).

Memory-safe: operates on the sparse contact matrix from hic_ingest. Insulation
score is computed per bin via a fixed-width window slide; A/B compartments via
the first eigenvector of the per-chromosome correlation matrix.

Scope of claims (validated 2026-05-28 against retained B35T1NC Micro-C author
TADs and GM12878 published subcompartments):
  - A/B compartments here are SUPPORTED on unbiased Hi-C (GM12878 chr19
    concordance 0.814 vs published subcompartments).
  - The insulation TAD boundaries here are VISUAL-QC / quick-look ONLY. This
    lightweight caller scores below a permuted null on the matched B35T1NC
    chr19 author TADs (best F1=0.049, recall_over_null=0.40), and observed/
    expected normalization does not lift it above null. Authoritative TAD /
    insulation biology routes to the cooltools-backed r_multiomics_factory/bulk
    lane (genome-wide F1=0.516, recall_over_null=6.84; chr19 matched-GT
    F1=0.25, recall_over_null=3.0). Do not present these boundaries as final
    TAD biology.

Outputs:
  adata.uns["hic_tad_boundaries"]  DataFrame: bin_id, chrom, position, insulation, is_boundary
  adata.uns["hic_compartments"]    DataFrame: bin_id, chrom, position, eigenvector_1, compartment ("A"|"B"|"low_information")
  adata.uns["hic_tad_metadata"]    Dict: per-chromosome A/B compartment claim-safety status
  runs/<run-id>/hic_tad/tad_summary.json
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..context import PipelineContext
from .._densify_policy import DensifyDecision, plan_densify

logger = logging.getLogger(__name__)


__references__ = {
    "Crane_insulation_2015": {
        "title": "Condensin-driven remodelling of X chromosome topology during dosage compensation",
        "authors": "Crane et al.",
        "journal": "Nature",
        "year": "2015",
        "doi": "10.1038/nature14450",
        "description": "Insulation-score TAD boundary detection method implemented here.",
    },
    "Lieberman_Aiden_HiC_2009": {
        "title": "Comprehensive Mapping of Long-Range Interactions Reveals Folding Principles of the Human Genome",
        "authors": "Lieberman-Aiden et al.",
        "journal": "Science",
        "year": "2009",
        "doi": "10.1126/science.1181369",
        "description": "Original A/B compartment definition via first eigenvector of the correlation matrix.",
    },
    "Nora_TADs_2012": {
        "title": "Spatial partitioning of the regulatory landscape of the X-inactivation centre",
        "authors": "Nora et al.",
        "journal": "Nature",
        "year": "2012",
        "doi": "10.1038/nature11049",
        "description": "Original TAD identification methodology.",
    },
}


def _insulation_score(mat: sp.csr_matrix, window_bins: int = 5) -> np.ndarray:
    """Crane-style insulation score: mean contacts crossing each candidate boundary.

    For bin ``i``, score the upstream window ``[i-window, i)`` against the
    downstream window ``[i, i+window)``. TAD boundaries should have low
    cross-window contact frequency. Empty windows are marked NaN so sparse
    unobserved chromosome tails are not called as artificial boundaries.
    """
    n = mat.shape[0]
    scores = np.zeros(n, dtype=np.float32)
    for i in range(n):
        left_lo = max(0, i - window_bins)
        left_hi = i
        right_lo = i
        right_hi = min(n, i + window_bins)
        if left_hi - left_lo < 1 or right_hi - right_lo < 1:
            scores[i] = np.nan
            continue
        decision = plan_densify(
            (left_hi - left_lo, right_hi - right_lo),
            mat.dtype,
            reason="Hi-C insulation score cross-boundary window",
        )
        if decision == DensifyDecision.ABORT:
            raise MemoryError(
                "hic_tad insulation window too large to densify safely; "
                f"window shape={(left_hi - left_lo, right_hi - right_lo)}"
            )
        if decision == DensifyDecision.CHUNK:
            logger.warning(
                "hic_tad insulation window is in CHUNK range; proceeding for bounded window shape=%s",
                (left_hi - left_lo, right_hi - right_lo),
            )
        # densify-allowed: bounded cross-boundary window guarded by plan_densify
        block = mat[left_lo:left_hi, right_lo:right_hi].toarray()
        # An unobserved window is NaN (skipped by _detect_boundaries), not a
        # boundary: this drops a true zero-cross-contact boundary but refuses to
        # fabricate boundaries on sparse unobserved tails (no silent compromise).
        if block.sum() <= 1e-8:
            scores[i] = np.nan
            continue
        scores[i] = float(block.mean())
    return scores


def _detect_boundaries(insulation: np.ndarray, k: float = 1.0) -> np.ndarray:
    """Local minima of insulation < (mean - k*std) are boundaries (boolean array)."""
    valid = ~np.isnan(insulation)
    if not valid.any():
        return np.zeros_like(insulation, dtype=bool)
    mean = float(np.nanmean(insulation))
    std = float(np.nanstd(insulation))
    threshold = mean - k * std
    boundaries = np.zeros_like(insulation, dtype=bool)
    for i in range(1, len(insulation) - 1):
        if not valid[i] or not valid[i - 1] or not valid[i + 1]:
            continue
        if insulation[i] < threshold and insulation[i] <= insulation[i - 1] and insulation[i] <= insulation[i + 1]:
            boundaries[i] = True
    return boundaries


def _zero_scores(bin_ids: np.ndarray, result: dict[int, float]) -> None:
    for b in bin_ids:
        result[int(b)] = 0.0


def _mark_low_information(
    chrom_name: str,
    bin_ids: np.ndarray,
    result: dict[int, float],
    status_by_chrom: dict[str, str],
    low_information_chromosomes: list[str],
) -> None:
    _zero_scores(bin_ids, result)
    status_by_chrom[chrom_name] = "low_information"
    low_information_chromosomes.append(chrom_name)


def _ab_compartments(mat: sp.csr_matrix, bins: pd.DataFrame) -> tuple[dict[int, float], dict[str, object]]:
    """Per-chromosome A/B compartment scores via first eigenvector of correlation matrix.

    Returns bin_id → eigenvector_1 value plus claim-safety metadata.
    """
    result: dict[int, float] = {}
    status_by_chrom: dict[str, str] = {}
    low_information_chromosomes: list[str] = []
    eps = 1e-8
    for chrom, group in bins.groupby("chrom"):
        chrom_name = str(chrom)
        bin_ids = group["bin_id"].to_numpy()
        if len(bin_ids) < 5:
            _mark_low_information(
                chrom_name, bin_ids, result, status_by_chrom, low_information_chromosomes
            )
            continue
        chrom_shape = (len(bin_ids), len(bin_ids))
        decision = plan_densify(
            chrom_shape,
            mat.dtype,
            reason="Hi-C chromosome contact submatrix for A/B compartment eigenvector",
        )
        if decision != DensifyDecision.GO:
            raise MemoryError(
                "hic_tad A/B compartment matrix too large to densify safely; "
                f"chrom={chrom!r} shape={chrom_shape} decision={decision}"
            )
        # densify-allowed: chromosome submatrix guarded by plan_densify before eigendecomposition
        chrom_mat = mat[bin_ids[:, None], bin_ids[None, :]].toarray()
        if chrom_mat.sum() <= eps:
            _mark_low_information(
                chrom_name, bin_ids, result, status_by_chrom, low_information_chromosomes
            )
            continue
        # Log-transform with pseudocount for numerical stability
        chrom_mat = np.log1p(chrom_mat)
        row_vars = chrom_mat.var(axis=1)
        variable_rows = row_vars > eps
        if int(variable_rows.sum()) < 2:
            _mark_low_information(
                chrom_name, bin_ids, result, status_by_chrom, low_information_chromosomes
            )
            continue
        variable_bin_ids = bin_ids[variable_rows]
        nonvariable_bin_ids = bin_ids[~variable_rows]
        for b in nonvariable_bin_ids:
            result[int(b)] = float("nan")
        chrom_mat = chrom_mat[variable_rows][:, variable_rows]
        # Build row correlations directly after screening constant rows; this
        # avoids np.corrcoef runtime warnings on low-information chromosomes.
        centered = chrom_mat - chrom_mat.mean(axis=1, keepdims=True)
        denom = np.sqrt(np.sum(centered * centered, axis=1, keepdims=True))
        denom = np.maximum(denom, eps)
        normalized = centered / denom
        corr = normalized @ normalized.T
        corr = np.clip(corr, -1.0, 1.0)
        try:
            _eigvals, eigvecs = np.linalg.eigh(corr)
            # Largest eigenvalue → last column
            ev1 = eigvecs[:, -1]
        except np.linalg.LinAlgError:
            _mark_low_information(
                chrom_name, bin_ids, result, status_by_chrom, low_information_chromosomes
            )
            continue
        for j, b in enumerate(variable_bin_ids):
            result[int(b)] = float(ev1[j])
        status_by_chrom[chrom_name] = "confident"

    if not status_by_chrom or len(low_information_chromosomes) == len(status_by_chrom):
        compartment_status = "low_information"
    elif low_information_chromosomes:
        compartment_status = "partial_low_information"
    else:
        compartment_status = "confident"
    metadata = {
        "compartment_status": compartment_status,
        "low_information_chromosomes": sorted(low_information_chromosomes),
        "compartment_status_by_chrom": status_by_chrom,
    }
    return result, metadata


class HiCTADModule:
    """TAD boundary + A/B compartment scoring on Hi-C contact matrix."""

    name = "hic_tad"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {
        "uns": ["hic_contact_matrix", "hic_bins"],
    }
    provides_keys: dict[str, list[str]] = {
        "uns": ["hic_tad_boundaries", "hic_compartments", "hic_tad_metadata"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata
        mat = adata.uns.get("hic_contact_matrix")
        bins = adata.uns.get("hic_bins")
        if mat is None or bins is None:
            ctx.status(self.name, "skipped", "hic_contact_matrix or hic_bins missing — run hic_ingest first")
            ctx.metadata["hic_tad_status"] = "skipped_no_input"
            return

        window_bins = int(getattr(ctx.cfg, "hic_tad_window_bins", 5))
        boundary_k = float(getattr(ctx.cfg, "hic_tad_boundary_k", 1.0))

        logger.info("%s: computing insulation (window=%d bins)", self.name, window_bins)
        insulation = _insulation_score(mat, window_bins=window_bins)
        is_boundary = _detect_boundaries(insulation, k=boundary_k)

        boundaries_df = bins.copy()
        boundaries_df["position"] = ((boundaries_df["start"].astype("int64") + boundaries_df["end"].astype("int64")) // 2)
        boundaries_df["insulation"] = insulation
        boundaries_df["is_boundary"] = is_boundary
        adata.uns["hic_tad_boundaries"] = boundaries_df

        logger.info("%s: computing A/B compartments via eigendecomposition", self.name)
        ab_scores, hic_tad_metadata = _ab_compartments(mat, bins)
        compart_df = bins.copy()
        compart_df["position"] = ((compart_df["start"].astype("int64") + compart_df["end"].astype("int64")) // 2)
        compart_df["eigenvector_1"] = compart_df["bin_id"].map(lambda b: ab_scores.get(int(b), 0.0))
        status_by_chrom = hic_tad_metadata["compartment_status_by_chrom"]
        compart_df["compartment"] = compart_df.apply(
            lambda row: (
                "low_information"
                if status_by_chrom.get(str(row["chrom"])) == "low_information"
                or pd.isna(row["eigenvector_1"])
                else ("A" if row["eigenvector_1"] >= 0 else "B")
            ),
            axis=1,
        )
        adata.uns["hic_compartments"] = compart_df
        adata.uns["hic_tad_metadata"] = hic_tad_metadata

        out_dir = ctx.run_dir / "hic_tad"
        out_dir.mkdir(parents=True, exist_ok=True)
        summary = {
            "n_bins": int(bins.shape[0]),
            "n_boundaries": int(is_boundary.sum()),
            "n_compartment_A": int((compart_df["compartment"] == "A").sum()),
            "n_compartment_B": int((compart_df["compartment"] == "B").sum()),
            "boundary_k": boundary_k,
            "window_bins": window_bins,
            "hic_tad_metadata": hic_tad_metadata,
        }
        (out_dir / "tad_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

        ctx.metadata["hic_tad_status"] = "ok"
        ctx.metadata["hic_n_boundaries"] = int(is_boundary.sum())
        ctx.metadata["hic_compartment_counts"] = {"A": summary["n_compartment_A"], "B": summary["n_compartment_B"]}
        ctx.metadata["hic_tad_metadata"] = hic_tad_metadata
        logger.info("%s: %d boundaries, %d A / %d B compartments",
                    self.name, summary["n_boundaries"], summary["n_compartment_A"], summary["n_compartment_B"])
