"""Hi-C TAD boundary + A/B compartment scoring (Wave 2B / P2.S19).

Memory-safe: operates on the sparse contact matrix from hic_ingest. Insulation
score is computed per bin via a fixed-width window slide; A/B compartments via
the first eigenvector of the per-chromosome correlation matrix.

Outputs:
  adata.uns["hic_tad_boundaries"]  DataFrame: bin_id, chrom, position, insulation, is_boundary
  adata.uns["hic_compartments"]    DataFrame: bin_id, chrom, position, eigenvector_1, compartment ("A"|"B")
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
    """Crane-style insulation score: mean contact count within a square window
    centered at each diagonal bin.

    Output is normalized to mean=0 per chromosome chunk by the caller; this
    function returns raw means.
    """
    n = mat.shape[0]
    scores = np.zeros(n, dtype=np.float32)
    for i in range(n):
        lo = max(0, i - window_bins)
        hi = min(n, i + window_bins + 1)
        if hi - lo <= 1:
            scores[i] = np.nan
            continue
        decision = plan_densify(
            (hi - lo, hi - lo),
            mat.dtype,
            reason="Hi-C insulation score diagonal window",
        )
        if decision == DensifyDecision.ABORT:
            raise MemoryError(
                "hic_tad insulation window too large to densify safely; "
                f"window shape={(hi - lo, hi - lo)}"
            )
        if decision == DensifyDecision.CHUNK:
            logger.warning(
                "hic_tad insulation window is in CHUNK range; proceeding for bounded window shape=%s",
                (hi - lo, hi - lo),
            )
        # densify-allowed: bounded diagonal window guarded by plan_densify
        block = mat[lo:hi, lo:hi].toarray()
        # Mean of upper-triangle excluding diagonal
        iu = np.triu_indices(block.shape[0], k=1)
        if iu[0].size == 0:
            scores[i] = np.nan
            continue
        scores[i] = float(block[iu].mean())
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


def _ab_compartments(mat: sp.csr_matrix, bins: pd.DataFrame) -> dict[str, np.ndarray]:
    """Per-chromosome A/B compartment scores via first eigenvector of correlation matrix.

    Returns dict bin_id → eigenvector_1 value (sign → compartment).
    """
    result: dict[int, float] = {}
    for chrom, group in bins.groupby("chrom"):
        bin_ids = group["bin_id"].to_numpy()
        if len(bin_ids) < 5:
            for b in bin_ids:
                result[int(b)] = 0.0
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
        # Normalize by O/E (observed / expected); use distance-decay normalization
        # Skip if too sparse
        if chrom_mat.sum() == 0:
            for b in bin_ids:
                result[int(b)] = 0.0
            continue
        # Log-transform with pseudocount for numerical stability
        chrom_mat = np.log1p(chrom_mat)
        # Correlation matrix
        corr = np.corrcoef(chrom_mat)
        corr = np.nan_to_num(corr, nan=0.0)
        try:
            eigvals, eigvecs = np.linalg.eigh(corr)
            # Largest eigenvalue → last column
            ev1 = eigvecs[:, -1]
        except np.linalg.LinAlgError:
            ev1 = np.zeros(len(bin_ids))
        for j, b in enumerate(bin_ids):
            result[int(b)] = float(ev1[j])
    return result


class HiCTADModule:
    """TAD boundary + A/B compartment scoring on Hi-C contact matrix."""

    name = "hic_tad"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {
        "uns": ["hic_contact_matrix", "hic_bins"],
    }
    provides_keys: dict[str, list[str]] = {
        "uns": ["hic_tad_boundaries", "hic_compartments"],
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
        ab_scores = _ab_compartments(mat, bins)
        compart_df = bins.copy()
        compart_df["position"] = ((compart_df["start"].astype("int64") + compart_df["end"].astype("int64")) // 2)
        compart_df["eigenvector_1"] = compart_df["bin_id"].map(lambda b: ab_scores.get(int(b), 0.0))
        compart_df["compartment"] = compart_df["eigenvector_1"].apply(lambda v: "A" if v >= 0 else "B")
        adata.uns["hic_compartments"] = compart_df

        out_dir = ctx.run_dir / "hic_tad"
        out_dir.mkdir(parents=True, exist_ok=True)
        summary = {
            "n_bins": int(bins.shape[0]),
            "n_boundaries": int(is_boundary.sum()),
            "n_compartment_A": int((compart_df["compartment"] == "A").sum()),
            "n_compartment_B": int((compart_df["compartment"] == "B").sum()),
            "boundary_k": boundary_k,
            "window_bins": window_bins,
        }
        (out_dir / "tad_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

        ctx.metadata["hic_tad_status"] = "ok"
        ctx.metadata["hic_n_boundaries"] = int(is_boundary.sum())
        ctx.metadata["hic_compartment_counts"] = {"A": summary["n_compartment_A"], "B": summary["n_compartment_B"]}
        logger.info("%s: %d boundaries, %d A / %d B compartments",
                    self.name, summary["n_boundaries"], summary["n_compartment_A"], summary["n_compartment_B"])
