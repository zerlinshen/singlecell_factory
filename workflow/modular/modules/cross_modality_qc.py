"""Cross-modality QC (Wave 2 / P1B.S9): barcode overlap across modalities.

Tri-state status (plan improvement #6):
  skipped_single_modality | ran | warn_divergence

Always writes runs/<run-id>/qc/cross_modality_status.json.
Also writes cross_modality_overlap.json when at least two modalities present.

Sparse-safe: uses obsm slice norms without densifying full X. Mirrors the
Wave 1 memory-safety lesson.
"""
from __future__ import annotations

import json
import logging

import numpy as np

from ..context import PipelineContext
from .modality_registry import OBSM_MODALITY_MAP

logger = logging.getLogger(__name__)

__references__: dict = {
    "project_local_utility": {
        "title": "Project-local utility — cross-modality barcode-overlap QC",
        "authors": "singlecell_factory contributors",
        "journal": "Internal documentation",
        "year": "2025",
        "doi": "PMID: 0",
        "description": "Jaccard-overlap QC across declared modalities; no external paper. Tier-3 project-local utility per docs/SCIENTIFIC_AUDIT_2026-05-15.md.",
    },
}

_DIVERGENCE_THRESHOLD = 0.05  # Jaccard < this -> warn_divergence


class CrossModalityQCModule:
    """Compute barcode overlap across modalities. Silent-skip when only RNA present."""

    name = "cross_modality_qc"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {"uns": ["modalities_present"]}
    provides_keys: dict[str, list[str]] = {}

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata
        modalities = list(adata.uns.get("modalities_present", ["singlecell_rna"]))
        qc_dir = ctx.run_dir / "qc"
        qc_dir.mkdir(parents=True, exist_ok=True)

        if len(modalities) < 2:
            status_doc = {
                "status": "skipped_single_modality",
                "modality_count": len(modalities),
                "modalities": modalities,
            }
            (qc_dir / "cross_modality_status.json").write_text(
                json.dumps(status_doc, indent=2), encoding="utf-8"
            )
            ctx.metadata["cross_modality_qc_status"] = "skipped_single_modality"
            logger.info("%s: single modality only — skip.", self.name)
            return

        # Build modality -> set-of-barcodes mapping by reading obsm norms.
        # No full-X densification (Wave 1 sparse-safety lesson).
        rna_barcodes = set(adata.obs_names)
        overlaps: dict[str, dict] = {}
        has_divergence = False

        # Track which modalities we've already counted (avoid double-counting from aliased keys)
        seen_modalities: set[str] = set()

        for obsm_key, modality in OBSM_MODALITY_MAP.items():
            if modality not in modalities or obsm_key not in adata.obsm or modality in seen_modalities:
                continue
            seen_modalities.add(modality)

            mat = np.asarray(adata.obsm[obsm_key])
            if mat.ndim == 1:
                present_mask = mat != 0
            else:
                present_mask = np.linalg.norm(mat, axis=1) > 0
            mod_barcodes = set(adata.obs_names[present_mask])

            intersection = len(rna_barcodes & mod_barcodes)
            union = len(rna_barcodes | mod_barcodes)
            jaccard = intersection / union if union > 0 else 0.0
            overlaps[modality] = {
                "rna_n": len(rna_barcodes),
                "modality_n": len(mod_barcodes),
                "intersection": intersection,
                "jaccard": round(jaccard, 4),
            }
            if jaccard < _DIVERGENCE_THRESHOLD:
                has_divergence = True
                logger.warning(
                    "%s: barcode overlap RNA vs %s is low (Jaccard=%.3f < %.2f)",
                    self.name, modality, jaccard, _DIVERGENCE_THRESHOLD,
                )

        final_status = "warn_divergence" if has_divergence else "ran"
        status_doc = {
            "status": final_status,
            "modality_count": len(modalities),
            "modalities": modalities,
        }
        (qc_dir / "cross_modality_status.json").write_text(
            json.dumps(status_doc, indent=2), encoding="utf-8"
        )
        (qc_dir / "cross_modality_overlap.json").write_text(
            json.dumps(overlaps, indent=2), encoding="utf-8"
        )
        ctx.metadata["cross_modality_qc_status"] = final_status
        ctx.metadata["cross_modality_overlaps"] = overlaps
        logger.info("%s: status=%s modalities=%s", self.name, final_status, modalities)
