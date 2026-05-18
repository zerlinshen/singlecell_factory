"""DRAFT — Wave 2 / Phase 1B P1B.S9: cross_modality_qc module.

DO NOT import or register until Wave 2 branch.

Tri-state status output per plan improvement #6:
  skipped_single_modality | ran | warn_divergence

Writes:
  runs/<run-id>/qc/cross_modality_status.json (always)
  runs/<run-id>/qc/cross_modality_overlap.json (when ran or warn_divergence)
"""
from __future__ import annotations

import json
import logging

import numpy as np

from ...context import PipelineContext

logger = logging.getLogger(__name__)

__references__: dict = {}

_DIVERGENCE_THRESHOLD = 0.05  # Jaccard < this = warn_divergence


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
            logger.info("%s: single modality — skip.", self.name)
            return

        # Compare barcode sets across modalities using obsm key indices
        rna_barcodes = set(adata.obs_names)
        overlaps: dict[str, dict] = {}
        has_divergence = False

        from ..._wave2_drafts.modality_registry import _OBSM_MODALITY_MAP
        for obsm_key, modality in _OBSM_MODALITY_MAP.items():
            if modality not in modalities or obsm_key not in adata.obsm:
                continue
            # All cells with non-zero norm in the obsm slot are "present" in this modality
            mat = np.asarray(adata.obsm[obsm_key])
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
