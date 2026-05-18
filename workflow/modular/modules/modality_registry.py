"""Modality registry: detect and inventory present omics modalities (Wave 2 / P1B.S8).

Idempotent and safe to re-run. Writes:
  adata.uns["modalities_present"]  -- list[str]
  runs/<run-id>/manifest/modalities.json
"""
from __future__ import annotations

import json
import logging

from ..context import PipelineContext

logger = logging.getLogger(__name__)

__references__: dict = {
    "project_local_utility": {
        "title": "Project-local utility — multi-omics modality inventory",
        "authors": "singlecell_factory contributors",
        "journal": "Internal documentation",
        "year": "2025",
        "doi": "PMID: 0",
        "description": "Detects present modalities by scanning obsm/obs keys; no external paper. Tier-3 project-local utility per docs/SCIENTIFIC_AUDIT_2026-05-15.md.",
    },
}

# Keys in obsm that indicate a non-RNA modality is present
OBSM_MODALITY_MAP = {
    "protein_adt": "protein_adt",
    "protein_clr": "protein_adt",
    "X_atac": "atac",
    "atac_peaks": "atac",
    "X_vdj": "vdj",
    "vdj_clonotypes": "vdj",
    "ribo_te": "ribo",
    "X_ribo": "ribo",
    "hic_compartment": "hic",
    "X_hic": "hic",
    "spatial": "spatial_transcriptomics",
    "X_spatial": "spatial_transcriptomics",
}

# Obs column flags that indicate a modality
OBS_FLAG_MAP = {
    "spatial_x": "spatial_transcriptomics",
    "spatial_y": "spatial_transcriptomics",
    "atac_qc_status": "atac",
    "clonotype_id": "vdj",
}


class ModalityRegistryModule:
    """Detect and inventory present modalities; idempotent and safe to run multiple times."""

    name = "modality_registry"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {"uns": ["modalities_present"]}

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata
        present: list[str] = ["singlecell_rna"]  # RNA always present

        for obsm_key, modality in OBSM_MODALITY_MAP.items():
            if obsm_key in adata.obsm and modality not in present:
                present.append(modality)

        for obs_key, modality in OBS_FLAG_MAP.items():
            if obs_key in adata.obs.columns and modality not in present:
                present.append(modality)

        adata.uns["modalities_present"] = present

        manifest_dir = ctx.run_dir / "manifest"
        manifest_dir.mkdir(parents=True, exist_ok=True)
        (manifest_dir / "modalities.json").write_text(
            json.dumps(
                {"modalities": present, "run_id": ctx.run_dir.name, "count": len(present)},
                indent=2,
            ),
            encoding="utf-8",
        )

        ctx.metadata["modalities_present"] = present
        ctx.metadata["modality_count"] = len(present)
        logger.info("%s: detected modalities: %s", self.name, present)
