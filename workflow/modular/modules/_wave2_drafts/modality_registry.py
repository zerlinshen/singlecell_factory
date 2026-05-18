"""DRAFT — Wave 2 / Phase 1B P1B.S8: modality_registry module.

DO NOT import or register until Wave 2 branch. Detection sources:
obsm keys, obs flags, manifest.extensions slots.

Writes:
  adata.uns["modalities_present"] (list[str])
  runs/<run-id>/manifest/modalities.json
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

from ...context import PipelineContext

logger = logging.getLogger(__name__)

__references__: dict = {}

# Keys in obsm that indicate a non-RNA modality is present
_OBSM_MODALITY_MAP = {
    "protein_adt": "protein_adt",
    "X_atac": "atac",
    "X_vdj": "vdj",
    "ribo_te": "ribo",
    "hic_compartment": "hic",
    "spatial": "spatial_transcriptomics",
}

# Obs column flags that indicate a modality
_OBS_FLAG_MAP = {
    "spatial_x": "spatial_transcriptomics",
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

        for obsm_key, modality in _OBSM_MODALITY_MAP.items():
            if obsm_key in adata.obsm and modality not in present:
                present.append(modality)

        for obs_key, modality in _OBS_FLAG_MAP.items():
            if obs_key in adata.obs.columns and modality not in present:
                present.append(modality)

        adata.uns["modalities_present"] = present

        manifest_dir = ctx.run_dir / "manifest"
        manifest_dir.mkdir(parents=True, exist_ok=True)
        (manifest_dir / "modalities.json").write_text(
            json.dumps({"modalities": present, "run_id": ctx.run_dir.name}, indent=2),
            encoding="utf-8",
        )

        ctx.metadata["modalities_present"] = present
        ctx.metadata["modality_count"] = len(present)
        logger.info("%s: detected modalities: %s", self.name, present)
