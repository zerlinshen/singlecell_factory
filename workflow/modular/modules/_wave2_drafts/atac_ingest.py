"""DRAFT — Wave 2 / Phase 2 P2.S11: atac_ingest module.

DO NOT import or register until Wave 2 branch.
Requires Signac (R-side) and/or SnapATAC2 (Python-side) in environment.

Inputs: fragment file or peak matrix
Outputs:
  adata.obsm["X_atac"] — LSI embedding
  adata.uns["atac_peaks"] — peak metadata DataFrame
  bundle extensions/atac/ — for R marker_db_module
"""
from __future__ import annotations

import logging

from ...context import PipelineContext

logger = logging.getLogger(__name__)

__references__ = {
    "SnapATAC2": {
        "title": "A fast, scalable and versatile tool for analysis of single-cell omics data",
        "authors": "Zhang et al.",
        "journal": "Nature Methods",
        "year": "2024",
        "doi": "10.1038/s41592-023-02139-9",
        "description": "Python ATAC-seq analysis framework",
    },
    "Signac": {
        "title": "Multimodal single-cell chromatin analysis with Signac",
        "authors": "Stuart et al.",
        "journal": "Nature Methods",
        "year": "2021",
        "doi": "10.1038/s41592-021-01282-5",
        "description": "R-side ATAC analysis and visualization",
    },
}


class ATACIngestModule:
    """ATAC-seq ingest: TF-IDF + LSI embedding from fragment file or peak matrix."""

    name = "atac_ingest"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {"obsm": ["X_atac"], "uns": ["atac_peaks"]}

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        atac_path = getattr(ctx.cfg, "atac_fragment_path", None)
        if atac_path is None:
            ctx.status(self.name, "skipped", "atac_fragment_path not set in config")
            ctx.metadata["atac_ingest_status"] = "skipped_no_input"
            return

        try:
            import snapatac2 as snap
        except ImportError:
            ctx.status(self.name, "skipped", "snapatac2 not installed")
            ctx.metadata["atac_ingest_status"] = "skipped_snapatac2_missing"
            return

        logger.info("%s: ingesting ATAC from %s", self.name, atac_path)
        # DRAFT: full implementation in Wave 2
        # snap.pp.make_peak_matrix(...)
        # snap.tl.spectral(...)
        # adata.obsm["X_atac"] = lsi_embedding
        raise NotImplementedError("atac_ingest is a Wave 2 draft — not yet implemented")
