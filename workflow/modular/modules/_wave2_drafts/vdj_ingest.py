"""DRAFT — Wave 2 / Phase 2 P2.S15: vdj_ingest module.

DO NOT import or register until Wave 2 branch.
Requires scirpy for VDJ receptor analysis.

Inputs: 10x cellranger VDJ output directory
Outputs:
  adata.obs["chain_pairing"] — scirpy chain pairing annotation
  adata.obs["clonotype_id"] — clonotype assignments
  adata.uns["vdj_clonotypes"] — clonotype DataFrame
  bundle extensions/vdj/
"""
from __future__ import annotations

import logging

from ...context import PipelineContext

logger = logging.getLogger(__name__)

__references__ = {
    "scirpy": {
        "title": "Scirpy: a Scanpy extension for analyzing single-cell T-cell receptor-sequencing data",
        "authors": "Sturm et al.",
        "journal": "Bioinformatics",
        "year": "2020",
        "doi": "10.1093/bioinformatics/btaa611",
        "description": "VDJ receptor analysis and clonotype calling",
    },
}


class VDJIngestModule:
    """VDJ receptor ingest: clonotype calling via scirpy."""

    name = "vdj_ingest"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {
        "obs": ["chain_pairing", "clonotype_id"],
        "uns": ["vdj_clonotypes"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        vdj_path = getattr(ctx.cfg, "vdj_dir", None)
        if vdj_path is None:
            ctx.status(self.name, "skipped", "vdj_dir not set in config")
            ctx.metadata["vdj_ingest_status"] = "skipped_no_input"
            return

        try:
            import scirpy as ir
        except ImportError:
            ctx.status(self.name, "skipped", "scirpy not installed")
            ctx.metadata["vdj_ingest_status"] = "skipped_scirpy_missing"
            return

        logger.info("%s: ingesting VDJ from %s", self.name, vdj_path)
        # DRAFT: full implementation in Wave 2
        # ir.io.read_10x_vdj(...)
        # ir.tl.chain_qc(...)
        # ir.tl.define_clonotypes(...)
        raise NotImplementedError("vdj_ingest is a Wave 2 draft — not yet implemented")
