"""DRAFT — Wave 2 / Phase 2 P2.S21: ribo_ingest module.

DO NOT import or register until Wave 2 branch.
Requires ribopy or riboseQC output for translation efficiency analysis.

Inputs: ribosome profiling alignment (BAM or pre-computed TE table)
Outputs:
  adata.uns["ribo_te"] — translation efficiency table (gene x sample)
  adata.obsm["ribo_te"] — per-cell TE embedding when single-cell ribo available
  runs/<run-id>/ribo/translation_efficiency.csv
  bundle extensions/ribo/
"""
from __future__ import annotations

import logging

from ...context import PipelineContext

logger = logging.getLogger(__name__)

__references__ = {
    "ribosomeProfilingQC": {
        "title": "ribosomeProfilingQC: quality control of ribosome profiling",
        "authors": "Ou et al.",
        "journal": "Bioconductor",
        "year": "2021",
        "doi": "10.18129/B9.bioc.ribosomeProfilingQC",
        "description": "R-side ribosome profiling QC and TE computation",
    },
}


class RiboIngestModule:
    """Ribosome profiling ingest: translation efficiency computation."""

    name = "ribo_ingest"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {"uns": ["ribo_te"]}

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        ribo_path = getattr(ctx.cfg, "ribo_bam_path", None)
        if ribo_path is None:
            ctx.status(self.name, "skipped", "ribo_bam_path not set in config")
            ctx.metadata["ribo_ingest_status"] = "skipped_no_input"
            return

        logger.info("%s: ingesting ribo profiling from %s", self.name, ribo_path)
        # DRAFT: full implementation in Wave 2
        # Compute translation efficiency (TE = ribo_counts / rna_counts per gene)
        # Write TE table, attach to uns
        raise NotImplementedError("ribo_ingest is a Wave 2 draft — not yet implemented")
