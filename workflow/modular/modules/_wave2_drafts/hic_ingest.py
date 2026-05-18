"""DRAFT — Wave 2 / Phase 2 P2.S18: hic_ingest module.

DO NOT import or register until Wave 2 branch.
Requires HiContacts (R) or cooltools (Python) for chromatin contact analysis.

Inputs: .cool/.mcool contact matrix
Outputs:
  adata.uns["hic_compartments"] — A/B compartment calls per chromosome
  adata.obsm["hic_compartment"] — per-cell compartment embedding (when single-cell Hi-C)
  runs/<run-id>/hic/compartments.csv
  bundle extensions/hic/
"""
from __future__ import annotations

import logging

from ...context import PipelineContext

logger = logging.getLogger(__name__)

__references__ = {
    "HiContacts": {
        "title": "HiContacts: an R/Bioconductor package for chromatin contacts",
        "authors": "Serizay et al.",
        "journal": "Bioconductor",
        "year": "2023",
        "doi": "10.18129/B9.bioc.HiContacts",
        "description": "R-side Hi-C contact analysis and compartment calling",
    },
    "cooltools": {
        "title": "Cooltools: enabling high-resolution Hi-C analysis in Python",
        "authors": "Open2C et al.",
        "journal": "PLOS Computational Biology",
        "year": "2022",
        "doi": "10.1371/journal.pcbi.1010802",
        "description": "Python-side Hi-C contact matrix analysis",
    },
}


class HiCIngestModule:
    """Hi-C ingest: A/B compartment calling from contact matrix."""

    name = "hic_ingest"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {"uns": ["hic_compartments"]}

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        hic_path = getattr(ctx.cfg, "hic_cool_path", None)
        if hic_path is None:
            ctx.status(self.name, "skipped", "hic_cool_path not set in config")
            ctx.metadata["hic_ingest_status"] = "skipped_no_input"
            return

        try:
            import cooler
        except ImportError:
            ctx.status(self.name, "skipped", "cooler not installed")
            ctx.metadata["hic_ingest_status"] = "skipped_cooler_missing"
            return

        logger.info("%s: ingesting Hi-C from %s", self.name, hic_path)
        # DRAFT: full implementation in Wave 2
        # cooler.Cooler(hic_path)
        # Compute eigenvectors for A/B compartments
        # Write compartments.csv, attach to uns
        raise NotImplementedError("hic_ingest is a Wave 2 draft — not yet implemented")
