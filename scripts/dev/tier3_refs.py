"""Tier-3 canonical citations (project-local utility modules).

These modules orchestrate or detect rather than apply a published algorithm,
so they get an empty __references__ dict with a documentation tag in the
description per the audit doc's Tier-3 convention.
"""
from __future__ import annotations


TIER_3_REFS: dict[str, dict[str, dict[str, str]]] = {
    "modality_registry.py": {
        "project_local_utility": {
            "title": "Project-local utility — multi-omics modality inventory",
            "authors": "singlecell_factory contributors",
            "journal": "Internal documentation",
            "year": "2025",
            "doi": "PMID: 0",
            "description": "Detects present modalities by scanning obsm/obs keys; no external paper. Tier-3 project-local utility per docs/SCIENTIFIC_AUDIT_2026-05-15.md.",
        },
    },
    "cross_modality_qc.py": {
        "project_local_utility": {
            "title": "Project-local utility — cross-modality barcode-overlap QC",
            "authors": "singlecell_factory contributors",
            "journal": "Internal documentation",
            "year": "2025",
            "doi": "PMID: 0",
            "description": "Jaccard-overlap QC across declared modalities; no external paper. Tier-3 project-local utility per docs/SCIENTIFIC_AUDIT_2026-05-15.md.",
        },
    },
}
