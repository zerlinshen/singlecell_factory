"""Canonical factory path resolution for the Bioinformatics Research Pipeline.

The three repositories are siblings under the physical suite root:
``/home/zerlinshen/Bioinformatics Research Pipeline``.  Do not depend on the
retired ``/home/zerlinshen/<repo>`` compatibility paths; they are intentionally
absent in the canonical-only layout.
"""

from __future__ import annotations

from pathlib import Path


SINGLECELL_FACTORY_ROOT = Path(__file__).resolve().parents[1]
SUITE_ROOT = SINGLECELL_FACTORY_ROOT.parent
R_MULTIOMICS_FACTORY_ROOT = SUITE_ROOT / "r_multiomics_factory"
PLOTTING_FACTORY_ROOT = SUITE_ROOT / "plotting_factory"
PROJECTS_ROOT = Path("/home/zerlinshen/projects")


def suite_child(name: str) -> Path:
    """Return a canonical sibling repository path under the suite root."""

    return SUITE_ROOT / name
