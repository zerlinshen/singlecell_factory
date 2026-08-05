"""Unit tests for LIANA claim-panel filters (Wave-4)."""
from __future__ import annotations

import pandas as pd

from scripts.filter_liana_claim_panels import filter_panels


def test_filter_panels_recovers_vegf_and_checkpoint_edges():
    df = pd.DataFrame(
        {
            "source": ["Plasma cell", "Myeloid/Macro", "Tumor epithelial"],
            "target": ["Tumor epithelial", "NK cell", "Fibroblast"],
            "ligand_complex": ["VEGFA", "HLA-DRA", "COL1A1"],
            "receptor_complex": ["EGFR", "LAG3", "ITGB1"],
            "magnitude_rank": [0.1, 0.05, 0.9],
            "specificity_rank": [0.2, 0.1, 0.8],
        }
    )
    panels = filter_panels(df)
    assert len(panels["vegf_egfr"]) == 1
    assert panels["vegf_egfr"].iloc[0]["ligand_complex"] == "VEGFA"
    assert len(panels["vegfa_egfr"]) == 1
    assert len(panels["classic_vegf_vegfr"]) == 0
    assert len(panels["checkpoint_immune"]) == 1
    assert panels["checkpoint_immune"].iloc[0]["receptor_complex"] == "LAG3"


def test_filter_panels_empty_when_no_hits():
    df = pd.DataFrame(
        {
            "source": ["A"],
            "target": ["B"],
            "ligand_complex": ["FOO"],
            "receptor_complex": ["BAR"],
            "magnitude_rank": [0.1],
        }
    )
    panels = filter_panels(df)
    assert panels["vegf_egfr"].empty
    assert panels["checkpoint_immune"].empty
