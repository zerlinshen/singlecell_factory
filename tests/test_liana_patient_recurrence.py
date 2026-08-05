"""Unit tests for patient-level LIANA recurrence helper (Wave-8)."""
from __future__ import annotations

import pandas as pd

from scripts.liana_patient_recurrence import (
    edge_keys,
    lr_keys,
    recurrence,
    summarize_histology,
)


def _liana_row(src, tgt, lig, rec, mag=0.1):
    return {
        "source": src,
        "target": tgt,
        "ligand_complex": lig,
        "receptor_complex": rec,
        "magnitude_rank": mag,
        "specificity_rank": mag,
    }


def test_lr_and_edge_keys():
    df = pd.DataFrame(
        [
            _liana_row("Macrophage", "T cell CD8", "CD274", "PDCD1"),
            _liana_row("Mast cell", "T cell CD4", "CD274", "PDCD1"),
        ]
    )
    assert lr_keys(df) == {"CD274|PDCD1"}
    assert len(edge_keys(df)) == 2


def test_recurrence_ranking():
    sets = {
        "d1": {"A|B", "C|D"},
        "d2": {"A|B"},
        "d3": {"A|B", "E|F"},
    }
    rows = recurrence(sets)
    assert rows[0]["key"] == "A|B"
    assert rows[0]["fraction"] == 1.0
    assert rows[0]["n_donors"] == 3


def test_summarize_histology_flagship():
    d1 = pd.DataFrame(
        [
            _liana_row("Macrophage", "T cell CD8", "CD274", "PDCD1", 0.01),
            _liana_row("Endothelial cell", "Epithelial cell malignant", "VEGFA", "KDR", 0.02),
            _liana_row("NK cell", "NK cell", "LGALS9", "HAVCR2", 0.03),
        ]
    )
    d2 = pd.DataFrame(
        [
            _liana_row("Macrophage", "T cell CD4", "CD80", "CTLA4", 0.01),
            _liana_row("Endothelial cell", "Epithelial cell malignant", "VEGFA", "EGFR", 0.02),
        ]
    )
    out = summarize_histology({"d1": d1, "d2": d2})
    assert out["n_donors"] == 2
    assert out["flagship_pair_recurrence"]["CD274|PDCD1"]["n_donors"] == 1
    assert out["flagship_pair_recurrence"]["VEGFA|KDR"]["n_donors"] == 1
    assert out["donors_with_any_classic_vegfr"] == 1
