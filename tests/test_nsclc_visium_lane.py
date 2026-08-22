"""Wave-9 NSCLC Visium lane regression tests (synthetic, deterministic).

Covers the lane wrapper (``analyze_adata`` panel scoring + Moran + panel
filtering + honesty stamps) and the paired tumour-vs-background contrast
builder. The LIANA engine itself is monkeypatched here — real-engine
coverage lives in ``scripts/test_joint_spatial_liana_smoke.py`` and the
Wave-9 evidence run.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import sys
from types import SimpleNamespace
from scipy import sparse

from scripts.nsclc_visium_spatial_liana import (
    BOUNDARY,
    SCORE_PANELS,
    analyze_adata,
    build_contrast,
)

PLANTED = [
    "STAB1", "FOLR2", "SLC40A1", "MERTK", "GPR34", "F13A1",
    "SELENOP", "C3", "IGF1",
    "SPP1", "APOE", "TREM2", "ABCA1", "APOC1", "FABP5",
    "LYVE1", "MRC1", "SIGLEC1", "C1QA", "CD163",
    "VEGFA", "KDR", "CD274", "PDCD1", "CTLA4", "CD80",
]


def _synthetic_visium(seed: int = 7) -> "ad.AnnData":
    import anndata as ad

    rng = np.random.default_rng(seed)
    n_spots, n_bg = 144, 300
    genes = [f"G{i}" for i in range(n_bg)] + PLANTED
    x = rng.poisson(0.5, size=(n_spots, len(genes))).astype(np.float32)
    # spatially patterned expression for planted genes (right half high)
    grid = np.array([(i // 12, i % 12) for i in range(n_spots)], dtype=float)
    right = grid[:, 1] >= 6
    for j, g in enumerate(genes):
        if g in PLANTED:
            x[right, j] += 20
    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"spot{i}" for i in range(n_spots)]),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm["spatial"] = grid
    return adata


def _fake_liana(adata, **kwargs) -> None:
    adata.uns["liana_res"] = pd.DataFrame(
        {
            "ligand_complex": ["VEGFA", "CD274", "LGALS9"],
            "receptor_complex": ["KDR", "PDCD1", "HAVCR2"],
            "source": ["cluster_0", "cluster_1", "cluster_0"],
            "target": ["cluster_1", "cluster_2", "cluster_2"],
            "magnitude_rank": [0.01, 0.02, 0.03],
            "specificity_rank": [0.01, 0.02, 0.03],
        }
    )


def _fake_spatial_neighbors(adata, **kwargs) -> None:
    coords = np.asarray(adata.obsm["spatial"])
    row: list[int] = []
    col: list[int] = []
    for i, source in enumerate(coords):
        distance = np.abs(coords - source).sum(axis=1)
        for j in np.flatnonzero(distance == 1):
            row.append(i)
            col.append(int(j))
    adata.obsp["spatial_connectivities"] = sparse.csr_matrix(
        (np.ones(len(row), dtype=np.float32), (row, col)),
        shape=(len(coords), len(coords)),
    )


def _force_weights_moran(*args, **kwargs) -> None:
    raise RuntimeError("synthetic test uses the deterministic weights fallback")


def test_analyze_adata_schema_and_honesty(monkeypatch):
    # The promoted RAPIDS environment intentionally omits LIANA because its
    # pandas constraint conflicts with RAPIDS 26.08. Supply only the mocked
    # interface this synthetic unit test consumes.
    fake_li = SimpleNamespace(mt=SimpleNamespace(rank_aggregate=_fake_liana))
    monkeypatch.setitem(sys.modules, "liana", fake_li)
    fake_sq = SimpleNamespace(
        gr=SimpleNamespace(
            spatial_neighbors=_fake_spatial_neighbors,
            spatial_autocorr=_force_weights_moran,
        )
    )
    monkeypatch.setitem(sys.modules, "squidpy", fake_sq)
    adata = _synthetic_visium()
    out = analyze_adata(adata, max_spots=144, seed=41)

    s = out["summary"]
    assert s["spatial_connectivities"] is True
    assert s["n_liana"] == 3
    assert s["n_vegf_egfr"] >= 1
    assert s["n_classic_vegf_vegfr"] >= 1
    assert s["n_checkpoint_any"] >= 1
    # Leiden cluster labels never match curated immune labels
    assert s["n_checkpoint_immune"] == 0
    assert "checkpoint_immune_note" in s
    assert s["panel_gene_coverage"]["de_zuani_stab1_foetal"] == "6/6"
    assert "moran_engine" in s

    moran = out["moran_table"]
    assert {"feature", "I"} <= set(moran.columns)
    assert len(moran) > 0
    # planted spatial pattern must be detectable (either moran engine)
    planted_i = moran.loc[moran["feature"] == "STAB1", "I"].iloc[0]
    assert planted_i > 0.2
    # panel score columns scored (not all-NaN)
    means = s["panel_score_means"]
    assert means["score_de_zuani_stab1_foetal"] is not None


def test_build_contrast_descriptive_stamps():
    section_map = pd.DataFrame(
        {
            "section": ["P24_T1", "P24_B1"],
            "patient": ["Patient 24", "Patient 24"],
            "histology": ["LUAD", "LUAD"],
            "condition": ["tumor", "background"],
        }
    )

    def _row(sec, ckpt, mean_shift):
        return {
            "section": sec,
            "summary": {
                "n_spots_analyzed": 1000,
                "n_clusters": 8,
                "n_liana": 5000,
                "n_vegf_egfr": 100,
                "n_classic_vegf_vegfr": 10,
                "n_checkpoint_any": ckpt,
                "n_checkpoint_immune": 0,
                "panel_score_means": {
                    f"score_{p}": mean_shift for p in SCORE_PANELS
                },
                "panel_score_moranI": {
                    f"score_{p}": 0.3 + mean_shift for p in SCORE_PANELS
                },
            },
        }

    rows = [_row("P24_T1", 40, 0.5), _row("P24_B1", 25, 0.1)]
    table, summary = build_contrast(rows, section_map)

    assert len(table) == 2
    assert summary["figure_parity"] is False
    assert summary["claimable_as_nsclc_tissue"] is True
    assert summary["scientific_boundary"] == BOUNDARY
    assert "not paper figure parity" in BOUNDARY
    assert len(summary["patients"]) == 1
    entry = summary["patients"][0]
    assert entry["patient"] == "Patient 24"
    assert entry["histology"] == "LUAD"
    assert entry["delta_n_checkpoint_any"] == pytest.approx(15.0)
    assert entry["panel_score_deltas"]["de_zuani_stab1_foetal"] == pytest.approx(0.4)
