"""Wave-9 LUCA patient expansion driver — pure-function tests.

Covers donor eligibility selection and stratified capping without touching
the atlas or LIANA (engine path is exercised by the Wave-9 evidence run).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import sparse

from scripts.liana_luca_patient_expansion import (
    MIN_CELLS_PER_TYPE,
    donor_table,
    load_existing_donor_stats,
    load_expression_subset,
    read_obs_columns,
    select_donors,
    stratified_cap,
)


def _obs() -> pd.DataFrame:
    rows = []
    # LUAD donors: D1 large, D2 medium, D3 too few cells, D4 low diversity
    rng = np.random.default_rng(0)
    types = [f"T{i}" for i in range(8)]
    for donor, n, n_types, disease in [
        ("D1", 2000, 8, "lung adenocarcinoma"),
        ("D2", 1500, 8, "lung adenocarcinoma"),
        ("D3", 500, 8, "lung adenocarcinoma"),
        ("D4", 1500, 3, "lung adenocarcinoma"),
        ("S1", 1800, 8, "squamous cell lung carcinoma"),
        ("S2", 1200, 8, "squamous cell lung carcinoma"),
    ]:
        ts = rng.choice(types[:n_types], size=n)
        rows.append(
            pd.DataFrame(
                {
                    "donor_id": donor,
                    "origin": "tumor_primary",
                    "disease": disease,
                    "cell_type_predicted": ts,
                }
            )
        )
    # non-tumour rows must be ignored
    rows.append(
        pd.DataFrame(
            {
                "donor_id": ["N1"] * 3000,
                "origin": "normal",
                "disease": "lung adenocarcinoma",
                "cell_type_predicted": ["T0"] * 3000,
            }
        )
    )
    return pd.concat(rows, ignore_index=True)


def test_donor_table_and_selection():
    dt = donor_table(_obs())
    assert set(dt["donor"]) == {"D1", "D2", "D3", "D4", "S1", "S2"}

    selected, ineligible, not_selected = select_donors(dt, "LUAD", 2)
    assert selected == ["D1", "D2"]  # top-2 by cells
    reasons = {e["donor"]: e["reason"] for e in ineligible}
    assert "D3" in reasons and "D4" in reasons
    assert str(MIN_CELLS_PER_TYPE) in reasons["D4"] or ">=" in reasons["D4"]
    assert not_selected == []  # D1/D2 are the only eligible LUAD donors

    selected_lusc, _, _ = select_donors(dt, "LUSC", 8)
    assert selected_lusc == ["S1", "S2"]  # fewer eligible than requested


def test_stratified_cap_properties():
    rng = np.random.default_rng(1)
    types = pd.Series(
        np.concatenate(
            [np.full(700, "a"), np.full(200, "b"), np.full(100, "c")]
        )
    )
    idx = np.arange(len(types))
    out = stratified_cap(idx, types, 300, "47:donorX")
    assert len(out) == 300
    kept = types.iloc[out].value_counts()
    assert set(kept.index) == {"a", "b", "c"}  # every type retained
    assert kept["a"] > kept["b"] > kept["c"]  # proportional ordering
    # determinism
    out2 = stratified_cap(idx, types, 300, "47:donorX")
    assert (out == out2).all()
    # below cap -> passthrough
    small = stratified_cap(np.arange(50), types.iloc[:50], 300, "47:donorX")
    assert len(small) == 50


def test_direct_h5_subset_ignores_large_slots(tmp_path):
    import anndata as ad

    matrix = sparse.csr_matrix(np.arange(30, dtype=np.float32).reshape(6, 5))
    atlas = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {
                "donor_id": pd.Categorical(["D1", "D1", "D2", "D2", "D3", "D3"]),
                "origin": pd.Categorical(["tumor_primary"] * 6),
                "disease": pd.Categorical(["lung adenocarcinoma"] * 6),
                "cell_type_predicted": pd.Categorical(["T", "B", "T", "B", "T", "B"]),
            },
            index=[f"cell{i}" for i in range(6)],
        ),
        var=pd.DataFrame(
            {"feature_name": pd.Categorical(["A", "B", "C", "D", "E"])},
            index=[f"ENSG{i}" for i in range(5)],
        ),
    )
    atlas.obsm["large_embedding"] = np.ones((6, 20), dtype=np.float32)
    atlas.obsp["large_graph"] = sparse.eye(6, format="csr")
    atlas.layers["counts"] = matrix.copy()
    atlas_path = tmp_path / "atlas.h5ad"
    atlas.write_h5ad(atlas_path)

    obs = read_obs_columns(
        atlas_path, ["donor_id", "origin", "disease", "cell_type_predicted"]
    )
    assert obs.loc[2, "donor_id"] == "D2"

    subset = load_expression_subset(atlas_path, np.array([1, 4]))
    np.testing.assert_array_equal(subset.X.toarray(), matrix[[1, 4]].toarray())
    assert subset.var_names.tolist() == ["A", "B", "C", "D", "E"]
    assert list(subset.obsm.keys()) == []
    assert list(subset.obsp.keys()) == []
    assert list(subset.layers.keys()) == []


def test_resume_checkpoint_requires_summary_and_csv(tmp_path):
    out_dir = tmp_path / "donor"
    table_dir = out_dir / "tables"
    table_dir.mkdir(parents=True)
    csv_path = table_dir / "cell_communication_liana.csv"
    csv_path.write_text(
        "ligand_complex,receptor_complex,source,target,magnitude_rank,specificity_rank\n"
        "VEGFA,KDR,T,B,0.1,0.2\n"
    )
    stats = {"n_liana": 1, "csv": str(csv_path), "n_cells": 1000, "n_types": 8}
    (out_dir / "donor_summary.json").write_text(json.dumps(stats))

    loaded = load_existing_donor_stats(out_dir)
    assert loaded is not None
    assert loaded["resumed"] is True

    csv_path.write_text(
        "ligand_complex,receptor_complex,source,target,magnitude_rank,specificity_rank\n"
    )
    assert load_existing_donor_stats(out_dir) is None

    csv_path.unlink()
    assert load_existing_donor_stats(out_dir) is None
