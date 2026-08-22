"""Network-free unit and negative tests for Reference Atlas OOD mapping and Census materialization."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from scripts.benchmark_reference_mapping import select_reference_hvg_genes
from scripts.materialize_cellxgene_reference import (
    materialize_reference,
    select_stratified_soma_joinids,
    validate_census_build_date,
)
from workflow.modular._reference_mapping import (
    ASSIGNMENT_STATUS_ACCEPTED,
    ASSIGNMENT_STATUS_REJECTED_BOTH,
    ASSIGNMENT_STATUS_REJECTED_LOW_CONF,
    ASSIGNMENT_STATUS_REJECTED_OOD_DIST,
    align_reference_genes,
    calibrate_reference_ood_threshold,
    compute_parity_metrics,
    map_knn_reference,
    normalize_l2,
    resolve_reference_device,
)
from workflow.modular.config import PipelineConfig
from workflow.modular.context import PipelineContext
from workflow.modular.modules.annotation import AnnotationModule


# ---------------------------------------------------------------------------
# 1. Census Materializer Unit & Negative Tests (Network-Free)
# ---------------------------------------------------------------------------

def test_census_alias_rejected_without_network():
    """Verify dynamic aliases and invalid formats fail immediately before network calls."""
    with pytest.raises(ValueError, match="Dynamic alias 'latest' is forbidden"):
        validate_census_build_date("latest")

    with pytest.raises(ValueError, match="Dynamic alias 'stable' is forbidden"):
        validate_census_build_date("stable")

    with pytest.raises(ValueError, match="Census build date cannot be empty"):
        validate_census_build_date("")

    with pytest.raises(ValueError, match="must match YYYY-MM-DD format"):
        validate_census_build_date("2025/11/08")


def test_census_stratified_selection_determinism():
    """Verify deterministic stratified soma_joinid selection."""
    obs_df = pd.DataFrame({
        "soma_joinid": np.arange(100),
        "cell_type": ["T cell"] * 40 + ["B cell"] * 40 + ["Monocyte"] * 20,
    })
    sel1, counts1 = select_stratified_soma_joinids(obs_df, label_key="cell_type", max_cells=30, max_cells_per_label=10, seed=42)
    sel2, counts2 = select_stratified_soma_joinids(obs_df, label_key="cell_type", max_cells=30, max_cells_per_label=10, seed=42)

    np.testing.assert_array_equal(sel1, sel2)
    assert counts1 == counts2
    assert len(sel1) == 30
    assert counts1["T cell"] == 10
    assert counts1["B cell"] == 10
    assert counts1["Monocyte"] == 10


def test_census_materialization_with_mock_client(tmp_path: Path):
    """Test materialization pipeline, atomic writing, and receipt integrity with fake client."""
    mock_client = MagicMock()
    mock_client.__version__ = "1.18.0"
    mock_client.get_census_version_description.return_value = {
        "release_build": "2025-11-08",
        "soma": {"uri": "s3://mock-census"},
    }

    n_cells = 50
    mock_obs = pd.DataFrame({
        "soma_joinid": np.arange(n_cells),
        "dataset_id": ["ds1"] * n_cells,
        "donor_id": ["donor1"] * 25 + ["donor2"] * 25,
        "cell_type": ["Neuron"] * 25 + ["Astrocyte"] * 25,
        "cell_type_ontology_term_id": ["CL:0000540"] * 25 + ["CL:0000127"] * 25,
        "assay": ["10x"] * n_cells,
        "tissue": ["brain"] * n_cells,
        "tissue_general": ["central nervous system"] * n_cells,
        "development_stage": ["adult"] * n_cells,
        "disease": ["normal"] * n_cells,
        "is_primary_data": [True] * n_cells,
    })
    mock_client.get_obs.return_value = mock_obs

    # Mock get_anndata using the exact requested SOMA coordinates, returned in
    # reverse order to exercise deterministic reordering in the materializer.
    def _fake_get_anndata(*_args, obs_coords, **_kwargs):
        selected = mock_obs.set_index("soma_joinid").loc[list(reversed(obs_coords))].reset_index()
        mock_x = sp.csr_matrix(np.ones((len(selected), 30), dtype=np.float32))
        selected.index = pd.Index([f"cell_{i}" for i in range(len(selected))])
        return ad.AnnData(
            X=mock_x,
            obs=selected,
            var=pd.DataFrame(index=[f"GENE_{j}" for j in range(30)]),
        )

    mock_client.get_anndata.side_effect = _fake_get_anndata

    receipt = materialize_reference(
        project_root=tmp_path,
        census_build="2025-11-08",
        max_cells=20,
        max_cells_per_label=10,
        seed=42,
        census_client=mock_client,
    )

    assert receipt["status"] == "completed"
    assert receipt["requested_census_build"] == "2025-11-08"
    assert receipt["n_cells_selected"] == 20
    assert len(receipt["software_provenance"]["materializer_script_sha256"]) == 64
    assert Path(receipt["artifacts"]["reference_h5ad_path"]).exists()
    get_anndata_kwargs = mock_client.get_anndata.call_args.kwargs
    assert "obs_column_names" in get_anndata_kwargs
    assert "column_names" not in get_anndata_kwargs

    # Reopen and check
    reopened = ad.read_h5ad(receipt["artifacts"]["reference_h5ad_path"])
    assert reopened.shape == (20, 30)
    assert reopened.obs["soma_joinid"].tolist() == sorted(reopened.obs["soma_joinid"].tolist())


def test_census_stratified_selection_handles_categorical_labels():
    obs_df = pd.DataFrame(
        {
            "soma_joinid": np.arange(12),
            "cell_type": pd.Categorical(
                ["Neuron"] * 6 + ["Astrocyte"] * 6,
                categories=["Neuron", "Astrocyte", "Unused category"],
            ),
        }
    )
    selected, counts = select_stratified_soma_joinids(
        obs_df, label_key="cell_type", max_cells=8, max_cells_per_label=4, seed=42
    )
    assert len(selected) == 8
    assert counts == {"Astrocyte": 4, "Neuron": 4}


def test_census_stratified_selection_rejects_casefolded_unknown_labels():
    obs_df = pd.DataFrame({"soma_joinid": [1, 2], "cell_type": ["Neuron", "unknown"]})
    with pytest.raises(ValueError, match="contains 1 invalid labels"):
        select_stratified_soma_joinids(obs_df, label_key="cell_type", max_cells=2)


def test_census_build_metadata_mismatch_fails_closed(tmp_path: Path):
    mock_client = MagicMock()
    mock_client.get_census_version_description.return_value = {"release_build": "2025-11-01"}
    with pytest.raises(RuntimeError, match="resolved '2025-11-01'"):
        materialize_reference(
            project_root=tmp_path,
            census_build="2025-11-08",
            census_client=mock_client,
        )
    receipts = list(tmp_path.glob("runs/*/python/reference/census_reference_receipt.json"))
    assert len(receipts) == 1
    payload = json.loads(receipts[0].read_text(encoding="utf-8"))
    assert payload["status"] == "failed_build_mismatch"


# ---------------------------------------------------------------------------
# 2. Gene Alignment & Matrix Normalization Tests
# ---------------------------------------------------------------------------

def test_gene_alignment_direct_var_names():
    q = ad.AnnData(X=np.zeros((5, 10)), var=pd.DataFrame(index=[f"G_{i}" for i in range(10)]))
    r = ad.AnnData(X=np.zeros((5, 10)), var=pd.DataFrame(index=[f"G_{i}" for i in range(5, 15)]))

    q_g, r_g, route, n_shared, q_f, r_f = align_reference_genes(q, r, min_shared=5)
    assert route == "exact_var_names"
    assert n_shared == 5
    assert q_g == [f"G_{i}" for i in range(5, 10)]


def test_gene_alignment_fails_when_below_min_shared():
    q = ad.AnnData(X=np.zeros((5, 10)), var=pd.DataFrame(index=[f"G_{i}" for i in range(10)]))
    r = ad.AnnData(X=np.zeros((5, 10)), var=pd.DataFrame(index=[f"OTHER_{i}" for i in range(10)]))

    with pytest.raises(ValueError, match="Insufficient unambiguous shared genes"):
        align_reference_genes(q, r, min_shared=5)


def test_gene_alignment_rejects_duplicate_var_names():
    q = ad.AnnData(X=np.zeros((2, 3)), var=pd.DataFrame(index=["G1", "G1", "G2"]))
    r = ad.AnnData(X=np.zeros((2, 3)), var=pd.DataFrame(index=["G1", "G2", "G3"]))
    with pytest.raises(ValueError, match="requires unique query and reference var_names"):
        align_reference_genes(q, r, min_shared=1)


def test_gene_alignment_rejects_ambiguous_stable_ids_without_symbol_fallback():
    query = ad.AnnData(
        X=np.ones((2, 3)),
        var=pd.DataFrame({"gene_id": ["ENSG1", "ENSG1", "ENSG3"]}, index=["qa", "qb", "qc"]),
    )
    reference = ad.AnnData(
        X=np.ones((2, 3)),
        var=pd.DataFrame({"feature_id": ["ENSG1", "ENSG2", "ENSG3"]}, index=["RA", "RB", "RC"]),
    )
    with pytest.raises(ValueError, match="Stable-ID alignment is ambiguous"):
        align_reference_genes(query, reference, min_shared=2)


def test_reference_hvg_feature_selection_excludes_query_signal():
    rng = np.random.default_rng(42)
    reference_x = np.log1p(rng.poisson(2.0, size=(40, 50))).astype(float)
    query_x = np.log1p(rng.poisson(2.0, size=(2, 50))).astype(float)
    obs_names = [f"ref_{i}" for i in range(40)] + ["query_1", "query_ood"]
    var_names = [f"G_{i}" for i in range(50)]
    adata = ad.AnnData(
        X=sp.csr_matrix(np.vstack([reference_x, query_x])),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )
    changed_query = ad.AnnData(
        X=sp.csr_matrix(np.vstack([reference_x, query_x * 1000.0])),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )
    reference_cells = [f"ref_{i}" for i in range(40)]
    selected = select_reference_hvg_genes(adata, reference_cells, max_genes=10)
    selected_after_query_perturbation = select_reference_hvg_genes(changed_query, reference_cells, max_genes=10)
    assert selected == selected_after_query_perturbation


def test_reference_device_auto_uses_cpu_without_validated_domain():
    device, receipt = resolve_reference_device("auto")
    assert device == "cpu"
    assert receipt["decision"] == "cpu_unvalidated_domain"


def test_reference_device_uses_promoted_domain_without_production_dual_run(tmp_path: Path, monkeypatch):
    policy_path = tmp_path / "policy.json"
    policy_path.write_text(
        json.dumps(
            {
                "modules": {
                    "reference_mapping_knn": {
                        "domains": {
                            "validated-domain": {
                                "status": "promoted",
                                "certificate_id": "cert-1",
                                "validated_backend": {"cuml_version": "26.08.00"},
                                "evidence_run": "/evidence/run",
                                "claim_scope": "bounded",
                            }
                        }
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "workflow.modular._reference_mapping.importlib.import_module",
        lambda name: type("FakeCuML", (), {"__version__": "26.08.00"})(),
    )
    device, receipt = resolve_reference_device(
        "auto", validation_domain="validated-domain", policy_path=policy_path
    )
    assert device == "gpu"
    assert receipt["decision"] == "gpu_promoted_for_validated_domain"
    assert receipt["certificate_id"] == "cert-1"


def test_reference_device_version_drift_routes_auto_to_cpu_and_blocks_explicit_gpu(tmp_path: Path, monkeypatch):
    policy_path = tmp_path / "policy.json"
    policy_path.write_text(
        json.dumps(
            {
                "modules": {
                    "reference_mapping_knn": {
                        "domains": {
                            "validated-domain": {
                                "status": "promoted",
                                "validated_backend": {"cuml_version": "26.08.00"},
                            }
                        }
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "workflow.modular._reference_mapping.importlib.import_module",
        lambda name: type("FakeCuML", (), {"__version__": "99.0"})(),
    )
    device, receipt = resolve_reference_device(
        "auto", validation_domain="validated-domain", policy_path=policy_path
    )
    assert device == "cpu"
    assert receipt["decision"] == "cpu_backend_version_drift"
    with pytest.raises(RuntimeError, match="revalidate before use"):
        resolve_reference_device("gpu", validation_domain="validated-domain", policy_path=policy_path)


def test_normalize_l2_preserves_sparsity():
    csr = sp.csr_matrix(np.array([[3.0, 4.0], [0.0, 5.0]], dtype=np.float32))
    normed = normalize_l2(csr)
    assert sp.issparse(normed)
    dense = normed.toarray()
    assert np.isclose(np.linalg.norm(dense[0]), 1.0)
    assert np.isclose(np.linalg.norm(dense[1]), 1.0)


# ---------------------------------------------------------------------------
# 3. OOD Calibration & KNN Mapping Tests
# ---------------------------------------------------------------------------

def test_ood_calibration_fixed_mode():
    ref = ad.AnnData(X=np.zeros((20, 10)), var=pd.DataFrame(index=[f"G_{i}" for i in range(10)]))
    thresh, receipt = calibrate_reference_ood_threshold(
        ref, ref_genes=[f"G_{i}" for i in range(10)], mode="fixed", fixed_threshold=0.45
    )
    assert thresh == 0.45
    assert receipt["ood_calibration_mode"] == "fixed"

    with pytest.raises(ValueError, match="requires an explicit float fixed_threshold"):
        calibrate_reference_ood_threshold(
            ref, ref_genes=[f"G_{i}" for i in range(10)], mode="fixed", fixed_threshold=None
        )


def test_ood_calibration_whole_group_split():
    """Verify group-level split avoids sample data leakage during calibration."""
    n = 60
    ref = ad.AnnData(
        X=np.random.normal(size=(n, 10)),
        obs=pd.DataFrame({
            "sample": ["S1"] * 20 + ["S2"] * 20 + ["S3"] * 20,
            "cell_type": ["A"] * 30 + ["B"] * 30,
        }),
        var=pd.DataFrame(index=[f"G_{i}" for i in range(10)]),
    )
    thresh, receipt = calibrate_reference_ood_threshold(
        ref,
        ref_genes=[f"G_{i}" for i in range(10)],
        mode="reference_quantile",
        quantile=0.90,
        group_key="sample",
        cal_fraction=0.33,
        seed=42,
    )
    assert 0.0 < thresh < 2.0
    assert "whole_group_split_on_sample" in receipt["split_method"]


def test_ood_calibration_fails_closed_for_missing_group_and_small_reference():
    ref = ad.AnnData(
        X=np.ones((30, 10)),
        obs=pd.DataFrame({"cell_type": ["A"] * 15 + ["B"] * 15}),
        var=pd.DataFrame(index=[f"G_{i}" for i in range(10)]),
    )
    with pytest.raises(ValueError, match="group key 'sample' is absent"):
        calibrate_reference_ood_threshold(ref, list(ref.var_names), group_key="sample", k=5)

    with pytest.raises(ValueError, match="requires an explicit whole-group obs key"):
        calibrate_reference_ood_threshold(ref, list(ref.var_names), k=5)

    small = ref[:10].copy()
    with pytest.raises(ValueError, match="requires at least 20 cells"):
        calibrate_reference_ood_threshold(small, list(small.var_names), group_key="cell_type", k=5)


def test_knn_mapping_never_reduces_k_or_accepts_unknown_reference_labels():
    q = ad.AnnData(X=np.ones((2, 5)), var=pd.DataFrame(index=[f"G{i}" for i in range(5)]))
    r = ad.AnnData(
        X=np.ones((4, 5)),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}),
        var=pd.DataFrame(index=[f"G{i}" for i in range(5)]),
    )
    with pytest.raises(ValueError, match="k=5 exceeds reference size"):
        map_knn_reference(q, r, "cell_type", k=5, min_shared_genes=5)

    r.obs.loc[r.obs.index[0], "cell_type"] = "unknown"
    with pytest.raises(ValueError, match="invalid labels"):
        map_knn_reference(q, r, "cell_type", k=3, min_shared_genes=5)


def test_rejected_cells_never_override_labels():
    """Verify rejection safety contract: rejected cells cannot override cell_type under any mode."""
    # Synthetic reference with 2 distinct clusters
    np.random.seed(42)
    ref_x = np.vstack([
        np.random.normal(loc=5.0, scale=0.5, size=(30, 20)),
        np.random.normal(loc=-5.0, scale=0.5, size=(30, 20)),
    ])
    ref = ad.AnnData(
        X=ref_x,
        obs=pd.DataFrame({"cell_type": ["Neuron"] * 30 + ["Astrocyte"] * 30}),
        var=pd.DataFrame(index=[f"G_{i}" for i in range(20)]),
    )

    # Query with 1 known (close to Neuron) and 1 far OOD cell (e.g. all orthogonal)
    query_x = np.vstack([
        np.random.normal(loc=5.0, scale=0.5, size=(1, 20)),   # known cell
        np.random.normal(loc=0.0, scale=0.1, size=(1, 20)),   # OOD noise cell
    ])
    query = ad.AnnData(
        X=query_x,
        obs=pd.DataFrame({
            "cell_type_marker": ["Marker_Neuron", "Marker_Microglia"],
            "annotation_confidence": [0.05, 0.05],  # low confidence -> candidates for override
        }, index=["known_cell", "ood_cell"]),
        var=pd.DataFrame(index=[f"G_{i}" for i in range(20)]),
    )
    query.obs["cell_type"] = query.obs["cell_type_marker"]

    # Run mapping with a conservative distance threshold
    map_df, meta = map_knn_reference(
        query_adata=query,
        ref_adata=ref,
        label_key="cell_type",
        k=5,
        min_confidence=0.6,
        distance_threshold=0.2,  # Strict threshold
        device="cpu",
        min_shared_genes=10,
    )

    # Known cell should be accepted, OOD cell should be rejected
    assert map_df.loc["known_cell", "reference_assignment_status"] == ASSIGNMENT_STATUS_ACCEPTED
    assert map_df.loc["known_cell", "reference_cell_type"] == "Neuron"

    assert map_df.loc["ood_cell", "reference_assignment_status"] in {
        ASSIGNMENT_STATUS_REJECTED_OOD_DIST,
        ASSIGNMENT_STATUS_REJECTED_BOTH,
    }
    assert map_df.loc["ood_cell", "reference_cell_type"] == "Unknown"
    assert bool(map_df.loc["ood_cell", "reference_ood"]) is True

    # Check override simulation
    accepted_mask = map_df["reference_assignment_status"] == ASSIGNMENT_STATUS_ACCEPTED
    query.obs.loc[accepted_mask, "cell_type"] = map_df.loc[accepted_mask, "reference_cell_type"]

    assert query.obs.loc["known_cell", "cell_type"] == "Neuron"
    assert query.obs.loc["ood_cell", "cell_type"] == "Marker_Microglia"  # UNTOUCHED!


def test_requested_reference_missing_file_fails_loud(tmp_path: Path):
    """Verify that a requested reference file that does not exist raises FileNotFoundError."""
    ctx = PipelineContext(
        cfg=PipelineConfig(
            project="test",
            output_dir=tmp_path,
            cellranger=None,
            reference_adata=tmp_path / "non_existent.h5ad",
        ),
        run_dir=tmp_path,
        figure_dir=tmp_path,
        table_dir=tmp_path,
    )
    adata = ad.AnnData(X=np.zeros((10, 10)), var=pd.DataFrame(index=[f"G_{i}" for i in range(10)]))
    adata.obs["cell_type"] = "Unknown"
    adata.obs["cell_type_marker"] = "Unknown"
    adata.obs["annotation_confidence"] = 0.0

    with pytest.raises(FileNotFoundError, match="Requested reference H5AD file does not exist"):
        AnnotationModule._try_reference_mapping(adata, ctx)


def test_invalid_programmatic_reference_override_mode_fails_before_io(tmp_path: Path):
    ctx = PipelineContext(
        cfg=PipelineConfig(
            project="test",
            output_dir=tmp_path,
            cellranger=None,
            reference_adata=tmp_path / "non_existent.h5ad",
            reference_override_mode="typo",
        ),
        run_dir=tmp_path,
        figure_dir=tmp_path,
        table_dir=tmp_path,
    )
    adata = ad.AnnData(X=np.zeros((1, 1)), var=pd.DataFrame(index=["G1"]))
    with pytest.raises(ValueError, match="Unsupported reference_override_mode"):
        AnnotationModule._try_reference_mapping(adata, ctx)


def test_explicit_gpu_request_fails_without_fallback(monkeypatch):
    """Verify that explicit GPU request raises RuntimeError without silent CPU fallback if cuML is missing."""
    import sys
    monkeypatch.setitem(sys.modules, "cuml", None)
    monkeypatch.setitem(sys.modules, "cuml.neighbors", None)

    q = ad.AnnData(X=np.zeros((5, 10)), var=pd.DataFrame(index=[f"G_{i}" for i in range(10)]))
    r = ad.AnnData(
        X=np.zeros((10, 10)),
        obs=pd.DataFrame({"cell_type": ["A"] * 5 + ["B"] * 5}),
        var=pd.DataFrame(index=[f"G_{i}" for i in range(10)]),
    )

    with pytest.raises(RuntimeError, match="GPU reference mapping requested but cuML/CuPy cannot be imported"):
        map_knn_reference(q, r, label_key="cell_type", k=5, device="gpu", min_shared_genes=5)


def test_cpu_gpu_parity_comparison_logic():
    """Verify parity metric evaluator fails on perturbed predictions."""
    idx = [f"cell_{i}" for i in range(100)]
    df_a = pd.DataFrame({
        "reference_predicted_label": ["A"] * 50 + ["B"] * 50,
        "reference_confidence": [0.9] * 100,
        "reference_distance": [0.1] * 100,
        "reference_assignment_status": [ASSIGNMENT_STATUS_ACCEPTED] * 100,
        "reference_cell_type": ["A"] * 50 + ["B"] * 50,
        "reference_ood": [False] * 100,
    }, index=idx)

    # Identical copy
    df_b = df_a.copy()
    metrics = compute_parity_metrics(df_a, df_b)
    assert metrics["parity_verdict"] == "PASS"
    assert metrics["candidate_label_agreement"] == 1.0
    assert metrics["assignment_status_agreement"] == 1.0

    # Perturb 20 cells in df_b
    df_b.loc[idx[:20], "reference_predicted_label"] = "C"
    df_b.loc[idx[:20], "reference_assignment_status"] = ASSIGNMENT_STATUS_REJECTED_LOW_CONF
    metrics_fail = compute_parity_metrics(df_a, df_b)
    assert metrics_fail["parity_verdict"] == "FAIL"
    assert metrics_fail["candidate_label_agreement"] == 0.8
