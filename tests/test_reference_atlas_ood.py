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

import scripts.benchmark_reference_mapping as benchmark_module
from scripts.benchmark_reference_mapping import (
    _json_sha256,
    _sha256_bytes,
    _write_frozen_input_manifest,
    compute_sha256,
    run_lane,
    select_reference_features,
    select_reference_hvg_genes,
    validate_feature_selection_receipt,
    verify_frozen_inputs,
    verify_lane_output,
)
from scripts.materialize_cellxgene_reference import (
    materialize_reference,
    required_observation_columns,
    select_stratified_soma_joinids,
    validate_census_build_date,
    verify_materialization_receipt,
)
from workflow.modular._reference_mapping import (
    ASSIGNMENT_STATUS_ACCEPTED,
    ASSIGNMENT_STATUS_REJECTED_BOTH,
    ASSIGNMENT_STATUS_REJECTED_LOW_CONF,
    ASSIGNMENT_STATUS_REJECTED_OOD_DIST,
    _configure_cuda_toolkit_path,
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
    assert get_anndata_kwargs["obs_column_names"] == required_observation_columns("cell_type", "donor_id")
    assert receipt["matrix_storage"] == "csr"
    assert receipt["environment_locks"]
    receipt_path = Path(receipt["artifacts"]["reference_h5ad_path"]).parent / "census_reference_receipt.json"
    assert verify_materialization_receipt(receipt_path)["status"] == "verified"

    # Reopen and check
    reopened = ad.read_h5ad(receipt["artifacts"]["reference_h5ad_path"])
    assert reopened.shape == (20, 30)
    assert reopened.obs["soma_joinid"].tolist() == sorted(reopened.obs["soma_joinid"].tolist())

    h5ad_path = Path(receipt["artifacts"]["reference_h5ad_path"])
    original_h5ad = h5ad_path.read_bytes()
    h5ad_path.write_bytes(original_h5ad + b"tamper")
    with pytest.raises(ValueError, match="Artifact SHA-256 mismatch"):
        verify_materialization_receipt(receipt_path)
    h5ad_path.write_bytes(original_h5ad)

    joinids_path = Path(receipt["artifacts"]["selected_joinids_path"])
    original_joinids = joinids_path.read_text(encoding="utf-8")
    joinids_path.write_text(original_joinids + "999999\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Artifact SHA-256 mismatch"):
        verify_materialization_receipt(receipt_path)
    joinids_path.write_text(original_joinids, encoding="utf-8")

    original_receipt = receipt_path.read_text(encoding="utf-8")
    receipt_payload = json.loads(original_receipt)
    receipt_payload["seed"] = 999
    receipt_path.write_text(json.dumps(receipt_payload), encoding="utf-8")
    with pytest.raises(ValueError, match="Receipt SHA-256"):
        verify_materialization_receipt(receipt_path)
    receipt_path.write_text(original_receipt, encoding="utf-8")


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
    assert payload["status"] == "failed_materialization"


def test_census_materializer_network_free_import_failure_writes_receipt(tmp_path: Path, monkeypatch):
    """Tests can import and exercise the materializer without Census installed."""
    import sys

    monkeypatch.setitem(sys.modules, "cellxgene_census", None)
    monkeypatch.setitem(sys.modules, "tiledbsoma", None)
    with pytest.raises(RuntimeError, match="dependencies not installed"):
        materialize_reference(project_root=tmp_path, census_build="2025-11-08")
    receipt_path = next(tmp_path.glob("runs/*/python/reference/census_reference_receipt.json"))
    assert json.loads(receipt_path.read_text(encoding="utf-8"))["status"] == "failed_materialization"


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


def test_reference_feature_selection_excludes_query_signal():
    """Query-only spikes cannot cross the explicit reference-selection boundary."""
    rng = np.random.default_rng(42)
    ref_x = np.zeros((40, 50), dtype=np.float32)
    for g in range(50):
        ref_x[:, g] = rng.normal(loc=1.0, scale=(g + 1) * 0.1, size=40)
    reference_cells = [f"ref_{i}" for i in range(40)]
    query_cells = [f"query_{i}" for i in range(20)]
    obs = pd.DataFrame(index=[*reference_cells, *query_cells])
    var = pd.DataFrame(index=[f"G_{i}" for i in range(50)])
    baseline = ad.AnnData(
        X=sp.vstack(
            [sp.csr_matrix(ref_x), sp.csr_matrix((len(query_cells), ref_x.shape[1]), dtype=np.float32)],
            format="csr",
        ),
        obs=obs,
        var=var,
    )

    spike_rows = np.arange(len(query_cells), dtype=np.int64)
    spiked_query = sp.csr_matrix(
        (
            np.full(len(query_cells), 1_000_000.0, dtype=np.float32),
            (spike_rows, np.zeros(len(query_cells), dtype=np.int64)),
        ),
        shape=(len(query_cells), ref_x.shape[1]),
    )
    spiked = ad.AnnData(
        X=sp.vstack([sp.csr_matrix(ref_x), spiked_query], format="csr"),
        obs=obs.copy(),
        var=var.copy(),
    )

    selected_1, receipt_1 = select_reference_hvg_genes(baseline, reference_cells, max_genes=15)
    selected_2, receipt_2 = select_reference_hvg_genes(spiked, reference_cells, max_genes=15)
    assert len(selected_1) == 15
    assert receipt_1["selection_route"] == "reference_only_sparse_variance_v1"
    assert receipt_1["independently_query_blind"] is True
    assert receipt_1["fit_population"] == "reference_cells_only"
    assert receipt_1["query_cells_consumed"] is False
    assert selected_1 == selected_2
    assert receipt_1["ordered_genes_sha256"] == receipt_2["ordered_genes_sha256"]
    assert _json_sha256(receipt_1) == _json_sha256(receipt_2)
    assert receipt_1["fit_cells_sha256"] == receipt_2["fit_cells_sha256"]
    assert receipt_1["n_fit_cells"] == receipt_2["n_fit_cells"] == len(reference_cells)
    assert "G_0" not in selected_1[:10]

    leaky_selected, leaky_receipt = select_reference_features(spiked, max_genes=15)
    assert leaky_selected[0] == "G_0"
    assert leaky_selected != selected_1
    assert leaky_receipt["n_fit_cells"] == spiked.n_obs

    with pytest.raises(TypeError, match="required positional argument"):
        select_reference_hvg_genes(baseline, max_genes=15)  # type: ignore[call-arg]


def test_reference_feature_selection_reference_sensitivity():
    """Verify that modifying only reference expression changes feature rankings in expected direction."""
    rng = np.random.default_rng(42)
    ref_x = np.zeros((40, 20), dtype=np.float32)
    # G_0 to G_18 have low variance, G_19 has moderate variance
    for g in range(19):
        ref_x[:, g] = rng.normal(loc=1.0, scale=0.1, size=40)
    ref_x[:, 19] = rng.normal(loc=1.0, scale=1.0, size=40)

    ref_adata = ad.AnnData(
        X=sp.csr_matrix(ref_x),
        obs=pd.DataFrame(index=[f"ref_{i}" for i in range(40)]),
        var=pd.DataFrame(index=[f"G_{i}" for i in range(20)]),
    )
    selected_initial, receipt_initial = select_reference_features(ref_adata, max_genes=5)
    assert selected_initial[0] == "G_19"

    # Now perturb reference expression so G_0 becomes the highest variance feature
    ref_x_modified = ref_x.copy()
    ref_x_modified[:, 0] = rng.normal(loc=1.0, scale=10.0, size=40)
    ref_adata_modified = ad.AnnData(
        X=sp.csr_matrix(ref_x_modified),
        obs=pd.DataFrame(index=[f"ref_{i}" for i in range(40)]),
        var=pd.DataFrame(index=[f"G_{i}" for i in range(20)]),
    )
    selected_modified, receipt_modified = select_reference_features(ref_adata_modified, max_genes=5)
    assert selected_modified[0] == "G_0"
    assert selected_modified != selected_initial
    assert receipt_modified["ordered_genes_sha256"] != receipt_initial["ordered_genes_sha256"]
    assert _json_sha256(receipt_modified) != _json_sha256(receipt_initial)


def test_reference_feature_selection_sparse_dense_determinism():
    """Verify deterministic output across CSR, CSC, and ndarray representations without densification."""
    rng = np.random.default_rng(42)
    dense_x = np.maximum(0.0, rng.normal(loc=1.0, scale=1.0, size=(50, 30))).astype(np.float32)
    csr_x = sp.csr_matrix(dense_x)
    csc_x = sp.csc_matrix(dense_x)

    obs = pd.DataFrame(index=[f"c_{i}" for i in range(50)])
    var = pd.DataFrame(index=[f"g_{i}" for i in range(30)])

    ad_csr = ad.AnnData(X=csr_x, obs=obs, var=var)
    ad_csc = ad.AnnData(X=csc_x, obs=obs, var=var)
    ad_dense = ad.AnnData(X=dense_x, obs=obs, var=var)

    genes_csr, receipt_csr = select_reference_features(ad_csr, max_genes=10)
    genes_csc, receipt_csc = select_reference_features(ad_csc, max_genes=10)
    genes_dense, receipt_dense = select_reference_features(ad_dense, max_genes=10)

    assert genes_csr == genes_csc == genes_dense
    assert receipt_csr["ordered_genes_sha256"] == receipt_csc["ordered_genes_sha256"] == receipt_dense["ordered_genes_sha256"]
    assert receipt_csr["matrix_storage_class"] == "csr"
    assert receipt_csc["matrix_storage_class"] == "csc"
    assert receipt_dense["matrix_storage_class"] == "ndarray"


def test_feature_selection_receipt_validation_fail_closed():
    """Verify validator fails closed on missing, false, old-route, or hash-mismatched receipts."""
    ref_cells = [f"r{i}" for i in range(10)]
    genes = [f"G{i}" for i in range(5)]
    fit_cells_sha = _sha256_bytes("\n".join(ref_cells).encode("utf-8") + b"\n")
    ordered_genes_sha = _sha256_bytes("\n".join(genes).encode("utf-8") + b"\n")

    valid_receipt = {
        "schema_version": "1.0",
        "selection_route": "reference_only_sparse_variance_v1",
        "independently_query_blind": True,
        "fit_population": "reference_cells_only",
        "query_cells_consumed": False,
        "n_fit_cells": 10,
        "fit_cells_sha256": fit_cells_sha,
        "n_input_features": 20,
        "matrix_storage_class": "csr",
        "score_method": "sparse_sample_variance",
        "tie_break_rule": "lexicographical_gene_symbol_ascending",
        "requested_cap": 5,
        "effective_cap": 5,
        "selected_count": 5,
        "ordered_genes_sha256": ordered_genes_sha,
    }
    # Valid receipt passes
    validate_feature_selection_receipt(valid_receipt, expected_reference_cells=ref_cells, expected_genes=genes)

    # 1. Non-dict receipt
    with pytest.raises(ValueError, match="must be a dictionary"):
        validate_feature_selection_receipt(None, expected_reference_cells=ref_cells, expected_genes=genes)  # type: ignore

    # 2. Old / historical route
    old_route = valid_receipt.copy()
    old_route["selection_route"] = "source_highly_variable_ordered_cap"
    with pytest.raises(ValueError, match="invalid or historical"):
        validate_feature_selection_receipt(old_route, expected_reference_cells=ref_cells, expected_genes=genes)

    # 3. independently_query_blind is False
    leakage = valid_receipt.copy()
    leakage["independently_query_blind"] = False
    with pytest.raises(ValueError, match="must declare independently_query_blind: true"):
        validate_feature_selection_receipt(leakage, expected_reference_cells=ref_cells, expected_genes=genes)

    # 4. fit_population not reference_cells_only
    pop = valid_receipt.copy()
    pop["fit_population"] = "reference_and_query"
    with pytest.raises(ValueError, match="fit_population must be 'reference_cells_only'"):
        validate_feature_selection_receipt(pop, expected_reference_cells=ref_cells, expected_genes=genes)

    # 5. query_cells_consumed is True
    consumed = valid_receipt.copy()
    consumed["query_cells_consumed"] = True
    with pytest.raises(ValueError, match="indicates query cells were consumed"):
        validate_feature_selection_receipt(consumed, expected_reference_cells=ref_cells, expected_genes=genes)

    # 6. Fit cell count / hash mismatch
    count_mismatch = valid_receipt.copy()
    count_mismatch["n_fit_cells"] = 9
    with pytest.raises(ValueError, match="fit cell count mismatch"):
        validate_feature_selection_receipt(count_mismatch, expected_reference_cells=ref_cells, expected_genes=genes)

    hash_mismatch = valid_receipt.copy()
    hash_mismatch["fit_cells_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="fit cells SHA-256 does not match"):
        validate_feature_selection_receipt(hash_mismatch, expected_reference_cells=ref_cells, expected_genes=genes)

    # 7. Ordered genes count / hash mismatch
    gene_count_mismatch = valid_receipt.copy()
    gene_count_mismatch["selected_count"] = 4
    with pytest.raises(ValueError, match="gene count mismatch"):
        validate_feature_selection_receipt(gene_count_mismatch, expected_reference_cells=ref_cells, expected_genes=genes)

    gene_hash_mismatch = valid_receipt.copy()
    gene_hash_mismatch["ordered_genes_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="ordered genes SHA-256 does not match"):
        validate_feature_selection_receipt(gene_hash_mismatch, expected_reference_cells=ref_cells, expected_genes=genes)


def test_reference_feature_selection_rejects_duplicate_names_and_zero_variance():
    """Verify input validation: duplicate var_names, duplicate obs_names, zero variance, max_genes < 1."""
    # max_genes < 1
    ad_valid = ad.AnnData(X=np.ones((5, 5)), var=pd.DataFrame(index=[f"g{i}" for i in range(5)]))
    with pytest.raises(ValueError, match="max_genes must be >= 1"):
        select_reference_features(ad_valid, max_genes=0)

    # Duplicate var_names
    ad_dup_var = ad.AnnData(X=np.ones((5, 3)), var=pd.DataFrame(index=["g1", "g1", "g2"]))
    with pytest.raises(ValueError, match="var_names must be unique"):
        select_reference_features(ad_dup_var, max_genes=2)

    # Duplicate obs_names
    ad_dup_obs = ad.AnnData(
        X=np.ones((3, 3)),
        obs=pd.DataFrame(index=["c1", "c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    with pytest.raises(ValueError, match="obs_names must be unique"):
        select_reference_features(ad_dup_obs, max_genes=2)

    # Zero variance for all features
    ad_const = ad.AnnData(
        X=np.full((10, 5), fill_value=3.0, dtype=np.float32),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(10)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(5)]),
    )
    with pytest.raises(ValueError, match="No usable variable features found"):
        select_reference_features(ad_const, max_genes=3)


def test_reference_device_auto_uses_cpu_without_validated_domain():
    assert PipelineConfig(project="test", output_dir=Path("/tmp"), cellranger=None).reference_device == "auto"
    device, receipt = resolve_reference_device("auto")
    assert device == "cpu"
    assert receipt["decision"] == "cpu_unvalidated_domain"
    device, receipt = resolve_reference_device("cpu")
    assert device == "cpu"
    assert receipt["decision"] == "explicit_cpu"


def test_reference_gpu_mode_off_routes_auto_cpu_and_blocks_explicit_gpu():
    device, receipt = resolve_reference_device(
        "auto", validation_domain="validated-domain", gpu_mode="off"
    )
    assert device == "cpu"
    assert receipt["decision"] == "cpu_gpu_mode_off"
    with pytest.raises(ValueError, match="gpu_mode='off'"):
        resolve_reference_device("gpu", validation_domain="validated-domain", gpu_mode="off")


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
        "workflow.modular._reference_mapping._preflight_reference_gpu",
        lambda: {"cuml_version": "26.08.00", "cuda_residency_verified": True},
    )
    device, receipt = resolve_reference_device(
        "auto", validation_domain="validated-domain", policy_path=policy_path
    )
    assert device == "gpu"
    assert receipt["decision"] == "gpu_promoted_for_validated_domain"
    assert receipt["certificate_id"] == "cert-1"
    assert receipt["production_dual_run"] is False


def test_reference_gpu_missing_backend_routes_auto_cpu_and_blocks_explicit(tmp_path: Path, monkeypatch):
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

    def _missing_backend():
        raise ImportError("missing test backend")

    monkeypatch.setattr("workflow.modular._reference_mapping._preflight_reference_gpu", _missing_backend)
    device, receipt = resolve_reference_device(
        "auto", validation_domain="validated-domain", policy_path=policy_path
    )
    assert device == "cpu"
    assert receipt["decision"] == "cpu_gpu_backend_unavailable"
    with pytest.raises(RuntimeError, match="No CPU fallback is allowed"):
        resolve_reference_device(
            "gpu", validation_domain="validated-domain", policy_path=policy_path
        )


def test_cuda_toolkit_resolution_replaces_invalid_inherited_path(tmp_path: Path, monkeypatch):
    fake_prefix = tmp_path / "fake-env"
    header = fake_prefix / "targets" / "x86_64-linux" / "include" / "cuda_runtime.h"
    header.parent.mkdir(parents=True)
    header.write_text("// test header\n", encoding="utf-8")
    monkeypatch.setattr("workflow.modular._reference_mapping.sys.prefix", str(fake_prefix))
    monkeypatch.setenv("CUDA_PATH", str(tmp_path / "missing-cuda"))
    receipt = _configure_cuda_toolkit_path()
    assert receipt["cuda_path"] == str(header.parents[1])
    assert receipt["cuda_path_source"] == "interpreter_target_toolkit"
    assert receipt["cuda_headers_verified"] is True


def test_reference_device_version_drift_routes_auto_cpu_and_blocks_explicit(tmp_path: Path, monkeypatch):
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
        "workflow.modular._reference_mapping._preflight_reference_gpu",
        lambda: {"cuml_version": "99.0", "cuda_residency_verified": True},
    )
    device, receipt = resolve_reference_device(
        "auto", validation_domain="validated-domain", policy_path=policy_path
    )
    assert device == "cpu"
    assert receipt["decision"] == "cpu_backend_version_drift"
    with pytest.raises(RuntimeError, match="revalidate before use"):
        resolve_reference_device(
            "gpu", validation_domain="validated-domain", policy_path=policy_path
        )


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
    assert receipt["distance_threshold"] == thresh
    assert receipt["fit_groups"] == sorted(receipt["fit_groups"])
    assert receipt["calibration_groups"] == sorted(receipt["calibration_groups"])
    assert not set(receipt["fit_groups"]) & set(receipt["calibration_groups"])
    assert len(receipt["fit_groups_sha256"]) == 64
    assert len(receipt["calibration_groups_sha256"]) == 64


def test_ood_calibration_is_query_independent_and_seed_repeatable():
    rng = np.random.default_rng(13)
    ref = ad.AnnData(
        X=rng.normal(size=(60, 8)),
        obs=pd.DataFrame({"sample": [f"S{idx // 10}" for idx in range(60)]}),
        var=pd.DataFrame(index=[f"G_{idx}" for idx in range(8)]),
    )
    query = ad.AnnData(X=rng.normal(size=(12, 8)), var=ref.var.copy())
    before = calibrate_reference_ood_threshold(ref, list(ref.var_names), group_key="sample", seed=42, k=5)
    query.X = query.X * 10_000.0  # Held-out query values never enter calibration.
    after = calibrate_reference_ood_threshold(ref, list(ref.var_names), group_key="sample", seed=42, k=5)
    assert before == after


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
            reference_ood_mode="fixed",
            reference_fixed_distance_threshold=0.5,
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


def test_requested_reference_quantile_requires_group_before_reference_io(tmp_path: Path):
    ctx = PipelineContext(
        cfg=PipelineConfig(
            project="test",
            output_dir=tmp_path,
            cellranger=None,
            reference_adata=tmp_path / "non_existent.h5ad",
            reference_ood_mode="reference_quantile",
            reference_calibration_group_key=None,
        ),
        run_dir=tmp_path,
        figure_dir=tmp_path,
        table_dir=tmp_path,
    )
    adata = ad.AnnData(X=np.zeros((1, 1)), var=pd.DataFrame(index=["G1"]))
    with pytest.raises(ValueError, match="requires reference_calibration_group_key before reference compute"):
        AnnotationModule._try_reference_mapping(adata, ctx)


def test_programmatic_gpu_mode_off_fails_before_reference_io(tmp_path: Path):
    ctx = PipelineContext(
        cfg=PipelineConfig(
            project="test",
            output_dir=tmp_path,
            cellranger=None,
            reference_adata=tmp_path / "non_existent.h5ad",
            reference_ood_mode="fixed",
            reference_fixed_distance_threshold=0.5,
            reference_device="gpu",
            gpu_mode="off",
        ),
        run_dir=tmp_path,
        figure_dir=tmp_path,
        table_dir=tmp_path,
    )
    adata = ad.AnnData(X=np.zeros((1, 1)), var=pd.DataFrame(index=["G1"]))
    with pytest.raises(ValueError, match="gpu_mode='off'"):
        AnnotationModule._try_reference_mapping(adata, ctx)


def test_no_reference_records_marker_only_skip(tmp_path: Path):
    ctx = PipelineContext(
        cfg=PipelineConfig(project="test", output_dir=tmp_path, cellranger=None),
        run_dir=tmp_path,
        figure_dir=tmp_path,
        table_dir=tmp_path,
    )
    adata = ad.AnnData(X=np.ones((2, 2)), var=pd.DataFrame(index=["G1", "G2"]))
    adata.obs["cell_type"] = ["MarkerA", "MarkerB"]
    AnnotationModule._try_reference_mapping(adata, ctx)
    assert adata.obs["cell_type"].tolist() == ["MarkerA", "MarkerB"]
    assert ctx.metadata["reference_mapping_status"] == "skipped_no_reference"
    assert adata.uns["annotation"]["reference_mapping"]["status"] == "skipped_no_reference"


def test_override_all_preserves_rejected_marker_labels(tmp_path: Path):
    reference = ad.AnnData(
        X=np.array([[1.0, 0.0], [0.9, 0.0], [0.0, 1.0], [0.0, 0.9]]),
        obs=pd.DataFrame({"cell_type": ["TypeA", "TypeA", "TypeB", "TypeB"]}),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    reference_path = tmp_path / "reference.h5ad"
    reference.write_h5ad(reference_path)
    query = ad.AnnData(
        X=np.array([[1.0, 0.0], [-1.0, 0.0]]),
        obs=pd.DataFrame(
            {
                "cell_type": ["MarkerA", "MarkerRejected"],
                "cell_type_marker": ["MarkerA", "MarkerRejected"],
                "annotation_confidence": [0.0, 0.0],
            },
            index=["accepted", "rejected"],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    ctx = PipelineContext(
        cfg=PipelineConfig(
            project="test",
            output_dir=tmp_path,
            cellranger=None,
            reference_adata=reference_path,
            reference_k=1,
            reference_min_shared_genes=2,
            reference_ood_mode="fixed",
            reference_fixed_distance_threshold=0.1,
            reference_override_mode="all",
        ),
        run_dir=tmp_path,
        figure_dir=tmp_path,
        table_dir=tmp_path,
    )
    AnnotationModule._try_reference_mapping(query, ctx)
    assert query.obs.loc["accepted", "cell_type"] == "TypeA"
    assert query.obs.loc["rejected", "reference_cell_type"] == "Unknown"
    assert query.obs.loc["rejected", "cell_type"] == "MarkerRejected"


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
    monkeypatch.setattr("workflow.modular._reference_mapping._configure_cuda_toolkit_path", lambda: {})

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

    label_only = df_a.copy()
    label_only.loc[idx[:20], "reference_predicted_label"] = "C"
    assert compute_parity_metrics(df_a, label_only)["parity_verdict"] == "FAIL"

    status_only = df_a.copy()
    status_only.loc[idx[:20], "reference_assignment_status"] = ASSIGNMENT_STATUS_REJECTED_LOW_CONF
    assert compute_parity_metrics(df_a, status_only)["parity_verdict"] == "FAIL"

    distance_only = df_a.copy()
    distance_only.loc[idx[0], "reference_distance"] = 0.5
    assert compute_parity_metrics(df_a, distance_only)["parity_verdict"] == "FAIL"


def _write_frozen_input_fixture(tmp_path: Path) -> Path:
    frozen = tmp_path / "frozen_inputs"
    frozen.mkdir()
    reference = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32))
    query = sp.csr_matrix(np.tile(np.array([[1.0, 0.0]], dtype=np.float32), (256, 1)))
    sp.save_npz(frozen / "reference_matrix.npz", reference)
    sp.save_npz(frozen / "query_matrix.npz", query)
    (frozen / "reference_cells.txt").write_text("r1\nr2\n", encoding="utf-8")
    (frozen / "query_cells.txt").write_text("".join(f"q{idx}\n" for idx in range(256)), encoding="utf-8")
    (frozen / "reference_labels.txt").write_text("A\nB\n", encoding="utf-8")
    (frozen / "query_proxy_labels.txt").write_text("".join("A\n" for _ in range(256)), encoding="utf-8")
    (frozen / "genes.txt").write_text("G1\nG2\n", encoding="utf-8")
    feature_receipt = {
        "schema_version": "1.0",
        "selection_route": "reference_only_sparse_variance_v1",
        "independently_query_blind": True,
        "fit_population": "reference_cells_only",
        "query_cells_consumed": False,
        "n_fit_cells": 2,
        "fit_cells_sha256": compute_sha256(frozen / "reference_cells.txt"),
        "n_input_features": 2,
        "matrix_storage_class": "csr",
        "score_method": "sparse_sample_variance",
        "tie_break_rule": "lexicographical_gene_symbol_ascending",
        "requested_cap": 2,
        "effective_cap": 2,
        "selected_count": 2,
        "ordered_genes_sha256": compute_sha256(frozen / "genes.txt"),
    }
    (frozen / "feature_selection_receipt.json").write_text(
        json.dumps(feature_receipt, indent=2) + "\n", encoding="utf-8"
    )
    _write_frozen_input_manifest(
        frozen,
        {"source_sha256": "a" * 64, "n_reference_cells": 2, "n_query_cells": 256, "n_genes": 2},
        {"distance_threshold": 0.123456789, "ood_calibration_mode": "reference_quantile"},
        {
            "k": 1,
            "min_confidence": 0.5,
            "seed": 42,
            "warmup_rows": 256,
            "measured_repetitions": 3,
            "consumed_output_columns": [],
        },
    )
    return frozen


@pytest.mark.parametrize(
    "artifact_name",
    [
        "reference_matrix.npz",
        "genes.txt",
        "ood_threshold_receipt.json",
        "parameters.json",
        "query_cells.txt",
        "feature_selection_receipt.json",
    ],
)
def test_process_lane_contract_rejects_hash_mismatch(tmp_path: Path, artifact_name: str):
    frozen = _write_frozen_input_fixture(tmp_path)
    target = frozen / artifact_name
    if artifact_name.endswith(".npz"):
        sp.save_npz(target, sp.csr_matrix(np.ones((2, 2), dtype=np.float32)))
    else:
        target.write_text(target.read_text(encoding="utf-8") + "tampered\n", encoding="utf-8")
    with pytest.raises(ValueError, match=f"Frozen input SHA-256 mismatch: {artifact_name}"):
        verify_frozen_inputs(frozen)


def test_schema3_frozen_inputs_require_feature_selection_proof(tmp_path: Path):
    frozen = _write_frozen_input_fixture(tmp_path)
    manifest_path = frozen / "input_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"].pop("feature_selection_receipt.json")
    manifest.pop("feature_selection_receipt_sha256")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Schema >=3 frozen input manifest"):
        verify_frozen_inputs(frozen)


def test_run_benchmark_requires_persisted_feature_selection_receipt_before_calibration(tmp_path: Path, monkeypatch):
    frozen = _write_frozen_input_fixture(tmp_path)
    feature_receipt = json.loads((frozen / "feature_selection_receipt.json").read_text(encoding="utf-8"))
    (frozen / "feature_selection_receipt.json").unlink()
    reference = ad.AnnData(
        X=sp.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame({"sample": ["S1", "S2"], "cell_type": ["A", "B"]}, index=["r1", "r2"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    query = ad.AnnData(
        X=sp.csr_matrix(np.ones((256, 2))),
        obs=pd.DataFrame({"cell_type": ["A"] * 256}, index=[f"q{idx}" for idx in range(256)]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )

    monkeypatch.setattr(
        benchmark_module,
        "prepare_trevino_split",
        lambda **_kwargs: {
            "reference_adata": reference,
            "query_adata": query,
            "genes": ["G1", "G2"],
            "frozen_dir": frozen,
            "split_receipt": {
                "source_sha256": "a" * 64,
                "held_out_ood_label": "Microglia",
                "n_reference_cells": 2,
                "n_query_cells": 256,
                "n_genes": 2,
                "feature_selection": feature_receipt,
            },
        },
    )
    monkeypatch.setattr(
        benchmark_module,
        "calibrate_reference_ood_threshold",
        lambda **_kwargs: pytest.fail("calibration must not run without persisted feature proof"),
    )
    monkeypatch.setattr(
        benchmark_module,
        "_run_lane_process",
        lambda **_kwargs: pytest.fail("child lanes must not launch without persisted feature proof"),
    )

    with pytest.raises(FileNotFoundError, match="missing before OOD calibration"):
        benchmark_module.run_benchmark(
            project_root=tmp_path,
            run_id="2026-08-23T2358Z-0000000",
            allow_dirty=True,
        )


def test_within_lane_scientific_nondeterminism_writes_failure_receipt(tmp_path: Path, monkeypatch):
    frozen = _write_frozen_input_fixture(tmp_path)
    calls = {"count": 0}

    def _nondeterministic_map(*, query_adata, **_kwargs):
        calls["count"] += 1
        label = "A" if calls["count"] % 2 else "B"
        frame = pd.DataFrame(
            {
                "reference_predicted_label": label,
                "reference_confidence": 1.0,
                "reference_distance": 0.1,
                "reference_assignment_status": ASSIGNMENT_STATUS_ACCEPTED,
                "reference_cell_type": label,
                "reference_ood": False,
            },
            index=query_adata.obs_names,
        )
        return frame, {"device_info": {"backend_residency": "host_cpu"}}

    monkeypatch.setattr(benchmark_module, "map_knn_reference", _nondeterministic_map)
    output_dir = tmp_path / "lane"
    with pytest.raises(RuntimeError, match="scientific output"):
        run_lane(
            device="cpu",
            input_dir=frozen,
            output_dir=output_dir,
            expected_input_manifest_sha256=compute_sha256(frozen / "input_manifest.json"),
        )
    assert json.loads((output_dir / "lane_receipt.json").read_text(encoding="utf-8"))["status"] == "failed"


@pytest.mark.parametrize(
    ("distance_step", "should_pass"),
    [
        (4e-7, True),
        (2e-6, False),
    ],
)
def test_within_lane_distance_repeatability_uses_predeclared_tolerance(
    tmp_path: Path,
    monkeypatch,
    distance_step: float,
    should_pass: bool,
):
    frozen = _write_frozen_input_fixture(tmp_path)
    calls = {"count": 0}

    def _distance_drift_map(*, query_adata, **_kwargs):
        calls["count"] += 1
        measured_repetition = max(calls["count"] - 2, 0)
        frame = pd.DataFrame(
            {
                "reference_predicted_label": "A",
                "reference_confidence": 1.0,
                "reference_distance": 0.1 + measured_repetition * distance_step,
                "reference_assignment_status": ASSIGNMENT_STATUS_ACCEPTED,
                "reference_cell_type": "A",
                "reference_ood": False,
            },
            index=query_adata.obs_names,
        )
        return frame, {"device_info": {"backend_residency": "host_cpu"}}

    monkeypatch.setattr(benchmark_module, "map_knn_reference", _distance_drift_map)
    output_dir = tmp_path / "lane"
    kwargs = {
        "device": "cpu",
        "input_dir": frozen,
        "output_dir": output_dir,
        "expected_input_manifest_sha256": compute_sha256(frozen / "input_manifest.json"),
    }
    if not should_pass:
        with pytest.raises(RuntimeError, match="exceed the predeclared"):
            run_lane(**kwargs)
        assert json.loads((output_dir / "lane_receipt.json").read_text(encoding="utf-8"))["status"] == "failed"
        return

    run_lane(**kwargs)
    _, metrics = verify_lane_output(
        output_dir,
        device="cpu",
        expected_input_manifest_sha256=kwargs["expected_input_manifest_sha256"],
    )
    repeatability = metrics["within_lane_repeatability"]
    assert repeatability["status"] == "pass"
    assert repeatability["contract"] == "exact_scientific_outputs_plus_bounded_float_distance"
    assert repeatability["raw_bitwise_identical"] is False
    assert repeatability["distance_max_abs_diff"] == pytest.approx(8e-7, abs=1e-12)


def test_failed_gpu_child_is_fail_not_promoted(tmp_path: Path, monkeypatch):
    frozen = _write_frozen_input_fixture(tmp_path)
    query_ids = [f"q{idx}" for idx in range(256)]
    query_labels = ["Microglia"] * 128 + ["A"] * 128
    query = ad.AnnData(
        X=sp.csr_matrix(np.ones((256, 2))),
        obs=pd.DataFrame({"cell_type": query_labels}, index=query_ids),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    reference = ad.AnnData(
        X=sp.csr_matrix(np.ones((2, 2))),
        obs=pd.DataFrame(
            {"sample": ["S1", "S2"], "cell_type": ["A", "B"]},
            index=["r1", "r2"],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    ref_cells = reference.obs_names.astype(str).tolist()
    ref_cells_sha = _sha256_bytes("\n".join(ref_cells).encode("utf-8") + b"\n")
    genes_sha = compute_sha256(frozen / "genes.txt")

    def _prepared(**_kwargs):
        return {
            "reference_adata": reference,
            "query_adata": query,
            "genes": ["G1", "G2"],
            "frozen_dir": frozen,
            "split_receipt": {
                "source_sha256": "a" * 64,
                "held_out_ood_label": "Microglia",
                "n_reference_cells": 2,
                "n_query_cells": 256,
                "n_genes": 2,
                "feature_selection": {
                    "schema_version": "1.0",
                    "selection_route": "reference_only_sparse_variance_v1",
                    "independently_query_blind": True,
                    "fit_population": "reference_cells_only",
                    "query_cells_consumed": False,
                    "n_fit_cells": 2,
                    "fit_cells_sha256": ref_cells_sha,
                    "n_input_features": 2,
                    "matrix_storage_class": "csr",
                    "score_method": "sparse_sample_variance",
                    "tie_break_rule": "lexicographical_gene_symbol_ascending",
                    "requested_cap": 2,
                    "effective_cap": 2,
                    "selected_count": 2,
                    "ordered_genes_sha256": genes_sha,
                },
            },
        }

    cpu_mapping = pd.DataFrame(
        {
            "reference_predicted_label": "A",
            "reference_confidence": 1.0,
            "reference_distance": 0.1,
            "reference_assignment_status": ASSIGNMENT_STATUS_ACCEPTED,
            "reference_cell_type": "A",
            "reference_ood": False,
        },
        index=query_ids,
    )
    cpu_metrics = {"median_wall_seconds": 0.1, "ru_maxrss_mb": 1.0, "metadata": {"device_info": {}}}

    monkeypatch.setattr(benchmark_module, "prepare_trevino_split", _prepared)
    monkeypatch.setattr(
        benchmark_module,
        "calibrate_reference_ood_threshold",
        lambda **_kwargs: (0.5, {"distance_threshold": 0.5}),
    )

    def _launch(*, device, output_dir, **_kwargs):
        output_dir.mkdir(parents=True, exist_ok=True)
        if device == "gpu":
            (output_dir / "lane_receipt.json").write_text(
                json.dumps({"status": "failed", "error_message": "forced GPU child failure"}), encoding="utf-8"
            )
            return {"exit_code": 1, "log_path": "forced-gpu.log", "vram_poll": {}}
        return {"exit_code": 0, "log_path": "forced-cpu.log", "vram_poll": {}}

    monkeypatch.setattr(benchmark_module, "_run_lane_process", _launch)
    monkeypatch.setattr(
        benchmark_module,
        "verify_lane_output",
        lambda _output_dir, *, device, **_kwargs: (cpu_mapping, cpu_metrics) if device == "cpu" else pytest.fail("GPU output must not be accepted"),
    )
    monkeypatch.setattr(benchmark_module, "_write_run_manifest", lambda **_kwargs: tmp_path / "manifest.json")
    result = benchmark_module.run_benchmark(
        project_root=tmp_path,
        run_id="2026-08-23T2359Z-0000000",
        allow_dirty=True,
    )
    assert result["technical_verdict"] == "FAIL_NOT_PROMOTED"
    assert "no CPU fallback" in result["failure_reason"]
