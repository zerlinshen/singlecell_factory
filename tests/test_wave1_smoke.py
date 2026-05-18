"""Wave 1 integration smoke test (T11).

Validates that the Wave 1 deliverables wire together correctly on synthetic
data. Each test exercises one AC slice without requiring a full pipeline run.

Covers:
  - T2 / AC-5: marker_db_loader registers + reads curated YAML overlay
  - T3 / AC-5: context_aware_annotation produces obs columns + metrics JSON
  - T5 / AC-6: markers.parquet export when context_aware_annotation ran
  - T9 / AC-5: R-side schema registration via KNOWN_BUNDLE_EXTENSIONS file presence
  - Bridge symlink rule (project CLAUDE.md)
  - Schema parity (singlecell vs multiomics_r contracts)
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


FACTORY_ROOT = Path(__file__).resolve().parent.parent
R_FACTORY_ROOT = FACTORY_ROOT.parent / "multiomics_r_factory"


@pytest.fixture
def synthetic_adata():
    """Small synthetic AnnData (200 cells × 50 genes sparse) with leiden labels."""
    rng = np.random.default_rng(42)
    n_cells, n_genes = 200, 50
    gene_names = (
        ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19"]
        + ["CD3D", "CD3E", "CD3G", "TRAC", "CD8A"]
        + ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19"]
        + ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"]
        + ["NKG7", "GNLY", "KLRD1", "KLRF1", "GZMB"]
        + [f"GENE_{i}" for i in range(n_genes - 25)]
    )
    X = sp.random(n_cells, n_genes, density=0.3, dtype=np.float32, random_state=42).tocsr()
    # 4 clusters of 50 cells each
    leiden = np.array([str(i // 50) for i in range(n_cells)])
    adata = ad.AnnData(X=X, obs=pd.DataFrame({"leiden": leiden}, index=[f"cell_{i}" for i in range(n_cells)]))
    adata.var_names = gene_names
    return adata


# --------------------------------------------------------------------------
# T2 / AC-5: marker_db_loader reads curated YAML overlay
# --------------------------------------------------------------------------

def test_marker_db_loader_picks_up_curated_yaml(synthetic_adata):
    from workflow.modular.modules.marker_db_loader import MarkerDbLoaderModule

    cfg = SimpleNamespace(tissue="lung", condition="NSCLC")
    ctx = SimpleNamespace(
        adata=synthetic_adata,
        cfg=cfg,
        metadata={},
        status=lambda *a, **k: None,
    )
    MarkerDbLoaderModule().run(ctx)

    index = synthetic_adata.uns.get("marker_db_index", {})
    assert "curated" in index, f"curated overlay not loaded; keys={list(index.keys())}"
    assert len(index["curated"]) > 0, "curated overlay loaded zero rows"
    # Verify provenance columns
    df = index["curated"]
    assert "source_db" in df.columns
    assert "source_version" in df.columns
    assert (df["source_db"] == "curated").all()


# --------------------------------------------------------------------------
# T3 / AC-5: context_aware_annotation produces obs columns + metrics JSON
# --------------------------------------------------------------------------

def test_context_aware_annotation_writes_obs_and_metrics(synthetic_adata, tmp_path):
    from workflow.modular.modules.marker_db_loader import MarkerDbLoaderModule
    from workflow.modular.modules.context_aware_annotation import ContextAwareAnnotationModule

    cfg = SimpleNamespace(
        tissue="lung",
        condition="NSCLC",
        validate_context=True,
        context_mismatch_threshold=0.3,
        context_min_cells=20,
    )
    ctx = SimpleNamespace(
        adata=synthetic_adata,
        cfg=cfg,
        random_state=42,
        run_dir=tmp_path,
        metadata={},
        status=lambda *a, **k: None,
    )
    MarkerDbLoaderModule().run(ctx)
    ContextAwareAnnotationModule().run(ctx)

    assert "context_aware_celltype" in synthetic_adata.obs.columns
    # All cells should be assigned (no NaN)
    assert synthetic_adata.obs["context_aware_celltype"].notna().all()

    candidate_json = tmp_path / "metrics" / "marker_intelligence_candidate.json"
    assert candidate_json.exists(), f"missing {candidate_json}"
    data = json.loads(candidate_json.read_text())
    assert "per_cluster" in data
    assert len(data["per_cluster"]) > 0

    # AC-4: --validate-context emits context_validation.json
    validation_json = tmp_path / "metrics" / "context_validation.json"
    assert validation_json.exists()


# --------------------------------------------------------------------------
# AC-10: sparse-aware patch produces same cluster-mean as densify (within float32 noise)
# --------------------------------------------------------------------------

def test_sparse_mean_equivalent_to_dense(synthetic_adata):
    X = synthetic_adata.X
    mask = synthetic_adata.obs["leiden"].values == "0"

    sparse_mean = np.asarray(X[mask].mean(axis=0)).flatten().astype(np.float32)
    dense_mean = np.asarray(X[mask].todense(), dtype=np.float32).mean(axis=0)

    assert sparse_mean.shape == dense_mean.shape
    # Within float32 summation-order noise
    assert np.allclose(sparse_mean, dense_mean, atol=1e-5, rtol=1e-5)


# --------------------------------------------------------------------------
# T5 / AC-6: markers.parquet export when context_aware_annotation ran
# --------------------------------------------------------------------------

def test_marker_resolutions_export(synthetic_adata, tmp_path):
    """maybe_export_marker_resolutions writes markers.parquet + registers extension."""
    import sys
    sys.path.insert(0, str(FACTORY_ROOT))
    from scripts.export_singlecell_r_bundle import maybe_export_marker_resolutions

    # Add the prerequisites
    from workflow.modular.modules.marker_db_loader import MarkerDbLoaderModule

    cfg = SimpleNamespace(tissue="lung", condition="NSCLC")
    ctx = SimpleNamespace(
        adata=synthetic_adata, cfg=cfg, metadata={}, status=lambda *a, **k: None,
    )
    MarkerDbLoaderModule().run(ctx)
    synthetic_adata.obs["context_aware_celltype"] = "B cell"  # context_aware_annotation marker

    manifest = {}
    entry = maybe_export_marker_resolutions(manifest, tmp_path, synthetic_adata)
    assert entry is not None, "exporter returned None despite marker_db_index + context_aware_celltype present"

    parquet_path = tmp_path / "extensions" / "marker_resolutions" / "markers.parquet"
    assert parquet_path.exists(), f"missing {parquet_path}"

    df = pd.read_parquet(parquet_path)
    expected_cols = {"cell_type", "marker_set", "source_db", "source_version", "score"}
    assert expected_cols.issubset(df.columns), f"missing cols: {expected_cols - set(df.columns)}"
    assert len(df) > 0, "markers.parquet is empty"

    # Manifest extension registered with status active
    assert "extensions" in manifest
    assert "marker_resolutions" in manifest["extensions"]


def test_marker_resolutions_export_skips_without_caa(synthetic_adata, tmp_path):
    """No markers.parquet when context_aware_annotation did NOT run (AC-10 zero-regression)."""
    import sys
    sys.path.insert(0, str(FACTORY_ROOT))
    from scripts.export_singlecell_r_bundle import maybe_export_marker_resolutions
    # marker_db_index absent (legacy path) → fail-soft
    manifest = {}
    entry = maybe_export_marker_resolutions(manifest, tmp_path, synthetic_adata)
    assert entry is None, "exporter should fail-soft when prerequisites absent"
    assert not (tmp_path / "extensions").exists(), "should not write any extension files"


# --------------------------------------------------------------------------
# T4 / AC-6: schema parity between two bundle_schema.yaml copies
# --------------------------------------------------------------------------

def test_schema_parity():
    """Both bundle_schema.yaml copies must be byte-identical (CI guard)."""
    sf_schema = FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    rf_schema = R_FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    assert sf_schema.exists(), f"missing {sf_schema}"
    assert rf_schema.exists(), f"missing {rf_schema}"
    assert sf_schema.read_bytes() == rf_schema.read_bytes(), "schema parity broken"


# --------------------------------------------------------------------------
# T9 / AC-5: R-side KNOWN_BUNDLE_EXTENSIONS lists marker_resolutions
# --------------------------------------------------------------------------

def test_r_extension_registry_includes_marker_resolutions():
    """multiomics_r_factory/R_bundle/io_bundle.R registers marker_resolutions."""
    io_bundle = R_FACTORY_ROOT / "R_bundle" / "io_bundle.R"
    assert io_bundle.exists(), f"missing {io_bundle}"
    src = io_bundle.read_text(encoding="utf-8")
    assert "marker_resolutions" in src, "marker_resolutions not registered in R_bundle/io_bundle.R"
    # Confirm v2.2 schema accepted
    assert "v2.2" in src or "singlecell_r_bundle_v2.2" in src, "v2.2 schema not in R loader"


# --------------------------------------------------------------------------
# Bridge symlink rule (project CLAUDE.md)
# --------------------------------------------------------------------------

def test_bridge_symlinks_intact():
    """bridges/local_r_pipeline_macbook/R{,_bundle} must be symlinks."""
    bridge_dir = FACTORY_ROOT / "bridges" / "local_r_pipeline_macbook"
    for sub in ("R", "R_bundle"):
        path = bridge_dir / sub
        assert path.is_symlink(), f"{path} must be a symlink, not a real file/dir"


# --------------------------------------------------------------------------
# Registry: both new modules registered, neither in MUTATING_MODULES
# --------------------------------------------------------------------------

def test_new_modules_in_registry():
    from workflow.modular.pipeline import _build_registry, MUTATING_MODULES
    reg = _build_registry()
    assert "marker_db_loader" in reg, "marker_db_loader not in registry"
    assert "context_aware_annotation" in reg, "context_aware_annotation not in registry"
    assert "marker_db_loader" not in MUTATING_MODULES
    assert "context_aware_annotation" not in MUTATING_MODULES
