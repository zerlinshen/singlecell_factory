"""Wave 2A integration smoke test.

Covers Phase 1B foundation (modality_registry, cross_modality_qc) + ATAC vertical slice
(atac_ingest, atac_module.R, bundle schema atac slot active).
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.io as sio
import scipy.sparse as sp


FACTORY_ROOT = Path(__file__).resolve().parent.parent
R_FACTORY_ROOT = FACTORY_ROOT.parent / "multiomics_r_factory"


@pytest.fixture
def synthetic_rna_adata():
    """RNA-only adata (no obsm modality slots)."""
    rng = np.random.default_rng(42)
    n_cells, n_genes = 200, 50
    X = sp.random(n_cells, n_genes, density=0.3, dtype=np.float32, random_state=42).tocsr()
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def synthetic_multimodal_adata(synthetic_rna_adata):
    """Synthetic RNA + ATAC + protein for cross-modality tests."""
    a = synthetic_rna_adata.copy()
    n = a.n_obs
    a.obsm["X_atac"] = np.random.default_rng(1).standard_normal((n, 10)).astype(np.float32)
    a.obsm["protein_clr"] = np.random.default_rng(2).standard_normal((n, 5)).astype(np.float32)
    return a


def _ctx_for(adata, tmp_path, **cfg_overrides):
    cfg = SimpleNamespace(**cfg_overrides)
    return SimpleNamespace(
        adata=adata, cfg=cfg, random_state=42,
        run_dir=tmp_path, metadata={}, status=lambda *a, **k: None,
    )


# --------------------------------------------------------------------------
# modality_registry — AC-9 prerequisite
# --------------------------------------------------------------------------

def test_modality_registry_rna_only(synthetic_rna_adata, tmp_path):
    from workflow.modular.modules.modality_registry import ModalityRegistryModule
    ctx = _ctx_for(synthetic_rna_adata, tmp_path)
    ModalityRegistryModule().run(ctx)
    assert synthetic_rna_adata.uns["modalities_present"] == ["singlecell_rna"]
    manifest = json.loads((tmp_path / "manifest" / "modalities.json").read_text())
    assert manifest["count"] == 1


def test_modality_registry_multimodal(synthetic_multimodal_adata, tmp_path):
    from workflow.modular.modules.modality_registry import ModalityRegistryModule
    ctx = _ctx_for(synthetic_multimodal_adata, tmp_path)
    ModalityRegistryModule().run(ctx)
    mods = synthetic_multimodal_adata.uns["modalities_present"]
    assert "singlecell_rna" in mods
    assert "atac" in mods
    assert "protein_adt" in mods
    assert len(mods) == 3


# --------------------------------------------------------------------------
# cross_modality_qc — tri-state status
# --------------------------------------------------------------------------

def test_cross_modality_qc_skip_single(synthetic_rna_adata, tmp_path):
    from workflow.modular.modules.modality_registry import ModalityRegistryModule
    from workflow.modular.modules.cross_modality_qc import CrossModalityQCModule
    ctx = _ctx_for(synthetic_rna_adata, tmp_path)
    ModalityRegistryModule().run(ctx)
    CrossModalityQCModule().run(ctx)
    status = json.loads((tmp_path / "qc" / "cross_modality_status.json").read_text())
    assert status["status"] == "skipped_single_modality"
    # No overlap file on skip
    assert not (tmp_path / "qc" / "cross_modality_overlap.json").exists()


def test_cross_modality_qc_runs_on_multimodal(synthetic_multimodal_adata, tmp_path):
    from workflow.modular.modules.modality_registry import ModalityRegistryModule
    from workflow.modular.modules.cross_modality_qc import CrossModalityQCModule
    ctx = _ctx_for(synthetic_multimodal_adata, tmp_path)
    ModalityRegistryModule().run(ctx)
    CrossModalityQCModule().run(ctx)
    status = json.loads((tmp_path / "qc" / "cross_modality_status.json").read_text())
    assert status["status"] in ("ran", "warn_divergence")
    overlap = json.loads((tmp_path / "qc" / "cross_modality_overlap.json").read_text())
    # Synthetic data has full barcode overlap (all cells in every modality) -> Jaccard ~1.0
    assert "atac" in overlap
    assert overlap["atac"]["jaccard"] > 0.5


# --------------------------------------------------------------------------
# atac_ingest — fail-soft + sparse-aware functional path
# --------------------------------------------------------------------------

def test_atac_ingest_skips_without_input(synthetic_rna_adata, tmp_path):
    from workflow.modular.modules.atac_ingest import ATACIngestModule
    ctx = _ctx_for(synthetic_rna_adata, tmp_path)
    ATACIngestModule().run(ctx)
    assert "X_atac" not in synthetic_rna_adata.obsm
    assert "atac_peaks" not in synthetic_rna_adata.obsm
    assert ctx.metadata["atac_ingest_status"].startswith("skipped")


def test_atac_ingest_sparse_lsi(synthetic_rna_adata, tmp_path):
    """End-to-end ATAC ingest on a synthetic sparse peak matrix."""
    from workflow.modular.modules.atac_ingest import ATACIngestModule

    n_cells = synthetic_rna_adata.n_obs
    n_peaks = 100
    rng = np.random.default_rng(42)
    # Binary-like sparse peak matrix (cells x peaks)
    peak_mat = sp.random(n_cells, n_peaks, density=0.1, format="csr", random_state=42)
    peak_mat.data = (peak_mat.data > 0).astype(np.float32)

    peak_mtx_path = tmp_path / "peaks.mtx"
    sio.mmwrite(str(peak_mtx_path), peak_mat)

    # Write peaks.bed
    peaks_bed = tmp_path / "peaks.bed"
    bed_rows = [f"chr1\t{i*100}\t{i*100+50}" for i in range(n_peaks)]
    peaks_bed.write_text("\n".join(bed_rows))

    ctx = _ctx_for(
        synthetic_rna_adata, tmp_path,
        atac_peak_matrix_path=str(peak_mtx_path),
        atac_peaks_bed_path=str(peaks_bed),
        atac_n_components=10,
    )
    ATACIngestModule().run(ctx)

    assert "atac_peaks" in synthetic_rna_adata.obsm
    assert sp.issparse(synthetic_rna_adata.obsm["atac_peaks"])
    assert synthetic_rna_adata.obsm["atac_peaks"].shape == (n_cells, n_peaks)
    assert "X_atac" in synthetic_rna_adata.obsm
    assert synthetic_rna_adata.obsm["X_atac"].shape == (n_cells, 10)
    assert "atac_peaks" in synthetic_rna_adata.uns
    assert "atac_var" in synthetic_rna_adata.uns
    assert ctx.metadata["atac_ingest_status"] == "ok"
    summary = json.loads((tmp_path / "atac_ingest" / "atac_summary.json").read_text())
    assert summary["n_peaks"] == n_peaks
    assert summary["peak_matrix_obsm_key"] == "atac_peaks"


# --------------------------------------------------------------------------
# Schema + R-side wiring
# --------------------------------------------------------------------------

def test_schema_parity_still_holds():
    sf = FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    rf = R_FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    assert sf.read_bytes() == rf.read_bytes(), "schema parity broken after Wave 2A edits"


def test_atac_slot_active():
    """Bundle schema's ATAC slot is active only because exporter+reader are wired."""
    import yaml
    sf = FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    doc = yaml.safe_load(sf.read_text())
    atac = doc["extensions_v22"]["atac"]
    assert atac["status"] == "active", f"atac status is {atac['status']}, expected active"
    assert "table" in atac and "lsi_path" in atac["table"]


def test_unfinished_modalities_reserved():
    """VDJ/Ribo/Hi-C remain reserved until they get exporter+reader vertical slices."""
    import yaml
    sf = FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    doc = yaml.safe_load(sf.read_text())
    for ext_name in ("vdj", "ribo", "hic"):
        assert doc["extensions_v22"][ext_name]["status"] == "reserved"


def test_r_extension_registry_includes_atac():
    io_bundle = R_FACTORY_ROOT / "R_bundle" / "io_bundle.R"
    src = io_bundle.read_text()
    assert '"atac"' in src
    # confirm it appears in the KNOWN_BUNDLE_EXTENSIONS line
    for line in src.splitlines():
        if line.strip().startswith("KNOWN_BUNDLE_EXTENSIONS"):
            assert '"atac"' in line, f"atac not in registry line: {line}"
            return
    raise AssertionError("KNOWN_BUNDLE_EXTENSIONS line not found")


def test_r_atac_module_present():
    atac_r = R_FACTORY_ROOT / "R" / "atac_module.R"
    assert atac_r.exists()
    src = atac_r.read_text()
    for fn in ("load_atac_extension", "plot_atac_lsi", "plot_atac_peak_count_dist"):
        assert fn in src, f"{fn} missing from atac_module.R"


# --------------------------------------------------------------------------
# Registry: all 3 new modules present, no MUTATING, no cycles
# --------------------------------------------------------------------------

def test_new_modules_in_registry():
    from workflow.modular.pipeline import _build_registry, MUTATING_MODULES, _resolve_execution_order
    from workflow.modular.module_catalog import MANDATORY_MODULES
    reg = _build_registry()
    for name in ("modality_registry", "cross_modality_qc", "atac_ingest"):
        assert name in reg, f"{name} not in registry"
        assert name not in MUTATING_MODULES

    # DAG resolves with all 3 new modules active alongside mandatories
    optional = ["modality_registry", "cross_modality_qc", "atac_ingest"]
    order = _resolve_execution_order(list(MANDATORY_MODULES), optional)
    for name in optional:
        assert name in order, f"{name} not in resolved order"
    # cross_modality_qc must come after modality_registry (declared dep)
    assert order.index("cross_modality_qc") > order.index("modality_registry")
