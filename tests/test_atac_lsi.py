"""Tests for workflow/modular/modules/atac_lsi.py (Wave 5 / US-W5-1)."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ctx(tmp_path: Path, adata):
    """Build a minimal PipelineContext-like mock for unit tests."""
    ctx = MagicMock()
    ctx.adata = adata
    ctx.figure_dir = tmp_path
    ctx.run_dir = tmp_path
    return ctx


def _make_adata(n_obs: int = 50, n_peaks: int = 200, seed: int = 42):
    """Build a minimal AnnData with sparse obsm['atac_peaks']."""
    import anndata

    rng = np.random.default_rng(seed)
    X = sp.random(n_obs, 100, density=0.1, format="csr", dtype=np.float32,
                  random_state=int(rng.integers(0, 2**31)))
    peak_mat = sp.random(n_obs, n_peaks, density=0.05, format="csr", dtype=np.float32,
                         random_state=int(rng.integers(0, 2**31)))
    adata = anndata.AnnData(X=X)
    adata.obsm["atac_peaks"] = peak_mat
    return adata


# ---------------------------------------------------------------------------
# Module-scoped fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def atac_lsi_mod():
    import importlib
    return importlib.import_module("workflow.modular.modules.atac_lsi")


@pytest.fixture(scope="module")
def ATACLSIModule(atac_lsi_mod):
    return atac_lsi_mod.ATACLSIModule


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_atac_lsi_reproduces_truncatedsvd_atol_1e6(tmp_path, ATACLSIModule):
    """Module X_lsi must match manual TruncatedSVD on same TF-IDF input to atol=1e-6."""
    from sklearn.decomposition import TruncatedSVD
    from workflow.modular.modules.atac_lsi import _tfidf_transform

    adata = _make_adata()
    peak_mat_copy = adata.obsm["atac_peaks"].copy()

    ctx = _make_ctx(tmp_path, adata)
    mod = ATACLSIModule()
    mod.run(ctx)

    # Reproduce manually using the same peak matrix copy
    tfidf = _tfidf_transform(peak_mat_copy)
    svd = TruncatedSVD(n_components=50, random_state=42)
    expected_full = svd.fit_transform(tfidf)
    expected = expected_full[:, 1:]  # drop first component

    np.testing.assert_allclose(
        adata.obsm["X_lsi"], expected, atol=1e-6,
        err_msg="X_lsi does not match manual TruncatedSVD to atol=1e-6",
    )


def test_atac_lsi_raises_on_missing_obsm(tmp_path, ATACLSIModule):
    """ModuleContractError raised when obsm['atac_peaks'] is absent."""
    from workflow.modular._contract_violation import ModuleContractError

    adata = _make_adata()
    del adata.obsm["atac_peaks"]
    ctx = _make_ctx(tmp_path, adata)

    mod = ATACLSIModule()
    with pytest.raises(ModuleContractError, match="atac_peaks"):
        mod.run(ctx)


def test_atac_lsi_raises_on_dense_obsm(tmp_path, ATACLSIModule):
    """ModuleContractError raised when obsm['atac_peaks'] is a dense ndarray."""
    from workflow.modular._contract_violation import ModuleContractError

    adata = _make_adata()
    adata.obsm["atac_peaks"] = adata.obsm["atac_peaks"].toarray()
    ctx = _make_ctx(tmp_path, adata)

    mod = ATACLSIModule()
    with pytest.raises(ModuleContractError, match="dense"):
        mod.run(ctx)


def test_atac_lsi_references_non_empty_high_if(atac_lsi_mod):
    """__references__ must contain at least one high-IF journal entry."""
    refs = atac_lsi_mod.__references__
    assert isinstance(refs, dict) and len(refs) >= 1, "__references__ is empty"
    high_if_journals = {"Cell", "Nat Methods", "Nat Genet"}
    found = any(
        entry.get("journal") in high_if_journals
        for entry in refs.values()
        if isinstance(entry, dict)
    )
    assert found, (
        f"No high-IF journal entry found in __references__. "
        f"Got: {[e.get('journal') for e in refs.values()]}"
    )


def test_atac_lsi_registry():
    """atac_lsi must be registered in _build_registry() after atac_qc and before multimodal_integration."""
    from workflow.modular.pipeline import _build_registry

    registry = _build_registry()
    assert "atac_lsi" in registry, "atac_lsi not found in pipeline registry"

    keys = list(registry.keys())
    idx_lsi = keys.index("atac_lsi")
    idx_qc = keys.index("atac_qc")
    idx_mm = keys.index("multimodal_integration")

    assert idx_qc < idx_lsi, "atac_lsi must appear after atac_qc in registry"
    assert idx_lsi < idx_mm, "atac_lsi must appear before multimodal_integration in registry"


def test_atac_lsi_variance_explained_plot(tmp_path, ATACLSIModule):
    """Variance-explained PNG must be written to figure_dir."""
    adata = _make_adata()
    ctx = _make_ctx(tmp_path, adata)

    mod = ATACLSIModule()
    mod.run(ctx)

    plot_path = tmp_path / "lsi_variance_explained.png"
    assert plot_path.exists(), f"PNG not found at {plot_path}"
    assert plot_path.stat().st_size > 0, "PNG is empty"
