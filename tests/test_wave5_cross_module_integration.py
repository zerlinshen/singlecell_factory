"""Wave 5 cross-module integration gate (US-W5-CI).

Exercises the full atac_lsi → multimodal_integration(WNN) → trajectory → peak_to_gene
sub-DAG on the Trevino synthetic fixture. All tests are marked ``wave5_blocking``
and MUST pass before any real Trevino run (US-W5-4).

v6.3 schema: ATAC peaks live in adata.obsm["atac_peaks"] (NOT layers).
"""
from __future__ import annotations

import tempfile
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

from workflow.modular._contract_violation import ModuleContractError

pytestmark = pytest.mark.wave5_blocking


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ctx(tmp_path: Path, adata: ad.AnnData):
    """Minimal PipelineContext-compatible namespace for module.run() calls."""

    class _Cfg:
        trajectory_root_cluster = None
        multimodal_integration = None

    figure_dir = tmp_path / "figures"
    table_dir = tmp_path / "tables"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    ctx = SimpleNamespace(
        adata=adata,
        cfg=_Cfg(),
        metadata={},
        figure_dir=figure_dir,
        table_dir=table_dir,
        run_dir=tmp_path,
        status=lambda *args, **kwargs: None,
        random_state=42,
    )
    return ctx


def _make_synthetic_adata(n_obs: int = 200, n_rna: int = 300, n_atac: int = 1000) -> ad.AnnData:
    """Small synthetic AnnData matching the v6.3 schema (obsm["atac_peaks"])."""
    from scripts.dev.load_trevino_2021 import load_synthetic
    adata = load_synthetic(n_obs=n_obs, n_rna=n_rna, n_atac=n_atac, seed=0)
    return adata


# ---------------------------------------------------------------------------
# Fixture: session-scoped synthetic adata (cheap to build once)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def synth_adata():
    return _make_synthetic_adata()


# ---------------------------------------------------------------------------
# Test 1: happy-path sub-DAG
# atac_lsi → multimodal_integration(WNN skip) → trajectory(X_pca fallback) → peak_to_gene
# ---------------------------------------------------------------------------

def test_happy_path_subdag(tmp_path):
    """Full sub-DAG runs without error on the synthetic fixture."""
    import scanpy as sc
    from workflow.modular.modules.atac_lsi import ATACLSIModule
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    from workflow.modular.modules.trajectory import TrajectoryModule
    from workflow.modular.modules.peak_to_gene import compute_peak_gene_linkages

    adata = _make_synthetic_adata()
    ctx = _make_ctx(tmp_path, adata)

    # Step 1: atac_lsi — requires obsm["atac_peaks"], provides obsm["X_lsi"]
    ctx.figure_dir = tmp_path / "atac_lsi"
    ctx.figure_dir.mkdir(parents=True, exist_ok=True)
    ATACLSIModule().run(ctx)
    assert "X_lsi" in adata.obsm, "atac_lsi must write obsm['X_lsi']"
    assert adata.obsm["X_lsi"].shape == (adata.n_obs, 49)

    # Step 2: inject X_pca placeholder (WNN integration needs it as RNA latent)
    rng = np.random.RandomState(0)
    adata.obsm["X_pca"] = rng.randn(adata.n_obs, 49).astype(np.float32)

    # Step 3: multimodal_integration with SC_MULTIMODAL_SECOND_OBSM=X_lsi
    # WNN engine will skip cleanly (no Rscript in CI); we validate the contract path.
    # We use engine="wnn" + env override so the module reaches the key-resolution logic.
    # Since Rscript is absent, it will skip with status "skipped_no_rscript".
    cfg = MultimodalIntegrationConfig(engine="wnn", _explicitly_set=False)
    mm_module = MultimodalIntegrationModule(cfg)
    # Inject X_lsi as second modality key via env (simulates SC_MULTIMODAL_SECOND_OBSM=X_lsi)
    # We do NOT monkeypatch here — we manually put X_lsi in obsm so the module can find it.
    # The explicit-set path is tested in the contract-raise tests below.
    mm_module.run(ctx)
    status = adata.uns.get("multimodal_status", {})
    assert status.get("engine") == "wnn"
    # WNN must have either skipped (no Rscript/no driver in CI) or succeeded
    assert status.get("status") in {"skipped", "ok"}, f"unexpected WNN status: {status}"

    # Step 4: inject a synthetic X_wnn so trajectory has a joint embedding
    adata.obsm["X_wnn"] = rng.randn(adata.n_obs, 2).astype(np.float32)

    # Step 5: trajectory — requires X_wnn or X_pca; computes diffmap + DPT
    sc.pp.neighbors(adata, use_rep="X_wnn", n_neighbors=5)
    # UMAP needed by trajectory's _plot_pseudotime_umap
    sc.tl.umap(adata)

    ctx.figure_dir = tmp_path / "trajectory"
    ctx.table_dir = tmp_path / "trajectory"
    ctx.figure_dir.mkdir(parents=True, exist_ok=True)

    TrajectoryModule().run(ctx)
    assert "dpt_pseudotime" in adata.obs, "trajectory must write obs['dpt_pseudotime']"

    # Step 6: peak_to_gene linkage (n_perms=10 for fast CI)
    compute_peak_gene_linkages(adata, n_perms=10)
    assert "peak_to_gene_top1000" in adata.uns
    assert "peak_to_gene_linkages" in adata.uns


# ---------------------------------------------------------------------------
# Test 2: atac_lsi contract raise — missing obsm["atac_peaks"]
# ---------------------------------------------------------------------------

def test_atac_lsi_raises_on_missing_atac_peaks(tmp_path):
    """atac_lsi must raise ModuleContractError when obsm['atac_peaks'] is absent."""
    from workflow.modular.modules.atac_lsi import ATACLSIModule
    from tests.fixtures.trevino_synthetic import make_missing_layer_variant

    adata = _make_synthetic_adata()
    adata_missing = make_missing_layer_variant(adata)
    ctx = _make_ctx(tmp_path, adata_missing)
    ctx.figure_dir = tmp_path

    with pytest.raises(ModuleContractError, match="atac_peaks"):
        ATACLSIModule().run(ctx)


# ---------------------------------------------------------------------------
# Test 3: atac_lsi contract raise — dense obsm["atac_peaks"]
# ---------------------------------------------------------------------------

def test_atac_lsi_raises_on_dense_atac_peaks(tmp_path):
    """atac_lsi must raise ModuleContractError when obsm['atac_peaks'] is dense."""
    from workflow.modular.modules.atac_lsi import ATACLSIModule
    from tests.fixtures.trevino_synthetic import make_atac_dense_variant

    adata = _make_synthetic_adata()
    adata_dense = make_atac_dense_variant(adata)
    ctx = _make_ctx(tmp_path, adata_dense)
    ctx.figure_dir = tmp_path

    with pytest.raises(ModuleContractError, match="dense"):
        ATACLSIModule().run(ctx)


# ---------------------------------------------------------------------------
# Test 4: multimodal_integration contract raise — explicit second key missing
# ---------------------------------------------------------------------------

def test_multimodal_raises_on_explicit_missing_second_key(tmp_path, monkeypatch):
    """multimodal_integration must raise ModuleContractError when SC_MULTIMODAL_SECOND_OBSM
    is set explicitly but the key is absent from adata.obsm."""
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )

    adata = _make_synthetic_adata()
    # Inject X_pca so RNA latent is satisfied; do NOT inject X_lsi
    adata.obsm["X_pca"] = np.random.RandomState(0).randn(adata.n_obs, 49).astype(np.float32)

    monkeypatch.setenv("SC_MULTIMODAL_ENGINE", "wnn")
    monkeypatch.setenv("SC_MULTIMODAL_SECOND_OBSM", "X_lsi")

    ctx = _make_ctx(tmp_path, adata)
    cfg = MultimodalIntegrationConfig()  # env vars drive engine + second key

    with pytest.raises(ModuleContractError, match="explicitly-configured second obsm key"):
        MultimodalIntegrationModule(cfg).run(ctx)


# ---------------------------------------------------------------------------
# Test 5: trajectory contract raise — neither X_wnn nor X_pca present
# ---------------------------------------------------------------------------

def test_trajectory_raises_on_missing_embeddings(tmp_path):
    """trajectory must raise ModuleContractError when both X_wnn and X_pca are absent."""
    from workflow.modular.modules.trajectory import TrajectoryModule

    adata = _make_synthetic_adata()
    # Confirm neither embedding is present
    assert "X_wnn" not in adata.obsm
    assert "X_pca" not in adata.obsm

    ctx = _make_ctx(tmp_path, adata)

    with pytest.raises(ModuleContractError, match="X_wnn"):
        TrajectoryModule().run(ctx)


# ---------------------------------------------------------------------------
# Test 6: peak_to_gene contract raise — dense adata.X
# ---------------------------------------------------------------------------

def test_peak_to_gene_raises_on_dense_x(tmp_path):
    """compute_peak_gene_linkages must raise ModuleContractError when adata.X is dense."""
    from workflow.modular.modules.peak_to_gene import compute_peak_gene_linkages

    adata = _make_synthetic_adata()
    # Densify adata.X — violates sparse-axis contract
    adata.X = np.asarray(adata.X.todense())

    with pytest.raises(ModuleContractError, match="sparse"):
        compute_peak_gene_linkages(adata, n_perms=10)


# ---------------------------------------------------------------------------
# Test 7: peak_to_gene contract raise — dense obsm["atac_peaks"]
# ---------------------------------------------------------------------------

def test_peak_to_gene_raises_on_dense_atac_peaks(tmp_path):
    """compute_peak_gene_linkages must raise ModuleContractError when obsm['atac_peaks'] is dense."""
    from workflow.modular.modules.peak_to_gene import compute_peak_gene_linkages
    from tests.fixtures.trevino_synthetic import make_atac_dense_variant

    adata = _make_synthetic_adata()
    adata_dense = make_atac_dense_variant(adata)

    with pytest.raises(ModuleContractError, match="sparse"):
        compute_peak_gene_linkages(adata_dense, n_perms=10)


# ---------------------------------------------------------------------------
# Test 8: obsm["X_lsi"] shape contract after atac_lsi run
# ---------------------------------------------------------------------------

def test_atac_lsi_output_shape(tmp_path, synth_adata):
    """atac_lsi output obsm['X_lsi'] must be (n_obs, 49) — first SVD component dropped."""
    from workflow.modular.modules.atac_lsi import ATACLSIModule

    adata = synth_adata.copy()
    ctx = _make_ctx(tmp_path, adata)
    ctx.figure_dir = tmp_path

    ATACLSIModule().run(ctx)

    assert adata.obsm["X_lsi"].shape == (adata.n_obs, 49), (
        f"Expected X_lsi shape ({adata.n_obs}, 49), got {adata.obsm['X_lsi'].shape}"
    )
    assert "lsi_variance_explained" in adata.uns
    assert len(adata.uns["lsi_variance_explained"]) == 50


# ---------------------------------------------------------------------------
# Test 9: multimodal env-override happy path (X_lsi present + engine=wnn)
# ---------------------------------------------------------------------------

def test_multimodal_env_override_with_xlsi_present(tmp_path, monkeypatch):
    """When SC_MULTIMODAL_SECOND_OBSM=X_lsi and X_lsi is present, no contract error raised."""
    from workflow.modular.modules.atac_lsi import ATACLSIModule
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )

    adata = _make_synthetic_adata()
    ctx = _make_ctx(tmp_path, adata)
    ctx.figure_dir = tmp_path

    # Run atac_lsi to produce X_lsi
    ATACLSIModule().run(ctx)
    assert "X_lsi" in adata.obsm

    # Inject X_pca placeholder
    adata.obsm["X_pca"] = np.random.RandomState(1).randn(adata.n_obs, 49).astype(np.float32)

    monkeypatch.setenv("SC_MULTIMODAL_ENGINE", "wnn")
    monkeypatch.setenv("SC_MULTIMODAL_SECOND_OBSM", "X_lsi")

    cfg = MultimodalIntegrationConfig()
    try:
        MultimodalIntegrationModule(cfg).run(ctx)
    except ModuleContractError:
        pytest.fail("ModuleContractError raised even though X_lsi is present in obsm")

    status = adata.uns.get("multimodal_status", {})
    assert status.get("engine") == "wnn", f"Expected engine=wnn, got {status}"
