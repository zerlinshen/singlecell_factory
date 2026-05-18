"""Tests for trajectory module contract gate (US-W5-8-CONTRACT)."""
import numpy as np
import pytest
import anndata as ad

from workflow.modular._contract_violation import ModuleContractError
from workflow.modular.modules.trajectory import TrajectoryModule


def _make_adata(**obsm_kwargs):
    adata = ad.AnnData(np.zeros((5, 3)))
    adata.obs["leiden"] = ["0", "0", "1", "1", "2"]
    for key, val in obsm_kwargs.items():
        adata.obsm[key] = val
    return adata


class _MinimalCtx:
    def __init__(self, adata):
        self.adata = adata
        self.metadata = {}

        class _Cfg:
            trajectory_root_cluster = None

        self.cfg = _Cfg()

        import tempfile, pathlib
        self._tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self._tmp.name)
        self.figure_dir = base / "figures"
        self.table_dir = base / "tables"
        self.figure_dir.mkdir()
        self.table_dir.mkdir()

    def __del__(self):
        try:
            self._tmp.cleanup()
        except Exception:
            pass


def test_requires_keys_uses_obsm_any_of():
    mod = TrajectoryModule()
    assert "obsm_any_of" in mod.requires_keys
    assert "X_wnn" in mod.requires_keys["obsm_any_of"]
    assert "X_pca" in mod.requires_keys["obsm_any_of"]
    assert "obsm" not in mod.requires_keys or "X_umap" not in mod.requires_keys.get("obsm", [])


def test_raises_module_contract_error_when_both_absent():
    adata = _make_adata()
    ctx = _MinimalCtx(adata)
    with pytest.raises(ModuleContractError, match="X_wnn"):
        TrajectoryModule().run(ctx)


def _patch_scanpy_tl(monkeypatch):
    import scanpy as sc

    def _stub_diffmap(a, **kw):
        a.obsm["X_diffmap"] = np.zeros((a.n_obs, 2), dtype=np.float32)

    def _stub_dpt(a, **kw):
        a.obs["dpt_pseudotime"] = np.linspace(0, 1, a.n_obs)

    monkeypatch.setattr(sc.tl, "diffmap", _stub_diffmap, raising=False)
    monkeypatch.setattr(sc.tl, "dpt", _stub_dpt, raising=False)
    monkeypatch.setattr(sc.tl, "paga", lambda a, **kw: None, raising=False)


def test_x_pca_used_when_only_pca_present(monkeypatch):
    pca_vals = np.random.default_rng(0).random((5, 2)).astype(np.float32)
    adata = _make_adata(X_pca=pca_vals)
    ctx = _MinimalCtx(adata)
    _patch_scanpy_tl(monkeypatch)
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_umap", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_paga", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_gene_trends", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_density", staticmethod(lambda a, c: None))

    TrajectoryModule().run(ctx)
    assert ctx.metadata["trajectory_embedding_key"] == "X_pca"


def test_x_wnn_preferred_over_x_pca(monkeypatch):
    wnn_vals = np.ones((5, 2), dtype=np.float32) * 99.0
    pca_vals = np.zeros((5, 2), dtype=np.float32)
    adata = _make_adata(X_wnn=wnn_vals, X_pca=pca_vals)
    ctx = _MinimalCtx(adata)
    _patch_scanpy_tl(monkeypatch)
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_umap", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_paga", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_gene_trends", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_density", staticmethod(lambda a, c: None))

    TrajectoryModule().run(ctx)
    assert ctx.metadata["trajectory_embedding_key"] == "X_wnn"
