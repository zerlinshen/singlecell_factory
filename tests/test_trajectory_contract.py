"""Tests for trajectory module contract gate (US-W5-8-CONTRACT)."""
import numpy as np
import pandas as pd
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
    def __init__(self, adata, *, root_cluster=None, root_justification=None):
        self.adata = adata
        self.metadata = {}

        class _Cfg:
            trajectory_root_cluster = root_cluster
            trajectory_root_justification = root_justification

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
    assert ctx.metadata["trajectory_claimable"] is False
    assert adata.uns["trajectory"]["inference_class"] == "exploratory_unrooted_dpt"


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


def test_explicit_biologically_justified_root_is_claimable(monkeypatch):
    pca_vals = np.random.default_rng(1).random((5, 2)).astype(np.float32)
    adata = _make_adata(X_pca=pca_vals)
    ctx = _MinimalCtx(
        adata,
        root_cluster="1",
        root_justification="Cluster 1 expresses validated early progenitor markers.",
    )
    _patch_scanpy_tl(monkeypatch)
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_umap", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_paga", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_gene_trends", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_density", staticmethod(lambda a, c: None))

    TrajectoryModule().run(ctx)

    provenance = adata.uns["trajectory"]
    assert provenance["claimable"] is True
    assert provenance["root_cluster"] == "1"
    assert provenance["root_cluster_found"] is True
    assert provenance["root_justification"].startswith("Cluster 1")
    assert provenance["inference_class"] == "biologically_rooted_dpt"


@pytest.mark.parametrize(
    ("root_cluster", "root_justification"),
    [
        (None, None),
        ("1", None),
        ("missing", "Expected early state"),
    ],
)
def test_missing_or_invalid_biological_root_is_non_claimable(
    monkeypatch, root_cluster, root_justification
):
    pca_vals = np.random.default_rng(2).random((5, 2)).astype(np.float32)
    adata = _make_adata(X_pca=pca_vals)
    ctx = _MinimalCtx(
        adata,
        root_cluster=root_cluster,
        root_justification=root_justification,
    )
    _patch_scanpy_tl(monkeypatch)
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_umap", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_paga", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_gene_trends", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_density", staticmethod(lambda a, c: None))

    TrajectoryModule().run(ctx)

    provenance = adata.uns["trajectory"]
    assert provenance["claimable"] is False
    assert provenance["claim_status"] == "non_claimable_missing_biologically_justified_root"
    assert provenance["root_selection_method"] == "exploratory_diffmap_minimum_fallback"
    pseudotime_table = pd.read_csv(ctx.table_dir / "dpt_pseudotime.csv")
    assert not pseudotime_table["trajectory_claimable"].any()
    assert set(pseudotime_table["trajectory_claim_status"]) == {
        "non_claimable_missing_biologically_justified_root"
    }
