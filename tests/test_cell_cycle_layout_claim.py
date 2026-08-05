"""Cell-cycle regress-after-layout honesty (Wave-2 design C)."""
from __future__ import annotations

import numpy as np
import anndata as ad
import pytest

from workflow.modular.modules.cell_cycle import CellCycleModule, S_GENES, G2M_GENES


class _Cfg:
    regress_cell_cycle = True


class _Ctx:
    def __init__(self, adata, regress=True):
        self.adata = adata
        self.cfg = _Cfg()
        self.cfg.regress_cell_cycle = regress
        self.metadata = {"clustering_claim_status": "claimable"}
        self.random_state = 0
        import tempfile
        from pathlib import Path

        self._tmp = tempfile.TemporaryDirectory()
        base = Path(self._tmp.name)
        self.figure_dir = base / "figures"
        self.table_dir = base / "tables"
        self.figure_dir.mkdir()
        self.table_dir.mkdir()
        self._status = []

    def status(self, *args, **kwargs):
        self._status.append((args, kwargs))

    def __del__(self):
        try:
            self._tmp.cleanup()
        except Exception:
            pass


def _adata_with_cc_genes():
    genes = list(dict.fromkeys(S_GENES[:15] + G2M_GENES[:15] + ["OTHER1", "OTHER2"]))
    n = 20
    rng = np.random.default_rng(0)
    a = ad.AnnData(rng.poisson(2.0, size=(n, len(genes))).astype(np.float32))
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obsm["X_umap"] = rng.normal(size=(n, 2)).astype(np.float32)
    return a


def test_regress_after_layout_stamps_exploratory(monkeypatch):
    import scanpy as sc

    def fake_score(adata, s_genes=None, g2m_genes=None, **kw):
        n = adata.n_obs
        adata.obs["S_score"] = np.linspace(0, 1, n)
        adata.obs["G2M_score"] = np.linspace(1, 0, n)
        adata.obs["phase"] = ["G1"] * n

    monkeypatch.setattr(sc.tl, "score_genes_cell_cycle", fake_score, raising=False)
    monkeypatch.setattr(sc.pp, "regress_out", lambda *a, **k: None, raising=False)
    monkeypatch.setattr(sc.pl, "umap", lambda *a, **k: None, raising=False)

    ctx = _Ctx(_adata_with_cc_genes(), regress=True)
    ctx.metadata["batch_risk"] = {"clustering_claim_status": "claimable"}
    CellCycleModule().run(ctx)
    assert ctx.metadata["cell_cycle_regressed_after_layout"] is True
    assert ctx.metadata["cell_cycle_layout_claimable"] is False
    assert ctx.metadata["clustering_claim_status"] == "exploratory"
    assert ctx.metadata["batch_risk"]["clustering_claim_status"] == "exploratory"
    assert (
        ctx.metadata["batch_risk"]["clustering_claim_reason"]
        == "exploratory_cell_cycle_regressed_after_layout"
    )
    assert ctx.adata.uns["cell_cycle"]["regressed_after_layout"] is True


def test_no_regress_does_not_force_layout_exploratory(monkeypatch):
    import scanpy as sc

    def fake_score(adata, s_genes=None, g2m_genes=None, **kw):
        n = adata.n_obs
        adata.obs["S_score"] = np.linspace(0, 1, n)
        adata.obs["G2M_score"] = np.linspace(1, 0, n)
        adata.obs["phase"] = ["G1"] * n

    monkeypatch.setattr(sc.tl, "score_genes_cell_cycle", fake_score, raising=False)
    monkeypatch.setattr(sc.pl, "umap", lambda *a, **k: None, raising=False)

    ctx = _Ctx(_adata_with_cc_genes(), regress=False)
    CellCycleModule().run(ctx)
    assert ctx.metadata["cell_cycle_regressed_after_layout"] is False
    assert ctx.metadata.get("clustering_claim_status") == "claimable"


def test_regress_before_layout_is_applied_before_clustering(monkeypatch):
    """Wave-3 option A: no X_pca/umap/leiden yet → embedding-affecting regress."""
    import scanpy as sc

    def fake_score(adata, s_genes=None, g2m_genes=None, **kw):
        n = adata.n_obs
        adata.obs["S_score"] = np.linspace(0, 1, n)
        adata.obs["G2M_score"] = np.linspace(1, 0, n)
        adata.obs["phase"] = ["G1"] * n

    called = {"regress": 0}

    def fake_regress(adata, keys, **kw):
        called["regress"] += 1

    monkeypatch.setattr(sc.tl, "score_genes_cell_cycle", fake_score, raising=False)
    monkeypatch.setattr(sc.pp, "regress_out", fake_regress, raising=False)

    adata = _adata_with_cc_genes()
    # Pre-layout: strip embeddings that _adata_with_cc_genes may have set
    if "X_umap" in adata.obsm:
        del adata.obsm["X_umap"]
    ctx = _Ctx(adata, regress=True)
    CellCycleModule().run(ctx)
    assert called["regress"] == 1
    assert ctx.metadata["cell_cycle_regressed_after_layout"] is False
    assert ctx.metadata["cell_cycle_regress_out_embedding_effect"] == (
        "applied_before_clustering"
    )
    assert ctx.metadata["cell_cycle_layout_claimable"] is True
    assert adata.uns["cell_cycle"]["embedding_effect"] == "applied_before_clustering"


def test_catalog_orders_cell_cycle_before_clustering_when_both_requested():
    from workflow.modular.pipeline import MANDATORY_MODULES, _resolve_execution_order

    order = _resolve_execution_order(
        list(MANDATORY_MODULES), ["cell_cycle", "clustering"]
    )
    assert "cell_cycle" in order and "clustering" in order
    assert order.index("cell_cycle") < order.index("clustering")
