"""Metacell engine/claim stamps (Wave-2 W2.2)."""
from __future__ import annotations

import numpy as np
import anndata as ad

from workflow.modular.modules.metacell import MetacellModule


class _Ctx:
    def __init__(self, adata):
        self.adata = adata
        self.metadata = {}
        import tempfile
        from pathlib import Path

        self._tmp = tempfile.TemporaryDirectory()
        base = Path(self._tmp.name)
        self.figure_dir = base / "figures"
        self.table_dir = base / "tables"
        self.figure_dir.mkdir()
        self.table_dir.mkdir()

    def __del__(self):
        try:
            self._tmp.cleanup()
        except Exception:
            pass


def _adata(n=80, p=20):
    rng = np.random.default_rng(0)
    a = ad.AnnData(rng.poisson(1.0, size=(n, p)).astype(np.float32))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var_names = [f"g{i}" for i in range(p)]
    a.obsm["X_pca"] = rng.normal(size=(n, 5)).astype(np.float32)
    a.obsm["X_umap"] = rng.normal(size=(n, 2)).astype(np.float32)
    a.obs["leiden"] = (["0"] * 40 + ["1"] * 40)
    return a


def test_seacells_failure_stamps_kmeans_nonclaimable(monkeypatch):
    mod = MetacellModule()

    def boom(adata, n_metacells):
        return None, "seacells_failed", "ImportError: no SEACells"

    monkeypatch.setattr(mod, "_try_seacells", staticmethod(boom))
    ctx = _Ctx(_adata())
    mod._run_impl(ctx)
    assert ctx.metadata["metacell_engine"] == "minibatch_kmeans_fallback"
    assert ctx.metadata["metacell_claimable"] is False
    assert "fallback" in ctx.metadata["metacell_claim_reason"]
    assert ctx.adata.uns["metacell"]["claimable"] is False
    assert "metacell" in ctx.adata.obs


def test_seacells_success_is_claimable(monkeypatch):
    mod = MetacellModule()
    n = 80
    labels = np.array([i % 10 for i in range(n)])

    def ok(adata, n_metacells):
        return labels, "seacells", "Persad_SEACells_2023"

    monkeypatch.setattr(mod, "_try_seacells", staticmethod(ok))
    ctx = _Ctx(_adata(n=n))
    mod._run_impl(ctx)
    assert ctx.metadata["metacell_engine"] == "seacells"
    assert ctx.metadata["metacell_claimable"] is True
