"""Regression tests for the cell_fate CellRank-fallback opt-in gate.

Audit 2026-05-22 finding C-2: the previous module silently swapped the cited
Lange 2022 CellRank method for a hand-rolled connectivity-diffusion proxy on
any non-ImportError exception, and never wrote engine provenance into
adata.uns. These tests pin the Principle 9 contract:

- ImportError remains a canonical missing-dependency fallback (always allowed),
  but adata.uns["cell_fate"]["engine"] must be the fallback engine string, not
  cellrank.
- Any other exception during CellRank now raises RuntimeError unless
  SC_ALLOW_CELLRANK_FALLBACK=1 is set explicitly.
- The unified adata.obsm["fate_probabilities"] shape (numpy float32 +
  adata.uns["cell_fate"]["fate_probability_columns"]) is written for both
  engines.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from workflow.modular.config import CellRangerConfig, PipelineConfig, QCConfig
from workflow.modular.context import PipelineContext
from workflow.modular.modules.cell_fate import (
    CellFateModule,
    _CELLRANK_FALLBACK_ENV,
    _FALLBACK_ENGINE,
    _PRIMARY_ENGINE,
)


def _adata_with_pseudotime(n_obs: int = 40, n_vars: int = 60, n_clusters: int = 3) -> ad.AnnData:
    rng = np.random.default_rng(11)
    x = sparse.csr_matrix(rng.poisson(1.0, size=(n_obs, n_vars)).astype("float32"))
    obs = pd.DataFrame(
        {
            "dpt_pseudotime": np.linspace(0.0, 1.0, n_obs, dtype=np.float32),
            "leiden": pd.Categorical(
                [str(i % n_clusters) for i in range(n_obs)]
            ),
        },
        index=[f"cell{i:03d}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene{i:03d}" for i in range(n_vars)])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    conn = sparse.eye(n_obs, format="csr", dtype=np.float32) + sparse.csr_matrix(
        (np.ones(n_obs - 1, dtype=np.float32),
         (np.arange(n_obs - 1), np.arange(1, n_obs))),
        shape=(n_obs, n_obs),
    )
    conn = conn + conn.T
    adata.obsp["connectivities"] = conn.tocsr()
    return adata


def _ctx(tmp_path: Path) -> PipelineContext:
    cfg = PipelineConfig(
        project="cell-fate-test",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        qc=QCConfig(min_genes=0, min_counts=0),
        random_state=11,
    )
    fig_dir = tmp_path / "figures"
    table_dir = tmp_path / "tables"
    fig_dir.mkdir()
    table_dir.mkdir()
    return PipelineContext(cfg=cfg, run_dir=tmp_path, figure_dir=fig_dir, table_dir=table_dir)


def test_import_error_falls_back_and_stamps_fallback_engine(monkeypatch, tmp_path):
    """ImportError -> canonical fallback path, but engine stamp must be 'fallback_*', not cellrank."""
    monkeypatch.delenv(_CELLRANK_FALLBACK_ENV, raising=False)
    adata = _adata_with_pseudotime()
    ctx = _ctx(tmp_path)
    ctx.adata = adata

    def fake_cellrank(adata_arg, ctx_arg):
        raise ImportError("cellrank not installed")

    monkeypatch.setattr(CellFateModule, "_run_cellrank", staticmethod(fake_cellrank))

    CellFateModule().run(ctx)

    assert adata.uns["cell_fate"]["engine"] == _FALLBACK_ENGINE
    assert adata.uns["cell_fate"]["primary_engine"] == _PRIMARY_ENGINE
    assert "cellrank_import_error" in adata.uns["cell_fate"]["fallback_reason"]
    assert ctx.metadata["cell_fate_engine"] == _FALLBACK_ENGINE
    assert ctx.metadata["cell_fate_method_actually_used"] == _FALLBACK_ENGINE
    assert isinstance(adata.obsm["fate_probabilities"], np.ndarray)
    assert adata.obsm["fate_probabilities"].dtype == np.float32
    assert "fate_probability_columns" in adata.uns["cell_fate"]


def test_runtime_error_raises_without_opt_in(monkeypatch, tmp_path):
    """Non-ImportError CellRank failure must raise unless the operator opts in."""
    monkeypatch.delenv(_CELLRANK_FALLBACK_ENV, raising=False)
    adata = _adata_with_pseudotime()
    ctx = _ctx(tmp_path)
    ctx.adata = adata

    def fake_cellrank(adata_arg, ctx_arg):
        raise RuntimeError("CellRank exploded mid-computation")

    monkeypatch.setattr(CellFateModule, "_run_cellrank", staticmethod(fake_cellrank))

    with pytest.raises(RuntimeError) as exc_info:
        CellFateModule().run(ctx)
    msg = str(exc_info.value)
    assert _CELLRANK_FALLBACK_ENV in msg
    assert "banned by default" in msg
    assert "Lange 2022" in msg


def test_runtime_error_with_opt_in_uses_fallback_and_acknowledges(monkeypatch, tmp_path):
    """SC_ALLOW_CELLRANK_FALLBACK=1 lets the fallback run and stamps an acknowledgement flag."""
    monkeypatch.setenv(_CELLRANK_FALLBACK_ENV, "1")
    adata = _adata_with_pseudotime()
    ctx = _ctx(tmp_path)
    ctx.adata = adata

    def fake_cellrank(adata_arg, ctx_arg):
        raise RuntimeError("CellRank exploded mid-computation")

    monkeypatch.setattr(CellFateModule, "_run_cellrank", staticmethod(fake_cellrank))

    CellFateModule().run(ctx)

    assert ctx.metadata["cell_fate_engine"] == _FALLBACK_ENGINE
    assert ctx.metadata["cell_fate_fallback_opt_in_acknowledged"] is True
    assert "cellrank_runtime_error_opt_in" in adata.uns["cell_fate"]["fallback_reason"]
