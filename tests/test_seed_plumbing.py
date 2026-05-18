"""Tests for P-CC.S26a: ctx.random_state plumbing through PipelineContext and stochastic modules."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from anndata import AnnData

from workflow.modular.config import CellRangerConfig, PipelineConfig
from workflow.modular.context import PipelineContext


def _make_ctx(random_state: int = 42, tmp_path: Path | None = None) -> PipelineContext:
    run_dir = tmp_path or Path("/tmp/seed_test_run")
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg = PipelineConfig(
        project="seed_test",
        output_dir=run_dir / "out",
        cellranger=CellRangerConfig(
            sample_root=run_dir,
            outs_dir=run_dir / "outs" / "filtered_feature_bc_matrix",
        ),
        random_state=random_state,
    )
    return PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
    )


def _synthetic_adata(n_obs: int = 20, n_vars: int = 30) -> AnnData:
    rng = np.random.default_rng(0)
    adata = AnnData(rng.integers(0, 10, size=(n_obs, n_vars)).astype(np.float32))
    adata.obs_names = [f"C{i}" for i in range(n_obs)]
    adata.var_names = [f"G{i}" for i in range(n_vars)]
    return adata


# ---------------------------------------------------------------------------
# ctx.random_state property
# ---------------------------------------------------------------------------

def test_ctx_random_state_default(tmp_path):
    ctx = _make_ctx(random_state=42, tmp_path=tmp_path)
    assert ctx.random_state == 42


def test_ctx_random_state_custom(tmp_path):
    ctx = _make_ctx(random_state=99, tmp_path=tmp_path)
    assert ctx.random_state == 99


def test_pipeline_config_random_state_field():
    cfg = PipelineConfig(
        project="x",
        output_dir=Path("/tmp/x"),
        cellranger=CellRangerConfig(
            sample_root=Path("/tmp/x"),
            outs_dir=Path("/tmp/x/outs/filtered_feature_bc_matrix"),
        ),
        random_state=77,
    )
    assert cfg.random_state == 77


# ---------------------------------------------------------------------------
# doublet_detection: metadata key written (fallback path — tiny dataset)
# ---------------------------------------------------------------------------

def test_random_state_propagates_to_doublet_detection(tmp_path):
    from workflow.modular.modules.doublet_detection import DoubletDetectionModule

    ctx = _make_ctx(random_state=42, tmp_path=tmp_path)
    # Tiny dataset triggers the fast fallback path (all singlets) — no scrublet needed
    tiny = _synthetic_adata(n_obs=5, n_vars=5)
    ctx.adata = tiny

    mod = DoubletDetectionModule()
    mod.run(ctx)

    assert ctx.metadata.get("doublet_random_state") == 42


# ---------------------------------------------------------------------------
# batch_correction: sync of clustering.random_state before post-correction steps
# ---------------------------------------------------------------------------

def test_batch_correction_syncs_random_state_to_clustering_cfg(tmp_path):
    from workflow.modular.modules.batch_correction import BatchCorrectionModule

    ctx = _make_ctx(random_state=42, tmp_path=tmp_path)
    ctx.cfg.clustering.random_state = 0  # intentionally wrong before sync
    ctx.adata = _synthetic_adata()
    ctx.adata.obs["leiden"] = "0"

    mod = BatchCorrectionModule()
    # Module skips because the batch key ("sample") is absent — that's fine;
    # the sync must happen before the skip guard.
    mod.run(ctx)

    assert ctx.cfg.clustering.random_state == 42


# ---------------------------------------------------------------------------
# clustering: random_state override and metadata recording
# ---------------------------------------------------------------------------

def test_clustering_random_state_synced_and_recorded(tmp_path, monkeypatch):
    from workflow.modular.modules import clustering as clustering_mod
    from workflow.modular.modules.clustering import ClusteringModule

    ctx = _make_ctx(random_state=42, tmp_path=tmp_path)
    adata = _synthetic_adata()
    adata.obs["leiden"] = "0"
    adata.obsm["X_pca"] = np.zeros((adata.n_obs, 2))
    adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2))
    ctx.adata = adata

    mod = ClusteringModule()

    # Stub out heavy compute while keeping the seed-sync and metadata lines
    def _fake_cpu(a, cfg, c):
        a.obs["leiden"] = "0"
        a.obsm.setdefault("X_pca", np.zeros((a.n_obs, 2)))
        a.obsm.setdefault("X_umap", np.zeros((a.n_obs, 2)))

    monkeypatch.setattr(mod, "_run_cpu", _fake_cpu)
    monkeypatch.setattr(mod, "_should_use_css", lambda *a, **kw: False)
    monkeypatch.setattr(mod, "_plot_umap_clusters", lambda *a, **kw: None)

    mod._run_impl(ctx)

    assert ctx.metadata.get("clustering_random_state") == 42
    assert ctx.cfg.clustering.random_state == 42


# ---------------------------------------------------------------------------
# differential_expression: metadata key written
# ---------------------------------------------------------------------------

def test_de_random_state_recorded(tmp_path, monkeypatch):
    from workflow.modular.modules.differential_expression import DifferentialExpressionModule

    ctx = _make_ctx(random_state=42, tmp_path=tmp_path)
    adata = _synthetic_adata()
    adata.obs["leiden"] = "0"
    ctx.adata = adata

    mod = DifferentialExpressionModule()

    # Stub out rank_genes_groups so we don't need a full pipeline setup
    import pandas as pd

    def _fake_run_cpu_rank(adata, method, rank_kwargs, n_genes, corr_method=None):
        cols = ["group", "names", "scores", "pvals_adj", "logfoldchanges", "pct_nz_group", "pct_nz_reference"]
        return pd.DataFrame(columns=cols)

    monkeypatch.setattr(mod, "_run_cpu_rank_genes_groups", _fake_run_cpu_rank)

    import unittest.mock as mock
    with mock.patch("workflow.modular.modules.differential_expression.gpu_available", return_value=False):
        mod._run_impl(ctx)

    assert ctx.metadata.get("de_random_state") == 42
