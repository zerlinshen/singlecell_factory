"""Clustering science gate (updated 2026-05-19 for F-1 CSS removal).

Per Plan F-1 (~/.omc/plans/nc-cell-clustering-final-strategy-plan.md,
Principle 2), CSS clustering is removed from the production science path.
This file's tests verify the F-1 removal contract:
  1. SC_CLUSTERING_ENGINE=sparse_exact disables CSS (still passes — sparse_exact
     remains a valid request and CSS never runs anyway post-F-1).
  2. Explicit SC_CLUSTERING_ENGINE=css raises RuntimeError loudly.
  3. scale_mode=massive no longer auto-activates CSS (F-4 remap + F-1 gate).
  4. scale_mode=large/standard still never activates CSS.

Historical: pre-2026-05-19 this file also tested CSS science correctness on
synthetic data; those tests are removed because the codepath is removed.
F-2 will add real-sparse_exact science gates here when that lane is implemented.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData


def _try_import_sparse_engine():
    try:
        from workflow.modular.modules.clustering import ClusteringModule  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.fixture(scope="module")
def synthetic_clustered_adata():
    """3 well-separated cell-type groups with planted marker genes."""
    rng = np.random.default_rng(13)
    n_per_group = 500
    n_genes = 800
    n_marker = 30  # markers per group

    # Each group has its own up-regulated genes
    base = rng.negative_binomial(2, 0.5, size=(n_per_group * 3, n_genes)).astype(np.float32)
    # Group A: genes [0:30] up
    base[:n_per_group, :n_marker] *= 6.0
    # Group B: genes [30:60] up
    base[n_per_group:2 * n_per_group, n_marker:2 * n_marker] *= 6.0
    # Group C: genes [60:90] up
    base[2 * n_per_group:, 2 * n_marker:3 * n_marker] *= 6.0

    ground_truth = np.array(["A"] * n_per_group + ["B"] * n_per_group + ["C"] * n_per_group)

    X_sparse = sp.csr_matrix(base)
    adata = AnnData(X_sparse)
    adata.var_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    adata.obs["ground_truth"] = pd.Categorical(ground_truth)
    adata.obs["sample"] = pd.Categorical(
        rng.choice(["s1", "s2", "s3", "s4"], size=adata.n_obs)
    )
    return adata


def _make_minimal_ctx(adata, scale_mode, tmp_path):
    """Construct a minimal PipelineContext stand-in for clustering tests."""
    from types import SimpleNamespace

    cfg = SimpleNamespace(
        scale_mode=scale_mode,
        gpu_mode="off",
        clustering=SimpleNamespace(
            target_sum=10_000,
            n_top_genes=200,
            n_pcs=20,
            n_neighbors=15,
            leiden_resolution=0.8,
            random_state=0,
            scale_data=False,
        ),
    )
    figure_dir = tmp_path / "figures"
    table_dir = tmp_path / "tables"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    ctx = SimpleNamespace(
        adata=adata,
        cfg=cfg,
        metadata={},
        figure_dir=figure_dir,
        table_dir=table_dir,
    )
    return ctx


@pytest.mark.skipif(not _try_import_sparse_engine(), reason="phase 2 not enabled")
def test_clustering_engine_flag_disables_css(synthetic_clustered_adata, tmp_path, monkeypatch):
    """SC_CLUSTERING_ENGINE=sparse_exact returns False from _should_use_css.

    Post-F-1: this still passes because _should_use_css now returns False for
    ALL inputs (CSS is removed). The sparse_exact-specific metadata
    annotation is preserved.
    """
    from workflow.modular.modules.clustering import ClusteringModule

    monkeypatch.setenv("SC_CLUSTERING_ENGINE", "sparse_exact")
    adata = synthetic_clustered_adata.copy()
    ctx = _make_minimal_ctx(adata, scale_mode="massive", tmp_path=tmp_path)
    module = ClusteringModule()
    use_css = module._should_use_css(ctx, adata)
    assert use_css is False
    assert ctx.metadata.get("css_status") == "removed_per_plan_F1"
    assert ctx.metadata.get("clustering_engine") == "sparse_exact"


@pytest.mark.skipif(not _try_import_sparse_engine(), reason="phase 2 not enabled")
def test_explicit_css_request_raises(synthetic_clustered_adata, tmp_path, monkeypatch):
    """F-1 invariant: SC_CLUSTERING_ENGINE=css must raise RuntimeError loudly.

    Pre-2026-05-19 this would activate CSS clustering. Post-F-1 (Plan
    ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md, Principle 2),
    CSS is removed from the production science path; an explicit request
    must fail at the gate rather than silently routing elsewhere.
    """
    from workflow.modular.modules.clustering import ClusteringModule

    monkeypatch.setenv("SC_CLUSTERING_ENGINE", "css")
    adata = synthetic_clustered_adata.copy()
    ctx = _make_minimal_ctx(adata, scale_mode="massive", tmp_path=tmp_path)
    module = ClusteringModule()
    with pytest.raises(RuntimeError, match="CSS clustering is removed"):
        module._should_use_css(ctx, adata)


@pytest.mark.skipif(not _try_import_sparse_engine(), reason="phase 2 not enabled")
def test_clustering_default_massive_no_longer_uses_css(synthetic_clustered_adata, tmp_path, monkeypatch):
    """F-1 + F-4 invariant: scale_mode=massive without explicit CSS request must NOT route to CSS.

    Pre-2026-05-19 this returned True (silent CSS activation via the massive
    preset). Post-F-1/F-4 the preset's clustering_engine slot is remapped to
    "auto" and the dispatcher's auto branch no longer activates CSS.
    """
    from workflow.modular.modules.clustering import ClusteringModule

    monkeypatch.delenv("SC_CLUSTERING_ENGINE", raising=False)
    adata = synthetic_clustered_adata.copy()
    ctx = _make_minimal_ctx(adata, scale_mode="massive", tmp_path=tmp_path)
    module = ClusteringModule()
    use_css = module._should_use_css(ctx, adata)
    assert use_css is False, (
        "F-1 / F-4 / Principle 2 violation: scale_mode=massive activated CSS without "
        "explicit request. CSS is removed from the production science path."
    )
    assert ctx.metadata.get("css_status") == "removed_per_plan_F1"


@pytest.mark.skipif(not _try_import_sparse_engine(), reason="phase 2 not enabled")
def test_clustering_default_large_uses_cpu_path(synthetic_clustered_adata, tmp_path, monkeypatch):
    """Default in large/standard mode never uses CSS (this was always true, F-1 strengthens it)."""
    from workflow.modular.modules.clustering import ClusteringModule

    monkeypatch.delenv("SC_CLUSTERING_ENGINE", raising=False)
    adata = synthetic_clustered_adata.copy()
    ctx = _make_minimal_ctx(adata, scale_mode="large", tmp_path=tmp_path)
    module = ClusteringModule()
    use_css = module._should_use_css(ctx, adata)
    assert use_css is False
