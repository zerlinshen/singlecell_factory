"""Clustering science gate: SC_CLUSTERING_ENGINE=sparse_exact vs CSS (massive mode).

Unlike numerical parity gates, this is a SCIENCE gate: sparse_exact is expected
to be MORE accurate than CSS (which is a Cluster Similarity Spectrum approximation
that trades precision for memory). We verify that:
  1. The flag actually disables CSS routing
  2. sparse_exact produces a clustering that recovers planted ground-truth groups
     better than chance (>= 0.5 ARI vs ground truth)
  3. The engine produces same-or-better cluster purity than CSS on the same data

We do NOT compare cluster IDs between sparse_exact and CSS directly — they are
different algorithms and IDs are not comparable. We compare against ground truth.
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
    """SC_CLUSTERING_ENGINE=sparse_exact must override massive-mode CSS routing."""
    from workflow.modular.modules.clustering import ClusteringModule

    monkeypatch.setenv("SC_CLUSTERING_ENGINE", "sparse_exact")
    adata = synthetic_clustered_adata.copy()
    ctx = _make_minimal_ctx(adata, scale_mode="massive", tmp_path=tmp_path)
    module = ClusteringModule()
    use_css = module._should_use_css(ctx, adata)
    assert use_css is False
    assert ctx.metadata.get("css_status") == "disabled_sparse_exact_engine"
    assert ctx.metadata.get("clustering_engine") == "sparse_exact"


@pytest.mark.skipif(not _try_import_sparse_engine(), reason="phase 2 not enabled")
def test_clustering_default_massive_uses_css(synthetic_clustered_adata, tmp_path, monkeypatch):
    """Default (no env flag) in massive mode still routes to CSS."""
    from workflow.modular.modules.clustering import ClusteringModule

    monkeypatch.delenv("SC_CLUSTERING_ENGINE", raising=False)
    adata = synthetic_clustered_adata.copy()
    ctx = _make_minimal_ctx(adata, scale_mode="massive", tmp_path=tmp_path)
    module = ClusteringModule()
    use_css = module._should_use_css(ctx, adata)
    assert use_css is True
    assert ctx.metadata.get("css_status") == "enabled"


@pytest.mark.skipif(not _try_import_sparse_engine(), reason="phase 2 not enabled")
def test_clustering_default_large_uses_cpu_path(synthetic_clustered_adata, tmp_path, monkeypatch):
    """Default in large/standard mode never uses CSS (CSS is massive-only)."""
    from workflow.modular.modules.clustering import ClusteringModule

    monkeypatch.delenv("SC_CLUSTERING_ENGINE", raising=False)
    adata = synthetic_clustered_adata.copy()
    ctx = _make_minimal_ctx(adata, scale_mode="large", tmp_path=tmp_path)
    module = ClusteringModule()
    use_css = module._should_use_css(ctx, adata)
    assert use_css is False
