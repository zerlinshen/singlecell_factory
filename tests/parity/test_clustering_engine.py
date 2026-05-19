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


# ----------------------------------------------------------------------------
# F-2 contract tests (Plan ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md,
# Principle 7): --clustering-engine=sparse_exact is now a dedicated CPU-only
# dispatch arm with strict input validation + forced arpack SVD + provenance.
# ----------------------------------------------------------------------------


@pytest.mark.skipif(not _try_import_sparse_engine(), reason="phase 2 not enabled")
def test_sparse_exact_raises_on_dense_input(synthetic_clustered_adata, tmp_path):
    """F-2: _run_sparse_exact must reject dense input (no silent randomized SVD).

    Pre-F-2 the flag fell through to _run_cpu which silently used
    svd_solver="randomized" for dense input — a mild compromise that the
    flag name "sparse_exact" did not advertise. Post-F-2 the dispatch arm
    raises loudly so the flag honors its name (Principle 7).
    """
    import numpy as np
    from anndata import AnnData
    from workflow.modular.modules.clustering import ClusteringModule

    # Dense input.
    adata_dense = AnnData(np.asarray(synthetic_clustered_adata.X.todense()))
    adata_dense.obs = synthetic_clustered_adata.obs.copy()
    adata_dense.var = synthetic_clustered_adata.var.copy()

    ctx = _make_minimal_ctx(adata_dense, scale_mode="standard", tmp_path=tmp_path)
    module = ClusteringModule()
    with pytest.raises(RuntimeError, match="sparse input matrix"):
        module._run_sparse_exact(adata_dense, ctx.cfg.clustering, ctx)


@pytest.mark.skipif(not _try_import_sparse_engine(), reason="phase 2 not enabled")
def test_sparse_exact_records_provenance_metadata(synthetic_clustered_adata, tmp_path, monkeypatch):
    """F-2: _run_sparse_exact must record full provenance for lane_manifest.json (G-C0).

    Heavy scanpy numerical steps (HVG, PCA, neighbors, UMAP, Leiden) are
    stubbed to no-ops because scanpy's HVG/Wilcoxon paths invoke numba JIT
    that does not compile on this tiny synthetic fixture in some host
    environments. The contract being tested is the F-2 provenance-recording
    contract, not the numerical correctness of scanpy on a 6-cell fixture.
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from workflow.modular.modules.clustering import ClusteringModule
    from workflow.modular.modules import clustering as clustering_mod

    adata = synthetic_clustered_adata.copy()
    ctx = _make_minimal_ctx(adata, scale_mode="standard", tmp_path=tmp_path)
    module = ClusteringModule()

    # Stub the heavy scanpy / scipy steps to no-ops while preserving the
    # contract that the test cares about: input validation + provenance.
    monkeypatch.setattr(sc.pp, "normalize_total", lambda a, **kw: None)
    monkeypatch.setattr(sc.pp, "log1p", lambda a, **kw: None)
    monkeypatch.setattr(sc.pp, "scale", lambda a, **kw: None)
    def _fake_hvg(a, **kw):
        a.var["highly_variable"] = True
    monkeypatch.setattr(sc.pp, "highly_variable_genes", _fake_hvg)
    def _fake_pca(a, **kw):
        a.obsm["X_pca"] = np.zeros((a.n_obs, kw.get("n_comps", 20)), dtype=np.float32)
    monkeypatch.setattr(sc.tl, "pca", _fake_pca)
    monkeypatch.setattr(module, "_plot_pca_variance", lambda a, c, n: None)
    def _fake_get_or_compute(*a, **kw):
        kw.get("compute_fn", lambda: None)()
    monkeypatch.setattr(clustering_mod, "get_or_compute_neighbors", _fake_get_or_compute)
    monkeypatch.setattr(sc.pp, "neighbors", lambda a, **kw: None)
    monkeypatch.setattr(sc.tl, "umap", lambda a, **kw: None)
    def _fake_leiden(a, **kw):
        a.obs["leiden"] = pd.Categorical(["0"] * a.n_obs)
    monkeypatch.setattr(sc.tl, "leiden", _fake_leiden)

    module._run_sparse_exact(adata, ctx.cfg.clustering, ctx)

    # F-2 provenance contract — lane_manifest.json reads these keys.
    assert ctx.metadata.get("clustering_engine") == "sparse_exact"
    assert ctx.metadata.get("clustering_lane") == "scanpy_cpu_sparse_exact"
    assert ctx.metadata.get("pca_solver") == "arpack", (
        "F-2 violation: pca_solver must be 'arpack' (forced); no randomized SVD silently."
    )
    assert ctx.metadata.get("pca_input_was_sparse") is True
    assert ctx.metadata.get("neighbors_backend") == "scanpy_umap_true_knn"
    assert ctx.metadata.get("leiden_backend") == "igraph_cpu"
    assert ctx.metadata.get("leiden_flavor") == "igraph"
    assert ctx.metadata.get("sparse_exact_principle_8_compatible") is True


@pytest.mark.skipif(not _try_import_sparse_engine(), reason="phase 2 not enabled")
def test_sparse_exact_dispatch_routes_through_run_impl(synthetic_clustered_adata, tmp_path, monkeypatch):
    """F-2: SC_CLUSTERING_ENGINE=sparse_exact end-to-end dispatch routes to _run_sparse_exact.

    Verifies the dispatch site in _run_impl correctly hands off to the new
    sparse_exact arm (not _run_cpu) and that the early-return path populates
    the standard downstream metadata (n_clusters, clustering_backend). The
    sparse_exact arm itself is fully stubbed to avoid the numba JIT hit on
    the synthetic fixture; what's tested here is the dispatch routing.
    """
    import pandas as pd
    from workflow.modular.modules.clustering import ClusteringModule

    monkeypatch.setenv("SC_CLUSTERING_ENGINE", "sparse_exact")
    adata = synthetic_clustered_adata.copy()
    ctx = _make_minimal_ctx(adata, scale_mode="standard", tmp_path=tmp_path)
    ctx.random_state = 0

    module = ClusteringModule()
    call_log = []

    def stub_sparse_exact(adata, cfg, ctx):
        call_log.append("sparse_exact")
        adata.obs["leiden"] = pd.Categorical(["0"] * adata.n_obs)
        ctx.metadata["clustering_engine"] = "sparse_exact"

    def stub_cpu(adata, cfg, ctx):
        call_log.append("cpu")
        adata.obs["leiden"] = pd.Categorical(["0"] * adata.n_obs)

    def stub_gpu(adata, cfg, ctx):
        call_log.append("gpu")
        adata.obs["leiden"] = pd.Categorical(["0"] * adata.n_obs)

    monkeypatch.setattr(module, "_run_sparse_exact", stub_sparse_exact)
    monkeypatch.setattr(module, "_run_cpu", stub_cpu)
    monkeypatch.setattr(module, "_run_gpu", stub_gpu)
    monkeypatch.setattr(ClusteringModule, "_plot_umap_clusters", lambda self, a, c: None)

    module._run_impl(ctx)

    assert call_log == ["sparse_exact"], (
        f"F-2 dispatch violation: call_log={call_log}. Expected sparse_exact, not CPU/GPU."
    )
    assert ctx.metadata.get("clustering_backend") == "cpu_sparse_exact"
    assert ctx.metadata.get("clustering_engine") == "sparse_exact"
