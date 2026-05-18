"""Wave 2B / US-006 — verify Hotspot 1 RSS instrumentation in clustering.py.

The 200k synthetic stress run (US-008 deliverable) will exercise the full
_run_gpu code path; this test only validates that the instrumentation
helper and the call sites are wired correctly, without requiring rapids.
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest


def test_log_rss_helper_writes_trace():
    """ClusteringModule._log_rss appends a {tag, peak_rss_mb} dict to ctx.metadata."""
    from workflow.modular.modules.clustering import ClusteringModule

    ctx = SimpleNamespace(metadata={})
    ClusteringModule._log_rss(ctx, "test:tag")
    trace = ctx.metadata.get("clustering_rss_trace")
    assert trace, "trace should be a non-empty list after first call"
    assert isinstance(trace, list)
    assert len(trace) == 1
    entry = trace[0]
    assert entry["tag"] == "test:tag"
    assert "peak_rss_mb" in entry
    assert isinstance(entry["peak_rss_mb"], (int, float))
    assert entry["peak_rss_mb"] > 0


def test_log_rss_accumulates_across_calls():
    """Multiple _log_rss calls append; preserve order."""
    from workflow.modular.modules.clustering import ClusteringModule

    ctx = SimpleNamespace(metadata={})
    for tag in ("a", "b", "c"):
        ClusteringModule._log_rss(ctx, tag)
    trace = ctx.metadata["clustering_rss_trace"]
    assert [e["tag"] for e in trace] == ["a", "b", "c"]


def test_log_rss_failure_does_not_raise():
    """If resource import fails, _log_rss must not break the caller."""
    from workflow.modular.modules.clustering import ClusteringModule

    # Pass an object that breaks metadata.setdefault to provoke the except branch.
    class BrokenCtx:
        @property
        def metadata(self):
            raise RuntimeError("intentional break")

    # Should NOT raise
    ClusteringModule._log_rss(BrokenCtx(), "should-not-raise")


def test_clone_for_gpu_lite_preserves_data_and_isolates_X():
    """US-007 mitigation: _clone_for_gpu_lite returns AnnData with deep-copied X.

    Validates the contract-preserving property: mutating the clone's X does
    NOT corrupt the original. This is the safety predicate that lets us
    skip adata.copy() while keeping the GPU-failure fallback semantically
    valid.
    """
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    import pandas as pd
    from workflow.modular.modules.clustering import ClusteringModule

    rng = np.random.default_rng(0)
    X = sp.random(100, 50, density=0.3, dtype=np.float32, random_state=0).tocsr()
    obs = pd.DataFrame({"sample": [f"S{i % 3}" for i in range(100)]},
                       index=[f"c{i}" for i in range(100)])
    var = pd.DataFrame(index=[f"g{i}" for i in range(50)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_pca"] = rng.normal(size=(100, 5)).astype(np.float32)

    clone = ClusteringModule._clone_for_gpu_lite(adata)

    assert clone.n_obs == adata.n_obs
    assert clone.n_vars == adata.n_vars
    assert sp.issparse(clone.X)
    assert clone.X is not adata.X, "X must be deep-copied (log1p/scale mutate in place)"
    assert clone.X.data.base is not adata.X.data.base or clone.X.data is not adata.X.data
    # In-place mutation on the clone must NOT affect the original.
    clone.X.data[:] = 0.0
    assert adata.X.data.sum() > 0, "original X was corrupted by clone mutation"
    # var must be a deep copy — scanpy.pp.highly_variable_genes adds columns.
    clone.var["highly_variable"] = True
    assert "highly_variable" not in adata.var.columns, "var leakage: clone modifications visible in original"


def test_clone_for_gpu_lite_allocates_less_python_heap_than_full_copy():
    """US-007 mitigation: tracemalloc-measured Python heap growth.

    adata.copy() deep-copies obs/var/obsm/varm/uns/layers; _clone_for_gpu_lite
    shares those by reference (deep-copying only X + var). Measured via
    tracemalloc, which captures actual Python heap allocation during the
    clone call — not just clone attribute nbytes (those are identical because
    both clones contain the same logical data).

    Budget: lite arm allocates LESS than copy arm. The absolute delta scales
    with the size of obs/obsm/uns; on bigger cohorts (800k cells + deep uns
    dicts) the savings are proportionally larger. Site-1 elimination (the X
    deep-copy itself) is deferred to Wave 3 per docs/HOTSPOT1_DIAGNOSIS.md.
    """
    import gc
    import tracemalloc
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    import pandas as pd
    from workflow.modular.modules.clustering import ClusteringModule

    def _build_adata():
        rng = np.random.default_rng(0)
        n, g = 20_000, 2_000
        X = sp.random(n, g, density=0.1, dtype=np.float32, random_state=0).tocsr()
        obs = pd.DataFrame({
            "sample": [f"S{i % 8}" for i in range(n)],
            "leiden": rng.integers(0, 6, size=n).astype(str),
            "cell_type": [f"T{i % 12}" for i in range(n)],
        }, index=[f"cell_{i:06d}" for i in range(n)])
        var = pd.DataFrame({"feature_type": ["gene"] * g}, index=[f"GENE_{i:05d}" for i in range(g)])
        ad_obj = ad.AnnData(X=X, obs=obs, var=var)
        ad_obj.obsm["X_pca"] = rng.normal(size=(n, 50)).astype(np.float32)
        # Deep uns blob — this is where adata.copy() pays the most overhead.
        ad_obj.uns["nc_repro_dict"] = {f"k{i}": rng.normal(size=500).tolist() for i in range(40)}
        return ad_obj

    def _alloc_for_clone(clone_fn):
        adata = _build_adata()
        gc.collect()
        tracemalloc.start()
        clone = clone_fn(adata)
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        del adata, clone
        gc.collect()
        return peak

    lite_alloc = _alloc_for_clone(ClusteringModule._clone_for_gpu_lite)
    copy_alloc = _alloc_for_clone(lambda a: a.copy())
    reduction_pct = 100.0 * (copy_alloc - lite_alloc) / copy_alloc if copy_alloc > 0 else 0.0
    print(f"\n  copy()-arm peak tracemalloc ≈ {copy_alloc/1e6:.2f} MB")
    print(f"  lite-arm  peak tracemalloc ≈ {lite_alloc/1e6:.2f} MB")
    print(f"  reduction ≈ {reduction_pct:.1f}%")
    assert lite_alloc < copy_alloc, (
        f"lite clone did not reduce heap allocation: lite={lite_alloc} copy={copy_alloc}"
    )


def test_instrumentation_call_sites_present():
    """Hotspot trace tags appear in clustering.py — covers Wave 3 M2 site names."""
    from pathlib import Path
    src = Path(__file__).resolve().parent.parent / "workflow" / "modular" / "modules" / "clustering.py"
    text = src.read_text(encoding="utf-8")
    expected_tags = [
        "hotspot1:before_inplace_gpu",
        "hotspot1:after_raw_preserve",
        "hotspot1:after_inplace_gpu_success",
        "hotspot1:_run_gpu:before_preserve_raw",
        "hotspot1:_run_gpu:after_preserve_raw",
        "hotspot1:_run_gpu:after_materialize_X",
    ]
    for tag in expected_tags:
        assert tag in text, f"missing instrumentation tag: {tag}"
