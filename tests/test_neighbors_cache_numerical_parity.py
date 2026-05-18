"""Parity test for the EFF-1 neighbors cache.

Runs a nano pipeline (neighbors -> umap -> leiden) twice:
  Run A: SC_CACHE_NEIGHBORS=0  (disabled, control)
  Run B: SC_CACHE_NEIGHBORS=1  (enabled, second call hits the cache)

Both runs use sc.settings.n_jobs = 1 and random_state=42 to make UMAP/Leiden
deterministic. With n_jobs=1 and a fixed random_state the exact floating-point
path is identical regardless of cache; parity must be bit-exact (no atol).

If this test fails, EFF-1 ships with SC_CACHE_NEIGHBORS inert (disabled by
default everywhere). Do NOT enable the cache in any code path without a passing
parity gate.
"""
from __future__ import annotations

import os

import numpy as np
import pytest
import scipy.sparse as sp


def _make_nano_adata():
    """Return a small synthetic AnnData suitable for a nano neighbors run.

    Bypasses scanpy's numba-backed normalize_total to avoid numba/Python 3.13
    JIT compilation issues inside pytest. Constructs a pre-normalized log1p
    dense matrix and manually adds a synthetic X_pca embedding, so the fixture
    exercises only the neighbors/umap/leiden path (which is what we're testing).
    """
    import anndata as ad
    import scanpy as sc

    rng = np.random.default_rng(0)
    n_obs, n_vars = 120, 80
    # Dense float64 avoids numba CSR path in normalize_total
    X_raw = rng.random((n_obs, n_vars)).astype(np.float64)
    # Manual log1p normalization (row-sum normalize then log1p)
    row_sums = X_raw.sum(axis=1, keepdims=True)
    X_norm = np.log1p(X_raw / row_sums * 1e4)
    adata = ad.AnnData(X=X_norm)
    # Synthetic PCA embedding — 10 dims, bypasses scanpy PCA (avoids numba)
    adata.obsm["X_pca"] = rng.standard_normal((n_obs, 10)).astype(np.float64)
    return adata


def _run_pipeline(adata_in, cache_env: str):
    """Run neighbors->umap->leiden on a copy of adata_in with the given env setting.

    Both cache-off (cache_env="0") and cache-on (cache_env="1") runs go through
    get_or_compute_neighbors so that the cache-on second call exercises the hit path.
    sc.settings.n_jobs=1 pins all downstream computations to a deterministic path.
    """
    import scanpy as sc
    from workflow.modular import _neighbors_cache
    from workflow.modular._neighbors_cache import get_or_compute_neighbors

    _neighbors_cache.clear()
    _neighbors_cache._stats["hits"] = 0
    _neighbors_cache._stats["misses"] = 0

    adata = adata_in.copy()
    sc.settings.n_jobs = 1

    old_env = os.environ.get("SC_CACHE_NEIGHBORS")
    os.environ["SC_CACHE_NEIGHBORS"] = cache_env
    try:
        def _compute():
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10, use_rep="X_pca", method="umap")

        # First call — always computes (cache miss or disabled)
        get_or_compute_neighbors(
            adata,
            backend_id="scanpy",
            method="umap",
            metric="euclidean",
            n_pcs=10,
            n_neighbors=10,
            use_rep="X_pca",
            knn=True,
            random_state=42,
            compute_fn=_compute,
        )

        if cache_env == "1":
            # Second call with identical params — must hit the cache
            get_or_compute_neighbors(
                adata,
                backend_id="scanpy",
                method="umap",
                metric="euclidean",
                n_pcs=10,
                n_neighbors=10,
                use_rep="X_pca",
                knn=True,
                random_state=42,
                compute_fn=_compute,
            )

        sc.tl.umap(adata, random_state=42)
        sc.tl.leiden(adata, resolution=0.5, flavor="igraph", directed=False, random_state=42)
    finally:
        if old_env is None:
            os.environ.pop("SC_CACHE_NEIGHBORS", None)
        else:
            os.environ["SC_CACHE_NEIGHBORS"] = old_env

    return adata


@pytest.fixture(scope="module")
def nano_base():
    return _make_nano_adata()


def test_cache_parity_connectivities(nano_base):
    """connectivities must be bit-exact between cache-off and cache-on runs."""
    adata_A = _run_pipeline(nano_base, "0")
    adata_B = _run_pipeline(nano_base, "1")
    assert np.array_equal(
        adata_A.obsp["connectivities"].data,
        adata_B.obsp["connectivities"].data,
    ), "connectivities.data differs between cache-off and cache-on"


def test_cache_parity_distances(nano_base):
    """distances must be bit-exact between cache-off and cache-on runs."""
    adata_A = _run_pipeline(nano_base, "0")
    adata_B = _run_pipeline(nano_base, "1")
    assert np.array_equal(
        adata_A.obsp["distances"].data,
        adata_B.obsp["distances"].data,
    ), "distances.data differs between cache-off and cache-on"


def test_cache_parity_umap(nano_base):
    """X_umap must be bit-exact (no atol) between cache-off and cache-on runs.

    Exact equality holds because n_jobs=1 and random_state=42 pin the UMAP
    random walk to a single deterministic path.
    """
    adata_A = _run_pipeline(nano_base, "0")
    adata_B = _run_pipeline(nano_base, "1")
    assert np.array_equal(
        adata_A.obsm["X_umap"],
        adata_B.obsm["X_umap"],
    ), "X_umap differs between cache-off and cache-on (no atol — must be exact)"


def test_cache_parity_leiden(nano_base):
    """leiden cluster labels must be identical between cache-off and cache-on runs."""
    adata_A = _run_pipeline(nano_base, "0")
    adata_B = _run_pipeline(nano_base, "1")
    assert np.array_equal(
        adata_A.obs["leiden"].values,
        adata_B.obs["leiden"].values,
    ), "leiden labels differ between cache-off and cache-on"


def test_cache_disabled_by_default(nano_base):
    """With SC_CACHE_NEIGHBORS unset the cache module must not alter behavior."""
    os.environ.pop("SC_CACHE_NEIGHBORS", None)
    from workflow.modular._neighbors_cache import _cache_enabled
    assert not _cache_enabled(), "cache must be disabled when SC_CACHE_NEIGHBORS is unset"


def test_backend_id_assertion():
    """get_or_compute_neighbors raises ValueError for unknown backend_id."""
    import anndata as ad
    import scipy.sparse as sp
    from workflow.modular._neighbors_cache import get_or_compute_neighbors

    rng = np.random.default_rng(1)
    X = sp.random(20, 10, density=0.5, format="csr", random_state=1).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obsm["X_pca"] = rng.random((20, 5)).astype(np.float32)

    old = os.environ.get("SC_CACHE_NEIGHBORS")
    os.environ["SC_CACHE_NEIGHBORS"] = "1"
    try:
        with pytest.raises(ValueError, match="backend_id"):
            get_or_compute_neighbors(
                adata,
                backend_id="unknown_backend",
                method="umap",
                metric="euclidean",
                n_pcs=5,
                n_neighbors=5,
                use_rep="X_pca",
                knn=True,
                random_state=0,
                compute_fn=lambda: None,
            )
    finally:
        if old is None:
            os.environ.pop("SC_CACHE_NEIGHBORS", None)
        else:
            os.environ["SC_CACHE_NEIGHBORS"] = old


def test_lru_cap():
    """Cache evicts oldest entry when 32-entry cap is exceeded."""
    import anndata as ad
    from workflow.modular import _neighbors_cache

    _neighbors_cache.clear()
    old = os.environ.get("SC_CACHE_NEIGHBORS")
    os.environ["SC_CACHE_NEIGHBORS"] = "1"
    try:
        # Insert 33 distinct entries by varying random_state
        for i in range(33):
            rng = np.random.default_rng(i)
            X = sp.random(20, 10, density=0.5, format="csr", random_state=i).astype(np.float32)
            adata = ad.AnnData(X=X)
            adata.obsm["X_pca"] = rng.random((20, 5)).astype(np.float32)
            import scipy.sparse as _sp

            conn = _sp.eye(20, format="csr")
            dist = _sp.eye(20, format="csr")
            adata.obsp["connectivities"] = conn
            adata.obsp["distances"] = dist

            from workflow.modular._neighbors_cache import get_or_compute_neighbors
            get_or_compute_neighbors(
                adata,
                backend_id="scanpy",
                method="umap",
                metric="euclidean",
                n_pcs=5,
                n_neighbors=5,
                use_rep="X_pca",
                knn=True,
                random_state=i,
                compute_fn=lambda: None,
            )
        assert len(_neighbors_cache._cache) <= 32, (
            f"LRU cap exceeded: {len(_neighbors_cache._cache)} entries"
        )
    finally:
        if old is None:
            os.environ.pop("SC_CACHE_NEIGHBORS", None)
        else:
            os.environ["SC_CACHE_NEIGHBORS"] = old
        _neighbors_cache.clear()
