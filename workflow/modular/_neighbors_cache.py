"""In-process neighbors/PCA result cache for the singlecell_factory pipeline.

Activation: set env var SC_CACHE_NEIGHBORS=1. Default (unset or 0) = disabled,
preserving current behavior bit-for-bit.

Cache key (tight 9-tuple):
  (backend_id, method, metric, n_pcs, n_neighbors, use_rep, knn, random_state,
   sha256(shape || X.indices || X.indptr || X.data))  # sparse
   sha256(shape || X.tobytes())                        # dense

backend_id must be one of {"scanpy", "rapids"}. Any other value raises ValueError.

Hash coverage: sha256 covers the matrix shape AND, for sparse inputs, the
coordinate arrays (indices/indptr) in addition to the nonzero values (X.data).
A coordinate-blind hash (values only) could return a stale/wrong cached graph
when two sparse matrices share nonzero values at different positions; including
indices/indptr/shape tightens correctness.

LRU cap: 32 entries. Beyond 32 unique keys the oldest entry is evicted (LRU
order). In typical pipelines <=5 unique tight-key combinations appear per run;
the cap exists to bound memory during parameter sweeps.

Non-cacheable modules: bbknn and any call site whose neighbors API cannot be
expressed in the tight key signature MUST bypass the cache entirely and call
sc.pp.neighbors / rsc.pp.neighbors directly. Bypassing is correct and safe.
List: ["bbknn"].

Eviction: _neighbors_cache.clear() is called by _shutdown.run_shutdown_cleanup()
at pipeline end (cli.py:731, top-level finally). Cache is in-process only and
is never persisted to disk.
"""
from __future__ import annotations

import hashlib
import logging
import os
import time
from collections import OrderedDict
from contextlib import contextmanager
from typing import Optional, Tuple

import numpy as np
from scipy import sparse

logger = logging.getLogger(__name__)

_VALID_BACKENDS = frozenset({"scanpy", "rapids"})
_MAX_ENTRIES = 32

_cache: OrderedDict[tuple, Tuple[object, object]] = OrderedDict()

_stats = {"hits": 0, "misses": 0, "hash_overhead_total": 0.0}


def _cache_enabled() -> bool:
    return os.environ.get("SC_CACHE_NEIGHBORS", "0").strip() == "1"


@contextmanager
def _timed():
    t0 = time.perf_counter()
    yield
    _stats["hash_overhead_total"] += time.perf_counter() - t0


def _compute_key(
    adata,
    *,
    backend_id: str,
    method: str,
    metric: str,
    n_pcs: Optional[int],
    n_neighbors: int,
    use_rep: str,
    knn: bool,
    random_state: int,
) -> tuple:
    if backend_id not in _VALID_BACKENDS:
        raise ValueError(
            f"neighbors cache: backend_id must be one of {_VALID_BACKENDS}, got {backend_id!r}"
        )
    X = adata.obsm.get(use_rep) if use_rep and use_rep in adata.obsm else adata.X
    with _timed():
        hasher = hashlib.sha256()
        if sparse.issparse(X):
            # Hash coordinates (indices/indptr) + shape alongside values, not just
            # X.data: a coordinate-blind hash can return a stale/wrong cached graph
            # when two sparse matrices share nonzero values at different positions.
            Xc = X.tocsr()
            hasher.update(np.asarray(Xc.shape, dtype=np.int64).tobytes())
            hasher.update(Xc.indices.tobytes())
            hasher.update(Xc.indptr.tobytes())
            hasher.update(Xc.data.tobytes())
        else:
            arr = np.asarray(X)
            hasher.update(np.asarray(arr.shape, dtype=np.int64).tobytes())
            hasher.update(arr.tobytes())
        data_hash = hasher.hexdigest()
    return (backend_id, method, metric, n_pcs, n_neighbors, use_rep, knn, random_state, data_hash)


def _looser_key_would_hit(key: tuple) -> bool:
    """Check if any cached entry differs only by backend_id/method/metric."""
    _, method, metric, n_pcs, n_neighbors, use_rep, knn, random_state, data_hash = key
    loose = (n_pcs, n_neighbors, use_rep, knn, random_state, data_hash)
    for cached_key in _cache:
        _, cm, cmetric, cn_pcs, cn_neighbors, cu_rep, cknn, crs, chash = cached_key
        if (cn_pcs, cn_neighbors, cu_rep, cknn, crs, chash) == loose:
            return True
    return False


def get_or_compute_neighbors(
    adata,
    *,
    backend_id: str,
    method: str,
    metric: str,
    n_pcs: Optional[int],
    n_neighbors: int,
    use_rep: str,
    knn: bool,
    random_state: int,
    compute_fn,
):
    """Look up (connectivities, distances) from cache or compute via compute_fn.

    compute_fn() must call sc.pp.neighbors / rsc.pp.neighbors and write results
    into adata.obsp. It receives no arguments. After compute_fn returns this
    function reads adata.obsp['connectivities'] and adata.obsp['distances'] and
    stores copies in the cache.

    Returns None; adata is mutated in place (same contract as bare neighbors call).
    """
    if not _cache_enabled():
        compute_fn()
        return

    key = _compute_key(
        adata,
        backend_id=backend_id,
        method=method,
        metric=metric,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        use_rep=use_rep,
        knn=knn,
        random_state=random_state,
    )

    if key in _cache:
        _stats["hits"] += 1
        _cache.move_to_end(key)
        conn, dist = _cache[key]
        adata.obsp["connectivities"] = conn.copy()
        adata.obsp["distances"] = dist.copy()
        logger.debug(
            "neighbors_cache hit: backend=%s method=%s metric=%s n_pcs=%s n_neighbors=%s use_rep=%s",
            backend_id, method, metric, n_pcs, n_neighbors, use_rep,
        )
        return

    _stats["misses"] += 1
    if _looser_key_would_hit(key):
        logger.warning(
            "neighbors_cache miss: key=%s backend=%s method=%s metric=%s; "
            "a looser key (same data/params, different backend/method/metric) would have hit — "
            "Round-1c telemetry: normalize call sites to raise hit rate.",
            key[:4], backend_id, method, metric,
        )
    else:
        logger.warning(
            "neighbors_cache miss: key=%s backend=%s method=%s metric=%s",
            key[:4], backend_id, method, metric,
        )

    compute_fn()

    conn = adata.obsp.get("connectivities")
    dist = adata.obsp.get("distances")
    if conn is not None and dist is not None:
        if len(_cache) >= _MAX_ENTRIES:
            _cache.popitem(last=False)
        _cache[key] = (conn.copy(), dist.copy())
        _cache.move_to_end(key)


def clear() -> None:
    """Evict all cached entries. Called by _shutdown.run_shutdown_cleanup()."""
    _cache.clear()
    logger.debug("neighbors_cache cleared")


def stats() -> dict:
    """Return a snapshot of cache hit/miss/hash_overhead stats."""
    return {
        "hits": _stats["hits"],
        "misses": _stats["misses"],
        "hash_overhead_total_s": _stats["hash_overhead_total"],
        "entries": len(_cache),
    }
