"""Benchmark script for EFF-1 neighbors cache.

Runs nano and small_real tiers (skips medium/full if not available locally).
Emits 3 columns per tier to <project-root>/runs/<run-id>/perf/eff1_bench.json:
  wall_no_cache        - full pipeline wall with SC_CACHE_NEIGHBORS=0
  wall_with_cache_total - full pipeline wall with SC_CACHE_NEIGHBORS=1 (includes hash cost)
  hash_overhead        - sum of all sha256(X.data.tobytes()) durations in the cache-on run

Also reports cache_hits, cache_misses, and per-miss reasons from WARN log lines.

Hash overhead is amortized into wall_with_cache_total (not subtracted) so the
headline number is honest. The isolated hash_overhead column lets observers see
how much of the with-cache wall is hashing vs. actual cache savings.

Acceptance gates:
  small_real: HARD — wall_with_cache_total < wall_no_cache (target >= 2x speedup)
  nano:       SOFT — any positive speedup is sufficient
  full_real:  OBSERVATIONAL — logs WARN if hash_overhead >= 0.5 * (savings - hash_overhead)
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Optional

# Ensure project root is on sys.path so `workflow` is importable when the
# script is invoked directly (e.g. `python scripts/bench_neighbors_cache.py`).
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def _make_nano_adata():
    """Synthetic 120-cell AnnData for the nano tier."""
    import anndata as ad
    import numpy as np
    import scanpy as sc
    import scipy.sparse as sp

    n_obs, n_vars = 120, 80
    X = sp.random(n_obs, n_vars, density=0.3, format="csr", random_state=0).astype("float32")
    adata = ad.AnnData(X=X)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=50, flavor="seurat")
    sc.tl.pca(adata, n_comps=10, mask_var="highly_variable")
    return adata


def _find_small_real_adata() -> Optional[object]:
    """Attempt to load a cached small_real (~40k cell) fixture h5ad."""
    import anndata as ad

    candidates = [
        "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/tests/data/small_real.h5ad",
        "/home/zerlinshen/projects/fixtures/small_real.h5ad",
    ]
    for p in candidates:
        if Path(p).exists():
            logger.info("Loading small_real fixture from %s", p)
            return ad.read_h5ad(p)
    return None


def _run_neighbors_pipeline(adata_in, cache_val: str, n_reps: int = 3):
    """Run neighbors repeatedly, returning wall time and cache stats."""
    import scanpy as sc
    from workflow.modular import _neighbors_cache
    from workflow.modular._neighbors_cache import get_or_compute_neighbors

    _neighbors_cache.clear()
    _neighbors_cache._stats["hits"] = 0
    _neighbors_cache._stats["misses"] = 0
    _neighbors_cache._stats["hash_overhead_total"] = 0.0

    sc.settings.n_jobs = 1
    old_env = os.environ.get("SC_CACHE_NEIGHBORS")
    os.environ["SC_CACHE_NEIGHBORS"] = cache_val

    try:
        adata = adata_in.copy()
        t0 = time.perf_counter()
        for _ in range(n_reps):
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
                compute_fn=lambda: sc.pp.neighbors(
                    adata, n_neighbors=10, n_pcs=10, use_rep="X_pca", method="umap"
                ),
            )
        wall = time.perf_counter() - t0
        snap = _neighbors_cache.stats()
    finally:
        if old_env is None:
            os.environ.pop("SC_CACHE_NEIGHBORS", None)
        else:
            os.environ["SC_CACHE_NEIGHBORS"] = old_env

    return wall, snap


def _bench_tier(tier_name: str, adata, n_reps: int = 3) -> dict:
    logger.info("Benchmarking tier=%s n_reps=%d ...", tier_name, n_reps)

    # Cache-off run: use direct scanpy call (cache disabled)
    import scanpy as sc
    from workflow.modular import _neighbors_cache

    sc.settings.n_jobs = 1
    adata_work = adata.copy()
    os.environ["SC_CACHE_NEIGHBORS"] = "0"
    t0 = time.perf_counter()
    for _ in range(n_reps):
        sc.pp.neighbors(adata_work, n_neighbors=10, n_pcs=10, use_rep="X_pca", method="umap")
    wall_no_cache = time.perf_counter() - t0
    os.environ.pop("SC_CACHE_NEIGHBORS", None)

    # Cache-on run
    wall_with_cache, snap = _run_neighbors_pipeline(adata, "1", n_reps=n_reps)

    hash_overhead = snap["hash_overhead_total_s"]
    hits = snap["hits"]
    misses = snap["misses"]

    speedup = wall_no_cache / wall_with_cache if wall_with_cache > 0 else float("inf")

    result = {
        "tier": tier_name,
        "n_reps": n_reps,
        "wall_no_cache": round(wall_no_cache, 4),
        "wall_with_cache_total": round(wall_with_cache, 4),
        "hash_overhead": round(hash_overhead, 6),
        "cache_hits": hits,
        "cache_misses": misses,
        "speedup": round(speedup, 3),
    }

    # Gate evaluation
    if tier_name == "small_real":
        if wall_with_cache < wall_no_cache:
            logger.info("HARD GATE PASS: small_real speedup=%.2fx", speedup)
        else:
            logger.error("HARD GATE FAIL: small_real cache is SLOWER (speedup=%.2fx)", speedup)
            result["gate_result"] = "FAIL"
        result["gate_result"] = "PASS" if wall_with_cache < wall_no_cache else "FAIL"

    elif tier_name == "full_real":
        savings = wall_no_cache - (wall_with_cache - hash_overhead)
        if savings > 0 and hash_overhead >= 0.5 * savings:
            logger.warning(
                "full_real: hash_overhead (%.3fs) >= 50%% of compute savings (%.3fs). "
                "EFF-1 does not pay off at this scale — Round-1c should evaluate persisted hash.",
                hash_overhead, savings,
            )
            result["full_real_hash_dominates"] = True
        result["gate_result"] = "OBSERVATIONAL"

    else:
        gate = "PASS" if wall_with_cache < wall_no_cache else "SOFT_FAIL"
        result["gate_result"] = gate
        logger.info("Tier %s gate=%s speedup=%.2fx", tier_name, gate, speedup)

    logger.info(
        "tier=%s no_cache=%.4fs with_cache=%.4fs hash=%.6fs hits=%d misses=%d speedup=%.2fx",
        tier_name, wall_no_cache, wall_with_cache, hash_overhead, hits, misses, speedup,
    )
    return result


def main():
    parser = argparse.ArgumentParser(description="EFF-1 neighbors cache benchmark")
    parser.add_argument(
        "--tier",
        choices=["nano", "small_real", "medium_real", "full_real", "all"],
        default="nano",
        help="Fixture tier to benchmark (default: nano)",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Directory to write eff1_bench.json (default: runs/bench_<timestamp>/perf/)",
    )
    parser.add_argument(
        "--n-reps",
        type=int,
        default=3,
        help="Number of repeated neighbors calls per tier (default: 3)",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    if args.out_dir:
        out_dir = Path(args.out_dir)
    else:
        import datetime
        run_id = f"bench_{datetime.datetime.now().strftime('%Y%m%dT%H%M%S')}"
        out_dir = project_root / "runs" / run_id / "perf"
    out_dir.mkdir(parents=True, exist_ok=True)

    tiers_to_run = (
        ["nano", "small_real", "medium_real", "full_real"]
        if args.tier == "all"
        else [args.tier]
    )

    results = []
    for tier in tiers_to_run:
        if tier == "nano":
            adata = _make_nano_adata()
        elif tier in ("small_real", "medium_real", "full_real"):
            adata = _find_small_real_adata()
            if adata is None:
                logger.warning(
                    "Fixture for tier=%s not found locally; running nano only. "
                    "Place fixture at tests/data/small_real.h5ad to enable this tier.",
                    tier,
                )
                if tier != "nano":
                    results.append({"tier": tier, "gate_result": "SKIPPED", "reason": "fixture_not_found"})
                    continue
                adata = _make_nano_adata()
        else:
            logger.warning("Unknown tier %s, skipping", tier)
            continue

        # Ensure X_pca is present
        import scanpy as sc
        if "X_pca" not in adata.obsm:
            sc.pp.highly_variable_genes(adata, n_top_genes=min(50, adata.n_vars), flavor="seurat")
            sc.tl.pca(adata, n_comps=min(10, adata.n_obs - 1, adata.n_vars - 1), mask_var=None)

        result = _bench_tier(tier, adata, n_reps=args.n_reps)
        results.append(result)

    out_path = out_dir / "eff1_bench.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info("Results written to %s", out_path)

    # Exit non-zero if any HARD gate failed
    hard_fails = [r for r in results if r.get("gate_result") == "FAIL"]
    if hard_fails:
        logger.error("HARD GATE FAILURE for tiers: %s", [r["tier"] for r in hard_fails])
        sys.exit(1)


if __name__ == "__main__":
    main()
