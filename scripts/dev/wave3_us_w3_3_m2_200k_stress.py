"""US-W3-3 — measure peak host RSS for M2 mechanism at 200k synthetic.

Compares:
  - pre-fix arm: adata.copy() + adata.raw = adata_copy (the historical pattern)
  - M2 arm: in-place mutation + adata.raw mandated (no .copy() of the main adata)

Result determines whether the xfail-strict gate on
tests/test_clustering_memory_regression.py::test_200k_post_fix_under_35GB_and_drops_30pct_vs_prefix
can be flipped (requires <35 GB AND >=30% reduction).
"""
from __future__ import annotations

import gc
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _read_vmrss_gb() -> float:
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) / (1024.0 ** 2)
    return float("nan")


def _build_synthetic(n_cells: int, n_genes: int = 20_000, density: float = 0.1, seed: int = 42):
    rng = np.random.default_rng(seed)
    X = sp.random(n_cells, n_genes, density=density, format="csr", dtype=np.float32, random_state=seed)
    X.data = np.ceil(X.data * 100).astype(np.float32)
    obs = pd.DataFrame(
        {"sample": [f"S{i % 8}" for i in range(n_cells)]},
        index=[f"cell_{i:06d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"GENE_{i:05d}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


def arm_pre_fix(n_cells: int) -> dict:
    """Historical pattern: adata.copy() + raw = adata_copy.

    Mirrors the pre-Wave-3 b' mitigation pattern. Used as the pre-fix baseline
    for the regression delta.
    """
    import scanpy as sc

    rss_pre = _read_vmrss_gb()
    adata = _build_synthetic(n_cells)
    t0 = time.time()
    # Hotspot 1 site 1: full host duplication.
    adata_gpu = adata.copy()
    rss_after_copy = _read_vmrss_gb()
    # Site 2: raw preservation on the clone.
    adata_gpu.raw = adata_gpu
    # Materialize on host.
    if sp.issparse(adata_gpu.X):
        adata_gpu.X = adata_gpu.X.tocsr().astype(np.float32)
    rss_after_materialize = _read_vmrss_gb()
    # Scanpy preprocessing (host-side, same as _run_gpu does).
    sc.pp.normalize_total(adata_gpu, target_sum=1e4)
    sc.pp.log1p(adata_gpu)
    rss_peak = _read_vmrss_gb()
    elapsed = time.time() - t0
    del adata_gpu, adata
    gc.collect()
    return {
        "arm": "pre_fix",
        "n_cells": n_cells,
        "rss_pre_gb": round(rss_pre, 2),
        "rss_after_copy_gb": round(rss_after_copy, 2),
        "rss_after_materialize_gb": round(rss_after_materialize, 2),
        "rss_peak_gb": round(rss_peak, 2),
        "elapsed_secs": round(elapsed, 1),
    }


def arm_m2(n_cells: int) -> dict:
    """Wave 3 M2 pattern: in-place mutation on original adata + raw mandated.

    Mirrors clustering.py `_run_impl` GPU branch under Wave 3 M2:
      - skip _clone_for_gpu_lite
      - adata.raw = adata.copy() (M2 invariant)
      - host-side scanpy preprocessing mutates adata.X in place
    """
    import scanpy as sc

    rss_pre = _read_vmrss_gb()
    adata = _build_synthetic(n_cells)
    t0 = time.time()
    # M2 invariant: mandate raw preservation (the only "copy" in this arm).
    adata.raw = adata.copy()
    rss_after_raw = _read_vmrss_gb()
    # No host-side _clone_for_gpu_lite — _run_gpu mutates adata in place.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    rss_peak = _read_vmrss_gb()
    elapsed = time.time() - t0
    del adata
    gc.collect()
    return {
        "arm": "m2",
        "n_cells": n_cells,
        "rss_pre_gb": round(rss_pre, 2),
        "rss_after_raw_gb": round(rss_after_raw, 2),
        "rss_peak_gb": round(rss_peak, 2),
        "elapsed_secs": round(elapsed, 1),
    }


def main() -> None:
    out: dict = {"us_w3_3_m2_200k_stress": {}}
    # Pre-fix arm
    pre = arm_pre_fix(200_000)
    out["us_w3_3_m2_200k_stress"]["pre_fix"] = pre
    gc.collect()
    # M2 arm
    m2 = arm_m2(200_000)
    out["us_w3_3_m2_200k_stress"]["m2"] = m2
    # Compute reduction
    if pre["rss_peak_gb"] > 0:
        reduction_pct = 100.0 * (pre["rss_peak_gb"] - m2["rss_peak_gb"]) / pre["rss_peak_gb"]
    else:
        reduction_pct = 0.0
    out["us_w3_3_m2_200k_stress"]["reduction_pct"] = round(reduction_pct, 1)
    out["us_w3_3_m2_200k_stress"]["budget_35GB_met"] = m2["rss_peak_gb"] < 35.0
    out["us_w3_3_m2_200k_stress"]["reduction_30pct_met"] = reduction_pct >= 30.0
    out["us_w3_3_m2_200k_stress"]["ac_met"] = m2["rss_peak_gb"] < 35.0 and reduction_pct >= 30.0
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
