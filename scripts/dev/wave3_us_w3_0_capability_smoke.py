"""US-W3-0 Stage 1 — capability smoke for M1 + M2 at 50k synthetic.

Run in the sc_gpu conda env (has rapids-singlecell + cupy + CUDA toolchain):
    /home/zerlinshen/conda/envs/sc_gpu/bin/python scripts/dev/wave3_us_w3_0_capability_smoke.py

Writes structured JSON to stdout. Designed to be invoked by Ralph iteration.
"""
from __future__ import annotations

import gc
import json
import os
import sys
import time
import traceback
from typing import Any

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


def _build_synthetic(n_cells: int, n_genes: int = 20_000, density: float = 0.1, seed: int = 42) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    X = sp.random(n_cells, n_genes, density=density, format="csr", dtype=np.float32, random_state=seed)
    # Cast to integer counts so normalize_total has something sane to scale.
    X.data = np.ceil(X.data * 100).astype(np.float32)
    obs = pd.DataFrame(
        {"sample": [f"S{i % 8}" for i in range(n_cells)]},
        index=[f"cell_{i:06d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"GENE_{i:05d}" for i in range(n_genes)])
    a = ad.AnnData(X=X, obs=obs, var=var)
    return a


def smoke_m1_gpu_preprocessing(n_cells: int = 50_000) -> dict[str, Any]:
    """M1 capability smoke: rsc.pp.normalize_total/log1p/scale on GPU sparse input."""
    result: dict[str, Any] = {"mechanism": "M1", "n_cells": n_cells}
    try:
        import rapids_singlecell as rsc
        import cupy as cp
        import cupyx.scipy.sparse as csp
    except Exception as exc:
        result["status"] = "import_failed"
        result["error"] = repr(exc)
        return result

    result["rsc_version"] = getattr(rsc, "__version__", "unknown")
    adata = _build_synthetic(n_cells)
    rss_pre = _read_vmrss_gb()
    result["rss_pre_upload_gb"] = round(rss_pre, 2)

    # GPU upload: convert host sparse CSR to cupy sparse CSR.
    t0 = time.time()
    try:
        X_gpu = csp.csr_matrix(adata.X)
        adata.X = X_gpu  # AnnData holds the GPU-resident matrix
        upload_secs = time.time() - t0
        result["upload_secs"] = round(upload_secs, 3)
    except Exception as exc:
        result["status"] = "gpu_upload_failed"
        result["error"] = repr(exc)
        result["traceback"] = traceback.format_exc()
        return result

    rss_after_upload = _read_vmrss_gb()
    result["rss_after_upload_gb"] = round(rss_after_upload, 2)

    # rsc.pp.normalize_total — verify it accepts GPU sparse and stays sparse.
    try:
        rsc.pp.normalize_total(adata, target_sum=1e4)
        result["normalize_total_ok"] = True
        result["normalize_total_X_type"] = type(adata.X).__name__
        result["normalize_total_X_is_gpu_sparse"] = csp.issparse(adata.X)
    except Exception as exc:
        result["status"] = "normalize_total_failed"
        result["error"] = repr(exc)
        result["traceback"] = traceback.format_exc()
        return result

    try:
        rsc.pp.log1p(adata)
        result["log1p_ok"] = True
        result["log1p_X_type"] = type(adata.X).__name__
        result["log1p_X_is_gpu_sparse"] = csp.issparse(adata.X)
    except Exception as exc:
        result["status"] = "log1p_failed"
        result["error"] = repr(exc)
        return result

    # HVG selects genes; needed before scale to keep dense conversion bounded.
    try:
        rsc.pp.highly_variable_genes(adata, n_top_genes=2_000, flavor="seurat")
        result["hvg_ok"] = True
        result["n_hvg"] = int(adata.var["highly_variable"].sum())
    except Exception as exc:
        result["status"] = "hvg_failed"
        result["error"] = repr(exc)
        return result

    # rsc.pp.scale on sparse with zero_center=False keeps it sparse.
    try:
        rsc.pp.scale(adata, max_value=10, zero_center=False)
        result["scale_ok"] = True
        result["scale_X_type"] = type(adata.X).__name__
        result["scale_X_is_gpu_sparse"] = csp.issparse(adata.X)
    except Exception as exc:
        result["status"] = "scale_failed"
        result["error"] = repr(exc)
        result["traceback"] = traceback.format_exc()
        return result

    try:
        rsc.pp.pca(adata, n_comps=50, mask_var="highly_variable")
        result["pca_ok"] = True
        # adata.obsm['X_pca'] may be cupy ndarray; convert to host shape check
        Xp = adata.obsm["X_pca"]
        result["X_pca_shape"] = list(getattr(Xp, "shape", (-1, -1)))
        result["X_pca_type"] = type(Xp).__name__
    except Exception as exc:
        result["status"] = "pca_failed"
        result["error"] = repr(exc)
        result["traceback"] = traceback.format_exc()
        return result

    rss_post = _read_vmrss_gb()
    result["rss_post_pca_gb"] = round(rss_post, 2)
    result["host_rss_growth_gb"] = round(rss_post - rss_pre, 3)
    result["status"] = "ok"
    return result


def smoke_m2_raw_preservation(n_cells: int = 50_000) -> dict[str, Any]:
    """M2 capability smoke: adata.raw preserved + restore-on-fallback dry run.

    No GPU required for this arm — exercises the host-side raw preservation
    + restore semantic that M2 depends on.
    """
    result: dict[str, Any] = {"mechanism": "M2", "n_cells": n_cells}
    try:
        import scanpy as sc
    except Exception as exc:
        result["status"] = "scanpy_import_failed"
        result["error"] = repr(exc)
        return result

    adata = _build_synthetic(n_cells)
    rss_pre = _read_vmrss_gb()
    result["rss_pre_gb"] = round(rss_pre, 2)

    # Mandate raw preservation (M2 invariant).
    adata.raw = adata.copy()
    result["raw_set_ok"] = True
    result["raw_X_shape"] = list(adata.raw.X.shape)

    # Mutate the main slot via the same scanpy pipeline _run_gpu uses.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    result["host_X_mutated"] = True

    # Now simulate GPU failure: restore from raw.
    try:
        restored = adata.raw.to_adata()
        result["restore_ok"] = True
        result["restore_X_shape"] = list(restored.X.shape)
        # The restored adata should carry the raw (un-normalized) counts; spot check.
        sample_data = restored.X[:5, :5].toarray() if sp.issparse(restored.X) else restored.X[:5, :5]
        result["restore_raw_max_value"] = float(sample_data.max())
    except Exception as exc:
        result["status"] = "restore_failed"
        result["error"] = repr(exc)
        return result

    rss_post = _read_vmrss_gb()
    result["rss_post_gb"] = round(rss_post, 2)
    result["host_rss_growth_gb"] = round(rss_post - rss_pre, 3)
    result["status"] = "ok"
    return result


def main() -> None:
    out = {
        "us_w3_0_stage1_capability_smoke": {
            "n_cells": 50_000,
            "host_total_ram_gb": int(
                int(os.popen("grep MemTotal /proc/meminfo | awk '{print $2}'").read().strip()) / (1024 ** 2)
            ) if os.path.exists("/proc/meminfo") else None,
        }
    }
    m1 = smoke_m1_gpu_preprocessing()
    out["us_w3_0_stage1_capability_smoke"]["M1"] = m1
    gc.collect()

    m2 = smoke_m2_raw_preservation()
    out["us_w3_0_stage1_capability_smoke"]["M2"] = m2

    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
