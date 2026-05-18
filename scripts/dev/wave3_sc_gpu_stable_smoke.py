"""Smoke test for sc_gpu_stable env — confirms rapids-singlecell PCA works
on RTX 5090 D v2 with the RAPIDS 25.10 / CUDA 12.9 stack.

Reproduces the 50k×20k×2000 HVG case that crashed under RAPIDS 26.02 / CUDA 13
with CUSOLVER_STATUS_INTERNAL_ERROR. If PCA succeeds here, the Wave 3 plan can
proceed to a real NC2024 GPU run.

Run:
  /home/zerlinshen/conda/envs/sc_gpu_stable/bin/python \
    scripts/dev/wave3_sc_gpu_stable_smoke.py
"""
from __future__ import annotations

import json
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def main() -> None:
    result: dict = {"env": "sc_gpu_stable"}
    try:
        import rapids_singlecell as rsc
        import cupy as cp
        import cupyx.scipy.sparse as csp
        result["rsc_version"] = rsc.__version__
        result["cupy_version"] = cp.__version__
        result["cuda_runtime"] = cp.cuda.runtime.runtimeGetVersion()
        result["compute_capability"] = cp.cuda.Device().compute_capability
        result["device_name"] = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    except Exception as exc:
        result["status"] = "import_failed"
        result["error"] = repr(exc)
        print(json.dumps(result, indent=2))
        return

    rng = np.random.default_rng(0)
    X = sp.random(50_000, 20_000, density=0.1, format="csr", dtype=np.float32, random_state=0)
    X.data = np.ceil(X.data * 100).astype(np.float32)
    obs = pd.DataFrame({"s": [f"S{i % 3}" for i in range(50_000)]}, index=[f"c{i}" for i in range(50_000)])
    var = pd.DataFrame(index=[f"g{i}" for i in range(20_000)])
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.X = csp.csr_matrix(a.X)

    t0 = time.time()
    try:
        rsc.pp.normalize_total(a, target_sum=1e4)
        rsc.pp.log1p(a)
        rsc.pp.highly_variable_genes(a, n_top_genes=2_000, flavor="seurat")
        rsc.pp.pca(a, n_comps=50, mask_var="highly_variable")
        rsc.pp.neighbors(a, n_neighbors=15, n_pcs=50, use_rep="X_pca")
        rsc.tl.umap(a)
        rsc.tl.leiden(a, resolution=1.0)
        elapsed = time.time() - t0
        result["status"] = "ok"
        result["wall_secs"] = round(elapsed, 1)
        result["n_clusters"] = int(pd.Series(a.obs["leiden"]).nunique())
    except Exception as exc:
        result["status"] = "pca_or_downstream_failed"
        result["error"] = repr(exc)[:200]

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
