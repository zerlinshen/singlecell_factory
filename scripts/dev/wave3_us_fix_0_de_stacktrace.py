"""US-FIX-0 — probe the Series.nonzero failure with a stack trace.

Builds a 1000-cell synthetic AnnData with 10 clusters and runs
sc.tl.rank_genes_groups (and rsc DE fallback if scanpy works). Captures
the exact top-frame where the AttributeError fires so US-FIX-3 can decide
whether to pin scanpy, pin pandas, or apply the monkeypatch.

Run:
  /home/zerlinshen/conda/envs/sc_gpu_stable/bin/python \
    scripts/dev/wave3_us_fix_0_de_stacktrace.py
"""
from __future__ import annotations

import json
import sys
import traceback

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def main() -> None:
    out: dict = {}
    rng = np.random.default_rng(0)
    n_cells, n_genes = 1000, 500
    X = sp.random(n_cells, n_genes, density=0.2, format="csr", dtype=np.float32, random_state=0)
    X.data = np.ceil(X.data * 100).astype(np.float32)
    obs = pd.DataFrame(
        {"leiden": rng.integers(0, 10, size=n_cells).astype(str)},
        index=[f"c{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    a = ad.AnnData(X=X, obs=obs, var=var)

    import scanpy as sc
    out["scanpy_version"] = sc.__version__
    out["pandas_version"] = pd.__version__
    out["numpy_version"] = np.__version__
    try:
        import rapids_singlecell as rsc
        out["rsc_version"] = rsc.__version__
    except Exception:
        out["rsc_version"] = None

    # Minimal normalize + log1p + HVG to make rank_genes_groups happy
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.obs["leiden"] = a.obs["leiden"].astype("category")

    # Probe 1: scanpy CPU rank_genes_groups
    try:
        sc.tl.rank_genes_groups(a, groupby="leiden", method="wilcoxon")
        out["scanpy_rank_genes_groups"] = "ok"
    except Exception as exc:
        tb = traceback.format_exc()
        out["scanpy_rank_genes_groups"] = "FAILED"
        out["scanpy_error"] = repr(exc)[:200]
        out["scanpy_top_frame"] = next(
            (line.strip() for line in tb.splitlines() if "_scanpy" in line or "scanpy/" in line),
            "no scanpy frame found"
        )
        out["scanpy_traceback_tail"] = tb.splitlines()[-8:]

    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
