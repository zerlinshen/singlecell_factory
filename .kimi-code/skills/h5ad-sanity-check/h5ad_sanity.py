#!/usr/bin/env python3
"""h5ad_sanity.py - structural sanity probe for AnnData outputs.

Prints JSON facts about an .h5ad file (backed read, safe for large files).
The calling agent compares the facts against the run contract and issues the verdict.

Usage: h5ad_sanity.py <file.h5ad>
Exit 0 = probe ran (facts printed); exit 1 = unreadable file; exit 2 = usage/deps error.
"""
import json
import sys


def main():
    if len(sys.argv) < 2:
        print(json.dumps({"error": "usage: h5ad_sanity.py <file.h5ad>"}))
        sys.exit(2)
    path = sys.argv[1]
    try:
        import anndata as ad
        import numpy as np
    except ImportError as e:
        print(json.dumps({"error": f"missing dependency: {e}. Run with a project env python."}))
        sys.exit(2)
    try:
        a = ad.read_h5ad(path, backed="r")
    except Exception as e:
        print(json.dumps({"error": f"cannot read {path}: {e}"}))
        sys.exit(1)

    out = {"path": path, "n_obs": int(a.n_obs), "n_vars": int(a.n_vars)}

    X = a.X
    out["x_type"] = type(X).__name__
    try:
        if hasattr(X, "nnz"):
            data = X.data
            out["x_nnz"] = int(X.nnz)
            out["x_density"] = round(float(X.nnz) / max(1, a.n_obs * a.n_vars), 6)
            if data is not None and len(data):
                out["x_dtype"] = str(data.dtype)
                out["x_nan"] = int(np.isnan(data).sum())
                out["x_negative"] = int((data < 0).sum())
        else:
            # Backed sparse (CSRDataset etc.) or dense arrays: probe a bounded
            # row slice so 900k-cell files stay memory-safe.
            rows = min(a.n_obs, 2000)
            sub = X[:rows, :]
            arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
            out["x_dtype"] = str(arr.dtype)
            out["x_sample_rows"] = int(rows)
            out["x_nan_sampled"] = int(np.isnan(arr).sum())
            out["x_negative_sampled"] = int((arr < 0).sum())
            out["x_density_sampled"] = round(float((arr != 0).mean()), 6)
    except Exception as e:
        out["x_probe_error"] = str(e)

    try:
        out["obs_columns"] = list(a.obs.columns)
        out["var_columns"] = list(a.var.columns)
    except Exception as e:
        out["obs_var_probe_error"] = str(e)
    try:
        out["layers"] = list(a.layers.keys())
        out["obsm"] = list(a.obsm.keys())
        out["uns_keys"] = list(a.uns.keys())
    except Exception as e:
        out["slots_probe_error"] = str(e)
    try:
        out["obs_name_duplicates"] = int(a.obs_names.duplicated().sum())
    except Exception:
        out["obs_name_duplicates"] = "probe-failed"
    try:
        names = [str(g) for g in a.var_names[:200000]]
        out["mt_gene_prefix_present"] = any(g.startswith(("MT-", "mt-")) for g in names)
    except Exception:
        pass

    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
