"""Check top expressed genes in acceptance gate clusters to verify their identity."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import anndata as ad
import numpy as np
import scipy.sparse as sp

print("Loading clustering checkpoint...", flush=True)
adata = ad.read_zarr(
    "/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426/"
    "nc2024_tumor_20260426_121139/.checkpoints/after_clustering.zarr"
)
print(f"n_obs={adata.n_obs}, n_vars={adata.n_vars}, clusters={adata.obs['leiden'].nunique()}", flush=True)

gate_clusters = {
    "7":  "expected B cell (MS4A1+CD79A+)",
    "10": "expected NK (KLRD1+GNLY+)",
    "13": "expected Plasma (IGLV3-1+MZB1+)",
    "14": "expected Myeloid (CD163+OLR1+)",
    "17": "expected Mast (TPSAB1+)",
}

for cluster_id, desc in gate_clusters.items():
    mask = (adata.obs["leiden"].astype(str) == cluster_id).values
    n_cells = int(mask.sum())
    x = adata.X[mask]
    if sp.issparse(x):
        mean_expr = np.asarray(x.mean(axis=0)).flatten()
    else:
        mean_expr = np.asarray(x).mean(axis=0).flatten()
    top10_idx = np.argsort(mean_expr)[-10:][::-1]
    top10 = [(adata.var_names[i], round(float(mean_expr[i]), 4)) for i in top10_idx]
    print(f"\nCluster {cluster_id} ({desc}, n={n_cells}):")
    print(f"  Top genes: {top10}")
    # Check canonical markers
    for g in ["MS4A1", "CD79A", "TPSAB1", "CPA3", "MZB1", "IGLV3-1", "CD163", "LYZ", "KLRD1", "GNLY", "KLRF1"]:
        if g in adata.var_names:
            idx = list(adata.var_names).index(g)
            print(f"  {g}: {mean_expr[idx]:.4f}")
