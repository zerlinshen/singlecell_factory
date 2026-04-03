from __future__ import annotations

import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors


def baseline_qc(
    adata,
    min_genes: int,
    max_genes: int,
    max_mito_pct: float,
    max_ribo_pct: float,
    min_cells: int,
):
    """Baseline QC implementation before memory optimization."""
    gene_upper = adata.var_names.astype(str).str.upper()
    adata.var["mt"] = gene_upper.str.startswith("MT-")
    adata.var["ribo"] = gene_upper.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)
    adata = adata[adata.obs["n_genes_by_counts"] >= min_genes, :].copy()
    adata = adata[adata.obs["n_genes_by_counts"] <= max_genes, :].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= max_mito_pct, :].copy()
    adata = adata[adata.obs["pct_counts_ribo"] <= max_ribo_pct, :].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata


def baseline_clustering(
    adata,
    target_sum: float,
    n_top_genes: int,
    n_neighbors: int,
    n_pcs: int,
    resolution: float,
    random_state: int,
):
    """Baseline clustering implementation before PCA memory optimization."""
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_top_genes)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns.get("pca", {})
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep="X_pca")
    sc.tl.umap(adata, random_state=random_state)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        flavor="igraph",
        directed=False,
        random_state=random_state,
    )
    return adata


def baseline_pseudo_velocity(adata, n_neighbors: int = 15):
    """Baseline pseudo-velocity with Python loop."""
    umap = adata.obsm["X_umap"]
    pseudotime = adata.obs["dpt_pseudotime"].to_numpy()
    knn = NearestNeighbors(n_neighbors=n_neighbors)
    knn.fit(umap)
    _, idx = knn.kneighbors(umap)
    velocity = np.zeros_like(umap)
    for i in range(umap.shape[0]):
        nbr = idx[i, 1:]
        dt = pseudotime[nbr] - pseudotime[i]
        direction = umap[nbr] - umap[i]
        if np.allclose(dt, 0.0):
            continue
        velocity[i] = (direction * dt[:, None]).mean(axis=0)
    return velocity
