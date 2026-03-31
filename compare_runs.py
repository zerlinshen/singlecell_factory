#!/usr/bin/env python
"""
Compare rerun outputs with existing outputs
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os

print("=" * 80)
print("Comparing Rerun Outputs with Existing Project Outputs")
print("=" * 80)

# Load both datasets
old_adata = sc.read('/home/zerlinshen/singlecell_factory/output/pbmc_1k_v3_scanpy_backup/data/processed.h5ad')
new_adata = sc.read('/home/zerlinshen/singlecell_factory/output/pbmc_1k_v3_scanpy/data/processed.h5ad')

print("\n--- 1. AnnData Dimensions ---")
print(f"Existing: Cells={old_adata.n_obs}, Genes={old_adata.n_vars}")
print(f"Rerun:    Cells={new_adata.n_obs}, Genes={new_adata.n_vars}")
if old_adata.shape == new_adata.shape:
    print("✓ EXACT MATCH")
else:
    print("✗ MISMATCH!")

print("\n--- 2. Cell and Gene Names ---")
cell_match = list(old_adata.obs_names) == list(new_adata.obs_names)
gene_match = list(old_adata.var_names) == list(new_adata.var_names)
print(f"Cell barcodes match: {cell_match}")
print(f"Gene names match:    {gene_match}")

print("\n--- 3. QC Metrics ---")
metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']
for metric in metrics:
    old_med = old_adata.obs[metric].median()
    new_med = new_adata.obs[metric].median()
    old_mean = old_adata.obs[metric].mean()
    new_mean = new_adata.obs[metric].mean()
    match_med = abs(old_med - new_med) < 0.1
    match_mean = abs(old_mean - new_mean) < 0.1
    print(f"{metric:20s} | Old med: {old_med:8.1f} | New med: {new_med:8.1f} | {'✓' if match_med else '✗'}")

print("\n--- 4. Leiden Clustering ---")
old_clusters = old_adata.obs['leiden'].value_counts().sort_index()
new_clusters = new_adata.obs['leiden'].value_counts().sort_index()
print(f"Existing cluster count: {len(old_clusters)}")
print(f"Rerun cluster count:    {len(new_clusters)}")
print("\nCluster cell counts:")
combined_clusters = pd.DataFrame({
    'Existing': old_clusters,
    'Rerun': new_clusters
}).fillna(0)
print(combined_clusters)

print("\n--- 5. Cluster Consistency Check ---")
# Check if clustering patterns are similar
if len(old_clusters) == len(new_clusters):
    # Check if the relative sizes are similar
    old_sizes = old_clusters.values / old_clusters.sum()
    new_sizes = new_clusters.values / new_clusters.sum()
    size_corr = np.corrcoef(old_sizes, new_sizes)[0, 1]
    print(f"Cluster size correlation: {size_corr:.3f}")
    if size_corr > 0.9:
        print("✓ HIGHLY CORRELATED - Biological structure preserved")
    elif size_corr > 0.7:
        print("⚠ MODERATELY CORRELATED - Acceptable variation")
    else:
        print("✗ LOW CORRELATION - Potential issue")

print("\n--- 6. Marker Genes ---")
# Check if both have marker gene results
old_has_markers = 'rank_genes_groups' in old_adata.uns
new_has_markers = 'rank_genes_groups' in new_adata.uns
print(f"Existing has marker results: {old_has_markers}")
print(f"Rerun has marker results:    {new_has_markers}")

if old_has_markers and new_has_markers:
    # Compare top markers for first cluster
    old_markers = old_adata.uns['rank_genes_groups']['names']['0'][:5]
    new_markers = new_adata.uns['rank_genes_groups']['names']['0'][:5]
    print(f"\nTop 5 markers for cluster 0:")
    print(f"  Existing: {list(old_markers)}")
    print(f"  Rerun:    {list(new_markers)}")
    overlap = len(set(old_markers) & set(new_markers))
    print(f"  Overlap: {overlap}/5")

print("\n--- 7. Dimensionality Reduction ---")
dims = ['pca', 'neighbors', 'umap']
for dim in dims:
    old_has = (dim in old_adata.obsm) or (dim in old_adata.uns) or (dim in old_adata.obsp)
    new_has = (dim in new_adata.obsm) or (dim in new_adata.uns) or (dim in new_adata.obsp)
    print(f"{dim:10s} | Existing: {old_has} | Rerun: {new_has} | {'✓' if old_has == new_has else '✗'}")

print("\n--- 8. Cell Type Annotation ---")
old_has_ct = 'cell_type' in old_adata.obs
new_has_ct = 'cell_type' in new_adata.obs
print(f"Existing has cell types: {old_has_ct}")
print(f"Rerun has cell types:    {new_has_ct}")

if old_has_ct and new_has_ct:
    old_ct_counts = old_adata.obs['cell_type'].value_counts()
    new_ct_counts = new_adata.obs['cell_type'].value_counts()
    print("\nCell type counts:")
    combined_ct = pd.DataFrame({
        'Existing': old_ct_counts,
        'Rerun': new_ct_counts
    }).fillna(0)
    print(combined_ct)

print("\n" + "=" * 80)
print("Comparison Complete")
print("=" * 80)
