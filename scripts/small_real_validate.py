"""A.2 small_real validation: preprocess, cluster, annotate with cluster_voting, write report."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("Loading small_real.h5ad...", flush=True)
adata = ad.read_h5ad(ROOT / "tests/data/staircase/small_real.h5ad")
print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

# Minimal preprocessing for clustering
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
sc.pp.pca(adata, n_comps=15)
sc.pp.neighbors(adata, n_pcs=15, n_neighbors=15)
sc.tl.leiden(adata, resolution=1.0)
sc.tl.umap(adata)
print(f"Leiden clusters: {adata.obs['leiden'].nunique()}", flush=True)

# Run AnnotationModule with cluster_voting
from workflow.modular.modules.annotation import AnnotationModule, DEFAULT_MARKERS
from workflow.modular.config import PipelineConfig, CellRangerConfig

out_dir = ROOT / "results" / "small_real_validate_20260426"
out_dir.mkdir(parents=True, exist_ok=True)
table_dir = out_dir / "tables"
table_dir.mkdir(exist_ok=True)
fig_dir = out_dir / "figures"
fig_dir.mkdir(exist_ok=True)

cfg = PipelineConfig(
    project="small_real_validate",
    output_dir=out_dir,
    cellranger=CellRangerConfig(
        sample_root=ROOT / "tests/data/staircase",
        outs_dir=ROOT / "tests/data/staircase",
    ),
    annotation_strategy="cluster_voting",
    annotation_confidence_threshold=0.1,
)


class _Ctx:
    def __init__(self):
        self.adata = adata
        self.cfg = cfg
        self.metadata = {}
        self.table_dir = table_dir
        self.figure_dir = fig_dir


ctx = _Ctx()
print("Running AnnotationModule (cluster_voting)...", flush=True)
AnnotationModule().run(ctx)
print("Annotation complete.", flush=True)
print(f"Unknown pct: {ctx.metadata.get('annotation_unknown_pct')}%", flush=True)

# Load outputs
majority = pd.read_csv(table_dir / "cluster_majority_cell_type.csv", index_col=0)
score_matrix = pd.read_csv(table_dir / "cluster_score_matrix.csv", index_col=0)
print("\n=== cluster_majority_cell_type ===")
print(majority.to_string())

# Verify canonical signal clusters
print("\n=== Marker signal verification ===", flush=True)
canonical_checks = [
    ("MS4A1", "B cell"),
    ("CD79A", "B cell"),
    ("CD3D", "T cell"),
    ("CD3E", "T cell"),
    ("TPSAB1", "Mast cell"),
    ("CPA3", "Mast cell"),
]

for gene, expected_ct in canonical_checks:
    if gene not in adata.var_names:
        print(f"  {gene}: NOT IN DATASET")
        continue
    # Find clusters with high expression of this gene
    gene_means = adata.obs.groupby("leiden", observed=True).apply(
        lambda grp, g=gene: float(
            adata[grp.index, g].X.toarray().mean()
            if hasattr(adata[grp.index, g].X, "toarray")
            else float(adata[grp.index, g].X.mean())
        )
    )
    top_cluster = gene_means.idxmax()
    assigned = majority.loc[top_cluster, "majority_cell_type"] if top_cluster in majority.index else "N/A"
    status = "PASS" if assigned == expected_ct else "FAIL"
    print(f"  [{status}] {gene} top cluster={top_cluster} expr={gene_means[top_cluster]:.3f} -> assigned='{assigned}' (expected '{expected_ct}')")

# Write validation report
report_dir = ROOT / "ops" / "nc2024_methodology_audit"
report_dir.mkdir(parents=True, exist_ok=True)
report_path = report_dir / "SMALL_REAL_VALIDATION_2026-04-26.md"

with open(report_path, "w") as f:
    f.write("# Small Real Validation Report — cluster_voting strategy\n\n")
    f.write(f"**Date:** 2026-04-26\n")
    f.write(f"**Dataset:** tests/data/staircase/small_real.h5ad ({adata.n_obs} cells)\n")
    f.write(f"**Strategy:** cluster_voting (aggregate score_genes per cluster then argmax)\n")
    f.write(f"**Leiden clusters:** {adata.obs['leiden'].nunique()}\n")
    f.write(f"**Unknown pct:** {ctx.metadata.get('annotation_unknown_pct')}%\n\n")
    f.write("## Cluster × Top Marker × Assigned Cell Type\n\n")
    f.write("| Cluster | Majority Cell Type | Top Score |\n")
    f.write("|---------|-------------------|----------|\n")
    for cluster in sorted(majority.index, key=lambda x: int(x) if str(x).isdigit() else x):
        ct = majority.loc[cluster, "majority_cell_type"]
        top_score = float(score_matrix.loc[cluster].max()) if cluster in score_matrix.index else float("nan")
        f.write(f"| {cluster} | {ct} | {top_score:.4f} |\n")
    f.write("\n## Marker Signal Verification\n\n")
    f.write("| Gene | Expected CT | Top Cluster | Assigned CT | Status |\n")
    f.write("|------|------------|-------------|-------------|--------|\n")
    for gene, expected_ct in canonical_checks:
        if gene not in adata.var_names:
            f.write(f"| {gene} | {expected_ct} | N/A | N/A | NOT IN DATASET |\n")
            continue
        gene_means = adata.obs.groupby("leiden", observed=True).apply(
            lambda grp, g=gene: float(
                adata[grp.index, g].X.toarray().mean()
                if hasattr(adata[grp.index, g].X, "toarray")
                else float(adata[grp.index, g].X.mean())
            )
        )
        top_cluster = gene_means.idxmax()
        assigned = majority.loc[top_cluster, "majority_cell_type"] if top_cluster in majority.index else "N/A"
        status = "PASS" if assigned == expected_ct else "FAIL"
        f.write(f"| {gene} | {expected_ct} | {top_cluster} | {assigned} | {status} |\n")

print(f"\nValidation report written: {report_path}", flush=True)
