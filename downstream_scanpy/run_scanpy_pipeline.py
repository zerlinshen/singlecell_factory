"""Simple Scanpy pipeline for PBMC analysis.

This is the original, simplified pipeline for quick analysis.
For full features (doublet detection, adaptive QC, comprehensive plots), 
use scanpy_pipeline_optimized.py instead.

Steps:
- Read input .h5ad or 10x MTX
- Basic QC (mitochondrial percent, gene counts)
- Normalize, log-transform, select highly variable genes (HVGs)
- Scale, PCA, neighbors, UMAP
- Leiden clustering
- Rank marker genes per cluster and save outputs

Output Structure:
    output/
    └── <project_name>/
        ├── data/           # Processed .h5ad files
        ├── figures/        # UMAP and marker plots
        └── tables/         # CSV results
"""
import argparse
import os
import re
from datetime import datetime
import scanpy as sc
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTPUT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "output"))
TIMESTAMP_PATTERN = re.compile(r'_\d{8}_\d{4,6}$')


def build_timestamped_project_name(project_name: str) -> str:
    """Append a timestamp unless the name already ends with one."""
    normalized = re.sub(r'[^A-Za-z0-9._-]+', '_', project_name.strip()).strip('_')
    if not normalized:
        normalized = 'unnamed_project'
    if TIMESTAMP_PATTERN.search(normalized):
        return normalized
    return f"{normalized}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

def setup_project_dirs(project_name: str, base_output_dir: str = DEFAULT_OUTPUT_ROOT) -> dict:
    """Create standardized project output directories."""
    project_dir = os.path.abspath(os.path.join(base_output_dir, project_name))
    dirs = {
        'root': project_dir,
        'data': os.path.join(project_dir, 'data'),
        'figures': os.path.join(project_dir, 'figures'),
        'tables': os.path.join(project_dir, 'tables'),
        'cache': os.path.join(project_dir, 'cache')  # Project-specific cache
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return dirs


def run(input_path, dirs):
    """Run the simplified analysis pipeline."""
    sc.settings.figdir = dirs['figures']
    sc.settings.autoshow = False
    
    print(f"Loading data from: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # QC
    print("Running QC...")
    if "gene_symbols" in adata.var.columns:
        gene_symbols = adata.var["gene_symbols"].astype(str).str.upper()
    else:
        gene_symbols = adata.var_names.astype(str).str.upper()
    adata.var["mt"] = gene_symbols.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    
    initial_cells = adata.n_obs
    adata = adata[adata.obs["n_genes_by_counts"] > 200, :].copy()
    adata = adata[adata.obs["pct_counts_mt"] < 5, :].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"QC complete: {initial_cells} → {adata.n_obs} cells retained")
    
    # Normalization
    print("Normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    
    # Dimensionality reduction and clustering
    print("Running dimensionality reduction and clustering...")
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, random_state=0)
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph", directed=False, random_state=0)
    
    # Marker genes
    print("Finding marker genes...")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    
    # Save results
    print(f"Saving results to: {dirs['root']}")
    adata.write(os.path.join(dirs['data'], 'processed.h5ad'))
    
    # Plots
    sc.pl.umap(adata, color=["leiden"], legend_loc="on data")
    plt.savefig(os.path.join(dirs['figures'], "umap_leiden.png"), 
                bbox_inches="tight", dpi=150)
    plt.close()
    
    # Marker dotplot
    markers = {
        "T": ["CD3D", "CD3E", "IL7R", "CCR7"],
        "CD8": ["CD8A", "CD8B", "GZMB", "NKG7"],
        "NK": ["GNLY", "NKG7", "FCGR3A"],
        "B": ["MS4A1", "CD79A", "CD74"],
        "Mono": ["LYZ", "LST1", "S100A8", "S100A9"],
        "DC": ["FCER1A", "CST3"],
        "Platelet": ["PPBP", "PF4"],
    }
    if "gene_symbols" in adata.raw.var.columns:
        symbol_to_var = dict(
            zip(
                adata.raw.var["gene_symbols"].astype(str).values,
                adata.raw.var_names.astype(str).values
            )
        )
        var_names = [symbol_to_var[g] for lst in markers.values() for g in lst if g in symbol_to_var]
    else:
        var_names = [g for lst in markers.values() for g in lst if g in adata.raw.var_names]
    if var_names:
        fig = sc.pl.dotplot(adata, var_names=var_names, 
                           groupby="leiden", return_fig=True)
        fig.savefig(os.path.join(dirs['figures'], 
                                "dotplot_markers.png"),
                   bbox_inches="tight", dpi=150)
        plt.close()
    
    # Save tables
    de_df = sc.get.rank_genes_groups_df(adata, None)
    de_df.to_csv(os.path.join(dirs['tables'], "marker_genes.csv"), index=False)
    
    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    cluster_counts.to_csv(os.path.join(dirs['tables'], "cluster_counts.csv"))
    
    print("Pipeline complete!")
    print(f"  Data: {dirs['data']}")
    print(f"  Figures: {dirs['figures']}")
    print(f"  Tables: {dirs['tables']}")


def main():
    p = argparse.ArgumentParser(
        description='Simple Scanpy scRNA-seq Pipeline',
        epilog="""
Examples:
  # From Cell Ranger output
  python run_scanpy_pipeline.py --tenx_dir=../pbmc_1k_v3_count/outs/filtered_feature_bc_matrix --project=pbmc1k_simple

  # From existing h5ad
  python run_scanpy_pipeline.py --input=data/processed.h5ad --project=pbmc1k_rerun

Output Structure:
  output/
  └── <project_name>/
      ├── data/processed.h5ad
      ├── figures/
      │   ├── umap_leiden.png
      │   └── dotplot_markers.png
      └── tables/
          ├── marker_genes.csv
          └── cluster_counts.csv
        """
    )
    
    p.add_argument("--input", 
                  help="Path to input .h5ad; use --tenx_dir for 10x directory")
    p.add_argument("--tenx_dir", 
                  help="Path to 10x filtered_feature_bc_matrix directory")
    p.add_argument("--tenx_var_names", default="gene_symbols", 
                  choices=["gene_symbols", "gene_ids"],
                  help="Gene name type when reading 10x MTX")
    p.add_argument("--project", required=True,
                  help="Base project name; timestamp is appended automatically")
    p.add_argument("--base_output_dir", default=DEFAULT_OUTPUT_ROOT,
                  help=f"Base output directory (default: {DEFAULT_OUTPUT_ROOT})")
    
    args = p.parse_args()
    
    # Setup directories
    project_name = build_timestamped_project_name(args.project)
    print(f"Using project name: {project_name}")
    dirs = setup_project_dirs(project_name, args.base_output_dir)
    
    if args.tenx_dir and os.path.isdir(args.tenx_dir):
        print(f"Reading 10x data from: {args.tenx_dir}")
        # Use project-specific cache directory
        import os as os_module
        original_dir = os_module.getcwd()
        os_module.chdir(dirs['cache'])
        try:
            adata = sc.read_10x_mtx(args.tenx_dir, 
                                   var_names=args.tenx_var_names, 
                                   cache=True)
        finally:
            os_module.chdir(original_dir)
        tmp_h5ad = os.path.join(dirs['data'], 'raw_input.h5ad')
        adata.write(tmp_h5ad)
        run(tmp_h5ad, dirs)
    elif args.input and os.path.isfile(args.input):
        run(args.input, dirs)
    else:
        raise SystemExit(
            "Provide --input .h5ad or --tenx_dir pointing to "
            "10x filtered_feature_bc_matrix"
        )


if __name__ == "__main__":
    main()
