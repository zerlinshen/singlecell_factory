# Project-Based Output Structure

## Overview

All pipeline outputs are now organized in a clean, project-dependent folder structure under `output/<project_name>/`.

## Directory Structure

```
singlecell_factory/
├── downstream_scanpy/          # Analysis scripts
│   ├── run_scanpy_pipeline.py
│   └── scanpy_pipeline_optimized.py
├── output/                     # **ALL RESULTS HERE**
│   ├── <project_1>/           # Each analysis run
│   ├── <project_2>/
│   └── ...
└── ...
```

## Project Directory Layout

Each project folder contains:

```
output/<project_name>/
├── data/                       # Processed data files
│   ├── after_qc.h5ad         # After QC filtering
│   └── processed.h5ad        # Final processed data
├── figures/                    # All visualizations
│   ├── qc/                    # Quality control plots
│   │   ├── violin_n_genes_by_counts.png
│   │   ├── violin_total_counts.png
│   │   ├── violin_pct_counts_mt.png
│   │   ├── violin_pct_counts_ribo.png
│   │   ├── scatter_counts_vs_mito.png
│   │   ├── scatter_counts_vs_genes.png
│   │   ├── scatter_genes_vs_mito.png
│   │   └── qc_histograms.png
│   ├── analysis/              # Dimensionality reduction
│   │   ├── pca_variance_ratio.png
│   │   ├── pca_qc.png
│   │   ├── pca_leiden.png
│   │   ├── umap_leiden.png
│   │   └── umap_qc.png
│   └── markers/               # Marker gene visualizations
│       ├── rank_genes_groups_leiden_top20.png
│       ├── dotplot__top5.png
│       ├── heatmap_top5.png
│       ├── dotplot_canonical_markers.png
│       └── violin_canonical_markers_*.png
├── tables/                     # CSV data tables
│   ├── filtering_stats.csv    # QC statistics
│   ├── marker_genes.csv       # Differential expression
│   ├── cluster_counts.csv     # Cells per cluster
│   ├── cluster_summary.csv    # Counts + percentages
│   └── pipeline_config.csv    # Parameters used
├── logs/                       # Pipeline logs
│   ├── pipeline.log           # Complete execution log
│   └── pipeline_config.json   # Configuration in JSON
└── cache/                      # ⭐ ISOLATED per-project cache
    └── *.h5ad                # Temporary scanpy cache files
```

## Usage Examples

### Create a New Project

```bash
cd /home/zerlinshen/singlecell_factory/downstream_scanpy

python scanpy_pipeline_optimized.py \
  --tenx_dir=/home/zerlinshen/singlecell_factory/pbmc_1k_v3_count/outs/filtered_feature_bc_matrix \
  --project=my_first_analysis
```

Result: `output/my_first_analysis/` is created with all results.

### Run Multiple Analyses

```bash
# Analysis 1: Default settings
python scanpy_pipeline_optimized.py --tenx_dir=... --project=exp1_default

# Analysis 2: Permissive QC
python scanpy_pipeline_optimized.py --tenx_dir=... --project=exp1_permissive \
  --min_genes=100 --max_mito_pct=15

# Analysis 3: High resolution
python scanpy_pipeline_optimized.py --tenx_dir=... --project=exp1_highres \
  --leiden_resolution=1.0
```

Result: Three separate folders for easy comparison:
```
output/
├── exp1_default/
├── exp1_permissive/
└── exp1_highres/
```

### Resume from Intermediate

```bash
# Use after_qc.h5ad from a previous project
python scanpy_pipeline_optimized.py \
  --input=output/exp1_default/data/after_qc.h5ad \
  --project=exp1_reanalysis \
  --leiden_resolution=0.8
```

## Isolation Guarantees

### Complete Separation of Concerns

| Component | Location | Isolation |
|-----------|----------|-----------|
| **Analysis Code** | `downstream_scanpy/` | Never modified by runs |
| **Project Results** | `output/<project>/` | Completely separate per project |
| **Cache Files** | `output/<project>/cache/` | Project-specific, no cross-contamination |
| **Temp Files** | `output/<project>/cache/` | Auto-cleaned per project |

### What This Means

✅ **Deleting a project** only removes that project's results  
✅ **Running multiple projects** in parallel won't interfere  
✅ **Cache files** are isolated - no data leakage between projects  
✅ **Code updates** don't affect existing project results  
✅ **Different users** can run their own projects safely  

### Example: Complete Independence

```bash
# Project A: High stringency
python scanpy_pipeline_optimized.py --tenx_dir=... --project=donorA_stringent

# Project B: Low stringency (same input)
python scanpy_pipeline_optimized.py --tenx_dir=... --project=donorA_permissive

# Both projects have their own:
# - cache/ (different .h5ad cache files)
# - data/ (different processed results)
# - figures/ (different plots)
# - logs/ (different execution logs)
```

**Deleting Project B does NOT affect Project A in any way.**

## Benefits

1. **Organization**: All related files in one place
2. **Reproducibility**: Config and logs preserved per project
3. **Comparison**: Easy to compare different parameter settings
4. **Sharing**: Simply zip the project folder to share
5. **Clean workspace**: No scattered result files
6. **Isolation**: Projects cannot contaminate each other

## Cleaning Up

To remove a specific project:
```bash
rm -rf output/<project_name>
```

To remove ALL results (start fresh):
```bash
rm -rf output/*
```

## Notes

- Project names should be descriptive and use underscores (e.g., `pbmc1k_v2`)
- If a project folder already exists, new results will be added/overwritten
- Each project is completely self-contained
- The `output/` directory is created automatically if it doesn't exist
