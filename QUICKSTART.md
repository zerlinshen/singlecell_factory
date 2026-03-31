# Quick Start Guide

## Clean Environment Setup

All previous results have been removed. Start fresh with the new project-based structure.

## Run the Pipeline

### 1. Optimized Pipeline (Recommended)

```bash
cd /home/zerlinshen/singlecell_factory/downstream_scanpy

# Basic run
/home/zerlinshen/conda/envs/sc10x/bin/python scanpy_pipeline_optimized.py \
  --tenx_dir=/home/zerlinshen/singlecell_factory/pbmc_1k_v3_count/outs/filtered_feature_bc_matrix \
  --project=pbmc1k_default
```

**Output location**: `output/pbmc1k_default/`

### 2. Review Adaptive QC Suggestions

The pipeline will analyze your data and suggest thresholds. Look for:
```
ADAPTIVE QC THRESHOLD SUGGESTIONS (MAD-based)
Current config -> min_genes: 200, max_mito: 10.0%
Suggested      -> min_genes: 156, max_genes: 5234, max_mito: 12.5%
Expected retention: 45.2% vs Current: 12.4%
```

### 3. Compare Different Settings

```bash
# Run with suggested (permissive) thresholds
/home/zerlinshen/conda/envs/sc10x/bin/python scanpy_pipeline_optimized.py \
  --tenx_dir=/home/zerlinshen/singlecell_factory/pbmc_1k_v3_count/outs/filtered_feature_bc_matrix \
  --project=pbmc1k_permissive \
  --min_genes=156 \
  --max_mito_pct=12.5

# Run with high resolution clustering
/home/zerlinshen/conda/envs/sc10x/bin/python scanpy_pipeline_optimized.py \
  --tenx_dir=/home/zerlinshen/singlecell_factory/pbmc_1k_v3_count/outs/filtered_feature_bc_matrix \
  --project=pbmc1k_highres \
  --leiden_resolution=1.0
```

Compare results in:
```
output/
├── pbmc1k_default/      # Default settings
├── pbmc1k_permissive/   # More cells retained
└── pbmc1k_highres/      # More clusters
```

### 4. Simple Pipeline (Fast)

For quick analysis without doublet detection:

```bash
/home/zerlinshen/conda/envs/sc10x/bin/python run_scanpy_pipeline.py \
  --tenx_dir=/home/zerlinshen/singlecell_factory/pbmc_1k_v3_count/outs/filtered_feature_bc_matrix \
  --project=pbmc1k_simple
```

## New Project Structure

All results are organized under `output/<project_name>/`:

```
output/
└── <project_name>/
    ├── data/           # .h5ad files
    ├── figures/        # Plots
    │   ├── qc/         # QC visualizations
    │   ├── analysis/   # PCA/UMAP
    │   └── markers/    # Marker genes
    ├── tables/         # CSV files
    └── logs/           # Logs and configs
```

## Key Improvements

| Feature | Before | After |
|---------|--------|-------|
| Output location | Scattered `results/` folders | `output/<project>/` |
| Multiple runs | Overwrite previous | Separate projects |
| QC thresholds | Fixed | Adaptive suggestions |
| Doublet detection | None | Scrublet integrated |
| Logs | Console only | File + console |

## Next Steps

1. Run with default settings
2. Check adaptive QC suggestions in the log
3. Re-run with adjusted thresholds if needed
4. Compare multiple projects side-by-side
