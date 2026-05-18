# Single-Cell RNA Sequencing Figure Drawing Pipeline

> **R sources are symlinked from `r_multiomics_factory/`.
> Do NOT edit files under `bridges/local_r_pipeline_macbook/R/` or `bridges/local_r_pipeline_macbook/R_bundle/` directly;
> edit the upstream files in `r_multiomics_factory/R/` or `r_multiomics_factory/R_bundle/` and the symlinks will reflect changes automatically.**

## Current Status: Remote R Pipeline

This folder name is historical. Despite `local_r_pipeline_macbook` in the path,
the current maintained plotting/reporting workflow runs on the remote server.

- Remote root: `/home/zerlinshen/singlecell_factory`
- Remote R runtime: `/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript`
- Legacy rollback runtime: `/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript`
- Remote Python runtime for bundle export: `/home/zerlinshen/conda/bin/conda run -n sc_gpu`
- Local Mac role: review/organization only, not the maintained R plotting runtime

## Remote R Environment

The maintained remote R environment for v2 parquet bundle plotting is:

- `/home/zerlinshen/conda/envs/r_multiomics_arrow`

The previous environment remains available for rollback:

- `/home/zerlinshen/conda/envs/r_multiomics`

As of `2026-04-30`, `r_multiomics_arrow` is the default wrapper runtime. It is a
clone of `r_multiomics` with conda-forge `r-arrow 24.0.0` / `libarrow 24.0.0`
added and was validated through `read_bundle()` plus
`scripts/plot_remote_bundle_large.R` on both NC2024 v2 bundles. Installed and
verified packages include:

- core object/plotting: `Seurat 5.4.0`, `SeuratObject 5.4.0`, `ggplot2`,
  `data.table`, `readr`, `patchwork`, `cowplot`, `ggrepel`, `viridis`
- large plotting: `ggrastr`, `scattermore`
- heatmaps/reporting: `pheatmap`, `ComplexHeatmap`, `circlize`
- bridge/integration helpers: `hdf5r`, `zellkonverter`, `harmony`,
  `BiocManager`, `R.utils`, `remotes`, `arrow`

Validation artifact:

- `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/r_env_dependency_smoke_20260424`

That smoke test used the existing manifest-backed NC2024 `r_bundle`, read all
`810218` rows of metadata/UMAP/marker bundle data, plotted a `100000`-cell
raster UMAP via `ggrastr`, generated `pheatmap` and `ComplexHeatmap` outputs,
and opened `final_adata.h5ad` via `hdf5r`.

`SeuratDisk` is not installed by default. The conda-forge build currently
conflicts with this R 4.5 environment through old `spatstat` requirements. For
small H5AD bridge checks, prefer `zellkonverter` or direct `hdf5r` inspection.
For NC2024-scale figure work, prefer the compact bundle path instead of direct
full-object Seurat conversion.

For large NC2024-style cohorts, use `scripts/run_remote_bundle_plot.sh`. The old
`scripts/pull_and_plot_remote_result.sh` name is kept only as a compatibility
wrapper and delegates to remote plotting.

## Factory-To-R Bridge Contract

The remote bridge is designed to let `singlecell_factory` do heavy computation
and sparse/full-object access, then let R focus on publication-style figures
from compact, manifest-backed tables.

Key controls:

- `R_PLOT_THREADS`: thread budget for R-side BLAS/data.table work. Default: `8`.
- `R_BUNDLE_MARKERS`: comma-separated marker genes to export from AnnData `X`
  for R plotting. The marker-dot plot uses all numeric marker columns present in
  `marker_expr.csv.gz`, so ROI-specific marker requests are plotted directly.
- `R_BUNDLE_OBS_COLS`: comma-separated `.obs` columns to add to the plotting
  bundle. The wrapper always keeps the selected `group_by` and `cluster_by`
  columns so plotting cannot drop its own grouping fields.
- `R_BUNDLE_OBSM`: comma-separated `.obsm` embeddings to add to the plotting
  bundle. The wrapper always keeps `X_umap` and `X_pca` so plotting and reuse
  validation keep their required stems.
- `PLOT_SCRIPT`: R plotting script used after bundle validation. Default:
  `/home/zerlinshen/r_multiomics_factory/scripts/plot_remote_bundle_large.R`,
  which is the v1/v2-aware upstream plotting entry point.
- `FORCE_R_BUNDLE_EXPORT=1`: force bundle regeneration even if the reuse guard
  passes.

Example:

```bash
cd /home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook
R_PLOT_THREADS=8 \
R_BUNDLE_MARKERS=ELF3,EPCAM,KRT8,KRT18,PTPRC,CD3E,LYZ,MS4A1,NKG7 \
bash scripts/run_remote_bundle_plot.sh \
  /home/zerlinshen/singlecell_factory/results/<run> \
  /home/zerlinshen/singlecell_factory/results/<run>/r_plots/main \
  cell_type leiden 200000
```

Reuse is intentionally conservative. The wrapper reuses an existing bundle only
when:

- manifest source paths match the requested run
- manifest `input_h5ad_bytes` and `input_h5ad_mtime_epoch` match the current
  `final_adata.h5ad`
- optional marker/obs/obsm requests match the manifest
- R-side validation passes required-file, byte-size, SHA256, schema, expression
  semantics, and cell-alignment checks

The optional obs/obsm controls are additive, not replacement controls. For
example, `R_BUNDLE_OBS_COLS=sample,patient` with `group_by=cell_type` and
`cluster_by=leiden` exports `cell_type,leiden,sample,patient`. Likewise,
`R_BUNDLE_OBSM=X_tsne` still exports `X_umap,X_pca,X_tsne`.

Plotting is optimized for large visual summaries: UMAP point layers use
`ggrastr` rasterization when available, while heatmaps and marker summaries stay
table-based. This keeps R focused on figure quality without forcing a full
Seurat conversion of the NC2024 object.

Latest bridge hardening validation:

- export + prettier plot smoke:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_pretty_export_20260424`
- validated reuse smoke:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_pretty_reuse_20260424`
- canonical default-bundle refresh:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_refresh_20260424`
- canonical default-bundle reuse:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_reuse_20260424`

This repository contains a modular R pipeline designed to accelerate and standardize figure drawing for single-cell RNA sequencing (scRNA-seq) analysis. The pipeline is built on top of the popular `Seurat` and `ggplot2` ecosystems, providing a consistent, publication-ready aesthetic across all plots.

## Architecture

The pipeline is organized into multiple functional modules located in the `R/` directory. These modules separate different stages of scRNA-seq analysis and define a unified style.

## Module provenance and references

This section tracks modules that are thin wrappers around public, recognized packages (including GitHub/CRAN packages) and the standard articles/tools where those methods were introduced.

| Module | Origin source | How it is used here | Reference |
|---|---|---|---|
| `R/io_bridge.R` | Seurat object API + AnnData bridge (`zellkonverter` / `SeuratDisk`) | Load/attach data for Linux-exported objects, `singlecell_factory` `.h5ad` handoff, and external metadata | Satija Lab (`Seurat`): https://satijalab.org/seurat/ ; Bioconductor `zellkonverter`: https://bioconductor.org/packages/zellkonverter/ |
| `R/preprocessing_module.R` | Seurat + sctransform | QC, normalization, PCA, clustering, UMAP/t-SNE | Satija et al., Seurat paper; Hafemeister & Satija, SCTransform |
| `R/qc_plots.R` | Seurat + ggplot2 | QC visual diagnostics for counts/features/mitochondrial reads | Wickham, ggplot2 |
| `R/dim_plots.R` | Seurat + ggplot2 | Dimensionality reduction plot wrappers | Satija Lab, Seurat |
| `R/expression_plots.R` | Seurat + ggplot2 + viridis | Feature, dot, and heatmap visualizations | Satija Lab, Seurat; Neuwirth (viridis palette) |
| `R/composition_plots.R` | dplyr/tidyr + ggplot2 | Cell composition summaries | Tidyverse ecosystem |
| `R/marker_module.R` | Seurat | Differential expression and marker summaries | Seurat marker workflows |
| `R/pipeline_steps.R` | Pipeline orchestration (internal) | Plot orchestration and reduction guards | Internal module |
| `R/batch_integration_module.R` | Seurat + Harmony | Batch correction via CCA/RPCA and Harmony correction | Stuart et al., 2019 (Seurat integration); Korsunsky et al., 2019 (Harmony) |
| `R/annotation_module.R` | Seurat + AddModuleScore | Lightweight module-score-based annotation | Seurat marker scoring workflow |
| `R/remote_bundle_manifest.R` | Internal bridge contract | Validate compact `singlecell_factory` R bundles, required files, byte sizes, SHA256 hashes, and provenance manifest before plotting | Internal module |

If you add or modify modules that come from external package functionality, please append a row to this table in the same format.
See [MODULE_UPDATE_SKILL.md](MODULE_UPDATE_SKILL.md) for the standard maintenance checklist before each module update.
You can also run `Rscript scripts/check_module_updates.R` to validate changed R modules are registered in README and dependencies are consistent.

## Methods references supporting this pipeline

These references document the methods behind the modules currently implemented in this pipeline:

### Core preprocessing and clustering modules
- Satija R, Farrell JA, et al. Seurat: Seurat toolkit for single-cell genomics. *Cell* (2015) and associated Seurat workflow documentation.
- Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology* 2019.
- Stuart T, Butler A, et al. Comprehensive integration of single-cell data. *Cell* 2019.
- Korsunsky I, et al. Fast, sensitive, and accurate integration of single-cell data with Harmony. *Nature Methods* 2019.

### Marker detection and scoring
- Butler A, et al. Differential expression and marker discovery as used in Seurat marker workflows (`FindAllMarkers`/`FindMarkers`).
- Seurat AddModuleScore-based scoring (module-score-based annotation logic).

### Trajectory/normalization and visual design
- Wickham H. ggplot2: Elegant Graphics for Data Analysis.
- Neuwirth E. viridis - Perceptually uniform colormaps for scientific visualizations.

### Reproducibility and reporting
- This pipeline records all module changes through `[MODULE_UPDATE_SKILL.md](MODULE_UPDATE_SKILL.md)` and verifies module-documentation synchronization via:
  - `Rscript scripts/check_module_updates.R`

### Skill: update module documentation on every module change

When adding/updating modules, we use the following mandatory checklist:

1. Add/modify a brief entry in `README` module list.
2. If external package functionality is used, add/update that row under "Module provenance and references".
3. Add/refresh a concise usage example in the "How to Run" section.
4. Record optional required package dependencies in `main.R` `required_pkgs`.

If you skip step 2, downstream reproducibility checks may fail in manuscript workflows.

### Modules Built:

1. **`R/theme_config.R`**: 
   - **What it does**: Establishes the core visual identity of your plots. It provides `theme_sc_publication()`, a clean `ggplot2` theme. It also defines centralized functions for fetching continuous (`viridis`) and discrete (`RColorBrewer`) color palettes.

2. **`R/qc_plots.R`**: 
   - **What it does**: Handles Quality Control visualizations. It contains functions like `plot_qc_violin()` to easily visualize `nFeature_RNA`, `nCount_RNA`, and `percent.mt` per sample, and `plot_qc_scatter()` to inspect correlations between these metrics.

3. **`R/dim_plots.R`**: 
   - **What it does**: Streamlines Dimensionality Reduction plots. It provides `plot_umap_custom()`, heavily customizing Seurat's default `DimPlot` to remove axis lines/ticks for a cleaner look and automatically apply the unified color schemes for your clusters.

4. **`R/expression_plots.R`**: 
   - **What it does**: Visualizes marker gene expression across cells and clusters. It includes `plot_feature_custom()` for specialized feature plots, `plot_dot_custom()` for structured dot plots of multi-gene markers, and `plot_heatmap_custom()` for global heatmaps.

5. **`R/composition_plots.R`**: 
   - **What it does**: Analyzes cellular composition. It includes `plot_cell_fractions()`, which automatically aggregates cell cluster counts per sample/condition and returns a 100% stacked bar chart, making it easy to see proportion shifts between biological groups.
6. **`R/preprocessing_module.R`**:
   - **What it does**: Provides optional end-to-end preprocessing (QC filtering, normalization, variable feature selection, PCA, neighbor graph, clustering, UMAP/t-SNE) for imported objects.
7. **`R/marker_module.R`**:
   - **What it does**: Provides marker analysis (`FindAllMarkers`/`FindMarkers`) and marker report outputs (`*_all.csv`, `*_top.csv`, volcano plots, optional heatmap).
8. **`R/batch_integration_module.R`**:
   - **What it does**: Adds optional batch/sample integration using Seurat anchor-based workflows (`cca`/`rpca`) or Harmony correction (`harmony`).
9. **`R/annotation_module.R`**:
   - **What it does**: Adds lightweight marker-set annotation by module scoring and writes per-cell annotation plus summary outputs.
10. **`R/io_bridge.R`** and **`R/cli_utils.R`**:
   - **What it does**: Handle Linux object import, metadata attachment, argument parsing, and marker-file parsing.
11. **`R/remote_bundle_manifest.R`**:
   - **What it does**: Validates compact remote handoff bundles before plotting. This is the preferred large-cohort bridge because it checks required files and schema provenance without loading the full expression matrix into R.

## Getting Started

To test this pipeline or adapt it for your own data, look at `main.R`. 

`main.R` is the executable pipeline entrypoint and now supports both the built-in demo and a Linux-handoff mode:
1. Loads all necessary dependencies and the custom modules.
2. Loads a Seurat object from Linux outputs (RDS / RDA / H5Seurat / H5AD) and optional metadata.
3. Iteratively calls the custom functions from each module.
4. Optionally runs preprocessing and marker analysis (off by default) before plotting.
5. Optionally applies QC filtering and saves:
   - one multi-page PDF (`results_figures.pdf` by default)
   - PNG previews of each panel.

If you run in demo mode (no `--input` argument), it will load `pbmc_small` and generate the same figure set.

### How to Run

1. Open your R console, terminal, or RStudio within this directory.
2. Install packages (if missing) and run:
 ```r
   source("main.R")
 ```
3. Or run from CLI:
```bash
 Rscript main.R --input /path/to/seurat_obj.rds --metadata /path/to/meta.csv --marker-file /path/to/markers.txt --out-dir figures
 ```
4. Recommended large-cohort handoff from `singlecell_factory`:
```bash
 cd /home/zerlinshen/singlecell_factory
 bash bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh \
   /home/zerlinshen/singlecell_factory/results/<run> \
   /home/zerlinshen/singlecell_factory/results/<run>/r_plots/main \
   cell_type leiden
```
The wrapper exports or reuses a compact bundle, validates either v1 CSV/TSV or
v2 parquet/JSON manifests, and then delegates plotting to the upstream
`r_multiomics_factory/scripts/plot_remote_bundle_large.R` entry point. Bundle
payloads copy metadata, embeddings, and selected marker genes only; they do not
convert the full expression matrix to dense data. The R validator checks
required file contracts, byte sizes, SHA256 hashes, schema, and expression
semantics before plotting. Marker values are exported from AnnData `X` as
stored (`source_X_as_stored`) and are intended for plotting/visual summaries
only, not new DE or quantitative expression claims without full-object
validation.

5. One-command remote plotting wrapper:
```bash
 cd /home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook
 bash scripts/run_remote_bundle_plot.sh \
   /home/zerlinshen/singlecell_factory/results/<run> \
   /home/zerlinshen/singlecell_factory/results/<run>/r_plots/main \
   cell_type leiden
```
This wrapper generates or refreshes the compact bundle and runs R on the remote
server. By default, it reuses an existing bundle only when the manifest source
paths match the requested run, the manifest is newer than `final_adata.h5ad`,
and the R validator passes byte-size/SHA256 integrity checks. Set
`FORCE_R_BUNDLE_EXPORT=1` to force a fresh export. The plotting scripts require
a manifest-backed bundle and fail early if required files, file integrity
checks, schema, expression semantics, or cell alignment are missing.

The legacy `scripts/pull_and_plot_remote_result.sh` entrypoint no longer pulls
data to a Mac R pipeline; it delegates to `run_remote_bundle_plot.sh` and writes
remote plots under the requested remote output directory.

6. Small/medium object import remains available when a full Seurat object is explicitly needed:
```bash
 Rscript main.R --input /path/to/final_adata.h5ad --out-dir figures --group-by cell_type --cluster-by leiden
```
This `.h5ad` path is typically the `final_adata.h5ad` produced under `singlecell_factory/results/<project_timestamp>/`. For large cohorts such as NC2024, prefer the compact bundle path above because direct `.h5ad` conversion can force large in-memory objects inside R/Seurat.
If you use `.h5ad`, install one of the optional bridges first:
```r
 BiocManager::install("zellkonverter")
 # or
 remotes::install_github("mojaveazure/seurat-disk")
```
7. Run preprocessing + clustering + markers:
 ```bash
 Rscript main.R --input /path/to/seurat_obj.rds --run-preprocess true --preprocess-cluster true \
   --run-markers true --marker-group-by seurat_clusters --out-dir figures --output-pdf comprehensive.pdf
 ```
8. Run batch integration and annotation:
 ```bash
 Rscript main.R --input /path/to/seurat_obj.rds --run-integration true --integration-batch-col batch \
   --integration-method rpca --run-annotation true --annotation-reference ./celltype_reference.csv \
   --annotation-assay RNA --out-dir figures --output-pdf comprehensive.pdf
 ```
9. Check the output directory for `results_figures.pdf` and individual PNG files.

### New Modules

The pipeline now includes:
- `R/cli_utils.R`  
  Parses command-line options and marker list files.
- `R/io_bridge.R`  
  Loads Seurat objects from Linux-compatible formats and merges metadata by cell barcode.
- `R/pipeline_steps.R`  
  Orchestrates QC filtering, checks/recomputes dimensionality reduction, and runs all plot generators.
