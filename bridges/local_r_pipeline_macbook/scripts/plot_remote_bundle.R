suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

source('R/theme_config.R')
source('R/qc_plots.R')
source('R/dim_plots.R')
source('R/expression_plots.R')
source('R/composition_plots.R')
source('R/pipeline_steps.R')
source('R_bundle/remote_bundle_manifest.R')

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop('Usage: Rscript scripts/plot_remote_bundle.R <bundle_dir> <out_dir> [group_by] [cluster_by]')
}

bundle_dir <- args[[1]]
out_dir <- args[[2]]
group_by <- if (length(args) >= 3) args[[3]] else 'cell_type'
cluster_by <- if (length(args) >= 4) args[[4]] else 'leiden'

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lock_path <- file.path(out_dir, ".plot_remote_bundle.lock")
if (file.exists(lock_path)) {
  old_pid <- suppressWarnings(as.integer(readLines(lock_path, warn = FALSE, n = 1)))
  if (!is.na(old_pid)) {
    alive <- suppressWarnings(system2("kill", c("-0", as.character(old_pid)), stdout = FALSE, stderr = FALSE))
    if (identical(alive, 0L)) {
      stop(sprintf("Another plot_remote_bundle.R process is already running for %s (pid=%s)", out_dir, old_pid))
    }
  }
}
writeLines(as.character(Sys.getpid()), lock_path)
on.exit({
  if (file.exists(lock_path)) unlink(lock_path, force = TRUE)
}, add = TRUE)

bundle <- validate_remote_bundle(
  bundle_dir = bundle_dir,
  required_stems = c('obs', 'marker_expr', 'X_pca', 'X_umap'),
  allow_missing_manifest = FALSE
)

obs <- read_bundle_csv(bundle$files[['obs']])
expr <- read_bundle_csv(bundle$files[['marker_expr']])
pca <- read_bundle_csv(bundle$files[['X_pca']])
umap <- read_bundle_csv(bundle$files[['X_umap']])

validate_bundle_table_dimensions(bundle, list(obs = obs, marker_expr = expr, X_pca = pca, X_umap = umap))
validate_bundle_cell_alignment(list(obs = obs, marker_expr = expr, X_pca = pca, X_umap = umap))

stop(
  "plot_remote_bundle.R: the singlecell_factory exporter writes bundles with ",
  "expression_value_scale='source_X_as_stored', not raw_counts. ",
  "Building a Seurat object and re-normalising from a compact bundle is not supported here. ",
  "For publication-quality figures from a compact bundle, use: ",
  "  Rscript scripts/plot_remote_bundle_large.R <bundle_dir> <out_dir> [group_by] [cluster_by] ",
  "For full Seurat analysis from scratch, load the .h5ad directly via main.R: ",
  "  Rscript multiomics_r_factory/main.R --input <final_adata.h5ad> --out-dir <out_dir>"
)
