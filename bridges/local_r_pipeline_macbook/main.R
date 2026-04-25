## =========================================
##  ScRNA-seq Figure Pipeline (Linux -> R)
## =========================================
##
## Run examples:
##   Rscript main.R --help
##   Rscript main.R --input data/seurat_obj.rds
##   Rscript main.R --input data/seurat_obj.rds --metadata data/meta.csv --marker-file data/markers.txt --out-dir figures

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "viridis",
                   "scales", "Seurat", "patchwork", "tools")
optional_input_pkgs <- c("zellkonverter", "SeuratDisk")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, repos = "http://cran.us.r-project.org")
}

source("R/cli_utils.R")
source("R/theme_config.R")
source("R/io_bridge.R")
source("R/preprocessing_module.R")
source("R/qc_plots.R")
source("R/dim_plots.R")
source("R/expression_plots.R")
source("R/composition_plots.R")
source("R/marker_module.R")
source("R/batch_integration_module.R")
source("R/annotation_module.R")
source("R/pipeline_steps.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
})

cat("Pipeline modules loaded successfully.\n")

args <- parse_pipeline_args()
if (isTRUE(args$help)) {
  print_pipeline_usage()
  quit(save = "no", status = 0)
}

if (!is.null(args$input)) {
  seurat_obj <- load_scRNA_object(
    path = args$input,
    metadata_path = args$metadata,
    metadata_cell_col = args$metadata_cell_col
  )
} else {
  cat("No --input provided; running built-in pbmc_small demo.\n")
  data("pbmc_small")
  seurat_obj <- pbmc_small
}

if (!is.null(args$marker_file)) {
  marker_genes <- read_marker_file(args$marker_file)
} else {
  marker_genes <- args$markers
}

if (isTRUE(args$run_preprocess)) {
  seurat_obj <- preprocess_scRNA_object(
    seurat_obj = seurat_obj,
    skip_existing = args$preprocess_skip_existing,
    min_features = args$qc_min_features,
    max_features = args$qc_max_features,
    max_mito = args$qc_max_mito,
    normalization = args$normalization,
    variable_features = args$preprocess_nfeatures,
    n_pcs = args$preprocess_pcs,
    run_clusters = args$preprocess_cluster,
    resolution = args$preprocess_resolution,
    run_umap = args$preprocess_umap,
    run_tsne = args$preprocess_tsne,
    dims = args$dims,
    seed = args$preprocess_seed
  )
} else {
  seurat_obj <- filter_sc_object_by_qc(
    seurat_obj = seurat_obj,
    min_features = args$qc_min_features,
    max_features = args$qc_max_features,
    max_mito = args$qc_max_mito
  )
}

if (isTRUE(args$run_integration)) {
  seurat_obj <- run_batch_integration(
    seurat_obj = seurat_obj,
    batch_col = args$integration_batch_col,
    method = args$integration_method,
    nfeatures = args$integration_nfeatures,
    dims = args$integration_dims,
    k_anchor = 5L
  )
}

if (isTRUE(args$run_annotation)) {
  if (is.null(args$annotation_reference)) {
    stop("--run-annotation true requires --annotation-reference.")
  }
  seurat_obj <- run_celltype_annotation(
    seurat_obj = seurat_obj,
    reference_path = args$annotation_reference,
    assay = args$annotation_assay,
    out_dir = args$out_dir,
    score_min = args$annotation_min_score
  )
}

if (isTRUE(args$run_markers)) {
  marker_result <- run_marker_analysis(
    seurat_obj = seurat_obj,
    out_dir = args$out_dir,
    group_by = args$marker_group_by,
    ident_1 = args$marker_ident_1,
    ident_2 = args$marker_ident_2,
    assay = args$marker_assay,
    min_pct = args$marker_min_pct,
    logfc_threshold = args$marker_logfc_threshold,
    only_pos = args$marker_only_pos,
    test_use = args$marker_test_use,
    top_n = args$marker_top_n,
    volcano_fdr = args$volcano_fdr,
    volcano_logfc = args$volcano_logfc,
    output_prefix = "marker_analysis"
  )

  if (!is.null(marker_result$top) && nrow(marker_result$top) > 0L) {
    marker_genes <- unique(marker_result$top$gene)
    marker_genes <- marker_genes[marker_genes %in% rownames(seurat_obj)]
  }
}

run_plot_suite(
  seurat_obj = seurat_obj,
  out_dir = args$out_dir,
  output_pdf = args$output_pdf,
  group_by = args$group_by,
  cluster_by = args$cluster_by,
  split_by = args$split_by,
  reduction = args$reduction,
  markers = marker_genes,
  marker_mode = if (isTRUE(args$run_markers)) "marker_analysis" else NULL,
  save_png = args$save_png,
  dims = args$dims
)

cat("Pipeline finished.\n")
