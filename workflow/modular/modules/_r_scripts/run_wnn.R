#!/usr/bin/env Rscript
# run_wnn.R -- standalone driver for Seurat WNN multimodal joint embedding.
#
# Called by the Python multimodal_integration module via subprocess. Writes a
# parquet (cell, WNN_1, WNN_2) at --out-parquet on success.
#
# CLI:
#   Rscript run_wnn.R \
#     --rna-parquet <path> \
#     --protein-parquet <path> \
#     --out-parquet <path>
#
# Exit codes:
#   0 -- success, out parquet written
#   2 -- Seurat unavailable, version too old, or any computation failure
#         (Python module catches non-zero and skips cleanly).

suppressWarnings(suppressMessages({
  parse_args <- function(argv) {
    out <- list()
    i <- 1L
    while (i <= length(argv)) {
      arg <- argv[[i]]
      if (arg == "--rna-parquet" && i < length(argv)) {
        out$rna_parquet <- argv[[i + 1L]]; i <- i + 2L
      } else if (arg == "--protein-parquet" && i < length(argv)) {
        out$protein_parquet <- argv[[i + 1L]]; i <- i + 2L
      } else if (arg == "--out-parquet" && i < length(argv)) {
        out$out_parquet <- argv[[i + 1L]]; i <- i + 2L
      } else {
        i <- i + 1L
      }
    }
    out
  }

  args <- parse_args(commandArgs(trailingOnly = TRUE))
  for (key in c("rna_parquet", "protein_parquet", "out_parquet")) {
    if (is.null(args[[key]]) || !nzchar(args[[key]])) {
      message(sprintf("[run_wnn.R] missing required arg --%s", gsub("_", "-", key)))
      quit(status = 2L, save = "no")
    }
  }
  if (!file.exists(args$rna_parquet)) {
    message(sprintf("[run_wnn.R] rna parquet not found: %s", args$rna_parquet))
    quit(status = 2L, save = "no")
  }
  if (!file.exists(args$protein_parquet)) {
    message(sprintf("[run_wnn.R] protein parquet not found: %s", args$protein_parquet))
    quit(status = 2L, save = "no")
  }

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    message("[run_wnn.R] Seurat is not installed; cannot run WNN.")
    quit(status = 2L, save = "no")
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    message("[run_wnn.R] arrow is not installed; cannot read/write parquet.")
    quit(status = 2L, save = "no")
  }

  # Seurat WNN (FindMultiModalNeighbors) requires Seurat >= 4.0.
  seurat_version <- tryCatch(utils::packageVersion("Seurat"),
                             error = function(e) NULL)
  if (is.null(seurat_version) ||
      utils::compareVersion(as.character(seurat_version), "4.0.0") < 0L) {
    message(sprintf(
      "[run_wnn.R] Seurat version %s is too old; need >= 4.0.0 for WNN.",
      as.character(seurat_version)
    ))
    quit(status = 2L, save = "no")
  }

  result <- tryCatch({
    rna_df <- as.data.frame(arrow::read_parquet(args$rna_parquet))
    prot_df <- as.data.frame(arrow::read_parquet(args$protein_parquet))

    cell_col <- function(df) {
      if ("cell" %in% colnames(df)) "cell"
      else if ("__index_level_0__" %in% colnames(df)) "__index_level_0__"
      else stop("input parquet missing 'cell' index column")
    }
    rna_cc <- cell_col(rna_df)
    prot_cc <- cell_col(prot_df)
    rna_cells <- as.character(rna_df[[rna_cc]])
    prot_cells <- as.character(prot_df[[prot_cc]])
    if (!identical(rna_cells, prot_cells)) {
      if (setequal(rna_cells, prot_cells)) {
        prot_df <- prot_df[match(rna_cells, prot_cells), , drop = FALSE]
        prot_cells <- as.character(prot_df[[prot_cc]])
      } else {
        stop("rna and protein parquets have different cell sets")
      }
    }

    rna_mat <- as.matrix(rna_df[, !colnames(rna_df) %in% c("cell", "__index_level_0__"), drop = FALSE])
    prot_mat <- as.matrix(prot_df[, !colnames(prot_df) %in% c("cell", "__index_level_0__"), drop = FALSE])
    storage.mode(rna_mat) <- "double"
    storage.mode(prot_mat) <- "double"
    rownames(rna_mat) <- rna_cells
    rownames(prot_mat) <- prot_cells

    # Build a minimal Seurat object: synthesize toy RNA counts so
    # CreateSeuratObject succeeds, then attach the user-provided RNA latent
    # and the second-modality latent as `DimReduc` slots so WNN can run on
    # them directly without re-running PCA / CLR.
    n_cells <- nrow(rna_mat)
    n_features_dummy <- 5L
    counts <- matrix(0, nrow = n_features_dummy, ncol = n_cells)
    rownames(counts) <- paste0("feat_", seq_len(n_features_dummy))
    colnames(counts) <- rna_cells
    seu <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA")

    # CreateDimReducObject requires a `key` ending in "_". We assign user
    # latents directly as embeddings without recomputing them.
    seu[["rna_latent"]] <- Seurat::CreateDimReducObject(
      embeddings = rna_mat, key = "rnaLatent_", assay = "RNA"
    )
    seu[["prot_latent"]] <- Seurat::CreateDimReducObject(
      embeddings = prot_mat, key = "protLatent_", assay = "RNA"
    )

    seu <- Seurat::FindMultiModalNeighbors(
      seu,
      reduction.list = list("rna_latent", "prot_latent"),
      dims.list = list(seq_len(ncol(rna_mat)), seq_len(ncol(prot_mat))),
      modality.weight.name = "RNA.weight",
      verbose = FALSE
    )
    seu <- Seurat::RunUMAP(
      seu, nn.name = "weighted.nn",
      reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
      verbose = FALSE
    )

    umap <- Seurat::Embeddings(seu, reduction = "wnn.umap")
    out_df <- data.frame(
      cell = rownames(umap),
      WNN_1 = umap[, 1],
      WNN_2 = umap[, 2],
      stringsAsFactors = FALSE
    )

    arrow::write_parquet(out_df, args$out_parquet)
    list(ok = TRUE)
  }, error = function(e) {
    message(sprintf("[run_wnn.R] WNN computation failed: %s", conditionMessage(e)))
    list(ok = FALSE)
  })

  if (!isTRUE(result$ok)) {
    quit(status = 2L, save = "no")
  }
  quit(status = 0L, save = "no")
}))
