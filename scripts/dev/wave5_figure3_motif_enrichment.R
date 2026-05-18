#!/usr/bin/env Rscript

options(warn = 1)

suppressPackageStartupMessages({
  library(Matrix)
})

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected positional argument: %s", key), call. = FALSE)
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    }
    out[[substring(key, 3)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

json_string <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\"', x)
  paste0('"', x, '"')
}

json_scalar <- function(x) {
  if (is.logical(x)) {
    return(ifelse(isTRUE(x), "true", "false"))
  }
  if (is.numeric(x)) {
    if (length(x) == 0 || is.na(x) || is.nan(x)) return("null")
    return(format(x, scientific = FALSE, trim = TRUE))
  }
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return("null")
  }
  json_string(as.character(x))
}

write_simple_json <- function(path, named_values) {
  lines <- c("{")
  n <- length(named_values)
  for (idx in seq_along(named_values)) {
    nm <- names(named_values)[[idx]]
    suffix <- if (idx < n) "," else ""
    lines <- c(lines, sprintf("  %s: %s%s", json_string(nm), json_scalar(named_values[[idx]]), suffix))
  }
  lines <- c(lines, "}")
  writeLines(lines, path)
}

coord_key <- function(chr, start, end) {
  paste(chr, as.character(start), as.character(end), sep = ":")
}

hypergeom_for_groups <- function(groupings, match_matrix, universe_label) {
  groups <- setdiff(sort(unique(groupings)), "background")
  result <- list()
  for (group in groups) {
    draw <- groupings == group
    bg <- !draw
    draws <- sum(draw)
    bg_draws <- sum(bg)
    if (draws < 10 || bg_draws < 10) {
      next
    }
    k <- Matrix::colSums(match_matrix[draw, , drop = FALSE])
    m <- Matrix::colSums(match_matrix[bg, , drop = FALSE])
    n <- bg_draws - m
    enrichment_p <- stats::phyper(k - 1, m, n, draws, lower.tail = FALSE)
    depletion_p <- stats::phyper(k, m, n, draws, lower.tail = TRUE)
    fg_rate <- k / draws
    bg_rate <- ifelse((m + n) > 0, m / (m + n), NA_real_)
    log2_fe <- log2((fg_rate + 1e-12) / (bg_rate + 1e-12))
    df <- data.frame(
      universe = universe_label,
      group = group,
      motif_id = colnames(match_matrix),
      motif_name = sub("^[^_]+_", "", colnames(match_matrix)),
      matches = as.numeric(k),
      draws = draws,
      bg_matches = as.numeric(m),
      bg_draws = bg_draws,
      enrichment_p = as.numeric(enrichment_p),
      depletion_p = as.numeric(depletion_p),
      p_value_two_sided = as.numeric(pmin(enrichment_p, depletion_p) * 2),
      log2_fold_enrichment = as.numeric(log2_fe),
      stringsAsFactors = FALSE
    )
    df$p_adj_bonferroni <- p.adjust(df$p_value_two_sided, method = "bonferroni")
    df$neg_log10_p_adj <- -log10(pmax(df$p_adj_bonferroni, .Machine$double.xmin))
    result[[length(result) + 1]] <- df
  }
  do.call(rbind, result)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("table-s3b", "bed", "motif-rds", "out-dir")
missing <- setdiff(required, names(args))
if (length(missing) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

out_dir <- args[["out-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cluster_col <- if ("cluster-col" %in% names(args)) args[["cluster-col"]] else "Link cluster"

table_s3b <- read.delim(args[["table-s3b"]], check.names = FALSE, stringsAsFactors = FALSE)
bed <- read.delim(gzfile(args[["bed"]]), header = FALSE, stringsAsFactors = FALSE)
motif <- readRDS(args[["motif-rds"]])

if (nrow(bed) != nrow(motif)) {
  stop(sprintf("BED rows (%d) != motif rows (%d)", nrow(bed), nrow(motif)), call. = FALSE)
}
if (!inherits(motif, "Matrix")) {
  stop(sprintf("Motif object must inherit Matrix, got: %s", paste(class(motif), collapse = ",")), call. = FALSE)
}

bed_key <- coord_key(bed[[1]], bed[[2]], bed[[3]])
names(bed_key) <- rownames(motif)
key_to_peak_id <- setNames(rownames(motif), bed_key)

table_key <- coord_key(table_s3b[["Peak chromosome"]], table_s3b[["Peak start"]], table_s3b[["Peak end"]])
table_s3b$peak_id <- unname(key_to_peak_id[table_key])
mapped_rows <- sum(!is.na(table_s3b$peak_id))
unique_mapped_peaks <- length(unique(table_s3b$peak_id[!is.na(table_s3b$peak_id)]))

if (!cluster_col %in% names(table_s3b)) {
  stop(sprintf("Table S3B does not contain cluster column: %s", cluster_col), call. = FALSE)
}
table_s3b[[cluster_col]] <- as.character(table_s3b[[cluster_col]])

cluster_peak_df <- unique(table_s3b[!is.na(table_s3b$peak_id), c(cluster_col, "peak_id")])
names(cluster_peak_df) <- c("cluster", "peak_id")
cluster_peak_df <- cluster_peak_df[order(cluster_peak_df$cluster, cluster_peak_df$peak_id), ]

cluster_counts <- as.data.frame(table(cluster_peak_df$cluster), stringsAsFactors = FALSE)
names(cluster_counts) <- c("cluster", "unique_peak_count")
cluster_counts <- cluster_counts[order(cluster_counts$cluster), ]
write.table(cluster_counts, file.path(out_dir, "figure3_table_s3b_cluster_peak_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(table_s3b, file.path(out_dir, "figure3_table_s3b_with_peak_ids.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Universe 1: the exact Table S3B linked-peak universe, one row per unique cluster/peak.
idx <- match(cluster_peak_df$peak_id, rownames(motif))
linked_match <- motif[idx, , drop = FALSE]
rownames(linked_match) <- cluster_peak_df$peak_id
linked_enr <- hypergeom_for_groups(cluster_peak_df$cluster, linked_match, "table_s3b_linked_peaks")

# Universe 2: each Table S3B cluster versus all other consensus peaks. This is a diagnostic
# background, not a claim that the paper used the full consensus universe.
all_rows <- seq_len(nrow(motif))
full_enr_parts <- list()
for (cluster in sort(unique(cluster_peak_df$cluster))) {
  peak_ids <- cluster_peak_df$peak_id[cluster_peak_df$cluster == cluster]
  foreground_idx <- match(unique(peak_ids), rownames(motif))
  groupings <- rep("background", nrow(motif))
  groupings[foreground_idx] <- cluster
  # Put foreground first only to keep output grouping deterministic; sparse row slicing is cheap enough here.
  full_enr_parts[[length(full_enr_parts) + 1]] <- hypergeom_for_groups(groupings, motif, "all_consensus_peaks")
}
full_enr <- do.call(rbind, full_enr_parts)
full_enr <- full_enr[full_enr$group != "background", ]

enrichment <- rbind(linked_enr, full_enr)
enrichment <- enrichment[order(enrichment$universe, enrichment$group, enrichment$p_adj_bonferroni, -enrichment$log2_fold_enrichment), ]
write.table(enrichment, file.path(out_dir, "figure3_table_s3b_motif_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

top <- do.call(rbind, lapply(split(enrichment, list(enrichment$universe, enrichment$group), drop = TRUE), function(df) {
  df <- df[order(df$p_adj_bonferroni, -df$log2_fold_enrichment), ]
  head(df, 10)
}))
write.table(top, file.path(out_dir, "figure3_table_s3b_top_motifs_by_cluster.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

top_enriched <- do.call(rbind, lapply(split(enrichment[enrichment$log2_fold_enrichment > 0, ],
                                            list(enrichment$universe[enrichment$log2_fold_enrichment > 0],
                                                 enrichment$group[enrichment$log2_fold_enrichment > 0]),
                                            drop = TRUE), function(df) {
  df <- df[order(df$enrichment_p, df$p_adj_bonferroni, -df$log2_fold_enrichment), ]
  head(df, 10)
}))
write.table(top_enriched, file.path(out_dir, "figure3_table_s3b_top_enriched_motifs_by_cluster.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write_simple_json(
  file.path(out_dir, "figure3_table_s3b_motif_readiness.json"),
  list(
    table_s3b_rows = nrow(table_s3b),
    table_s3b_rows_mapped_to_peak_id = mapped_rows,
    table_s3b_unique_mapped_peaks = unique_mapped_peaks,
    motif_matrix_rows = nrow(motif),
    motif_matrix_cols = ncol(motif),
    bed_rows = nrow(bed),
    all_table_rows_mapped = mapped_rows == nrow(table_s3b),
    cluster_count = length(unique(cluster_peak_df$cluster)),
    status = ifelse(mapped_rows == nrow(table_s3b), "READY_MOTIF_ENRICHMENT_DIAGNOSTIC", "MAPPING_GAP")
  )
)

cat(sprintf("mapped_rows=%d/%d unique_peaks=%d clusters=%d motifs=%d\n",
            mapped_rows, nrow(table_s3b), unique_mapped_peaks,
            length(unique(cluster_peak_df$cluster)), ncol(motif)))
cat(sprintf("wrote=%s\n", normalizePath(out_dir)))
