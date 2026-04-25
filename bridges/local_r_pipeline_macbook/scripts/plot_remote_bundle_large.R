suppressPackageStartupMessages({
  library(ggplot2)
})

source("R_bundle/remote_bundle_manifest.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/plot_remote_bundle_large.R <bundle_dir> <out_dir> [group_by] [cluster_by] [max_cells]")
}

bundle_dir <- args[[1]]
out_dir <- args[[2]]
group_by <- if (length(args) >= 3) args[[3]] else "cell_type"
cluster_by <- if (length(args) >= 4) args[[4]] else "leiden"
max_cells <- if (length(args) >= 5) as.integer(args[[5]]) else 200000L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (is.na(max_cells) || max_cells <= 0L) {
  stop("max_cells must be a positive integer")
}

if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(as.integer(Sys.getenv("R_DATATABLE_NUM_THREADS", "8")))
}

bundle <- validate_remote_bundle(
  bundle_dir = bundle_dir,
  required_stems = c("obs", "X_umap", "marker_expr"),
  allow_missing_manifest = FALSE
)

obs <- read_bundle_csv(bundle$files[["obs"]])
umap <- read_bundle_csv(bundle$files[["X_umap"]])
expr <- read_bundle_csv(bundle$files[["marker_expr"]])

validate_bundle_table_dimensions(bundle, list(obs = obs, X_umap = umap, marker_expr = expr))
validate_bundle_cell_alignment(list(obs = obs, X_umap = umap, marker_expr = expr))

if (!group_by %in% colnames(obs)) {
  if ("cell_type_hint" %in% colnames(obs)) {
    obs[[group_by]] <- obs[["cell_type_hint"]]
  } else if ("CellName" %in% colnames(obs)) {
    obs[[group_by]] <- obs[["CellName"]]
  } else {
    stop(sprintf("Missing group_by column: %s", group_by))
  }
}

if (!cluster_by %in% colnames(obs)) {
  stop(sprintf("Missing cluster_by column: %s", cluster_by))
}

obs[[group_by]] <- as.factor(obs[[group_by]])
obs[[cluster_by]] <- as.factor(obs[[cluster_by]])

theme_nc2024 <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#e8e2d7", linewidth = 0.2),
      plot.title = element_text(face = "bold", color = "#1f2933"),
      plot.subtitle = element_text(color = "#52616b"),
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(0.65, "lines"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#1f2933")
    )
}

point_layer <- function(size = 0.12, alpha = 0.55) {
  layer <- geom_point(size = size, alpha = alpha, stroke = 0)
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    return(ggrastr::rasterise(layer, dpi = 300))
  }
  layer
}

df <- cbind(
  data.frame(cell = rownames(obs), obs, check.names = FALSE),
  data.frame(UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], check.names = FALSE)
)

set.seed(1)
total_cells <- nrow(df)
sampled_for_plotting <- total_cells > max_cells
sampling_note <- if (sampled_for_plotting) {
  sprintf(
    "Sampled %s of %s cells for plotting; expression scale: %s",
    format(max_cells, big.mark = ",", scientific = FALSE),
    format(total_cells, big.mark = ",", scientific = FALSE),
    bundle$manifest[["expression_value_scale"]]
  )
} else {
  sprintf(
    "All %s cells shown; expression scale: %s",
    format(total_cells, big.mark = ",", scientific = FALSE),
    bundle$manifest[["expression_value_scale"]]
  )
}
if (nrow(df) > max_cells) {
  keep <- sample.int(nrow(df), max_cells)
  plot_df <- df[keep, , drop = FALSE]
  expr_df <- expr[keep, , drop = FALSE]
} else {
  plot_df <- df
  expr_df <- expr
}

umap_group <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[group_by]])) +
  point_layer() +
  coord_equal() +
  theme_nc2024(base_size = 18) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(
    title = sprintf("NC2024 UMAP by %s", group_by),
    subtitle = sampling_note,
    x = "UMAP 1",
    y = "UMAP 2",
    color = group_by
  )

ggsave(file.path(out_dir, "umap_by_group.png"), umap_group, width = 11.5, height = 8.8, dpi = 240)

umap_cluster <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[cluster_by]])) +
  point_layer() +
  coord_equal() +
  theme_nc2024(base_size = 18) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 15)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
  labs(
    title = sprintf("NC2024 UMAP by %s", cluster_by),
    subtitle = sampling_note,
    x = "UMAP 1",
    y = "UMAP 2",
    color = cluster_by
  )

ggsave(file.path(out_dir, "umap_by_cluster.png"), umap_cluster, width = 11.5, height = 8.8, dpi = 240)

frac_tab <- as.data.frame(prop.table(table(obs[[group_by]])))
colnames(frac_tab) <- c("group", "fraction")
cell_fraction <- ggplot(frac_tab, aes(x = reorder(group, fraction), y = fraction)) +
  geom_col(fill = "#2f6f73", width = 0.78) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_nc2024(base_size = 16) +
  labs(title = "NC2024 cell fraction", x = NULL, y = "Fraction of cells")

ggsave(file.path(out_dir, "cell_fraction.png"), cell_fraction, width = 9, height = 6.5, dpi = 220)

marker_names <- colnames(expr_df)
marker_names <- marker_names[vapply(expr_df[, marker_names, drop = FALSE], is.numeric, logical(1))]
if (!length(marker_names)) {
  stop("No numeric marker genes found in marker_expr.")
}
marker_plot_df <- cbind(
  plot_df[, c(cluster_by), drop = FALSE],
  expr_df[, marker_names, drop = FALSE]
)
colnames(marker_plot_df)[1] <- "cluster"

marker_long <- do.call(
  rbind,
  lapply(marker_names, function(g) {
    data.frame(cluster = marker_plot_df$cluster, gene = g, expr = marker_plot_df[[g]], check.names = FALSE)
  })
)

marker_dot <- aggregate(expr ~ cluster + gene, data = marker_long, FUN = mean)
marker_dot$cluster <- factor(marker_dot$cluster, levels = rev(sort(unique(as.character(marker_dot$cluster)))))
marker_dot_plot <- ggplot(marker_dot, aes(x = gene, y = cluster, color = expr, size = expr)) +
  geom_point(alpha = 0.9) +
  scale_color_viridis_c(option = "magma") +
  scale_size(range = c(2.5, 12), guide = "none") +
  theme_nc2024(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  labs(
    title = "NC2024 marker mean exported expression by cluster",
    subtitle = sampling_note,
    x = NULL,
    y = cluster_by,
    color = "exported expression"
  )

ggsave(file.path(out_dir, "marker_dot.png"), marker_dot_plot, width = 12.5, height = 8.5, dpi = 240)

chunk_size <- 8L
marker_chunks <- split(marker_names, ceiling(seq_along(marker_names) / chunk_size))
marker_panel_plots <- lapply(seq_along(marker_chunks), function(i) {
  genes <- marker_chunks[[i]]
  panel_df <- marker_dot[marker_dot$gene %in% genes, , drop = FALSE]
  panel_df$gene <- factor(panel_df$gene, levels = genes)
  p <- ggplot(panel_df, aes(x = gene, y = cluster, color = expr, size = expr)) +
    geom_point(alpha = 0.92) +
    scale_color_viridis_c(option = "magma") +
    scale_size(range = c(3, 13), guide = "none") +
    theme_nc2024(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1),
      legend.position = "right",
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 15)
    ) +
    labs(
      title = sprintf("NC2024 marker panel %d/%d", i, length(marker_chunks)),
      subtitle = sampling_note,
      x = NULL,
      y = cluster_by,
      color = "exported expression"
    )
  ggsave(file.path(out_dir, sprintf("marker_dot_panel_%02d.png", i)), p, width = 10.8, height = 8.2, dpi = 240)
  p
})

pdf(file.path(out_dir, paste0(basename(normalizePath(out_dir, mustWork = FALSE)), ".pdf")), width = 8, height = 6)
print(umap_group)
print(umap_cluster)
print(cell_fraction)
print(marker_dot_plot)
for (panel_plot in marker_panel_plots) {
  print(panel_plot)
}
dev.off()

message("Saved figures in '", out_dir, "'.")
