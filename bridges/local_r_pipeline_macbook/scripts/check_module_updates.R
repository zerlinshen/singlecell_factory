#!/usr/bin/env Rscript

read_lines_safe <- function(path) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  readLines(path, warn = FALSE)
}

extract_module_names <- function(lines) {
  m <- regmatches(lines, gregexpr("`R/[A-Za-z0-9_.-]+\\.R`", lines))
  modules <- unique(unlist(m, use.names = FALSE))
  modules <- modules[nzchar(modules)]
  gsub("`", "", modules)
}

extract_required_pkgs <- function(main_lines) {
  pkg_line <- grep("required_pkgs <- c\\(", main_lines, fixed = TRUE)
  if (length(pkg_line) == 0L) return(character(0))
  start <- pkg_line[[1L]]
  txt <- main_lines[start]
  i <- start
  while (i < length(main_lines) && !grepl("\\)", main_lines[i], fixed = TRUE)) {
    i <- i + 1L
    txt <- c(txt, main_lines[i])
  }
  txt <- paste(txt, collapse = " ")
  txt <- gsub("required_pkgs <- c\\(|\\)", "", txt, fixed = FALSE)
  pkgs <- trimws(unlist(strsplit(txt, ",")))
  pkgs <- pkgs[nzchar(pkgs)]
  pkgs <- gsub('["]', "", pkgs)
  trimws(pkgs)
}

extract_imported_pkgs <- function(path) {
  lines <- read_lines_safe(path)
  hits <- gregexpr("(?<![A-Za-z0-9_])(library|requireNamespace|require)\\(([^),\"]+|\"[^\"]+\"|'[^']+')", lines, perl = TRUE)
  out <- character(0)
  for (i in seq_along(lines)) {
    if (hits[[i]][[1L]] == -1L) next
    m <- regmatches(lines[i], hits[[i]])
    for (x in m) {
      lib <- sub(".*\\((.+)", "\\1", x, perl = TRUE)
      lib <- sub("(^[[:space:]]*['\"]?)|(['\"].*$)", "", lib)
      lib <- trimws(lib)
      lib <- gsub("[[:space:]].*$", "", lib, perl = TRUE)
      if (nzchar(lib) && lib != "TRUE" && lib != "FALSE") {
        out <- c(out, lib)
      }
    }
  }
  unique(out)
}

normalize_paths <- function(paths) {
  normalizePath(paths, mustWork = TRUE)
}

if (Sys.getenv("R_LIBS_USER") == "") {
  Sys.setenv(R_LIBS_USER = tempdir())
}

git_root <- function() {
  out <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = TRUE)
  if (length(out) == 0L || any(grepl("fatal:", out, fixed = TRUE))) {
    stop("Not inside a git repository. Please run this script from project directory.")
  }
  out[[1L]]
}

root <- git_root()
setwd(root)

readme <- read_lines_safe("README.md")
readme_modules <- extract_module_names(readme)
main_lines <- read_lines_safe("main.R")
required_pkgs <- extract_required_pkgs(main_lines)

changed <- system2(
  "git",
  c("diff", "--name-only", "--diff-filter=AM", "HEAD", "--", "R"),
  stdout = TRUE,
  stderr = TRUE
)

if (length(changed) == 0L) {
  cat("No added/modified R files detected by git diff --name-only --diff-filter=AM HEAD -- R\n")
  quit(status = 0)
}

changed <- changed[endsWith(changed, ".R")]
if (length(changed) == 0L) {
  cat("No changed R scripts detected in this module scan.\n")
  quit(status = 0)
}

changed <- normalizePath(changed)
changed_names <- basename(changed)
changed_mods <- file.path("R", changed_names)
changed_base <- paste0("`", changed_mods, "`")

cat("Checking module updates for:\n")
cat(paste0("  - ", changed_mods, collapse = "\n"), "\n\n")

missing_readme <- character(0)
missing_provenance <- character(0)
missing_modules <- intersect(changed_mods, c("R/pipeline_steps.R", "R/pipeline_steps.r")) # placeholder

for (m in changed_mods) {
  in_list <- any(grepl(gsub("[.]", "\\.", sprintf("`%s`", m)), readme, fixed = FALSE))
  if (!in_list) missing_readme <- c(missing_readme, m)
}

provenance_rows <- grep("^\\| `R/", readme, fixed = FALSE)
for (m in changed_mods) {
  has_row <- any(grepl(sprintf("^\\| `%s`", m), readme[provenance_rows], fixed = TRUE))
  if (!has_row) missing_provenance <- c(missing_provenance, m)
}

missing_pkg_rows <- character(0)
missing_pkg <- c()
for (path in changed) {
  imports <- setdiff(extract_imported_pkgs(path), c("Seurat", "ggplot2", "dplyr", "tidyr", "RColorBrewer", "viridis", "scales", "tools", "patchwork"))
  if (length(imports) > 0L) {
    not_listed <- setdiff(imports, required_pkgs)
    if (length(not_listed) > 0L) missing_pkg <- c(missing_pkg, paste0(basename(path), ": ", paste(not_listed, collapse = ", ")))
  }
}

issues <- FALSE
if (length(missing_readme) > 0L) {
  issues <- TRUE
  cat("Missing module list entries:\n")
  cat(paste0("  - ", missing_readme, collapse = "\n"), "\n\n")
}
if (length(missing_provenance) > 0L) {
  issues <- TRUE
  cat("Missing provenance table rows in README module table:\n")
  cat(paste0("  - ", missing_provenance, collapse = "\n"), "\n\n")
}
if (length(missing_pkg) > 0L) {
  issues <- TRUE
  cat("Packages used by changed modules not listed in required_pkgs:\n")
  cat(paste0("  - ", unique(missing_pkg), collapse = "\n"), "\n\n")
}

if (issues) {
  cat("Module update check failed. Please update README/module references according to MODULE_UPDATE_SKILL.md.\n")
  quit(status = 1)
}

cat("Module update check passed.\n")
quit(status = 0)
