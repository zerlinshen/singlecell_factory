# US-W4-13 — R testthat seed for the bridges/ R scripts.
#
# Minimal regression coverage for scripts/plot_remote_bundle.R — the entry
# point used by the remote-to-local report bridge. Full integration testing
# (real Seurat + bundle fixture) is Wave 5 follow-up.
#
# Run via: Rscript -e 'testthat::test_dir("bridges/local_r_pipeline_macbook/tests/testthat")'
# requires r_multiomics conda env: /home/zerlinshen/conda/envs/r_multiomics/bin/Rscript

library(testthat)

# Resolve paths from THIS file's location (tests/testthat/) up to bridges/.
# testthat sets the cwd to the test file's dir, so paths must be relative
# to bridges/local_r_pipeline_macbook/tests/testthat/.
script_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
plot_script <- file.path(script_root, "scripts", "plot_remote_bundle.R")

test_that("plot_remote_bundle.R exists and is non-empty", {
  expect_true(file.exists(plot_script),
              info = "plot_remote_bundle.R must exist at bridges/local_r_pipeline_macbook/scripts/")
  info <- file.info(plot_script)
  expect_gt(info$size, 0)
})

test_that("plot_remote_bundle.R declares CLI usage message", {
  contents <- paste(readLines(plot_script), collapse = "\n")
  expect_match(contents, "Rscript scripts/plot_remote_bundle.R")
  expect_match(contents, "bundle_dir.*out_dir",
               info = "usage line should reference the two positional args")
})

test_that("plot_remote_bundle.R sources the expected R modules", {
  contents <- paste(readLines(plot_script), collapse = "\n")
  expected_sources <- c(
    "R/theme_config.R",
    "R/qc_plots.R",
    "R/dim_plots.R",
    "R_bundle/remote_bundle_manifest.R"
  )
  for (src in expected_sources) {
    expect_match(contents, fixed = TRUE, regexp = src,
                 info = sprintf("script must source %s", src))
  }
})

test_that("R/ and R_bundle/ are symlinks per CLAUDE.md bridge contract", {
  r_dir <- file.path(script_root, "R")
  r_bundle <- file.path(script_root, "R_bundle")
  # These paths exist as symlinks pointing into ../../multiomics_r_factory/.
  # Use Sys.readlink to confirm both are symlinks (file.exists fails when
  # cwd doesn't help and a relative symlink target can't resolve).
  expect_true(nzchar(Sys.readlink(r_dir)),
              info = "bridges R must be a symlink (CLAUDE.md bridge rule)")
  expect_true(nzchar(Sys.readlink(r_bundle)),
              info = "bridges R_bundle must be a symlink (CLAUDE.md bridge rule)")
})
