# 2026-04-24 NC2024 Bridge Code Review Hardening

## Task

Review and improve the remote cooperation between `singlecell_factory` and the
remote R plotting pipeline.

## Objective

- Let `singlecell_factory` remain the high-throughput sparse/full-object engine.
- Let R focus on publication-style figures from compact bundle tables.
- Improve stale-output safety, ROI configurability, and large-cohort plotting
  quality.
- Update remote and local README/run-memory surfaces.

## What Changed

- `scripts/export_singlecell_r_bundle.py`
  - Adds `input_h5ad_bytes` and `input_h5ad_mtime_epoch` to
    `bundle_manifest.tsv`.
- `R/remote_bundle_manifest.R`
  - Reads bundle CSV/GZ files through `data.table::fread()` when available while
    preserving row names.
- `scripts/run_remote_bundle_plot.sh`
  - Adds `R_PLOT_THREADS`, `R_BUNDLE_MARKERS`, `R_BUNDLE_OBS_COLS`, and
    `R_BUNDLE_OBSM`.
  - Treats `R_BUNDLE_OBS_COLS` and `R_BUNDLE_OBSM` as additive controls, not
    replacements, so selected `group_by`, selected `cluster_by`, `X_umap`, and
    `X_pca` stay available for plotting and validation.
  - Reuse now requires source path, source size, source mtime epoch, optional
    request matching, and R-side SHA256 validation.
  - R validation receives the bundle path via environment variable instead of
    interpolating the path into an R expression.
- `scripts/plot_remote_bundle_large.R`
  - Uses `ggrastr` rasterization when available.
  - Applies a cleaner NC2024 plotting theme.
  - Marker-dot plots now use every numeric marker column in `marker_expr.csv.gz`
    instead of a hard-coded marker list.

## Validation

- Python compile passed for the exporter.
- `bash -n` passed for the wrapper.
- R parse passed for the manifest validator and large plotter.
- Targeted exporter pytest passed with explicit conda prefix:
  `CONDA_PREFIX=/home/zerlinshen/conda/envs/sc_gpu PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 /home/zerlinshen/conda/envs/sc_gpu/bin/python -m pytest -o addopts="" tests/test_singlecell_r_bundle_export.py -q`
  returned `4 passed`.
- Running pytest without `CONDA_PREFIX` failed during collection through
  `anndata -> zarr -> cupy`; this was an environment activation issue, not a
  code assertion failure.
- Export + prettier plot smoke:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_pretty_export_20260424`
- Reuse smoke:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_pretty_reuse_20260424`
- The reuse smoke printed `Reusing validated R bundle`.
- Ralph follow-up refreshed the canonical default `r_bundle` itself with the
  stricter manifest keys, then verified default-path reuse:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_refresh_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_reuse_20260424`
- The canonical reuse smoke printed `Reusing validated R bundle`.
- Post-deslop canonical reuse smoke also printed `Reusing validated R bundle`:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_reuse_post_deslop_20260424`
- Post-contract-fix canonical reuse smoke again printed
  `Reusing validated R bundle`:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_reuse_after_contract_fix_20260424`
- Custom-contract smoke requested `PDCD1,LAG3,CTLA4`, extra obs columns
  `sample,patient`, and `X_umap`; the wrapper exported the effective contract
  `cell_type,leiden,sample,patient` and `X_umap,X_pca`, generated plots, then
  reused that same bundle:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_custom_contract_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_custom_contract_reuse_20260424`
- No prepare or full-cohort stage-1 rerun was launched.

## Artifact Classification

- Canonical:
  - `/home/zerlinshen/singlecell_factory/scripts/export_singlecell_r_bundle.py`
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/R/remote_bundle_manifest.R`
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh`
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/plot_remote_bundle_large.R`
- Evidence-only:
  - `/tmp/nc2024_bridge_review_bundle_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_pretty_export_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_pretty_reuse_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_refresh_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_reuse_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_reuse_post_deslop_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_reuse_after_contract_fix_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_custom_contract_20260424`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_custom_contract_reuse_20260424`

## Remaining Risks

- The temporary `/tmp/nc2024_bridge_review_bundle_20260424` is evidence-only and
  can be deleted later if space pressure appears.
- Full Seurat conversion is still intentionally discouraged for NC2024-scale
  objects; this lane improves compact-bundle figures, not R-side full-object
  analysis.

## Next Operator Notes

- For ROI plot requests, set `R_BUNDLE_MARKERS` explicitly and use the wrapper.
- For dense UMAP-like figures, keep `R_PLOT_THREADS` reasonable to avoid
  oversubscribing the remote server.
- Treat bridge smoke outputs as operational evidence, not new biological claim
  evidence.
