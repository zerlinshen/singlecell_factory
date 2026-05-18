# 2026-04-30 - NC2024 R Arrow Environment Validation

## Objective

- Provide an `arrow`-enabled R runtime for current NC2024 v2 parquet bundles.
- Validate the runtime through the real `multiomics_r_factory` bundle reader and
  plotting wrapper without rerunning the cohort or opening large H5AD objects.

## Starting Context

- `r_multiomics` used R `4.5.3` and did not have R package `arrow`.
- NC2024 v2 manifests already validated, but R-side v2 parquet plotting was
  blocked by missing `arrow`.
- Canonical v2 bundle targets:
  - `/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/r_bundle`
  - `/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2/r_bundle`

## What Was Run

- Host: `ubuntu-tail`
- Execution mode: `controller_validation`
- No cohort rerun happened.
- Environment provenance:
  `/home/zerlinshen/singlecell_factory/ops/env_records/r_arrow_20260430T095459+0800`

Commands included:

```bash
/home/zerlinshen/conda/bin/conda list -n r_multiomics --explicit
/home/zerlinshen/conda/bin/conda create -y -n r_multiomics_arrow --clone r_multiomics
/home/zerlinshen/conda/bin/conda install -n r_multiomics_arrow -c conda-forge --strict-channel-priority r-arrow --dry-run
/home/zerlinshen/conda/bin/conda install -y -n r_multiomics_arrow -c conda-forge --strict-channel-priority r-arrow
/home/zerlinshen/conda/bin/conda run -n r_multiomics_arrow Rscript -e 'library(arrow); packageVersion("arrow")'
/home/zerlinshen/conda/bin/conda run -n r_multiomics_arrow Rscript scripts/plot_remote_bundle_large.R <bundle> <out> cell_type leiden <max_cells>
```

## Outcome

- `success`: `r_multiomics_arrow` was created and validated.
- `success`: conda solver kept `r-base 4.5.3` and installed conda-forge
  `r-arrow 24.0.0` / `libarrow 24.0.0`.
- `success`: B/H v2 bundle read and plotted.
- `success`: tumor v2 bundle read and plotted with capped plotting sample.
- `success`: `run_remote_bundle_plot.sh` now defaults to
  `/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript`.
- `success`: `run_remote_bundle_plot.sh` now reuses v2
  `bundle_manifest.json` parquet bundles and delegates plotting to the upstream
  v1/v2-aware `multiomics_r_factory/scripts/plot_remote_bundle_large.R`.

## Evidence

- B/H validation:
  - `read_bundle_ok`
  - `obs_n=6426`
  - `obsm_keys=X_umap,X_pca`
  - `expr_dim=6426x9`
  - output:
    `/home/zerlinshen/multiomics_r_factory/output/nc2024_bh_arrow_smoke_20260430T095628+0800`
- Tumor validation:
  - `read_bundle_ok`
  - `obs_n=803784`
  - `obsm_keys=X_umap,X_pca`
  - `expr_dim=803784x9`
  - output:
    `/home/zerlinshen/multiomics_r_factory/output/nc2024_tumor_arrow_smoke_20260430T095652+0800`
- Wrapper validation:
  - B/H wrapper output:
    `/home/zerlinshen/multiomics_r_factory/output/nc2024_bh_arrow_wrapper_smoke_20260430T100416+0800`
  - tumor wrapper output:
    `/home/zerlinshen/multiomics_r_factory/output/nc2024_tumor_arrow_wrapper_smoke_20260430T100432+0800`
  - both wrapper logs reported `Reusing validated R bundle`.

## Problems Encountered

- None during install or validation.
- The original `r_multiomics` environment remains without `arrow`; it is now a
  rollback runtime rather than the default v2 parquet plotting runtime.

## What Was Changed

- Added environment:
  `/home/zerlinshen/conda/envs/r_multiomics_arrow`
- Updated default remote R wrapper runtime:
  `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh`
- Hardened the same wrapper so v2 JSON manifests are reusable and old bridge
  v1-only plotting scripts are bypassed by default.
- Updated relevant README/profile/run-memory documentation in
  `singlecell_factory` and `multiomics_r_factory`.

## Cautions for the Next Run

- Use `r_multiomics_arrow` for v2 parquet plotting.
- Use `r_multiomics` only as rollback unless it is intentionally upgraded.
- Compact R bundle values remain plotting/reporting handoff data, not a basis
  for new quantitative DE claims without full-object validation.

## Classification

- canonical:
  - `/home/zerlinshen/conda/envs/r_multiomics_arrow`
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh`
- evidence-only:
  - `/home/zerlinshen/singlecell_factory/ops/env_records/r_arrow_20260430T095459+0800`
  - `/home/zerlinshen/multiomics_r_factory/output/nc2024_bh_arrow_smoke_20260430T095628+0800`
  - `/home/zerlinshen/multiomics_r_factory/output/nc2024_tumor_arrow_smoke_20260430T095652+0800`
  - `/home/zerlinshen/multiomics_r_factory/output/nc2024_bh_arrow_wrapper_smoke_20260430T100416+0800`
  - `/home/zerlinshen/multiomics_r_factory/output/nc2024_tumor_arrow_wrapper_smoke_20260430T100432+0800`
- superseded:
  - treating `r_multiomics` as the default runtime for current v2 parquet plots
- failed exploratory:
  - none
