# 2026-04-24 NC2024 R Environment Hardening

## Task

Install and verify remote R packages recommended for NC2024 reporting,
visualization, and small/medium R-side follow-up while preserving the
remote-first plotting contract.

## Objective

- Equip `/home/zerlinshen/conda/envs/r_multiomics` with useful reporting and
  bridge packages.
- Validate the installed packages on the existing NC2024 `r_bundle`.
- Avoid full-cohort reruns and avoid direct full-object Seurat conversion.
- Update remote and local README surfaces in the same workstream.

## What Was Attempted

- Read current local and remote `before-every-run/LATEST.md`.
- Inspected installed packages in `r_multiomics`.
- Installed missing conda packages:
  `r-readr`, `r-ggrastr`, `r-pheatmap`,
  `bioconductor-complexheatmap`, `r-circlize`, `r-hdf5r`,
  `r-biocmanager`, `r-harmony`, `r-r.utils`, and `r-remotes`.
- Used the existing fresh full-cohort `r_bundle` for smoke validation.
- Checked whether `r-seuratdisk` could be installed reproducibly from conda.

## What Succeeded

- Seurat was already present as `Seurat 5.4.0` and `SeuratObject 5.4.0`.
- Installed/verified:
  `readr 2.2.0`, `ggrastr 1.0.2`, `scattermore 1.2`,
  `pheatmap 1.0.13`, `ComplexHeatmap 2.26.1`, `circlize 0.4.18`,
  `hdf5r 1.3.12`, `harmony 1.2.4`, `BiocManager 1.30.27`,
  `R.utils 2.13.0`, `zellkonverter 1.20.1`, and `remotes 2.5.0`.
- Validation output:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/r_env_dependency_smoke_20260424`
- Validation produced:
  `NC2024_R_ENV_VALIDATION_SUMMARY.txt`,
  `nc2024_r_env_smoke_raster_umap.png`,
  `nc2024_r_env_smoke_raster_umap.pdf`,
  `nc2024_r_env_smoke_pheatmap.png`, and
  `nc2024_r_env_smoke_complexheatmap.pdf`.
- Validation read all `810218` rows from the compact bundle, plotted `100000`
  sampled cells with `ggrastr`, generated heatmaps with `pheatmap` and
  `ComplexHeatmap`, and opened `final_adata.h5ad` via `hdf5r`.

## What Failed Or Was Rejected

- Initial validation showed that `data.table::fread()` requires `R.utils` for
  direct `.gz` reads. This was fixed by installing `r-r.utils`.
- `r-seuratdisk` was rejected for the current environment because the
  conda-forge build requires old `spatstat` versions incompatible with R 4.5 and
  the current `zellkonverter` stack.
- GitHub/source installation of `SeuratDisk` was not used because it would be
  less reproducible and could alter dependencies outside conda's solved state.

## Artifact Classification

- Canonical environment:
  `/home/zerlinshen/conda/envs/r_multiomics`
- Evidence-only validation output:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/r_env_dependency_smoke_20260424`
- Canonical docs updated:
  `/home/zerlinshen/singlecell_factory/README.md`
  `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/README.md`
  local reproduction `README.md`
  local `agent_runs/2026-04-24-r-environment-hardening/RUN.md`

## Remaining Risks

- `SeuratDisk` is unavailable in this environment through a clean conda solve.
  Prefer `zellkonverter` or direct `hdf5r` checks for small H5AD work.
- For NC2024-scale objects, full Seurat conversion remains discouraged because
  it can materialize large in-memory objects. Use compact bundles for plotting.

## Next Operator Notes

- Use the compact bundle path for large-cohort R plots.
- Use `ggrastr`/`scattermore` for dense UMAP-style plotting.
- Use `pheatmap`/`ComplexHeatmap` for summary heatmaps from compact tables.
- Use `harmony` only for explicit small/medium R-side integration checks, not
  as an implicit replacement for the canonical Python pipeline outputs.
