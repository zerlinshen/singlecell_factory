# 2026-04-24 - NC2024 fresh rerun H5AD schema inspection

## Objective

- Close the residual risk from the fresh direct `massive` rerun by opening
  `final_adata.h5ad`, not merely checking that the file exists.
- Equip the inspection branch with AnnData/HDF5 tooling if needed.
- Preserve concrete schema evidence before PDF packaging and downstream
  discovery analysis.

## Starting Context

- Fresh direct `massive` rerun:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
- Prior evidence:
  - `module_status.csv` reported all eight stage-1 modules as `ok`
  - `run_manifest.json` existed
  - `final_adata.h5ad` existed
- Residual risk before this inspection:
  - the final object had not yet been opened for AnnData/HDF5 schema validation

## What Was Run

- Host:
  - `ubuntu-tail`
- Working directory:
  - `/home/zerlinshen/singlecell_factory`
- Execution mode:
  - `debug_massive-inspection`
- Toolchain check:
  - `conda run -n sc_gpu python -c "import h5py; import anndata"`
- Inspection approach:
  - `anndata.read_h5ad(..., backed="r")`
  - `h5py.File(..., "r")`
- Schema artifact:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/schema_inspection/final_adata_schema_inspection.json`
- Cross-artifact validation artifact:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/schema_inspection/final_adata_schema_validation.json`

## Outcome

- `success`
- No dependency installation was needed.
- Existing remote versions were sufficient:
  - `anndata 0.12.10`
  - `h5py 3.16.0`

## Evidence

- `final_adata.h5ad` opened successfully in backed mode.
- AnnData shape:
  - `810218 x 30374`
- Cross-artifact validation agreement with `run_manifest.json`:
  - `cells_after_doublet_removal = 810218`
  - `genes_after_qc = 30374`
- Cross-artifact validation agreement with `module_status.csv`:
  - all eight stage-1 modules are `ok`
- HDF5 `X`:
  - encoding: `csr_matrix`
  - shape: `810218 x 30374`
  - data dtype: `float32`
- Embeddings:
  - `X_css = 810218 x 223`
  - `X_pca = 810218 x 20`
  - `X_umap = 810218 x 2`
- `obs`:
  - `61` columns
  - includes sample/patient metadata, QC fields, doublet fields, Leiden clusters,
    `cell_type`, immune subtype/signature fields, and TME score fields
- Metadata breadth:
  - `81` samples
  - `26` patients/donors
  - `11` cell types
  - `2` conditions
  - `4` disease labels

## Problems Encountered

- Long SSH here-doc commands with `conda run` produced no useful stdout in this
  runtime, despite returning success.
- Resolution was to stage a temporary local inspector script, `scp` it to
  `/tmp/inspect_h5ad_schema.py`, and execute that script directly under
  `conda run -n sc_gpu`.

## What Was Changed

- No pipeline code was changed.
- No remote dependency installation was required.
- Added schema evidence under the fresh rerun directory:
  - `schema_inspection/final_adata_schema_inspection.json`
- Added explicit manifest/module-status validation evidence under the fresh rerun directory:
  - `schema_inspection/final_adata_schema_validation.json`
- Synced a local mirror for report/package use:
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/figure_compare/assets/fresh_rerun_final_adata_schema_inspection.json`
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/figure_compare/assets/fresh_rerun_final_adata_schema_validation.json`

## Resolution

- The previous residual risk, "final_adata.h5ad exists but has not been opened",
  is closed.
- The fresh rerun now has object-level AnnData/HDF5 validation evidence.

## Cautions for the Next Run

- Continue using backed-mode AnnData or direct HDF5 inspection for this 12G file.
- Avoid long SSH here-doc inspection commands with `conda run`; prefer a staged
  script or one-line `python -c` checks.
- Do not trigger another full-cohort stage-1 run unless a new module-correctness
  question appears.

## Improvement Ideas

- The next productive branch remains the previously identified sequence:
  - PDF reproduction package from the fresh rerun
  - module-surface review for missed results
  - ELF3-focused analysis

## Classification

- `canonical`
