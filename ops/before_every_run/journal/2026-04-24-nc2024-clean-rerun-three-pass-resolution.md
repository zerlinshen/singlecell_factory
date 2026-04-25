# 2026-04-24 - NC2024 clean rerun three-pass resolution

## Objective

- Execute one clean full-cohort rerun with all eligible modules after the
  pseudobulk contract repair.
- Verify that `pseudobulk_de` can succeed pipeline-native, not only via a
  separate recovery directory.

## Starting Context

- Current canonical input remained:
  `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- Current canonical extended run remained protected during pre-run cleanup:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003`
- Frozen rerun contract:
  - `execution_mode=debug_massive`
  - `SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only`
  - eligible modules only
  - explicit `--pseudobulk-exploratory-group-vs-rest`
  - exclude `validate_cbioportal`, `rna_velocity`, `cnv_inference`, `evolution`

## What Was Run

- Pass 1 rerun:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_180203`
- Pass 2 rerun:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_190353`
- Pass 3 rerun:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652`

## Outcome

- `success after three-pass diagnosis/resolution`

## Evidence

- Final successful rerun:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652`
- Final successful status file:
  `module_status.csv`
- Final successful manifest:
  `run_manifest.json`
- Row-level contract audit:
  `module_reconciliation.tsv`
- Readable remote R rerender:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652/r_plots/phase5_readable_20260424`
- Local Phase-5 package:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase5_clean_full_cohort_rerun_20260424`

## Problems Encountered

- Pass 1 completed but `pseudobulk_de` still failed in-pipeline with:
  `_cs_matrix.sum() got an unexpected keyword argument 'keepdims'`
- Root cause analysis showed that on the lazy prepared-Zarr path,
  `cellranger` copied lazy `adata.X` into `adata.layers["counts"]`, so
  pseudobulk later received a lazy proxy rather than a concrete CSR matrix.
- Pass 2 changed `cellranger` to materialize the full counts layer early.
  This removed the lazy-proxy bug but pushed too much memory into the whole
  pipeline and the process was killed in a later stage.

## What Was Changed

- `cellranger.py` was first patched to materialize counts eagerly, then refined.
- Final winning fix:
  - keep the counts layer logically attached in `cellranger`
  - delay materialization until `pseudobulk_de._aggregate()`
  - when `pseudobulk_de` actually runs, convert lazy/eager backends into a
    stable CSR matrix on demand
- Added targeted regression coverage for lazy-count materialization.

## Resolution

- Pass 3 used the deferred-materialization fix and completed successfully.
- `module_status.csv` now shows:
  - `pseudobulk_de = ok`
  - all planned modules present and `ok`
  - all excluded modules absent by contract
- `run_manifest.json` now records:
  - `pseudobulk_de_status = completed`
  - `pseudobulk_de_mode = exploratory`
  - `pseudobulk_de_significant_genes = 153348`
- `module_reconciliation.tsv` explicitly proves:
  - planned modules all present
  - excluded modules all absent with reasons

## Cautions for the Next Run

- Do not reintroduce eager full-cohort counts-layer materialization at
  `cellranger` load time for massive runs.
- The winning pattern is:
  - keep lazy inputs through the general pipeline
  - materialize counts only inside pseudobulk aggregation
- Continue to treat `validate_cbioportal` as a separate deterministic sidecar,
  not part of the clean rerun contract.

## Classification

- `canonical`
  - successful third-pass rerun:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652`
- `evidence-only`
  - first-pass rerun:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_180203`
  - second-pass rerun:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_190353`
- `superseded`
  - none yet deleted; retirement decision should wait until all local/remote
    source-of-truth surfaces are updated
