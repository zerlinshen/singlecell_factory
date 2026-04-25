# 2026-04-24 - NC2024 fresh massive rerun

## Objective

- Execute a fresh direct `massive` stage-1 rerun for the NC2024 full cohort.
- Re-establish live end-to-end reproducibility evidence on top of the canonical prepared input.
- Keep this round scope-tight and do not yet mix in PDF packaging, missed-result discovery, or ELF3 analysis.

## Starting Context

- Canonical prepared input:
  - `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- Canonical prior direct stage-1 successful run:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222`
- Canonical prior controller fallback successful run:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552`
- Carry-forward rule before this run:
  - use direct `massive`
  - do not use `large` as the main debug lane
  - do not rerun prepare if canonical prepared input is still valid

## What Was Run

- Host:
  - `ubuntu-tail`
- Working directory:
  - `/home/zerlinshen/singlecell_factory`
- Execution mode:
  - `debug_massive`
- Command:
  - `bash scripts/run_emtab13526_full_cohort_stage1.sh massive NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO`
- Run directory:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`

## Outcome

- `success`

## Evidence

- `module_status.csv`:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/module_status.csv`
- `run_manifest.json`:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/run_manifest.json`
- Final object:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/final_adata.h5ad`
- Required downstream stage-1 artifacts confirmed:
  - `annotation/cell_type_annotation.csv`
  - `composition/composition_proportions.csv`
  - `immune_phenotyping/immune_phenotyping.csv`
  - `tumor_microenvironment/tme_scores_per_cell.csv`
- Checkpoints confirmed:
  - `after_cellranger.json`
  - `after_qc.json`
  - `after_doublet_detection.json`
  - `after_clustering.json`
  - `after_annotation.json`
  - `after_composition.json`

## Problems Encountered

- None that blocked completion.
- `doublet_detection` and `clustering` were the longest stages, but both completed successfully.
- The stage launcher log file itself stayed mostly empty because the child pipeline emitted useful outputs into the run directory rather than into the outer shell log.

## What Was Changed

- No repo code changes were required for this fresh rerun.
- This round was an execution-only reproducibility check.

## Resolution

- Fresh direct `massive` stage-1 rerun completed successfully on `2026-04-24`.
- All eight validated stage-1 modules were `ok` again:
  - `cellranger`
  - `qc`
  - `doublet_detection`
  - `clustering`
  - `annotation`
  - `composition`
  - `immune_phenotyping`
  - `tumor_microenvironment`
- This re-establishes live reproducibility evidence beyond the earlier 2026-04-23 successful baseline.

## Cautions for the Next Run

- Keep the next round focused on PDF packaging / missed-result discovery / ELF3 analysis rather than repeating stage-1 again without a new reason.
- Preserve the distinction between:
  - baseline reproducibility proof
  - expanded discovery analysis
- Continue using direct `massive` if another module-correctness question arises.

## Improvement Ideas

- Next round should harvest this fresh rerun into a dedicated reproduction package and PDF.
- Then expand into module-surface review for missed results and an ELF3-focused analysis thread.

## Classification

- `canonical`
