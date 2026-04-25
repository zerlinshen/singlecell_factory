# 2026-04-23 - NC2024 stage-1 recovery baseline

## Objective

- Stabilize the NC2024 full-cohort execution path.
- Prove that the full-cohort stage-1 lane can complete end to end.
- Prove that the controller fallback path `large -> massive` is truthful.

## Starting Context

- Canonical full-cohort input eventually validated at:
  - `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- Key scale facts:
  - `81` samples
  - `884050` retained cells
- Early confusion around `positive_barcodes_total = 96117359` was misleading and should not be treated as the final analysis object.

## What Was Run

- direct `massive` stage-1 recovery runs
- controller fallback validation run
- follow-up article-level single-cell analyses on top of the canonical successful run

## Outcome

- `success`

## Evidence

- direct successful stage-1:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222/module_status.csv`
- controller successful fallback:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552/module_status.csv`
- fallback proof log:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_WITH_FALLBACK.launch.log`

## Problems Encountered

- `large` repeatedly died by capacity after early modules
- `qc.py` failed on `Dataset2D.drop`
- `doublet_detection.py` grouped Scrublet failed on lazy `subX`
- `clustering.py` still had unsafe lazy matrix handoffs
- controller had truthfulness and `prepared_input.zarr` existence-check issues

## What Was Changed

- launcher/bootstrap fixes
- controller truthfulness and zarr existence-check fixes
- `qc.py` axis-table materialization
- `doublet_detection.py` matrix materialization before Scrublet
- `clustering.py` matrix materialization in CSS / CPU / GPU paths

## Resolution

- direct `massive` stage-1 completed successfully
- controller `large -> massive` completed successfully with `FINAL_STATUS=SUCCESS_MASSIVE`
- article-level follow-up analyses were able to proceed from the stable baseline

## Cautions for the Next Run

- Do not use `large` as the main debug lane for this cohort.
- Before launching any new remote run, read this note and `LATEST.md`.
- Prefer direct `massive` when the question is module correctness rather than orchestration truth.

## Improvement Ideas

- make controller smarter about reusing frozen prepared input
- keep updating the run memory after every major remote execution
- maintain canonical/superseded labeling on run directories to simplify cleanup

## Classification

- `canonical`
