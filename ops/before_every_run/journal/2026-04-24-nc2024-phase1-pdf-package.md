# 2026-04-24 - NC2024 phase-1 PDF reproduction package

## Objective

- Build the first formal PDF reproduction package after the fresh direct
  `massive` full-cohort rerun.
- Keep phase 1 focused on execution evidence, final-object schema validation,
  current claim support, and figure-comparison status.
- Do not mix in missed-result discovery or ELF3-focused analysis yet.

## Starting Context

- Fresh direct `massive` rerun:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
- Schema inspection:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/schema_inspection/final_adata_schema_inspection.json`
- Cross-artifact validation:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/schema_inspection/final_adata_schema_validation.json`

## What Was Built

- Package directory:
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase1_fresh_massive_reproduction_20260424`
- PDF:
  - `NC2024_phase1_fresh_massive_reproduction_20260424.pdf`
- Markdown source:
  - `NC2024_phase1_fresh_massive_reproduction_20260424.md`
- Evidence manifest:
  - `evidence_manifest.json`
- Included lightweight assets:
  - Fig. 1, Fig. 2, and Fig. 10 comparison boards
  - fresh rerun schema inspection / validation JSON mirrors
  - checkpoint hierarchy report artifacts
  - macrophage correlation / programme / cluster summaries
  - source markdown documents used to build the package

## Outcome

- `success`

## Evidence

- PDF generated with bundled Python `reportlab 4.4.9`.
- PDF size:
  - `2389997` bytes after the evidence-contract normalization rebuild
- PDF page count:
  - `4`
- Manifest validation:
  - `schema_validation_summary.pass = true`
  - all schema validation checks are true
  - shape remains `810218 x 30374`
  - `X_encoding = csr_matrix`
- Package manifest includes:
  - `13` local assets
  - `10` source docs

## What Changed

- Added a local package builder:
  - `build_nc2024_phase1_package.py`
- Added package output under:
  - `reproduction_packages/NC2024_phase1_fresh_massive_reproduction_20260424/`
- Updated reproduction workspace README to point to the phase-1 package.
- Updated this run-memory summary so the next branch starts from:
  - missed-result discovery
  - ELF3-focused analysis

## Cautions for the Next Branch

- The phase-1 PDF is an evidence package, not a new biological discovery run.
- The current checkpoint hierarchy remains `partial`.
- The foetal-like macrophage programme remains `partial`.
- Spatial claims remain unsupported in this scRNA-only lane.

## Improvement Ideas

- Next branch should inventory module surfaces for missed results.
- After that, run the ELF3-focused analysis thread using the fresh rerun as the
  validated baseline.

## Classification

- `canonical`
