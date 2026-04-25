# 2026-04-24 - NC2024 Phase-4 report and pseudobulk contract hardening

## Objective

- Regenerate the Phase-4 PDF into a readable report without launching a new
  full-cohort run.
- Turn `pseudobulk_de` into a contrast-aware eligible module rather than an
  always-on exploratory lane.
- Make `rna_velocity` skip honestly when the dataset lacks true splicing
  modality.
- Update remote/local guidance so Python and R plotting are documented as
  complementary, not mutually exclusive.
- Add a reproduce-stage run-retention skill so bulky old reproduce/test outputs
  can be cleaned after evidence capture.

## Starting Context

- Current source of truth run:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003`
- Current local package:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase4_extended_full_cohort_real_run_20260424`
- Carry-forward constraints:
  - no blind full-cohort rerun
  - keep original `module_status.csv` failure rows honest
  - preserve sparse/float32 behavior unless a method explicitly requires more

## What Was Run

- Local Phase-4 report builder rewrite only; no new full-cohort pipeline run.
- Remote code changes plus targeted tests on `ubuntu-tail`.
- Remote verification:
  - `python -m py_compile` for modified pipeline files
  - `pytest -q -o addopts='' tests/test_reliability_fixes.py`
  - `pytest -q -o addopts='' tests/test_rna_velocity_module.py -k "not_mutating or integration"`

## Outcome

- `success`
- Phase-4 PDF is now regenerated as a readable deck-style report.
- A first-pass unreadable backup PDF is preserved inside the same package.
- `pseudobulk_de` now has explicit confirmatory vs exploratory semantics in the
  code/config surface.
- `rna_velocity` now records a clear `skipped` status when true splicing inputs
  are absent.
- A new `reproduce-run-retention` skill now exists locally.

## Evidence

- New readable PDF:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase4_extended_full_cohort_real_run_20260424/NC2024_phase4_extended_full_cohort_real_run_20260424.pdf`
- First-pass backup:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase4_extended_full_cohort_real_run_20260424/NC2024_phase4_extended_full_cohort_real_run_20260424.first_pass_unreadable.pdf`
- Readable asset cards:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase4_extended_full_cohort_real_run_20260424/report_assets/figures_readable/`
- Visual sanity previews:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase4_extended_full_cohort_real_run_20260424/report_assets/page_previews_readable/`
- Remote modified files:
  - `/home/zerlinshen/singlecell_factory/workflow/modular/config.py`
  - `/home/zerlinshen/singlecell_factory/workflow/modular/cli.py`
  - `/home/zerlinshen/singlecell_factory/workflow/modular/modules/pseudobulk_de.py`
  - `/home/zerlinshen/singlecell_factory/workflow/modular/modules/rna_velocity.py`
  - `/home/zerlinshen/singlecell_factory/tests/test_reliability_fixes.py`
  - `/home/zerlinshen/singlecell_factory/README.md`
- New local skill:
  - `/Users/zerlinshen/.codex/skills/reproduce-run-retention/SKILL.md`

## Problems Encountered

- The original Phase-4 builder compressed tables and stacked dense figures,
  making the PDF hard to read.
- `pseudobulk_de` defaulted to exploratory group-vs-rest behavior even when no
  confirmatory contrast contract existed.
- `rna_velocity` could be requested on a dataset without spliced/unspliced
  modality, which turned an expected ineligibility into a misleading failure.
- One pseudobulk compatibility pass briefly broke a legacy test helper
  signature; this was fixed before final verification.

## What Was Changed

- Replaced the Phase-4 builder with a deck-style readable report generator.
- Added a `PseudobulkConfig` to the remote modular config/CLI surface.
- Added explicit pseudobulk flags:
  - `--pseudobulk-sample-col`
  - `--pseudobulk-group-col`
  - `--pseudobulk-contrast-col`
  - `--pseudobulk-contrast-a`
  - `--pseudobulk-contrast-b`
  - `--pseudobulk-contrast-json`
  - `--pseudobulk-exploratory-group-vs-rest`
- Changed pseudobulk semantics:
  - explicit contrast contract -> confirmatory mode
  - exploratory output only when explicitly enabled
  - missing counts or missing contract -> `skipped`, not silent pseudo-result
- Added RNA velocity eligibility skip behavior for missing splicing modality.
- Updated remote/local README and profile/playbook guidance to state:
  - all modules = all eligible modules
  - Python plots still matter
  - R is preferred for prettier large-cohort figures
- Added `reproduce-run-retention` skill for reproduce-stage storage governance.

## Resolution

- Remote regression tests passed after the compatibility fix:
  - `tests/test_reliability_fixes.py`: `9 passed`
  - `tests/test_rna_velocity_module.py -k "not_mutating or integration"`:
    `2 passed`
- Remote `py_compile` passed for all modified code files.
- Local PDF now has readable one-figure-per-page sections and readable status
  tables; the old PDF is preserved only as a backup.
- A second-pass remote R rerender now exists and is preferred for report use:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/r_plots/extended_real_run_readable_20260424`
  - larger `umap_by_group.png`
  - larger `umap_by_cluster.png`
  - split `marker_dot_panel_01.png` to `marker_dot_panel_03.png`
  - local mirror:
    `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase4_extended_full_cohort_real_run_20260424/figures/r_plots_readable`

## Cautions for the Next Run

- Do not treat `pseudobulk_de` exploratory output as a confirmatory paper-grade
  result unless an explicit contrast contract was supplied.
- Do not request `rna_velocity` on NC2024 full cohort unless true splicing
  inputs are available.
- Do not describe the new plotting contract as “R replaces Python”. The correct
  framing is Python module-native + R publication-style when better.
- The new reproduce-stage retention skill governs old reproduce/test outputs,
  not raw data or current source-of-truth final artifacts.

## Improvement Ideas

- Add a first-class recovery manifest surface so pipeline failure and recovery
  can be linked more formally than a separate recovery folder.
- Add provenance hashes for locally mirrored remote figures/tables.
- Add direct wrapper tests for pseudobulk CLI flags, not only module tests.

## Classification

- `canonical`
  - current extended full-cohort real run
  - regenerated readable Phase-4 package
- `evidence-only`
  - unreadable first-pass backup PDF
  - preview renders under `report_assets/page_previews_readable/`
- `superseded`
  - the old in-place Phase-4 unreadable layout as the preferred reading surface
- `failed exploratory`
  - none newly created in this round
