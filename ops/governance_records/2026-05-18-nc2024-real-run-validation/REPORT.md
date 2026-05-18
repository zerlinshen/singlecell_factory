# NC2024 Real Run Validation After Structure Change

Control-plane validation report only; not a scientific output.

- Generated UTC: `2026-05-18T09:03:03Z`
- Verdict: `REAL_RUN_PASS_PROJECT_ARCHITECTURE_CONFIRMED`
- Run ID: `2026-05-18T0900Z-13c2c88`
- Project run dir: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`
- Factory repo: `/home/zerlinshen/singlecell_factory`

## Cleanup Before Run

- Deleted legacy factory NC2024 entries: `48`
- Reclaimed target bytes: `30667532081`
- `results/` before cleanup: `29G	results`
- `results/` after cleanup: `27M	results`

## True Run Outcome

- First launch exposed a real architecture mismatch: old launcher generated an invalid non-hex run id.
- Recovery run used a valid `<timestamp>-<7hex>` run id and completed successfully.
- `launch_wave3_nc2024.sh` was patched to generate future run ids from the current factory git short SHA.
- Final AnnData shape: `5281 x 19504`
- Checked obs columns: `['leiden', 'cell_type', 'context_aware_celltype', 'context_aware_substate', 'sample']`

## Module Status

- `cellranger`: `ok` — completed
- `marker_db_loader`: `ok` — completed
- `qc`: `ok` — completed
- `doublet_detection`: `ok` — completed
- `clustering`: `ok` — completed
- `annotation`: `ok` — completed
- `batch_correction`: `skipped` — Only one batch found in 'sample'.
- `context_aware_annotation`: `ok` — completed
- `differential_expression`: `ok` — completed

Interpretation: `batch_correction` is an expected skip for this subsample because only one batch was present; all other modules completed.

## Governance Validation

- Single new run: `pass`, severity `{}`
- All NC project runs: `pass`, severity `{'info': 2, 'warning': 1}`
- New run has root `manifest.json`, Python `run_manifest.json`, Python `module_status.csv`, and project-owned figures.

## Factory Purity

- Factory NC2024 leftovers after cleanup: `0`
- Factory `results/` size after true run: `27M	/home/zerlinshen/singlecell_factory/results`
- No new NC2024 scientific output was written under `singlecell_factory/results/`.

## Artifacts

- JSON: `ops/governance_records/2026-05-18-nc2024-real-run-validation/real_run_validation.json`
- Single-run governance report: `ops/governance_records/2026-05-18-nc2024-real-run-validation/project_validations/nc-reproduction_2026-05-18T0900Z-13c2c88.governance_validation.json`
- Launch log: `/home/zerlinshen/projects/nc-reproduction/runs/launch_wave3_subsample_structure_validation_2026-05-18T0900Z-13c2c88.log`
- Initial failed-launch log: `/home/zerlinshen/projects/nc-reproduction/runs/launch_wave3_subsample_structure_validation_20260518T090019Z.log`
