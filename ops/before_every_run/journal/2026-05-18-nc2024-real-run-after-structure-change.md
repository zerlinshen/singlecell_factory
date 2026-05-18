# NC2024 real run after structure change — 2026-05-18

Task: clean legacy NC2024 factory records/results, then run one true NC2024 pipeline pass under the new project-root architecture.
Date/time UTC: `2026-05-18T09:03:38Z`
Execution mode: `controller_validation` / structure validation via real subsample run

## Objective

Prove the new architecture with an actual NC2024 pipeline execution, not only a governance validator pass:

- old factory `results/NC2024*` / `results/nc2024*` records are removed after user approval
- raw/reference data and project runs are preserved
- new scientific outputs land under `/home/zerlinshen/projects/nc-reproduction/runs/`
- factory receives only governance/control-plane records

## What was attempted

1. Cleaned legacy factory NC2024 result paths.
2. Tried existing `launch_wave3_nc2024.sh` with `STAGE_A=subsample`.
3. First launch failed early because the old launcher generated invalid run id `YYYY-MM-DDTHHMMZ-wave3-subsample`, while new `project_paths` requires `YYYY-MM-DDTHHMMZ-<7hex>`.
4. Recovered by launching the same real NC2024 subsample command with valid run id `2026-05-18T0900Z-13c2c88`.
5. Patched `/home/zerlinshen/projects/nc-reproduction/runs/launch_wave3_nc2024.sh` so future run ids use the current factory git short SHA.
6. Validated run artifacts and governance reports.

## What succeeded

- Cleanup record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-cleanup-legacy-nc2024-results/`.
- Deleted 48 legacy factory NC2024 entries, `30667532081` bytes; `results/` dropped from `29G` to `27M`.
- True run completed:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/`
- Final AnnData: `5281 x 19504`.
- Module status: 8 `ok`, 1 expected `skipped` (`batch_correction`, one batch only).
- New run governance validation: `pass` with empty severity counts.
- All NC project run validation remains `pass` with legacy warnings only for older runs.
- Factory `results/` has zero NC2024/nc2024 leftovers after the true run.
- Report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-real-run-validation/REPORT.md`.

## What failed or required recovery

- First launcher invocation failed before pipeline work because of invalid run-id format.
- Recovery was direct and successful; launcher was patched so this should not recur.

## Artifact classification

canonical:
- `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/manifest.json`
- `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/python/NC2024_WAVE3_SUBSAMPLE_P15_T1_20260518_170048/run_manifest.json`
- `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/python/NC2024_WAVE3_SUBSAMPLE_P15_T1_20260518_170048/module_status.csv`
- `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/python/NC2024_WAVE3_SUBSAMPLE_P15_T1_20260518_170048/final_adata.h5ad`

evidence-only:
- `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-real-run-validation/`
- launch logs under `/home/zerlinshen/projects/nc-reproduction/runs/launch_wave3_subsample_structure_validation_*.log`

superseded:
- deleted factory legacy result paths listed in `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-cleanup-legacy-nc2024-results/cleanup_manifest.json`

failed exploratory:
- initial failed launcher log for invalid run id only; no scientific partial output was produced from that attempt

## Next operator should remember

Use the new run id grammar. Future NC2024 outputs must be project-owned. Do not recreate `singlecell_factory/results/NC2024*`; factory should hold only code, contracts, validators, docs, and governance records.
