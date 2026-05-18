# NC2024 legacy artifact cleanup — 2026-05-18

- Date/time local: `2026-05-18T18:22:42+08:00`
- Objective: remove current off-mainline validation/test artifacts after the
  module-flow run had been recorded as evidence-only.
- Cleanup record:
  `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-legacy-artifact-cleanup/`.

## Preserved

- Current NC2024 source of truth:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.
- Project policy:
  `/home/zerlinshen/projects/nc-reproduction/ledger/project_retention_policy.yaml`.
- Module-flow validation report:
  `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-module-flow-validation/REPORT.md`.
- Raw/prepared inputs, configs, notebooks, launch wrappers, source code,
  environments, and governance records.
- Historical lightweight validation evidence:
  `/home/zerlinshen/singlecell_factory/results/small_real_validate_20260426`.

## Deleted

- Evidence-only validation run:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88`.
- Factory smoke result directories:
  `/home/zerlinshen/singlecell_factory/results/test_20260516_*`.
- Legacy result-side ops residue:
  `/home/zerlinshen/singlecell_factory/results/ops`.
- Empty legacy factory output directory:
  `/home/zerlinshen/singlecell_factory/output`.

## Verification

- `canonical_exists=yes`
- `policy_exists=yes`
- `validation_report_exists=yes`
- `deleted_run_exists=no`
- `factory_test_dirs=0`
- `factory_results_ops_exists=no`
- `factory_output_exists=no`
- Reclaimed bytes: `1614159948` (~1.50 GiB).

## Artifact Classification

- `canonical`: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.
- `canonical control-plane`: cleanup record and module-flow validation report.
- `deleted evidence-only`: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88`.
- `deleted legacy/smoke`: factory `test_20260516_*`, `results/ops`, and `output`.

## Next Operator Notes

- Do not cite `2026-05-18T1013Z-13c2c88` as an existing run directory; cite the
  preserved module-flow validation report instead.
- Factory `results/` now intentionally retains only
  `small_real_validate_20260426` as lightweight historical validation evidence.
