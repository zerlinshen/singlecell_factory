# NC2024 Legacy Artifact Cleanup — 2026-05-18

Control-plane cleanup record only. This is not scientific evidence.

## Objective

Remove current off-mainline validation/test artifacts so they do not clutter the
code-facing factory tree or project-root run surface.

## Preserved

- Current NC2024 source of truth:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`
- Project policy:
  `/home/zerlinshen/projects/nc-reproduction/ledger/project_retention_policy.yaml`
- Module-flow validation report:
  `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-module-flow-validation/REPORT.md`
- Raw/prepared inputs, configs, notebooks, launch wrappers, source code,
  environments, and governance records.
- Historical lightweight validation evidence:
  `/home/zerlinshen/singlecell_factory/results/small_real_validate_20260426`

## Deleted

- Evidence-only validation run:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88`
- Factory smoke result directories:
  `/home/zerlinshen/singlecell_factory/results/test_20260516_*`
- Legacy result-side ops residue:
  `/home/zerlinshen/singlecell_factory/results/ops`
- Empty legacy factory output directory:
  `/home/zerlinshen/singlecell_factory/output`

## Why Deletion Was Safe

- The deleted run had already been classified as `evidence-only`.
- Its module status, run path, final shape, fallback behavior, test-contract
  drift, and governance validation result are preserved in:
  `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-module-flow-validation/REPORT.md`
- The remote and local before-every-run memory already record that the canonical
  source of truth remains `2026-05-18T0900Z-13c2c88`.
- Smoke/test directories under factory `results/` were not canonical project
  science and had no current policy pointer.

## Verification

Post-delete check:

- `canonical_exists=yes`
- `policy_exists=yes`
- `validation_report_exists=yes`
- `deleted_run_exists=no`
- `factory_test_dirs=0`
- `factory_results_ops_exists=no`
- `factory_output_exists=no`

Current factory `results/` now contains only:

- `/home/zerlinshen/singlecell_factory/results/small_real_validate_20260426`

Reclaimed bytes: `1614159948` (~1.50 GiB).

Inventory files:

- `pre_delete_inventory.txt`
- `bytes_before.tsv`
- `bytes_after.tsv`
- `post_delete_check.txt`
- `reclaimed_bytes.txt`

## Verdict

`LEGACY_ARTIFACT_CLEANUP_COMPLETE`

The cleanup preserved canonical NC2024 evidence and removed off-mainline
validation/test artifacts that could confuse future code and project navigation.
