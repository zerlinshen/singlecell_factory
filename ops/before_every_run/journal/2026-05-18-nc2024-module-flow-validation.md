# NC2024 module flow validation — 2026-05-18

- Date/time UTC: `2026-05-18T10:13:38Z`
- Objective: Run a controlled NC2024 project-root validation to verify selected
  modules and execution flow after Linux file-governance cleanup.
- Execution mode: `controller_validation`.
- Remote repo: `/home/zerlinshen/singlecell_factory`.
- Project root: `/home/zerlinshen/projects/nc-reproduction`.
- Run id: `2026-05-18T1013Z-13c2c88`.
- Stage: `subsample` / `P15_T1`.
- Launch wrapper:
  `/home/zerlinshen/projects/nc-reproduction/runs/launch_wave3_nc2024.sh`.
- Launch log:
  `/home/zerlinshen/projects/nc-reproduction/runs/launch_wave3_subsample_validation_20260518T101338Z.log`.
- Output run:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88`.
- Governance validation record:
  `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-module-flow-validation/`.

## Outcome

- Run completed successfully.
- Final AnnData: `5281 x 19504`.
- Module status: 8 `ok`, 1 expected `skipped`.
- Expected skip: `batch_correction` skipped because only one `sample` batch was
  present.
- Project governance validator: `pass`, `control_plane_only = true`.
- Current canonical source of truth remains
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.

## Module Status

| Module | Status | Note |
| --- | --- | --- |
| `cellranger` | `ok` | existing prepared input used |
| `marker_db_loader` | `ok` | 5 DBs / 198 marker entries loaded |
| `qc` | `ok` | completed |
| `doublet_detection` | `ok` | fallback all-singlets after RAPIDS scrublet dtype issue |
| `clustering` | `ok` | GPU path completed |
| `annotation` | `ok` | completed |
| `batch_correction` | `skipped` | expected single-batch skip |
| `context_aware_annotation` | `ok` | context mismatches 0 |
| `differential_expression` | `ok` | GPU DE attempted, CPU fallback used |

## Checks

Passed:

- `scripts/ci/verify_module_references.py`: `45/45 modules verified; 2 project-local`.
- Targeted DAG/project-root/governance tests: `29 passed in 1.08s`.
- Minimal modular pipeline test: `1 passed`.
- Source manifests, module status, checkpoint JSON/H5AD files, figures, DE
  outputs, context metrics, and run ledger exist.

Known test-contract drift:

- Broader targeted pytest including all `tests/test_modular.py` produced
  `78 passed`, `5 failed`.
- Failures reflect stale test expectations:
  - DAG-size assertion expects `31`, but current dependencies are `42`;
  - GPU fallback tests expect old CPU fallback while current default policy is
    `raise`;
  - copy-object expectation conflicts with current M2 in-place raw-preservation
    strategy.

## Artifact Classification

- Current source of truth: existing run
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.
- New run: `evidence-only` module-flow validation.
- Governance report: `canonical` control-plane evidence for this validation
  round.
- Cleanup candidates: the new `1.6G` evidence-only run can be considered in a
  future retention cleanup after this journal, governance report, and validator
  outputs are preserved. Do not delete it in this round.

## Next Operator Notes

- Do not promote `2026-05-18T1013Z-13c2c88` to project source of truth without
  a separate scientific review and policy update.
- The run validates module/process flow, not full article-scale reproduction.
- The Plasma-cell label drift (`65` cells) versus the current source-of-truth
  run is worth noting but is not a process failure.
- If test health becomes the next task, update `tests/test_modular.py` to match
  the current expanded module DAG and explicit GPU failure policy.
