# NC2024 Module Flow Validation — 2026-05-18

Control-plane validation record only. This is not a new scientific source of
truth.

## Objective

Run a controlled NC2024 project-root validation to check whether the selected
modules and execution flow still behave as expected after the Linux file
governance work.

## Execution Mode

- `execution_mode = controller_validation`
- Remote repo: `/home/zerlinshen/singlecell_factory`
- Project root: `/home/zerlinshen/projects/nc-reproduction`
- Run id: `2026-05-18T1013Z-13c2c88`
- Stage: `subsample` / `P15_T1`
- Launch wrapper:
  `/home/zerlinshen/projects/nc-reproduction/runs/launch_wave3_nc2024.sh`
- Launch log:
  `/home/zerlinshen/projects/nc-reproduction/runs/launch_wave3_subsample_validation_20260518T101338Z.log`
- Output run:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88`

## Module Status

`module_status.csv`:

| Module | Status | Message |
| --- | --- | --- |
| `cellranger` | `ok` | completed |
| `marker_db_loader` | `ok` | completed |
| `qc` | `ok` | completed |
| `doublet_detection` | `ok` | completed |
| `clustering` | `ok` | completed |
| `annotation` | `ok` | completed |
| `batch_correction` | `skipped` | Only one batch found in `sample`. |
| `context_aware_annotation` | `ok` | completed |
| `differential_expression` | `ok` | completed |

This matches the current single-sample validation expectation: all selected
modules pass, and `batch_correction` is skipped because the run contains only
one batch/sample.

## Output Checks

- Final AnnData:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88/python/NC2024_WAVE3_SUBSAMPLE_P15_T1_20260518_181339/final_adata.h5ad`
- Final shape: `5281 x 19504`
- `X_pca` and `X_umap` are present.
- Root `manifest.json`, producer `run_manifest.json`, `module_status.csv`,
  checkpoints, QC figures, clustering figures, annotation outputs,
  context-validation metrics, DE tables, and run ledger are present.
- Project governance validator:
  `overall_status = pass`, `control_plane_only = true`, 1 run inspected.
- Validator outputs:
  - `nc-reproduction_2026-05-18T1013Z-13c2c88.governance_validation.json`
  - `nc-reproduction_2026-05-18T1013Z-13c2c88.governance_validation.md`

## Runtime Fallbacks And Expected Skips

- `batch_correction`: expected skip because only one batch was present.
- `doublet_detection`: `rapids_singlecell` scrublet rejected the data type and
  the pipeline used `fallback_all_singlets`; this is recorded as `ok`.
- `clustering`: GPU backend ran successfully (`clustering_backend = gpu`).
- `differential_expression`: attempted GPU DE, then fell back to CPU because
  current `rapids_singlecell.tl` lacks `rank_genes_groups`; this fallback is
  recorded in the manifest.

## Comparison To Current Source Of Truth

Current canonical NC2024 project source of truth remains:

`/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`

The new validation run has matching technical dimensions and module statuses:

- `cells_after_qc = 5281`
- `genes_after_qc = 19504`
- `context_validation_mismatches = 0`
- `batch_correction_status = skipped_single_batch`
- `de_backend = cpu`

Observed biological-label drift:

- Current source of truth had context-aware cell types:
  `T cell`, `Myeloid cell`, `B cell`.
- This validation run had:
  `T cell`, `Myeloid cell`, `B cell`, plus `Plasma cell` (`65` cells).

This is not a module-flow failure, but it means this run should remain
`evidence-only` unless a later review intentionally promotes it and updates
`ledger/project_retention_policy.yaml`.

## Module And Flow Tests

Passed:

- `python scripts/ci/verify_module_references.py`
  - `verify_module_references: OK (45/45 modules verified; 2 project-local)`
- `pytest --no-cov -q tests/test_full_dag_topo.py tests/test_module_catalog_ext.py tests/test_project_paths_full.py tests/test_project_governance_validator.py`
  - `29 passed in 1.08s`
- `pytest --no-cov -q tests/test_modular.py::test_modular_pipeline_minimal`
  - `1 passed`, with expected QC plotting warnings.

Known test-contract drift:

- `pytest --no-cov -q tests/test_full_dag_topo.py tests/test_module_catalog_ext.py tests/test_project_paths_full.py tests/test_project_governance_validator.py tests/test_modular.py`
  produced `78 passed`, `5 failed`.
- Failures are in `tests/test_modular.py`:
  - stale DAG-size assertion expects `31`, current `MODULE_DEPENDENCIES` has
    `42`;
  - GPU fallback tests expect old CPU-fallback behavior, while current default
    `SC_GPU_FAILURE_POLICY=raise` poisons adata and raises
    `ClusteringContractViolation`;
  - one test expects GPU path to receive a copied AnnData object, while current
    M2 behavior intentionally runs in-place with raw preservation.

These are test-contract drift findings, not failures of the project-root
validation run.

## Verdict

`MODULE_FLOW_VALIDATION_PASS_WITH_TEST_CONTRACT_DRIFT`

The selected NC2024 validation modules and project-root flow behave as expected
for the controlled `P15_T1` run. Do not promote this run over the current
canonical source of truth without a separate scientific review and policy
update.
