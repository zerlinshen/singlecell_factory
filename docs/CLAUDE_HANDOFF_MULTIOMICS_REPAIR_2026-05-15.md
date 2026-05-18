# Claude Handoff: Single Cell / Multi-Omics Repair

Generated: 2026-05-15T22:23:41Z

Audience: Claude or the next repository agent.

Repos covered:
- `/home/zerlinshen/singlecell_factory`
- `/home/zerlinshen/multiomics_r_factory`

## Executive State

The requested Ralph repair loop for the Single Cell, Factory, and Multi-OMX cross-repo multi-omics work is complete and was explicitly terminalized.

Current OMX state:
- `ralph` is inactive.
- `omx cancel` reported `Cancelled: ralph`.
- `omx state list-active --json` returned no active modes.
- `omx status --json` reported `ralph: inactive (phase: cancelled)`.
- Session state file: `/home/zerlinshen/.omx/state/sessions/019e2d69-ff01-7131-9b59-548f8134293d/ralph-state.json`.
- That state file now has `active: false`, `current_phase: cancelled`, and `completed_at: 2026-05-15T22:16:55Z`.

Main outcome:
- Project-root run-id handling is guarded before directory creation.
- Python-to-R bundle contract is synced at schema hash `6c40efe4a4df8ea8069bc1ab63761f51d96e58a0246416da2782640fa503eb62`.
- Bundle schema `v2.2` now has an active ATAC extension and reserved VDJ/Ribo/Hi-C extension slots.
- ATAC ingestion, LSI, peak-to-gene, R export, and R loading/plotting have targeted smoke and contract evidence.
- Densify-sensitive code paths are audited with explicit bounded markers.
- Documentation now records module rationale and supporting literature.

Important caveat:
- Real RNA bundle smoke passed on existing local data.
- No real ATAC h5ad was found locally, so ATAC was validated with synthetic Python-to-R round-trip tests rather than a real ATAC dataset.

## What Was Changed

### `singlecell_factory`

Project-root and run-id governance:
- `workflow/modular/project_paths.py`
  - Uses fallback run id sha `0000000`.
  - Validates `run_id` before creating directories.
- `workflow/modular/cli.py`
  - Puts `RunLedger` under `_run_dir` when project-root mode is active.
- `scripts/pack_run_for_mac.sh`
  - Project-root mode now packages from the actual bundle directory with `-C "$(dirname "$BUNDLE_DIR")" "$(basename "$BUNDLE_DIR")"`.
  - This prevents packaging the legacy `r_bundle` path by mistake.

ATAC module chain:
- `workflow/modular/modules/atac_ingest.py`
  - Writes canonical sparse counts to `obsm["atac_peaks"]`.
  - Keeps compatibility `obsm["X_atac"]`.
  - Stores peak metadata in `uns["atac_peaks"]` and `uns["atac_var"]`.
- `workflow/modular/modules/atac_lsi.py`
  - Fixes TF-IDF/SVD component bound to `min(50, n_peaks)`.
  - Raises for fewer than 2 peaks.
- `workflow/modular/modules/peak_to_gene.py`
  - Requires sparse `obsm["atac_peaks"]` and peak metadata.
  - Marks bounded densify use explicitly.
- `workflow/modular/module_catalog.py`
  - Records ATAC dependencies/order and reserved extension descriptions.

Bundle export and contract:
- `scripts/export_singlecell_r_bundle.py`
  - Adds schema `v2.2` support.
  - Adds `--include-atac`, `atac_lsi_obsm_key`, and `atac_peaks_uns_key`.
  - Exports `extensions/atac/lsi.parquet` with a cell index.
  - Exports `extensions/atac/peaks.parquet` with `peak_id`, `chrom`, `start`, and `end`.
  - Marks ATAC extension active in the manifest.
  - Requires `--schema-version v2.2` when `--include-atac` is used.
- `contracts/bundle_schema.yaml`
  - Marks ATAC active.
  - Keeps VDJ/Ribo/Hi-C reserved.
  - Clarifies table fields and sync hash.
- `contracts/.expected_sha256`
  - Synced to hash `6c40efe4a4df8ea8069bc1ab63761f51d96e58a0246416da2782640fa503eb62`.

Densify governance:
- `workflow/modular/modules/differential_expression.py`
- `workflow/modular/modules/hic_tad.py`
- `tests/test_densify_audit.py`
  - Adds or verifies tokenized audit coverage for executable `toarray` / `todense` usage and `densify-allowed` markers.

Docs and tests:
- `docs/MULTIOMICS_MODULE_RATIONALE.md`
  - Documents module rationale and literature support for ATAC, LSI, peak-to-gene, reserved VDJ/Ribo/Hi-C slots, and densify policy.
- `README.md`
  - Documents bundle `v2.1` / `v2.2`, ATAC, and reserved extension behavior.
- Added or updated targeted tests:
  - `tests/test_project_paths_full.py`
  - `tests/test_cli_project_root_smoke.py`
  - `tests/test_r_bundle_contract.py`
  - `tests/test_wave2a_smoke.py`
  - `tests/test_wave2b_atac_extensions_smoke.py`
  - `tests/test_atac_lsi.py`
  - `tests/test_wave2b_vdj_smoke.py`
  - `tests/test_wave2b_ribo_smoke.py`
  - `tests/test_wave2b_hic_smoke.py`
  - `tests/test_densify_audit.py`

### `multiomics_r_factory`

Run-id parsing:
- `R/cli_utils.R`
  - Adds `RUN_ID_PATTERN`.
  - Adds `validate_run_id`.
  - Rejects unsafe run ids such as `../../bad`.

Bundle reader and ATAC support:
- `R_bundle/io_bundle.R`
  - Accepts `singlecell_r_bundle_v2.2`.
  - Knows extension keys for `marker_resolutions`, `atac`, `vdj`, `hic`, and `ribo`.
  - Returns reserved extension metadata without pretending data exists.
  - Loads active ATAC through `ATAC_MODULE_R_PATH` or repo fallback.
  - Attaches loaded ATAC data under `bundle$extensions$atac$data`.
- `R/atac_module.R`
  - Loads ATAC extension metadata, LSI, and peak tables.
  - Handles reserved wrappers.
  - Preserves bundle root context.
  - Checks LSI rownames and order against bundle observations.
  - `plot_atac_lsi(group_by = ...)` now resolves grouping variables from `bundle$obs`.
  - Adds peak-count plotting support.

Contract and docs:
- `contracts/bundle_schema.yaml`
- `contracts/.expected_sha256`
  - Synced to hash `6c40efe4a4df8ea8069bc1ab63761f51d96e58a0246416da2782640fa503eb62`.
- `README.md`
  - Documents schema `v2.2`, ATAC support, schema hash, and `ATAC_MODULE_R_PATH`.

## Fresh Verification Evidence

Fresh verification was rerun after the Stop-hook reminder and after the architect review fixes.

Python project-root and CLI smoke:

```bash
cd /home/zerlinshen/singlecell_factory
python -m pytest -q tests/test_cli_project_root_smoke.py tests/test_project_paths_full.py --no-cov
```

Result: `19 passed`.

Python-to-R ATAC contract:

```bash
cd /home/zerlinshen/singlecell_factory
python -m pytest -q tests/test_r_bundle_contract.py::test_atac_extension_round_trip -m r_contract -o addopts=''
```

Result: `1 passed`.

ATAC module and wave smoke:

```bash
cd /home/zerlinshen/singlecell_factory
python -m pytest -q tests/test_wave2a_smoke.py tests/test_wave2b_atac_extensions_smoke.py tests/test_atac_lsi.py --no-cov
```

Result: `24 passed`.

Cross-repo contract parity:

```bash
cd /home/zerlinshen/multiomics_r_factory
bash tools/check_contracts_cross_repo.sh
```

Result: contract parity OK at sha `6c40efe4a4df8ea8069bc1ab63761f51d96e58a0246416da2782640fa503eb62`.

Reserved extension smoke:

```bash
cd /home/zerlinshen/singlecell_factory
python -m pytest -q tests/test_wave2b_vdj_smoke.py tests/test_wave2b_ribo_smoke.py tests/test_wave2b_hic_smoke.py --no-cov
```

Result: `25 passed`.

Known note: two NumPy runtime warnings are emitted by the Hi-C fixture and were not blockers.

Densify audit and policy:

```bash
cd /home/zerlinshen/singlecell_factory
python -m pytest -q tests/test_densify_audit.py tests/test_densify_policy.py --no-cov
```

Result: `10 passed`.

R run-id validation:

```bash
/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript -e "source('/home/zerlinshen/multiomics_r_factory/R/cli_utils.R'); stopifnot(validate_run_id('2026-05-14T0001Z-aaaaaaa')); stopifnot(!validate_run_id('../../bad')); cat('RUN_ID_VALIDATION_OK\n')"
```

Result: `RUN_ID_VALIDATION_OK`.

Real local RNA bundle smoke:

Input:
- `/home/zerlinshen/singlecell_factory/rna_velocity_pseudotime_analysis/runtime/results/lung_with_loom_v14/integrated.h5ad`

Result:
- `PY_SCHEMA=singlecell_r_bundle_v2.1`
- `REAL_SCHEMA=singlecell_r_bundle_v2.1`
- `REAL_OBS=500`
- `REAL_UMAP=500,2`
- `REAL_EXPR=500,3`

## Architect Review Status

A native architect subagent reviewed the Ralph result.

First pass found two blockers:
- `pack_run_for_mac.sh` still packaged the old `r_bundle` basename in project-root mode.
- `plot_atac_lsi(group_by = ...)` did not resolve grouping metadata from `bundle$obs`.

Both were fixed and retested.

Follow-up architect review reported:
- No blockers.
- Packaging fix verified.
- ATAC grouping fix verified.
- Targeted tests passed.

## Rationale and Literature Location

The main rationale document is:

`/home/zerlinshen/singlecell_factory/docs/MULTIOMICS_MODULE_RATIONALE.md`

It records the rationale for:
- ATAC count storage as sparse peak-by-cell data.
- TF-IDF plus truncated SVD/LSI for scATAC.
- Peak-to-gene linkage as a bounded bridge layer.
- Reserved slots for VDJ, Ribo, and Hi-C until mature implementations exist.
- Explicit densify markers for sparse matrix safety.

Claude should treat that document as the local rationale source before editing the modules again.

## Current Risks and Cautions

- The repos were already dirty. Do not assume every dirty file was created in this turn.
- Do not revert unrelated user or generated changes.
- Inspect diffs carefully before committing.
- ATAC has synthetic end-to-end evidence but not real local ATAC data evidence.
- The real-data smoke only proves RNA bundle compatibility for a 500-cell subset from an existing h5ad.
- Hi-C reserved smoke emits known fixture warnings; the reserved contract itself passes.

## Suggested Next Steps for Claude

1. If asked to continue implementation, start by checking scoped diffs in both repos rather than resetting or cleaning.
2. If asked to commit, split logically:
   - project-root/run-id governance
   - schema/contract sync
   - ATAC export/R import
   - reserved extensions
   - densify audit
   - docs/tests
3. If real ATAC validation becomes available, run a small h5ad through `--schema-version v2.2 --include-atac` and verify R import plus `plot_atac_lsi(group_by = ...)`.
4. Avoid launching huge NC2024 or remote pipelines unless explicitly requested.

