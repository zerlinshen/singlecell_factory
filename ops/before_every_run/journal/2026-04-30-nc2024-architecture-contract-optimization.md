# 2026-04-30 - NC2024 Architecture Contract Optimization

## Objective

- Review and optimize architecture boundaries for `singlecell_factory` and
  `multiomics_r_factory`: single-cell module hierarchy, Python-to-R bundle
  handoff, and no-rerun remote controller validation.

## Starting Context

- Canonical NC2024 v2 outputs:
  - `/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2`
  - `/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2`
- Carry-forward lesson: do not launch another full-cohort run for governance or
  handoff validation when existing v2 artifacts can prove the contract.
- Existing dirty repo state was preserved:
  - `/home/zerlinshen/singlecell_factory/.githooks/pre-commit`
  - `/home/zerlinshen/singlecell_factory/README.md`
  - `/home/zerlinshen/multiomics_r_factory/README.md`
  - untracked notebooks/reports and `.githooks/` / `.omc/`

## What Was Run

- Host: `ubuntu-tail`
- Working directories:
  - `/home/zerlinshen/singlecell_factory`
  - `/home/zerlinshen/multiomics_r_factory`
- Execution mode: `controller_validation`
- No pipeline execution happened. No H5AD was opened for full-object analysis.

Commands:

```bash
/home/zerlinshen/conda/bin/conda run -n sc_gpu python -m pytest -q -o addopts="" \
  tests/test_module_catalog.py \
  tests/test_singlecell_r_bundle_export.py \
  tests/test_modular.py::test_new_modules_in_dag \
  tests/test_modular.py::test_dependency_resolution_auto_includes \
  tests/test_scale_mode_preset_compat.py

python3 scripts/validate_nc2024_architecture_contract.py

bash scripts/ci/check_bridge_symlink.sh

/home/zerlinshen/conda/bin/conda run -n r_multiomics Rscript -e \
  "source('R_bundle/io_bundle.R'); cat('io_bundle_source_ok\n')"

/home/zerlinshen/conda/bin/conda run -n r_multiomics Rscript -e \
  "source('R_bundle/io_bundle.R'); m <- jsonlite::fromJSON('/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/r_bundle/bundle_manifest.json', simplifyVector=TRUE, simplifyDataFrame=FALSE); validate_bundle_v2_manifest(m); cat('nc2024_bh_v2_manifest_ok\n')"
```

## Outcome

- `success` for Python module hierarchy/bundle tests and NC2024 controller
  validation.
- `partial` for R-side v2 plotting: manifest validation passes, but full parquet
  plotting cannot run until R package `arrow` is installed in `r_multiomics`.

## Evidence

- Focused Python tests: `16 passed, 1 warning`.
- Architecture validation:
  - `status = ok`
  - tumor v2 exported cells: `803784`
  - B/H v2 exported cells: `6426`
  - latest ledgers:
    - `/home/zerlinshen/singlecell_factory/ops/run_ledger/nc2024_tumor_20260426_v2_20260426_202319.json`
    - `/home/zerlinshen/singlecell_factory/ops/run_ledger/nc2024_bh_20260426_v2_20260426_202320.json`
- Bridge symlink check:
  - `bridges/local_r_pipeline_macbook/R -> ../../../multiomics_r_factory/R`
  - `bridges/local_r_pipeline_macbook/R_bundle -> ../../../multiomics_r_factory/R_bundle`
- R manifest validation:
  - `io_bundle_source_ok`
  - `nc2024_bh_v2_manifest_ok`
  - `requireNamespace("arrow", quietly=TRUE)` returned `FALSE`

## Problems Encountered

- `pytest` default project addopts include coverage flags, but the active
  environment does not provide pytest-cov. The focused test run used
  `-o addopts=""`.
- `r_multiomics` lacks R package `arrow`, so `scripts/plot_remote_bundle_large.R`
  cannot plot current v2 parquet bundles yet. No dependency was installed
  because the workspace rule says no new dependencies without explicit request.

## What Was Changed

- `singlecell_factory`:
  - added `workflow/modular/module_catalog.py`
  - derived `MODULE_DEPENDENCIES` and CLI optional-module help from the catalog
  - added v2 mtx sidecar records to compact bundle manifests
  - added `scripts/validate_nc2024_architecture_contract.py`
  - added/updated focused tests
  - updated `README.md`, `AI_AGENT_PROTOCOL.md`, and `PROTOCOL.md`
- `multiomics_r_factory`:
  - fixed streaming read of `marker_expr.mtx.gz`
  - added v2 manifest validation to `R_bundle/io_bundle.R`
  - routed both bundle plotting scripts through `read_bundle()`
  - documented v1/v2 bundle contract and the current `arrow` environment gap

## Resolution

- The architecture contract is now executable and validated against the current
  NC2024 v2 artifacts without rerunning the cohort.
- The upstream/downstream bundle boundary is stricter: v2 mtx sidecars are now
  manifest-backed, and R refuses malformed v2 manifests before plotting.

## Cautions for the Next Run

- Do not cite this as a new biological run; it is a controller/contract
  validation round.
- Install or activate an R environment with `arrow` before using the current v2
  parquet bundles in `multiomics_r_factory/scripts/plot_remote_bundle*.R`.
- Keep `bridges/local_r_pipeline_macbook/R` and `R_bundle` as symlinks only.

## Improvement Ideas

- Extend the architecture validator to check `mac_assets/*.tgz.sha256`.
- Add a tiny local v2 parquet fixture for R-side CI once `arrow` is available.
- Consider promoting recovery manifests into a first-class pipeline artifact.

## Classification

- canonical:
  - `/home/zerlinshen/singlecell_factory/workflow/modular/module_catalog.py`
  - `/home/zerlinshen/singlecell_factory/scripts/validate_nc2024_architecture_contract.py`
  - `/home/zerlinshen/multiomics_r_factory/R_bundle/io_bundle.R`
- evidence-only:
  - current NC2024 v2 run dirs and run ledgers listed above
- superseded:
  - scattered module-list strings in CLI/docs as the primary source of hierarchy
- failed exploratory:
  - none
