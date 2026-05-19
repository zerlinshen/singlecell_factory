# 2026-05-18 - two real-dataset final factory validation

## Objective

- Prove the rebuilt three-factory management model with two real datasets: NC2024 NSCLC and Cell/Trevino.
- Validate `singlecell_factory -> r_multiomics_factory -> plotting_factory` project-root handoff without rerunning heavy upstream compute.
- Align produced bundles/R outputs with paper/process-data evidence and record remaining scientific boundaries honestly.

## Starting Context

- NC2024 canonical project-root run: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.
- Cell/Trevino human-facing evidence run: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88`.
- Cell/Trevino linked pipeline run for bridge validation: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88`.
- Carry-forward lessons: project roots are the scientific output surface; factory `results/` must not regain canonical NC2024/Cell outputs; Cell/Trevino remains a conditional public-resource reproduction, not FASTQ/fragments/BPNet exact parity.

## What Was Run

- Host: remote server direct shell under `/home/zerlinshen`.
- Working directories:
  - `/home/zerlinshen/singlecell_factory`
  - `/home/zerlinshen/r_multiomics_factory`
  - `/home/zerlinshen/plotting_factory`
- Controlling actions:
  - Re-exported NC2024 and Cell/Trevino compact R bundles into project-root `python/bundle/` directories.
  - Ran `r_multiomics_factory/scripts/plot_remote_bundle.R` for both project roots.
  - Added a two-dataset final validator and refreshed NC2024 architecture validation for the current project-root run.
- Run directories:
  - NC2024: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`
  - Cell/Trevino bridge: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88`
  - Cell/Trevino paper evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88`

## Outcome

- `success` for the three-factory project-root bridge and governance validation.
- Scientific interpretation remains bounded:
  - NC2024 validation is retained P15_T1 project/module/bridge evidence, not full article-scale multi-cohort annotation.
  - Cell/Trevino validation is conditional public-resource reproduction plus bridge validation, not raw FASTQ/fragments/BPNet parity.

## Evidence

- Governance report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-two-real-dataset-final-validation/REPORT.md`.
- Machine JSON: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-two-real-dataset-final-validation/two_real_dataset_final_validation.json`.
- NC2024 final shape: `5281 x 19504`; module status: 8 `ok`, `batch_correction` expected `skipped`.
- NC2024 bundle: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/python/bundle/`.
- NC2024 R outputs: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/r/` (8 files).
- Cell/Trevino final shape: `55653 x 25519`; all linked pipeline modules `ok`.
- Cell/Trevino bundle: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/bundle/`.
- Cell/Trevino R outputs: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/r/` (8 files).
- Cell/Trevino paper-evidence quality gate: `conditional`, with expected resource gaps only; no unexpected false checks, missing paths, or hash mismatches.

## Problems Encountered

- Old NC2024 architecture validator still targeted deleted legacy `singlecell_factory/results/nc2024_*_v2` paths instead of current project-root truth.
- `scripts/export_singlecell_r_bundle.py` documentation said `--output` is ignored with `--project-root`, but argparse still required `--output`.
- `r_multiomics_factory/scripts/plot_remote_bundle.R` used lung-only marker defaults, producing weak marker choices for Cell/Trevino neural data.
- `R_bundle/io_bundle.R` looked for the root run manifest in the wrong location for project-root bundles under `<run>/python/bundle`.
- `plotting_factory` R helpers used deprecated `element_rect(size=...)` with current ggplot2.

## What Was Changed

- `singlecell_factory/scripts/validate_nc2024_architecture_contract.py`: retargeted to current NC2024 project-root run, bundle, R outputs, and bridge symlinks.
- `singlecell_factory/scripts/validate_two_realdata_final.py`: added two-real-dataset no-rerun gate and report writer.
- `singlecell_factory/scripts/export_singlecell_r_bundle.py`: made `--output` required only for legacy no-project-root exports.
- `r_multiomics_factory/R_bundle/io_bundle.R`: fixed project-root manifest lookup and preserved SHA mismatch warning semantics.
- `r_multiomics_factory/scripts/plot_remote_bundle.R`: added dataset-aware marker selection for NC2024 and Cell/Trevino.
- `plotting_factory/r/general/dim_plots.R` and `r/general/expression_plots.R`: switched border styling from deprecated `size` to `linewidth`.
- Project retention policies for NC2024 and Cell/Trevino now preserve the final bridge-validation bundles, R outputs, and governance report.

## Resolution

- Both real datasets exported project-root bundles and produced project-root R outputs.
- Final validator passed with `PASS_TWO_REAL_DATASET_FACTORY_BRIDGE_VALIDATED`.
- NC2024-only architecture validator passed with `--verify-sha`.
- Focused Python tests, R source smoke, and compilation checks passed before this note; the full verification matrix is tracked in the Autopilot handoff.

## Cautions for the Next Run

- R-factory SHA warnings can be expected when old run manifests record `r_factory_sha_at_manifest_write` before a later bundle export records a newer `r_factory_sha_at_export`; treat as provenance unless the project policy requires a clean SHA match.
- Do not cite NC2024 P15_T1 as full article-scale multi-cohort truth.
- Do not upgrade Cell/Trevino conditional public-resource evidence to exact raw-data parity without staging missing FASTQ/fragments/BPNet resources.
- Keep `singlecell_factory/results/` out of the scientific source-of-truth path.
- Local Mac mirror path was not present in this Linux session, so this note updates the remote canonical journal only.

## Improvement Ideas

- Add a small CI job that runs `validate_two_realdata_final.py` in no-output mode against pinned project roots when those roots are mounted.
- Add a stable marker-preference config for future non-lung datasets instead of growing hard-coded marker lists.
- Consider recording project-root bundle/R-output preservation in a shared schema block for all reproduction projects.

## Classification

- `canonical`: source-of-truth project policies and retained source runs named above.
- `canonical`: final governance record `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-two-real-dataset-final-validation/`.
- `evidence-only`: project-root `python/bundle/` and `r/` bridge outputs generated for this final validation.
- `superseded`: legacy factory-root NC2024 result assumptions inside older validator/docs.
- `failed exploratory`: none retained from this round.
