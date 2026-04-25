# 2026-04-24 - NC2024 checkpoint audit formalized

## Objective

- Formalize the pair-level subtype checkpoint audit into a repo-side rerunnable script.
- Verify that the formalized repo-side script reproduces the evidence-only audit outputs.
- Refresh the local Fig. 2 comparison face so it consumes pair-level checkpoint audit outputs instead of only CTLA4-dominant top20 tables.

## Starting Context

- Prior evidence-only checkpoint audit run existed at:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105`
- Existing raw subtype LIANA inputs remained:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/luad/cell_communication/cell_communication_liana.csv`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/lusc/cell_communication/cell_communication_liana.csv`
- Existing current top20 checkpoint summaries remained:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_adenocarcinoma_checkpoint_top20.csv`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_squamous_cell_carcinoma_checkpoint_top20.csv`
- Carry-forward rule retained:
  - do not rerun full-cohort stage-1 baseline for this checkpoint-hierarchy question

## What Was Run

- Host:
  - `ubuntu-tail`
- Working directory:
  - `/home/zerlinshen/singlecell_factory`
- Repo-side script added:
  - `/home/zerlinshen/singlecell_factory/scripts/audit_nc2024_subtype_checkpoint_pairs.py`
- Focused tests added:
  - `/home/zerlinshen/singlecell_factory/tests/test_nc2024_subtype_checkpoint_audit_script.py`
- Documentation updated:
  - `/home/zerlinshen/singlecell_factory/README.md`
  - `/home/zerlinshen/singlecell_factory/PROTOCOL.md`
- Formalized run directory:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_FORMALIZED_AUTO_20260423_165106`
- Verification commands:
  - `/home/zerlinshen/conda/bin/conda run -n sc_gpu python -m pytest -q -o addopts='' tests/test_nc2024_subtype_checkpoint_audit_script.py`
  - `/home/zerlinshen/conda/bin/conda run -n sc_gpu python scripts/audit_nc2024_subtype_checkpoint_pairs.py --luad-raw-csv ... --lusc-raw-csv ... --luad-top20-csv ... --lusc-top20-csv ... --output-dir /home/zerlinshen/singlecell_factory/results --project NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_FORMALIZED_AUTO_20260423_165106`

## Outcome

- `success`

## Evidence

- Script test pass:
  - `4 passed in 0.05s`
- Formalized `module_status.csv`:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_FORMALIZED_AUTO_20260423_165106/module_status.csv`
- Formalized `run_manifest.json`:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_FORMALIZED_AUTO_20260423_165106/run_manifest.json`
- Formalized outputs:
  - `checkpoint_pair_summary.csv`
  - `checkpoint_pair_best_rows.csv`
  - `checkpoint_pair_visibility.csv`
  - `checkpoint_pair_presence_matrix.csv`
- Equivalence check:
  - summary / best_rows / visibility CSV diffs against the prior evidence-only audit were empty

## Problems Encountered

- The first pytest attempt inherited repo-level `--cov` addopts and failed because the environment did not have the expected coverage plugin active in that invocation path.
- This was resolved by running focused pytest with `-o addopts=''`.
- The first remote script version used a dataclass shape that was awkward under the importlib-based test loader and was simplified to plain mapping constants.

## What Was Changed

- Added repo-side script:
  - `/home/zerlinshen/singlecell_factory/scripts/audit_nc2024_subtype_checkpoint_pairs.py`
- Added focused tests:
  - `/home/zerlinshen/singlecell_factory/tests/test_nc2024_subtype_checkpoint_audit_script.py`
- Updated remote docs:
  - `/home/zerlinshen/singlecell_factory/README.md`
  - `/home/zerlinshen/singlecell_factory/PROTOCOL.md`
- Refreshed the local Fig. 2 board build surface and notes on the Mac coordination side.

## Resolution

- The subtype checkpoint audit is now a repo-side rerunnable surface rather than a one-off script inside a result directory.
- The formalized run reproduced the prior evidence-only checkpoint conclusions.
- The checkpoint-hierarchy conclusion remains:
  - raw subtype LIANA outputs contain paper-relevant checkpoint pairs beyond CTLA4
  - the current paper-facing checkpoint summary surface is still partial
- This is now backed by a stable repo-side script, focused tests, and updated docs.

## Cautions for the Next Run

- Keep `raw pair present` distinct from `paper-faithful checkpoint hierarchy reproduced`.
- Continue to avoid full-cohort stage-1 reruns for this question unless the subtype LIANA inputs themselves become invalid.
- The next likely lane is not another checkpoint audit rerun; it is either:
  - a more paper-faithful checkpoint summary surface, or
  - the separate foetal/onco-foetal comparator lane.

## Improvement Ideas

- If checkpoint hierarchy needs to be tightened further, formalize a second-stage summary script that ranks and presents paper-highlighted pairs more directly for LUAD/LUSC.
- Keep the local board assets synced from the formalized audit surface rather than from ad hoc result-dir scripts.

## Classification

- `canonical`
