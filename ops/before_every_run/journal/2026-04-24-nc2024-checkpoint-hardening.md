# 2026-04-24 - NC2024 checkpoint hardening

## Objective

- Implement the targeted solutions for the remaining checkpoint-audit risks.
- Add a stable repo-side verification wrapper.
- Add a higher-level checkpoint hierarchy summary surface.
- Reduce script-only testing fragility by moving core audit logic into an importable module.

## Starting Context

- Canonical prepared input remained:
  - `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- Canonical full-cohort stage-1 runs remained unchanged.
- Existing formalized audit run existed at:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_FORMALIZED_AUTO_20260423_165106`
- Current carry-forward lesson before this work:
  - raw subtype LIANA pair presence and paper-facing checkpoint hierarchy are separate questions

## What Was Run

- Host:
  - `ubuntu-tail`
- Working directory:
  - `/home/zerlinshen/singlecell_factory`
- New repo-side surfaces added:
  - `/home/zerlinshen/singlecell_factory/workflow/modular/reporting/nc2024_checkpoint_audit.py`
  - `/home/zerlinshen/singlecell_factory/scripts/audit_nc2024_subtype_checkpoint_pairs.py`
  - `/home/zerlinshen/singlecell_factory/scripts/summarize_nc2024_checkpoint_hierarchy.py`
  - `/home/zerlinshen/singlecell_factory/scripts/verify_nc2024_subtype_checkpoint_audit.sh`
- Focused tests:
  - `/home/zerlinshen/singlecell_factory/tests/test_nc2024_subtype_checkpoint_audit_script.py`
- Verification commands:
  - `/home/zerlinshen/conda/bin/conda run -n sc_gpu python -m pytest -q -o addopts='' tests/test_nc2024_subtype_checkpoint_audit_script.py`
  - `bash scripts/verify_nc2024_subtype_checkpoint_audit.sh`
- Verified run directory:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_VERIFIED_AUTO_20260423_172758`

## Outcome

- `success`

## Evidence

- Focused pytest:
  - `6 passed in 0.06s`
- Verified audit run artifact:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_VERIFIED_AUTO_20260423_172758`
- Verified hierarchy summary outputs:
  - `checkpoint_hierarchy_report.csv`
  - `checkpoint_hierarchy_report.json`
  - `checkpoint_hierarchy_memo.md`
- Wrapper proof:
  - `bash scripts/verify_nc2024_subtype_checkpoint_audit.sh` completed and returned `RUN_DIR=/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_VERIFIED_AUTO_20260423_172758`
- Fresh local Fig. 2 refresh consumed the new hierarchy report and regenerated the board on the Mac coordination side.

## Problems Encountered

- Bare remote `pytest` inherited repo-level `--cov` addopts that were not suitable for this focused path.
- This was resolved by explicitly using:
  - `python -m pytest -q -o addopts='' tests/test_nc2024_subtype_checkpoint_audit_script.py`
- A brittle markdown-string assertion in the hierarchy-summary test caused a false-negative and was replaced with a more stable assertion.

## What Was Changed

- Moved core audit logic into an importable module:
  - `/home/zerlinshen/singlecell_factory/workflow/modular/reporting/nc2024_checkpoint_audit.py`
- Kept CLI scripts thin:
  - audit CLI
  - hierarchy-summary CLI
  - verification wrapper
- Updated remote docs:
  - `/home/zerlinshen/singlecell_factory/README.md`
  - `/home/zerlinshen/singlecell_factory/PROTOCOL.md`
- Updated local Fig. 2 builder and notes to consume the verified hierarchy report.

## Resolution

- The subtype checkpoint lane now has:
  - stable core logic in an importable module
  - a thin audit CLI
  - a thin hierarchy-summary CLI
  - a one-command verification wrapper
- The hierarchy surface now goes beyond raw pair presence and produces a reusable claim/report layer.
- The scientific conclusion remains unchanged:
  - raw subtype LIANA outputs contain paper-relevant checkpoint pairs beyond CTLA4
  - checkpoint hierarchy still remains `partial`

## Cautions for the Next Run

- Do not confuse `supported_raw_only` with full paper-facing hierarchy reproduction.
- Continue to use the wrapper or the same `-o addopts=''` pytest pattern for focused verification.
- This work still does not touch the foetal/onco-foetal comparator gap.

## Improvement Ideas

- If checkpoint closure needs another step, the next surface should be a more paper-faithful ranking/presentation layer rather than another raw pair audit.
- If the hierarchy report becomes widely reused, consider promoting its schema into a dedicated documented reporting contract.

## Classification

- `canonical`
