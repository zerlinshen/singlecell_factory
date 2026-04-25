# 2026-04-24 - NC2024 subtype checkpoint audit

## Objective

- Execute the first targeted post-baseline lane for the NC2024 single-cell claim-closure work.
- Determine whether the paper-relevant subtype checkpoint pairs are absent from the data or merely hidden by the current summary surface.
- Avoid re-running the full-cohort stage-1 baseline unless a concrete blocker forces that branch.

## Starting Context

- Canonical prepared input remains:
  - `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- Canonical direct successful stage-1 run remains:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222`
- Canonical controller fallback successful run remains:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552`
- Carry-forward lesson used before this run:
  - treat `large` as a capacity probe only
  - do not rerun full-cohort stage-1 baseline without a concrete blocker
  - use current canonical artifacts as the source of truth
- Existing subtype communication inputs were already present:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900`

## What Was Run

- Host:
  - `ubuntu-tail`
- Working directory:
  - `/home/zerlinshen/singlecell_factory/results`
- Execution mode:
  - `debug_massive` carry-forward default for NC2024, but no new stage-1 compute was launched in this round
- Command or controlling script:
  - ad hoc pair-level checkpoint audit script written into the result directory and executed with:
    - `/home/zerlinshen/conda/bin/conda run -n sc_gpu python /home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105/audit_checkpoint_pairs.py`
- Run directory:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105`

## Outcome

- `success`

## Evidence

- `module_status.csv`:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105/module_status.csv`
- `run_manifest.json`:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105/run_manifest.json`
- pair-level summary:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105/checkpoint_pair_summary.csv`
- pair-level best rows:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105/checkpoint_pair_best_rows.csv`
- raw-vs-current-top20 visibility table:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105/checkpoint_pair_visibility.csv`
- presence matrix:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105/checkpoint_pair_presence_matrix.csv`

## Problems Encountered

- The first attempt used bare remote `python3` and failed because `pandas` was unavailable in that interpreter.
- This was an environment problem, not a logic problem in the audit itself.

## What Was Changed

- No repo code was changed in `/home/zerlinshen/singlecell_factory` during this round.
- A reproducible audit script was written into the result directory itself:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_AUTO_20260423_162105/audit_checkpoint_pairs.py`
- Local reproduction docs were updated to cite the new pair-level checkpoint evidence.

## Resolution

- Re-ran the audit under `sc_gpu`, which has the required Python packages.
- The audit showed that the current checkpoint top20 tables are not the whole story.
- Key carry-forward result:
  - LUAD raw subtype LIANA outputs already contain `LGALS9-HAVCR2`, `NECTIN2-CD226`, `NECTIN2-TIGIT`, and `HLA-F-LILRB1/2`
  - LUSC raw subtype LIANA outputs contain those same axes plus `NECTIN3-TIGIT` and `CD96-NECTIN1`
  - current checkpoint top20 tables still surface only the CTLA4 axis
- This means the immediate checkpoint gap is now a summary-surface / paper-facing postprocess problem, not simple absence of the raw checkpoint pairs.

## Cautions for the Next Run

- Do not infer checkpoint-pair absence from the current subtype top20 tables alone.
- Do not rerun the full-cohort stage-1 baseline merely to chase subtype checkpoint hierarchy.
- Preserve the distinction between:
  - raw subtype LIANA pair presence
  - current paper-facing checkpoint summary surface
- The new audit artifact is evidence-rich, but it is not a replacement for a formal reusable repo-side script or module yet.

## Improvement Ideas

- Formalize the pair-level checkpoint audit as a repo-side rerunnable script or supported postprocess surface.
- Refresh the local Fig. 2 comparison logic so it can cite pair-level checkpoint evidence rather than only CTLA4-dominant top20 tables.
- Keep foetal/onco-foetal comparator work as a separate next lane rather than mixing it into subtype checkpoint closure.

## Classification

- `evidence-only`
