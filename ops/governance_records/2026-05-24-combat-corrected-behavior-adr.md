# ADR: ComBat batch correction corrected-behavior (WS1 #4)

Date: 2026-05-24
Plan: `.omc/plans/2026-05-24-pipeline-governance-remediation.md` (WS1 #4)
Scope: `workflow/modular/modules/batch_correction.py`

## Defect

`sc.pp.combat(adata, key=batch_key)` rewrites only `adata.X`. The
post-correction `use_rep` selector had no `combat` branch, so it left
`use_rep="X_pca"`. The downstream neighbors/UMAP/Leiden recompute and the
`'after'` batch-mixing metric therefore ran on the STALE pre-correction
`X_pca`, while `batch_correction_status` was still recorded as `completed`.
Net effect: ComBat was a silent no-op for embedding/clustering — the corrected
counts never propagated to the latent structure that every downstream module
consumes.

## Fix

After `sc.pp.combat`, recompute PCA on the corrected `adata.X` into
`obsm['X_pca_combat']` (non-destructively — the pre-combat `X_pca` is restored
for the `'before'` baseline) and set `use_rep='X_pca_combat'`. `X_pca_combat`
is added to the CPU-coercion rep-key list. `batch_correction_use_rep` is now
recorded in metadata for every run.

## Combat-usage inventory

Searched the factory tree for `run_manifest.json` recording
`batch_method=combat` / `method=combat`: ZERO matches. Consistent with the
2026-05 factory-project separation (no scientific run outputs land inside the
factory repo). No prior in-tree combat outputs exist to invalidate; any future
audit of external project runs that selected combat before 2026-05-24 should
treat their embedding/clustering as having been computed on the UNCORRECTED
`X_pca` (the no-op), not on combat-corrected counts.

## Consequences

- A combat run now changes neighbors/UMAP/Leiden relative to the pre-combat
  baseline (intended — this is the correction actually taking effect).
- `obsm['X_pca_combat']` is a new representation key (non-breaking addition).
- Default `batch_method` remains `harmony`; combat is opt-in, so production
  NC/NG runs are unaffected.
