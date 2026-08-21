# LUSC held-out reference projection real smoke — 2026-08-13

## Objective and route

- Validate the existing `r_multiomics_factory` Seurat reference-projection
  module on retained real LUSC/LuCA data without adding a new module or making
  a biological claim.
- Dev OS task routing resolved the unchanged bounded request to
  `domain_execution` → `skill:execute-and-recover-pipeline` with
  `controller_validation`; the multi-model Dev OS chain was not used for the
  scientific run.

## Inputs and frozen design

- Project: `/home/zerlinshen/projects/lusc-gt-concordance-20260728`.
- Profile/tier: `lusc_zilionis_heldout_v1` / `smoke`.
- Held out: `Zilionis_Klein_2019` query study.
- Materialized split: 1,000 query cells, 4,781 reference cells, 24 labels.
- Scientific parameters stayed frozen: LogNormalize, 2,000 HVGs, 30 PCs,
  dimensions 1:30, seed 20260813, 8 thread caps.

## Recovery

- Retained failed run:
  `/home/zerlinshen/projects/lusc-gt-concordance-20260728/runs/2026-08-13T1321Z-4a6c1e2/`.
- Deterministic root cause: Seurat `FindTransferAnchors` captured 1.28 GiB of
  globals through `future`, exceeding its default 500 MiB protection limit.
- Minimal recovery: the governed R driver now scopes the documented
  `future.globals.maxSize` option to 2 GiB around projection only, restores the
  previous option immediately, records the limit in output parameters, and
  carries a synchronized entrypoint SHA-256 in the maturity registry.

## Successful evidence

- Canonical smoke run:
  `/home/zerlinshen/projects/lusc-gt-concordance-20260728/runs/2026-08-13T1328Z-7c9d3f1/`.
- Materialization: exit 0 in ~5 s; projection: exit 0 in ~19 s; sealed scoring:
  exit 0 in ~2 s.
- Namespace attestation: project root, ground truth, and user home were not
  mounted; measured namespace and mount-plan checks passed.
- Coverage 1.0; macro-F1 0.813375; balanced accuracy 0.850754; exact match
  0.826; weighted F1 0.817991. Bootstrap 95% macro-F1 interval:
  [0.742187, 0.854262].
- Regression evidence: 54 Python contract tests PASS; R entrypoint checks PASS;
  sandbox integration/negative-path suite PASS; R parse and `git diff --check`
  PASS.

## Carry-forward

- This is real smoke evidence, not the frozen full real tier. Do not populate
  maturity evidence IDs, promote `reference_projection/seurat`, or make a
  biological/publication claim yet.
- Before the full 5,736-query/24,000-reference tier, size the `future` globals
  ceiling against the larger reference object or prove a lower-memory Seurat
  path; do not assume the smoke's 2 GiB ceiling scales to the full tier.
- Preserve the three-process leakage boundary: trusted materializer →
  GT-unmounted sandboxed projection → independent sealed scorer.
