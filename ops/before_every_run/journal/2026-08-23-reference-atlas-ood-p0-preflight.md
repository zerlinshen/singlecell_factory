# 2026-08-23 - Reference Atlas OOD P0 preflight

## Objective

- Execute the approved bounded Reference Atlas OOD P0 amendment without
  changing the promoted RAPIDS environment or making a biological claim.
- Materialize a 600-cell Census reference and run a fresh process-isolated
  Trevino CPU/GPU technical benchmark from a clean implementation SHA.

## Clean implementation and fixed inputs

- Factory worktree: `dev-os/20260823-reference-atlas-ood-p0` at
  `2fd81ea7edc8f8660b881345b92150f507f7aa1e` (`2fd81ea`), clean before launch.
- Trevino source:
  `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/final_adata.h5ad`.
- Frozen source SHA-256:
  `44b3cf6f98ef890ad9a2723f275f8ac63b51148f5a685e38cfe802d561282ad7`.
- Census lock SHA-256 values: Conda/base
  `0dcdc3a831b7dcddb64d5382735300bd0024342d42960a5f165dc8e369c4e809`;
  pip `817ad3f2d7bb0394c167f21e42ff47452c0536f292fb7ad6b368daea66577ea3`.
- GPU explicit-package inventory SHA-256 before launch:
  `ac6b64abcfd68ae4a7293a12db27b0f63225c83e1758d7b024326c1880cae21b`.

## Environment and execution boundary

- Census lane: `sc_census_io`, Python 3.12.14; recreation must use both exact
  lock layers in a temporary prefix.
- Benchmark lane: `sc_gpu_rapids2608`, Python 3.13.15; NVIDIA GeForce RTX
  5090 D v2, driver 595.84, 24,455 MiB VRAM.
- Project-owned output root:
  `/home/zerlinshen/projects/reference-atlas-ood-validation/runs/`.
- No output may be written under this factory worktree. `SC_REQUIRE_PROJECT_ROOT=1`
  will be set for materialization and benchmark launches.

## Preconditions and stop conditions

- Reference mapping defaults to CPU. Explicit GPU requires cuML/CUDA residency
  and must fail rather than falling back to sklearn.
- Do not alter the source hash, held-out Microglia negative control, upstream
  retained feature set, threshold contract, labels, or caps after observing a
  result.
- Retain a complete failure receipt and stop promotion if Census recreation,
  materialization verification, frozen-input identity, GPU residency, OOD
  gates, or parity gates fails. One new Census run is allowed only for a
  demonstrably transient network failure.

## Preflight classification

- This is a bounded technical validation. Trevino labels are pipeline-derived
  proxy labels, not biological ground truth.
- The 2026-08-22 cumulative-process 9.42x result is historical and
  non-routing; new attributable timing must come only from the isolated run.
