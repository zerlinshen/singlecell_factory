# 2026-08-23 - Reference Atlas OOD P0 results

## Bounded outcome

- The amendment is complete as a technical P0 implementation and evidence pass.
- It is not a biological validation, a performance claim, or a production-GPU
  promotion.

## Census materialization

- Run: `/home/zerlinshen/projects/reference-atlas-ood-validation/runs/2026-08-22T1510Z-7a2a0e6/`.
- Status: `completed`; `--verify-only` passed.
- CSR reference: `600 x 61,497`; H5AD SHA-256
  `f3ae9fe89f93726e9ce431150f29226beabb4c8c3eb55c734bb330468e9f4307`.
- Selected join-ID SHA-256:
  `50c7fdcceccd2c29cd672183d13dace5aabcec725892657e656a92fd0c92c93e`.
- Census environment recreation used both lock layers. The RAPIDS explicit
  inventory SHA-256 remained
  `ac6b64abcfd68ae4a7293a12db27b0f63225c83e1758d7b024326c1880cae21b`.

## Isolated benchmark

- Initial fresh child run `/runs/2026-08-22T1527Z-7a2a0e6/` failed when the
  cuML child could not discover CUDA headers; it produced a failure receipt
  without fallback.
- After committing the child CUDA-toolkit bootstrap, clean SHA `3749a2d` ran
  `/runs/2026-08-22T1534Z-3749a2d/` from source SHA-256
  `44b3cf6f98ef890ad9a2723f275f8ac63b51148f5a685e38cfe802d561282ad7` and
  frozen-input SHA-256
  `d5a8d21f872a734ab4d2821e79809cb14ba4355e7b8a6055598a1d8212220d89`.
- CPU completed with median `0.9071675369923469` seconds and child peak RSS
  `3380.578125` MiB. The fixed held-out Microglia OOD gates passed.
- GPU loaded the same frozen input, verified the child CUDA header path, and
  was observed at 522 MiB VRAM. Its three raw consumed-output hashes differed,
  so it wrote a failure receipt and the benchmark verdict is `FAIL_NOT_PROMOTED`.
  No CPU fallback, parity result, GPU timing claim, or promotion is allowed.

## Non-negotiable carry-forward

- Do not change the source, split, labels, held-out negative control, retained
  feature set, threshold contract, or caps to repair this result.
- Any future GPU work must start from this receipt and prove raw within-lane
  repeatability as well as the predeclared integrity, OOD, residency, and
  parity gates.
