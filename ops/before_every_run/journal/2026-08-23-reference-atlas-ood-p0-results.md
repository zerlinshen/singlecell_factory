# 2026-08-23 - Reference Atlas OOD P0 results

## Bounded outcome

- The amendment is complete as a technical P0 implementation and evidence pass.
- The clean process-isolated replacement run passes the declared scientific,
  repeatability, residency, and parity gates and promotes a routing certificate
  only for `trevino-fetal-cortex-v1` with cuML `26.08.00`.
- It is not biological ground truth or a general GPU-module promotion.

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

### Scientific repeatability amendment and clean rerun

- The failed raw-hash rule above conflated byte identity with scientific
  identity. Before rerunning, the contract was amended to require exact
  candidate label, confidence, assignment status, final reference type, and OOD
  values, while applying a predeclared `rtol=atol=1e-6` only to floating-point
  mean-neighbor distance. Every repetition is still persisted and SHA-verified.
- A verifier bug in `2026-08-22T1548Z-ea425e8` compared parsed CSV
  reserialization with original file bytes; that failed run is retained. The
  verifier was corrected before the final clean commit and run.
- Clean run:
  `/home/zerlinshen/projects/reference-atlas-ood-validation/runs/2026-08-22T1558Z-c38053d/`.
  Factory SHA `c38053d`, frozen-input SHA-256
  `d5a8d21f872a734ab4d2821e79809cb14ba4355e7b8a6055598a1d8212220d89`.
- Exact scientific-output SHA-256 across CPU/GPU and all repetitions:
  `a24437b8da0c99b0af2179c9caf49ada80de01f29c744864a93b6ce4d5202a83`.
  Maximum GPU within-lane distance drift was `4.7684e-7`; maximum CPU/GPU
  distance difference was `9.2e-7`.
- Held-out Microglia rejection recall `0.9564`, known acceptance coverage
  `0.9319`, and known all-cell macro-F1 with rejected cells treated as unknown
  `0.8243`; every predeclared gate passed.
- CPU/GPU medians were `0.9080`/`0.08227` seconds (`11.04x`). The GPU child was
  CUDA resident and observed at 522 MiB peak PID-scoped VRAM. Production runs
  do not dual-run: exact promoted domain/version selects GPU; unknown or drifted
  cases select CPU before mapping; explicit GPU failure remains fail-loud.

## Non-negotiable carry-forward

- Do not change the source, split, labels, held-out negative control, retained
  feature set, threshold contract, or caps to repair this result.
- Any future certificate or backend-version change must start from these
  receipts and rerun the predeclared integrity, exact-scientific-output,
  bounded-distance, OOD, residency, and parity gates.
