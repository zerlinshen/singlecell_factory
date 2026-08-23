# 2026-08-23 - Reference Atlas query-blind certificate revocation

## Run intent and preflight

- Execution mode: bounded real-data benchmark, not production analysis.
- Canonical factory SHA at launch: `0b66858`; tracked tree clean.
- Source H5AD:
  `/home/zerlinshen/projects/pipeline-scientific-audit-20260805/runs/2026-08-05T1345Z-4e5392e/python/final_adata.h5ad`.
- Verified source SHA-256:
  `44b3cf6f98ef890ad9a2723f275f8ac63b51148f5a685e38cfe802d561282ad7`.
- Preflight capacity: 221 GB free disk, 83 GiB available RAM, swap unused,
  and 22.9 GB free GPU memory.
- Environment: `/home/zerlinshen/conda/envs/sc_gpu_rapids2608/bin/python`.

## Canonical run

- Run root:
  `/home/zerlinshen/projects/reference-atlas-ood-validation/runs/2026-08-23T1054Z-0b66858/`.
- Factory command: `scripts/benchmark_reference_mapping.py` with the verified
  source hash and project-owned run ID.
- Technical verdict: `FAIL_NOT_PROMOTED`; no production CPU/GPU dual run and
  no fallback inside the failed GPU benchmark lane.
- Frozen-input manifest SHA-256:
  `31e30d69d46b6d1deb262899f2508e942542d3f35035d629334bfdfc29e0984c`.

## Scientific-method result

- Feature fitting route: `reference_only_sparse_variance_v1`.
- Receipt proves 9,496 reference cells, 25,519 input genes, 3,000 selected
  genes, `query_cells_consumed: false`, and
  `independently_query_blind: true`.
- Query size: 4,143 cells; labels are pipeline-derived proxies.
- All four predeclared OOD gates failed on the CPU lane:
  - known acceptance coverage `0.5127118644` (required at least `0.80`)
  - known all-cell macro-F1 with rejected cells as unknown `0.5748506830`
    (required at least `0.70`)
  - held-out Microglia rejection recall `0.4305177112` (required at least
    `0.80`)
  - OOD-minus-known rejection-rate delta `-0.0567704244` (required at least
    `0.30`)
- CPU repetitions were exact. The GPU lane failed because
  `reference_confidence` changed for one cell in repetition 2 and all three
  pairwise distance comparisons exceeded `rtol=atol=1e-6`; the true maximum
  drift was `2.1905e-6` (repetition 1 versus 3).
- Read-only CPU/GPU evidence comparison found one candidate-label difference,
  five confidence differences, and zero differences in assignment status,
  final cell type, or OOD flag. For the one candidate difference, the CPU
  candidate matched the known proxy label, but both lanes rejected the cell at
  confidence `0.4`; this is not biological proof that CPU is correct.

## Routing decision

- Revoke certificate `reference-knn-trevino-20260822-v2` because its source-wide
  retained-feature mask was not independently query-blind.
- Preserve the historical `11.04x` measurement as superseded performance
  evidence; do not use it for current routing or a biological claim.
- Current policy status:
  `revoked_after_query_blind_validation_failure`.
- Production `auto` routes `trevino-fetal-cortex-v1` directly to CPU before
  mapping. Explicit GPU fails closed because no promoted certificate exists.
- CPU is the conservative lane, not biological ground truth. A future GPU
  promotion requires a referenced methodological correction plus a new clean
  real-data run passing query-blind OOD, repeatability, parity, and residency
  gates.

## Evidence and retention

- Preserve the source H5AD, frozen inputs, receipts, three CPU repetitions,
  three GPU repetition artifacts, logs, and failure manifest.
- Current run is negative canonical evidence and must not be cleaned as a
  redundant failure.
- The local Mac mirror was not available on this Ubuntu host; this journal and
  project-owned run remain the authoritative handoff surfaces.
