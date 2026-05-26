# Integration Calibration Bench — Trevino RNA Atlas — 2026-05-26

## Scope

- Task: T4 of the integration-bench team run (worker-2). CALIBRATION / RELATIVE
  lane ONLY — relative scoreboard of batch-integration embeddings on the real
  Trevino GSE162170 RNA atlas. NEVER emits a Harmony-vs-scVI winner verdict.
- Execution mode: controller_validation, bounded real-data run under explicit
  timeouts. Host: `/home/zerlinshen` on `ubuntu-tail`, env `sc_gpu`.
- Outputs (NEVER inside factory tree):
  `/home/zerlinshen/projects/integration-bench-20260526/`.

## What ran

- Assembled the real RNA atlas from the PINNED G5 path
  `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_geo/`
  (`GSE162170_rna_counts.tsv.gz` + `GSE162170_rna_cell_metadata.txt.gz`):
  57,868 cells x 33,355 genes, raw integer counts -> `layers['counts']` for scVI.
- Methods (6 embeddings, all PERSISTED to `embeddings/`): baseline `X_pca`;
  Harmony via harmonypy-direct; scVI seed sweep (seeds 0/1/2); shuffle-label
  negative control. (Harmony-extreme-theta=100 neg control attempted but went
  SINGULAR on full data — recorded, not hidden.)
- Scored with the pinned scib-metrics engine (worker-1 T1 harness); audit built
  via T3 verdict_tier + over-correcting machinery.

## Gates (all satisfied)

- G5: two-part entry guard PASS — batch_key=`sample` 8 levels spanning 3 technical
  batches. The single-batch multiome subset (`data/raw/trevino_2021_brain`,
  Sample.Batch=b2020_11 single level despite 3 replicate-library Sample.IDs)
  HARD-FAILS part (b); proof in `g5_single_batch_hardfail_proof.json`.
- G2: bio scored vs AUTHOR `seurat_clusters` (23 levels c0-c22), source stamped
  `author_seurat_partition_23cl_NOT_celltype`; cluster_names->celltype join NOT
  attempted (symlink confirmed DEAD).
- G4: within-age replicate-pair conditioned concordance reported PER AGE
  (w16/w20 clean technical-batch; w21+w24 CONFOUNDED, share b2020_02). Never pooled.
- G6: output stamped `calibration_relative_only`, verdict_tier
  `relative_ranking_only`, no winner claim.

## Key findings

- CRITICAL production finding (logged to `.omc/plans/open-questions.md`): the
  DEFAULT production Harmony backend (`BatchCorrectionModule._run_harmony`) is
  BROKEN in sc_gpu on BOTH paths — GPU rapids/cupy `CUBLAS_STATUS_NOT_INITIALIZED`
  (+ CUDA-context corruption cascading to scVI) and CPU scanpy1.12<->harmonypy0.2.0
  wrapper shape bug. Bench used harmonypy.run_harmony DIRECT (canonical
  Korsunsky-2019 algorithm; scanpy wrapper bypassed) — documented as a deviation,
  not relabeled as production output. Fix = separate production PR.
- Negative-control falsifiability: the shuffle-label control is the robust anchor
  (ARI 0.486 -> 0.0002 = total cluster-structure collapse). Harmony-extreme-theta
  is NOT a reliable over-corrector (singular at theta=100 on full data; ~real
  harmony at lower theta). Aggregate over-correcting flag uses a DATA-DERIVED
  T=0.346 (baseline_bio - shuffle_bio) -> 0 flags on real methods-under-test
  (scVI/Harmony do not shred biology). A structure-sensitive detector (ARI/NMI
  vs author partition, pre-registered drop>50%) reported SIDE-BY-SIDE; control
  fires (PASS), no method-under-test flagged.
- scVI seed band tight: bio mean 0.612 +/- 0.013, mixing 0.437 +/- 0.003 (n=3).

## Cautions / lessons for future runs

- The counts TSV is an R `write.table` RAGGED header: header has N cell-ID fields,
  each DATA row has N+1 (unlabeled leading gene-id). Take cell-id order from
  pandas' post-`index_col=0` column names; a naive `header[1:]` causes an
  OBS OFF-BY-ONE (caught + fixed here). Verified alignment with a known-cell probe.
- Do NOT dense-read the genes x cells TSV (>50 GB RSS OOM) and do NOT use
  `np.fromstring` at full width (heap corruption / SIGABRT). Use a pandas CHUNKED
  sparse reader (peak RSS ~11 GB).
- Long foreground shells get SIGTERM'd by the harness (~100s); launch heavy runs
  via `setsid` (fully detached) and avoid zsh globs in the launch line (nomatch
  aborts the command). conda-run buffers child stdout until exit.
- Embeddings are persisted (`embeddings/method_embeddings.h5ad` + per-method
  `.npy` + `obs.csv`) so the next phase can apply LABEL-AGNOSTIC cluster-
  separability / anti-over-merging metrics (the discovery requirement) WITHOUT
  re-running scVI. A `rescore_from_embeddings.py` helper re-scores in minutes.

## Verification

- Full run `FULL_EXIT=0`; scoreboard + audit md/JSON written; rescore (data-derived
  T) `RESCORE_EXIT=0`. native-vs-scib agreement near-perfect (ASW abs_diff ~0).
- All compute under explicit timeouts in sc_gpu (scib-metrics 0.5.9 / chex 0.1.91
  / plottable 0.1.5 pins asserted; jax CPU-only — LISI/kBET CPU fallback correct).
