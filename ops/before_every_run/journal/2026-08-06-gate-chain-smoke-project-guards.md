# Gate-chain smoke run: project-guards dual-gate first full-chain exercise — 2026-08-06

## Objective

Validate the newly installed project-guards chain (2026-08-06) end-to-end on a real
but minimal pipeline run: raw-data immutability respected, structural + vision QC
gates applied to run outputs before acceptance. This was a tooling-validation run,
not science.

## What was attempted

- Input: 3,000-cell random subsample (seed 41) of the canonical
  `data/raw/wch_lung_cancer_atlas/prepared_input.h5ad`, written to
  `/tmp/gate-chain-smoke/input_3k.h5ad` (raw tree untouched — immutability guard
  respected by construction).
- Command (execution_mode = debug_massive, minimal standard lane):
  `sc_gpu/bin/python -m workflow.modular.cli --project gate-chain-smoke
  --project-root /tmp/gate-chain-smoke --input-h5ad .../input_3k.h5ad
  --scale-mode standard --localcores 8 --localmem 16 --allow-dirty`
- First attempt refused by the repo dirty-guard (uncommitted hook/skill tooling
  changes from the same day's work) — provenance guard works as designed;
  `--allow-dirty` used to record the diff instead of committing.

## What succeeded

- Run `runs/2026-08-06T0831Z-d2fc4a3` (project_root /tmp/gate-chain-smoke):
  7/7 modules ok (cellranger, qc, ambient_correction [skipped, no triggers],
  doublet_detection, clustering, annotation, differential_expression).
- GPU clustering failed with CUSOLVER_STATUS_INTERNAL_ERROR and the policy
  fallback (restore-cpu) engaged cleanly — truthful fallback recorded in log.
- **Gate 1 (h5ad-sanity-check) on final_adata.h5ad: PASS** — 2,909 × 16,808;
  0 duplicate barcodes; 0 NaN / 0 negatives in 2,000-row sample; density 8.2%;
  `layers=[counts]`, `obsm=[X_pca, X_umap]`, full uns provenance chain present.
- **Gate 2 (plot-qc-review, vision rubric) on 14 figures: PASS with NOTES** —
  12 PASS with per-panel evidence (QC → clusters → annotation → markers
  cross-consistent); 2 NOTES recorded as caveats.

## What failed / remains risky

- Interpreter segfaulted during shutdown garbage collection AFTER run ledger +
  manifest were written (torch.distributed/zarr teardown path). Artifacts intact;
  cosmetic in outcome but alarming in logs — watch for it being misread as a run
  failure in future triage.
- NOTES-1 (from vision gate): doublet zero-calls (0/2,909) vs 0.11% expected
  prior; threshold 0.013 sits just above observed score max (~0.01). Matches the
  pipeline's own DOUBLET_UNDERCALL_DIAGNOSTIC. Suspect prior/threshold
  parameterization, not biology — verify before ever citing doublet-freeness.
- NOTES-2: de_volcano y=20 drawing-floor pileup (2,806/4,500 rows) + extreme
  log2FC tail; spot-check `marker_genes_all.csv` ≥10% detection filter before
  citing specific markers.
- Batch risk flags (32 batches, no integration) are expected on this smoke input;
  all clustering/annotation outputs are EXPLORATORY by manifest.

## Lessons (promote if they recur)

1. Probe h5ad with the WRITER env: sc_gpu-written final_adata.h5ad is unreadable
   by older anndata (IOSpec encoding_type='null'); sc_rna_velocity_pseudotime
   env failed, sc_gpu env read fine. Cross-env readability is a portability NOTE.
   (Folded into h5ad-sanity-check SKILL.md same day.)
2. The degenerate-spike rubric gap (F1=1.0 ×10 repeats: determinism vs silent
   scorer) was found in the 2026-08-06 validation-figures review and the
   disambiguation rule (NOTES → check per-repeat metrics table → then PASS/FAIL)
   was added to plot-qc-review the same day.
3. `--allow-dirty` is the sanctioned path for tooling-only dirty trees; do not
   commit just to satisfy the guard.

## Artifact classification

- `/tmp/gate-chain-smoke/runs/2026-08-06T0831Z-d2fc4a3/`: **evidence-only**
  (gate-chain validation; not a canonical science run). Cleanup candidate after
  harvest; /tmp auto-cleans. This journal entry is the durable record.
- `data/raw/wch_lung_cancer_atlas/prepared_input.h5ad`: canonical (untouched).

## What the next operator should remember

- Before accepting any run's figures: `/skill:plot-qc-review`; before accepting
  any .h5ad: `/skill:h5ad-sanity-check` (writer env). Both gates are project
  skills under `.kimi-code/skills/`; hooks enforce raw immutability at the shell.
