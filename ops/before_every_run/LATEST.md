# CURRENT STATUS BANNER — 2026-08-05 (read first)

**Yost / NG2025 / LKB1 research line is RETIRED.** Do not execute historical
Yost renderers, reopen deleted project roots, or treat entries below as live
work instructions.

- Retired project roots (absent by design):
  - `/home/zerlinshen/projects/yost-stk11-lkb1-storyline`  <!-- [retired-path-doc-ok] -->
  - `/home/zerlinshen/projects/ng2025-3d-genome`  <!-- [retired-path-doc-ok] -->
  - `/home/zerlinshen/projects/lkb1-frontend-rehearsal`  <!-- [retired-path-doc-ok] -->
- Evidence capsule (historical only):
  `/home/zerlinshen/projects/pipeline-publication-baselines/runs/retired_yost_ng2025_lkb1_v1/`
  governed by suite `contracts/retired_research_lines.yaml` (`sealed_verified`).
- Active figure-parity registry (Trevino-only, 8 entries):
  `governance/figure_parity_references_2026-08-04.json` under the suite root.
- Forbidden plotting route IDs must not reappear in the live plotting registry
  (see `forbidden_active_route_ids` in the retirement contract).

Sections below this banner are **append-only historical journal**. Paths that
mention deleted figure packages or `render_yost_*` scripts are archival and
must not be re-run.

## Active Wave-9 memory — LUCA patient LIANA / E-MTAB-13530

- Canonical W9.2 run:
  `/home/zerlinshen/projects/pipeline-scientific-audit-20260805/runs/2026-08-05T1345Z-4e5392e/`.
  It completed 8 LUAD + 8 LUSC donors (16/16, 0 errors) in 69.33 s with
  1.79 GB maximum RSS; strict `--resume` replay completed in 1.51 s.
- **Regression of previous issue:** `scanpy.read_h5ad(..., backed="r")`
  materialized the 892k-cell LUCA atlas-wide `obsp` graphs. The kernel killed
  the first W9.2 process at 91.7 GB RSS. The fixed driver reads only required
  HDF5 observation columns and selected CSR expression rows, processes donors
  sequentially, and checkpoints each donor. Preserve this sparse/staged route.
- W9.1 is complete (unblocked the same day). The zero-payload block lifted:
  a content-first HTTPS re-probe returned real bytes (HDF5 magic at offset 0,
  real tar header), and all 32 E-MTAB-13530 payloads were staged and passed
  content validation (size + magic/members + non-zero byte; SHA-256 in
  `DOWNLOAD_VALIDATION_20260805_HTTPS_STAGED.json`). Canonical W9.1 run:
  `/home/zerlinshen/projects/pipeline-scientific-audit-20260805/runs/2026-08-05T1459Z-4e5392e/`
  — 16/16 sections, 0 errors, 75.81 s wall, ~2.3 GB max RSS; paired
  tumour-vs-background contrast on 4 patients. F4-01/F4-02 stay `partial`;
  `figure_parity=false`. **Treat the EBI route as flaky: always validate
  content, never size.** Earlier zero-filled placeholders remain quarantined
  under `downloads/invalid_zero_placeholders_20260805/`.
- Final validation passed all 34/34 suite gates at 2026-08-05T15:21Z with the
  W9.1-complete record set (focused tests 6/6). The raw-data root records
  `staged_valid` with 16/16 sections. No commit or push was performed.
- Full entries:
  `journal/2026-08-05-wave9-luca-liana-oom-recovery.md`,
  `journal/2026-08-05-wave9-nsclc-visium-w91-unblocked.md`.

---

# Suite Optimization (env + figure-parity gate, no run) - 2026-05-29

- Doc/code optimization round; NO pipeline run executed. Evidence:
  `/home/zerlinshen/projects/pipeline-optimization-20260529/`
  (`ASSESSMENT.md`, `ROUND3_BLOCKER_SHARPENING.md`,
  `preflight/SUNDARAM_143GB_PREFLIGHT.md`, human conclusion under `ledger/`).
- ENV CHANGE (affects future remote runs): the default `r_multiomics` conda env
  now has `r-arrow=24.0.0` installed (renv-pin match) and it is declared in
  `r_multiomics_factory/envs/r_multiomics.yml`. Parquet roundtrip verified TRUE.
  Either `r_multiomics` or `r_multiomics_arrow` can now read/plot v2/v2.1/v2.2
  parquet bundles. Doc: `r_multiomics_factory/docs/ENVIRONMENT_REPRODUCIBILITY.md`.
- Suite gate is now 11 gates: figure-parity gate re-wired (gate 10,
  conditional-tolerant; output->/tmp). Active registry
  `governance/figure_parity_references_2026-05-29.json` (22 produced figures
  header/existence-validated; reference parity pending curation - NOT parity proof).
- Did NOT start the 143.6 GB Sundaram download (preflight plan only) or any
  >30 GB staging. Round-3 external blockers remain open with sharpened data asks.
- Verification: `bash scripts/run_all_gates.sh` -> 11/11 ALL GATES PASS;
  repo-doc-sync strict NO DRIFT (suite root + r_multiomics_factory);
  independent code-reviewer APPROVE + verifier VERIFIED.

---

# Pipeline Validation Real-Science Follow-Up - 2026-05-28

- Objective: finish the remaining real-data scientific validation and cleanup
  reproducible bulky intermediates.
- Evidence root:
  `/home/zerlinshen/projects/pipeline-validation-20260528/`.
- Current verdict: `PASS_SUPPORTED_NOT_FINAL_WITH_REVIEW_FLAGS`.
- Supported:
  - targeted sensitivity passes for LUSC AT1, LUSC cDC2, Trevino inhibitory
    interneuron, and Trevino intermediate progenitor;
  - LUSC origin DE has a true-count-compatible matched/dataset-aware sanity pass
    on 43 samples across 5 datasets, with permutation q<0.1 = 0 and 12/12
    sentinel directions correct;
  - GM12878 unbiased Hi-C chr19 A/B compartments pass published-subcompartment
    concordance (`0.814`, sign-fixed).
- Still review/conditional:
  - LUSC DC mature targeted sensitivity;
  - LUSC tumor-stage DE, because the usable contrast is single-dataset only;
  - factory TAD boundaries, because B35T1NC Micro-C author-TAD concordance is
    below null (`F1=0.049`, `recall_over_null=0.40` best tested).
- Code/docs updated:
  `workflow/modular/modules/hic_tad.py`,
  `tests/test_wave2b_hic_smoke.py`, `README.md`,
  `docs/MULTIOMICS_MODULE_RATIONALE.md`, and `AGENTS.md`.
- Verification:
  `python3 -m py_compile` pass; `jq empty` pass on all validation summary JSONs;
  no `NaN`/`Infinity` in summary JSON; `pytest --no-cov -q
  tests/test_wave2b_hic_smoke.py` = 14 passed; `git diff --check` pass.
- Cleanup:
  generated contact TSVs and Trevino `.checkpoints/*.h5ad` were deleted after
  compact evidence was retained. Validation root is 1.8M and `/home/zerlinshen`
  has 168G available.
- Carry-forward:
  do not claim final TAD biology until a stronger TAD caller is integrated and
  validated on the retained B35T1NC/cooltools-style ground truth.
- Full journal:
  `journal/2026-05-28-pipeline-validation-real-science-followup.md`.

---

# Pipeline Validation Post-Review Reconciliation - 2026-05-28

- Objective: fix the G008 review blockers in the integration-biology/multiomics
  validation round without weakening the real-data requirement.
- Evidence root:
  `/home/zerlinshen/projects/pipeline-validation-20260527/`.
- Scripts rerun:
  - `scripts/run_lusc_g002_g003_marker_mixing.py`
  - `scripts/run_trevino_g002_g003_g004_marker_mixing.py`
  - `scripts/run_lusc_g005_annotation_de_sanity.py`
  - `scripts/run_g006_hicstraw_factory_bridge.py`
- Current results:
  - G002 marker retention: LUSC 23/23 and Trevino 9/9 marker panels pass; all
    marker label-shuffle negative controls pass.
  - G003/G004: embedding negative controls fire on LUSC 24/24 and Trevino 9/9
    labels, but purity/rare-population flags remain and must not be hidden.
  - G005: pseudobulk DE no longer falls back to normalized/log expression. It
    uses only 71/87 raw-count-compatible samples after excluding 16 fractional
    count samples from `Guo_Zhang_2018` and `Maynard_Bivona_2020`; origin and
    tumor-stage contrasts remain dataset-confounded.
  - G006: the real LUSC H3K27ac HiChIP `.hic` bridge now requires and passes all
    7 Python/R technical gates; R-side `load_hic_extension` is `ok`.
  - G008: independent code review approved the fixes; architect/science status is
    `WATCH`; final registered verdict is `PASS_SUPPORTED_NOT_FINAL`.
- Carry-forward:
  - Do not round fractional-count samples for DE. Exclude or recover true raw
    counts.
  - Do not upgrade G003/G004 to clean passes; preserve flagged labels.
  - Do not claim unbiased 3D-genome biology from the H3K27ac HiChIP chr21 bridge.
  - Downstream summaries must use `PASS_SUPPORTED_NOT_FINAL`, not a clean global
    pass.

---

# Q22 v2 — Production change re-run (extreme-theta in `_compute_embeddings`) — 2026-05-27

- Trigger: lead added the extreme-theta over-corrector control to
  `integration_select.py::_compute_embeddings` (REPORTED-NOT-GATING per §D3).
  Re-ran the gate to verify no regression + produce a fresh design-consistent artifact.
- STEP 1 tests: `tests/test_integration_select_wiring.py` +
  `tests/test_discovery_integration_gate.py` -> **45 passed, 0 failed**. No test
  updates needed (wiring tests monkeypatch `_compute_embeddings`; gate tests cover
  REPORTED-NOT-GATING semantics).
- STEP 2 fresh artifact: `runs/lusc_dataset_axis_gate_v2/` (new cache key
  `1529e328…` != prior `b26ad3ac…` => fresh recompute, no replay). 92,430 cells,
  scVI 100 ep / seeds (0,1,2) / GPU. exit 0, ~22 min.
- VERDICT (design-consistent, task STEP 3): **Q22 PASSES.** (a) recommends REAL
  integration `harmony` (floors mix 0.1038>0.0520, dist 0.5884>=0.5410, iso
  0.4827>=0.4785; scVI band fails mixing+isolation). (b) the REQUIRED SHUFFLE
  anchor is disqualified by BOTH floors (dist 0.4782<0.5410 AND iso 0.4693<0.4785),
  fired=True — falsifiability proof holds.
- extreme-theta (REPORTED-NOT-GATING) — TWO calls DISAGREE this run:
  PRODUCTION `_compute_embeddings` call went SINGULAR at theta=100
  (`_LinAlgError`, logged in `integration_audit.json.candidates_failed`) =>
  production JSON `extreme_theta_control_present:FALSE` (expected flip to TRUE did
  NOT materialize, but the change DID run — failure is the known theta=100
  singularity at scale, 3rd confirmation). The driver's SEPARATE re-score
  SUCCEEDED => `overcorrector_audit.json` records it present (mix 0.0288, dist
  0.5659, iso 0.4890; survives dist+iso, fails mixing; fired=False). Either way
  NON-GATING; the shuffle anchor carries the Q22 (b) proof.
- Driver's OWN printed "Q22 PASS: False" keys (b) on the extreme-theta, not the
  shuffle — that is the wrong anchor per the design-consistent criterion;
  superseded by the shuffle-anchor verdict above. NO code committed.

---

# Q22 — LUSC Integration-Selection Gate (REAL, GPU) — 2026-05-27

- Task Q22: per-run discovery integration gate on a REAL strong-batch LUSC cohort,
  GPU path (P1 `bind_cuda_context` fix). Branch `wave6-trevino-v5.1`, env `sc_gpu`.
  Full journal: `journal/2026-05-27-q22-lusc-integration-gate-gpu.md`. NO code committed.
- Cohort: LuCA core 892k -> squamous subset (MONDO:0005097) 92,430 cells, batch_key
  `dataset` (9 datasets, cross-study + cross-platform). NOT subsampled. Outputs under
  `/home/zerlinshen/projects/lusc-integration-gate-20260527/runs/lusc_dataset_axis_gate/`.
- VERDICT: **Q22 PASSES** (both conditions). (a) gate recommends a REAL integration
  = `harmony` (cleared all 3 floors: mixing 0.104>0.052, distinctness 0.591>=0.542,
  isolation 0.489>=0.480). scVI did NOT survive (under-mixed at 100 epochs:
  mixing 0.020<floor, isolation 0.474<floor). (b) the SHUFFLE over-merge anchor
  FIRED — distinctness 0.478<0.542 AND isolation 0.469<0.480 (disqualified);
  scoreboard shows its over-merger signature (mixing 0.530/kBET 0.921 but LOWEST
  bio). chosen=harmony, backend=`direct`, X_pca_harmony [92430,15], converged 11/50.
- P1 GPU FIX CONFIRMED on real 92k data: pre-rapids cuBLAS/cuSOLVER warm held; GPU
  rsc PCA (6.1s, no CUSOLVER error) + scVI 3x100ep on GPU torch (no CUBLAS error)
  + downstream all in ONE process. Total wall ~22 min.
- CAVEAT (regression, same family as 2026-05-26): extreme-theta=100 over-corrector
  candidate went SINGULAR (`linalg.inv` zero diagonal) — could not be materialized.
  Known numerical limit at scale; gate treats extreme-theta as REPORTED-NOT-GATING,
  shuffle is the robust anchor and fired. Next time use theta 20-50 or ridge.
- Lesson: the gate's over-correction detector is BIO-CONSERVATION-based
  (distinctness/isolation floors), not mixing — shuffle has max mixing yet is
  disqualified. Do NOT chase extreme-theta to "pass" Q22.

---

# Integration Calibration Bench — Trevino RNA Atlas — 2026-05-26

- Task T4 (integration-bench team, worker-2): CALIBRATION/RELATIVE lane on the real
  Trevino GSE162170 RNA atlas (57,868 cells x 33,355 genes). Outputs under
  `/home/zerlinshen/projects/integration-bench-20260526/` (never in factory tree).
  Full journal: `journal/2026-05-26-integration-bench-trevino-calibration.md`.
- 6 embeddings persisted to `embeddings/` (baseline, scVI seed 0/1/2, harmonypy-direct
  Harmony, shuffle-label control) for next-phase label-agnostic separability metrics
  WITHOUT re-running scVI. Scoreboard + audit md/JSON written; `FULL_EXIT=0`.
- Gates: G5 two-part guard PASS (8 samples / 3 technical batches; single-batch multiome
  subset hard-fails part (b)); G2 vs author seurat_clusters 23cl (cluster_names join
  refused — dead symlink); G4 per-age (w16/w20 clean, w21+w24 confounded on b2020_02);
  G6 calibration_relative_only + relative_ranking_only, no winner.
- CRITICAL finding (logged to `.omc/plans/open-questions.md`): production DEFAULT Harmony
  backend `BatchCorrectionModule._run_harmony` is BROKEN in sc_gpu on BOTH paths (rapids
  GPU cupy CUBLAS_NOT_INITIALIZED + context corruption; scanpy1.12<->harmonypy0.2.0
  wrapper shape bug). Bench used harmonypy.run_harmony DIRECT (canonical algo; scanpy
  wrapper bypassed). Fix = separate production PR.
- Neg-control: shuffle-label is the robust over-correction anchor (ARI 0.486->0.0002);
  Harmony-extreme-theta singular at theta=100 on full data. Data-derived aggregate
  T=0.346 -> 0 flags on real methods; structure-sensitive detector (ARI/NMI) control
  fires (PASS). Over-correction metric SHOULD be cluster-separability based (discovery
  requirement) — re-scoped next phase.
- Cautions: ragged-header off-by-one trap (use pandas post-index_col columns, not
  header[1:]); never dense-read or np.fromstring the genes x cells TSV (OOM / SIGABRT —
  use pandas chunked sparse, ~11 GB); launch heavy runs via setsid (harness SIGTERMs
  long foreground shells ~100s; no zsh globs in launch line).

---

# Factory Governance Real-Data Validation - 2026-05-25

## Scope

- Task: G004 of the factory-governance ultragoal, responding to the critique
  that current gates were mostly structural and should not be treated as
  scientific validation.
- Execution mode: controller_validation with bounded real-data proof.
- Host: `/home/zerlinshen` on `ubuntu-tail`.

## What Changed

- Recorded failed HGMM exploratory rerun scaffold:
  `/home/zerlinshen/projects/hgmm-smoke/runs/2026-05-24T1930Z-0321773`.
  The process hung after scaffold creation and was terminated; the run now has a
  README marking it as failed exploratory evidence only.
- Created governed Round9 validation run:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-24T2347Z-0321773`.
- Exported a v2.1 compact R bundle from retained real Round9 LUSC consensus
  AnnData using current `singlecell_factory` exporter.
- Rendered `7` PNG figures plus `1` PDF with `r_multiomics_factory` in the
  `r_multiomics_arrow` environment through the plotting-factory bridge.
- Added project evidence and human-review records for the bounded validation
  claim.

## Verification

- Export and R plotting both completed under explicit timeouts.
- Bundle manifest records `5000` real cells, `X_umap`, `X_pca`, and `17`
  marker genes.
- PNG QA passed: `7/7` nonblank rendered PNGs.
- `python scripts/validate_run_bundle.py` passed with warnings only.
- `bash scripts/run_all_gates.sh` passed `8/8`.

## Caveats

- This is cross-factory real-data handoff validation, not a full raw-input
  pipeline reproduction and not NG2025 end-to-end proof.
- Default `r_multiomics` lacks R package `arrow`; `r_multiomics_arrow` works.
  Canonical environment repair remains open.
- New journal:
  `before-every-run/journal/2026-05-25-factory-governance-realdata-validation.md`.

---

# Remote Project Retention Cleanup And Suite Gitification - 2026-05-24

## Scope

- Task: execute the approved Ralph plan combining project retention cleanup and
  pipeline/suite gitification.
- Execution mode: controller_validation for governance/retention only.
- Remote host: `/home/zerlinshen` on `ubuntu-tail`.

## What Changed

- Created project cleanup manifest:
  `/home/zerlinshen/projects/ledger/cleanup_candidates/2026-05-24-project-retention-gitification.tsv`.
- Created missing retention policies for:
  - `/home/zerlinshen/projects/round9-singlecell-comparison/ledger/project_retention_policy.yaml`
  - `/home/zerlinshen/projects/hgmm-smoke/ledger/project_retention_policy.yaml`
- Preserved metadata/evidence records under per-project
  `ledger/cleanup_records/2026-05-24-project-retention-gitification/`.
- Deleted only explicit superseded candidates with preservation records:
  - `wave5-trevino/runs/20260516T0931Z-d192836f1bb0`
  - `wave5-trevino/runs/20260517T1005Z-7b539a5`
  - `hgmm-smoke/runs/2026-05-21T1625Z-9907a7d`
  - `hgmm-smoke/runs/2026-05-21T1626Z-9907a7d`
  - `hgmm-smoke/runs/2026-05-21T1630Z-9907a7d`
- Reclaimed approximately 6.07 GB; `/home/zerlinshen/projects` dropped from
  about 104G to about 98G.
- Updated suite root docs/ignore for a lightweight governance git repository.
- Remediated retired paths in `/home/zerlinshen/projects-bootstrap` and smoke
  tested `omc-new-project`.

## Preserved / Deferred

- Preserved Wave5 source-of-truth run `20260517T1436Z-13c2c88` and linked
  pipeline run `2026-05-17T2004Z-13c2c88`.
- Preserved Round9 comparison lanes and initial evidence-only baseline.
- Preserved NG2025 data/open/staging and source-of-truth run family.
- Deferred large Wave5 2026-05-19 runs pending a dedicated Cell cleanup
  inventory.

## Verification

- Retention checks passed with manifest coverage, source-of-truth existence,
  preservation records, and non-git project/data root checks.
- Suite gate and git commit/push status are recorded in the active Codex
  session final report.

---

# Latest Yost STK11/LKB1 figure cleanup - 2026-05-23

- User requested preserving only the final 14 Yost/STK11 figures and deleting
  outdated images locally and remotely.
- Retained canonical remote package:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.
- Retained Mac archive:
  `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.
- Deleted remote generated figure packages:
  `3d_lung_author_style_single_figures_20260522`,
  `3d_lung_stk11_ppt_20260521`,
  `3d_lung_stk11_single_figures_20260522`,
  `3d_lung_yost_author_code_figures_20260522`,
  `debug_individual_01`,
  `yost_2025_exact_panel_heatmaps_20260523T000054`,
  `yost_2025_exact_panel_heatmaps_20260523T000529`,
  `yost_2025_exact_panel_heatmaps_20260523T001058`, and
  `yost_stk11_storyline_single_figures_20260523T010004`.
- Deleted local generated figure packages under the Mac handoff `figures/`
  directory except the final individual package.
- Verification after cleanup: both local and remote figure-package listings show
  only the final individual package; the final package still contains 14 PNG,
  14 SVG, 14 PDF, and 14 TIFF outputs.
- Preserved boundary: raw data, prepared inputs, source code, original
  references, project ledgers, README files, and final package metadata were not
  deleted.
- Updated project policy cleanup record and local handoff `FIGURE_INDEX.md` /
  `CONCLUSION.md`.
- New journal:
  `before-every-run/journal/2026-05-23-yost-stk11-outdated-figure-cleanup.md`.

---

# Latest Yost STK11/LKB1 individual single-script figure update - 2026-05-23

- Corrected the storyline rendering rule after user review: final handoff figures
  must be generated one by one from real data/run tables, not copied from old
  rendered packages and not cropped from a larger generated figure.
- Current canonical remote storyline package:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.
- Superseded for final handoff:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_single_figures_20260523T010004`.
- Current renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_stk11_storyline_individual_figures.R`.
- Individual figure scripts:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/yost_stk11_individual/render_XX_*.R`.
- Compatibility wrapper:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_stk11_storyline_single_figures.R`
  now forwards to the individual renderer and must not reintroduce copy/crop
  behavior.
- Output: 14 standalone figures, each exported as PNG/SVG/PDF/TIFF.
  Figures 01-03 are direct Yost Fig. 1f context correspondences; 04-10 are
  STK11/LKB1 real-run extension figures; 11-13 add real NSCLC single-cell
  context from WCH and Round9; 14 is a metadata audit showing that current
  remote NSCLC resources lack audited STK11 genotype, therapy-response, and
  survival labels.
- Verification: renderer completed on `ubuntu-tail`; output counts are 14 PNG,
  14 SVG, 14 PDF, and 14 TIFF; `SCRIPT_MANIFEST.tsv` has 14/14 exit status 0;
  `PNG_QA.tsv` has 14/14 nonblank high-resolution PNGs; `FIGURE_STATUS.tsv`
  marks all rows `independent_data_render_no_copy_no_crop`; grep found no
  `file.copy` or `copy_exact_bundle` in the new individual scripts/driver.
- Current Mac handoff copy:
  `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.
- Caveats: CN and methylation remain proxy-only; STK11-mut versus WT LUAD is
  exploratory (`n=2` versus `n=2`); WCH and Round9 support cell-state/checkpoint
  context but not STK11-genotype or clinical-efficacy claims.
- New journal:
  `before-every-run/journal/2026-05-23-yost-stk11-individual-single-script-per-figure.md`.
- New human review log:
  `/home/zerlinshen/projects/ng2025-3d-genome/ledger/human_review/2026-05-23-yost-stk11-individual-single-script-per-figure-human-review.md`.

---

# Latest Yost STK11/LKB1 storyline single-figure package update - 2026-05-23

- Completed a remote-rendered STK11/LKB1 storyline package using only real
  NG2025/Yost project data and real run outputs.
- Current canonical remote storyline package:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_single_figures_20260523T010004`.
- Added R renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_stk11_storyline_single_figures.R`.
- Supporting exact-panel package:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T001058`.
- Output: 10 standalone single figures, each exported as PNG/SVG/PDF/TIFF.
  The first three figures are direct Yost Fig. 1f context correspondences; the
  remaining STK11/LKB1 storyline figures are extension outputs and are not
  labeled as original paper panels.
- Verification: renderer completed on `ubuntu-tail`; output counts are 10
  PNG, 10 SVG, 10 PDF, and 10 TIFF; R PNG QA passed for all 10 PNGs; `git diff
  --check` passed for the storyline renderer, exact-panel renderer title
  correction, and `r/report/README.md`; visual spot-check confirmed the CALDER
  title is no longer clipped and the STK11 extension panel has no non-original
  title.
- Caveats: CN and methylation remain proxy-only; STK11-mut versus WT LUAD is
  exploratory (`n=2` versus `n=2`); bulk ATAC motif context is not
  STK11-genotype-specific; external single-cell support remains a metadata-audit
  gap.
- New journal:
  `before-every-run/journal/2026-05-23-yost-stk11-storyline-single-figures.md`.
- New human review log:
  `/home/zerlinshen/projects/ng2025-3d-genome/ledger/human_review/2026-05-23-yost-stk11-storyline-single-figures-human-review.md`.

---

# Latest Yost exact-panel heatmap/contact package update - 2026-05-23

- Completed the Yost-only exact-panel heatmap/contact/track rerun requested by
  the handoff
  `/home/zerlinshen/projects/ng2025-3d-genome/NEXT_AGENT_YOST_EXACT_PANEL_REPRODUCTION.md`.
- Current canonical remote package:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T001058`.
- Added R renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_exact_panel_heatmaps.R`.
- Output status summary from `PANEL_STATUS.tsv`:
  - `EXACT_AUTHOR_CODE_DATA_AVAILABLE`: 2 panels, Fig. 1f H3K27ac 1D-signal
    and FitHiChIP loop sample-correlation heatmaps from full staged matrices.
  - `LUNG_SUBSET_EXACT`: 2 panels, LUAD/LUSC CALDER comp-rank correlation
    heatmap and chr8 HiChIP contact maps.
  - `AUTHOR_CODE_ADAPTED`: 1 MYC/PVT1 track/loop panel.
  - `STK11_EXTENSION_NOT_ORIGINAL_PANEL`: 1 separated STK11/LKB1 extension.
  - `UNSUPPORTED_INPUT_GAP`: full-cohort CALDER Fig. 1f remains blocked until
    all-cancer CALDER `all_sub_compartments.tsv` files are staged.
- Comparison sheet:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T001058/comparisons/Yost2025_NatGenet_Fig1_original_vs_exact_panel_outputs.png`.
- Verification: renderer completed without contact-map coordinate warnings after
  switching contact maps to `geom_tile`; package has 18 PNG files including
  original references, 7 SVG, 7 PDF, and 7 TIFF outputs; R PNG QA passed for all
  18 PNGs; filename check found `0` Miao/Qiaowei outputs; STK11 files are absent
  from reproduced/lung-subset panel directories and present only under
  `04_STK11_LKB1_extension`; `git diff --check` passed for the exact-panel
  renderer.
- New journal:
  `before-every-run/journal/2026-05-23-yost-exact-panel-heatmaps.md`.
- New human review log:
  `/home/zerlinshen/projects/reproductions/ledger/human_review/2026-05-23-yost-exact-panel-heatmaps-human-review.md`.

---

# Latest Yost author-code-informed figure package update - 2026-05-22

- Replaced the mixed lung 3D-genome figure-generation lane with a Yost-only
  author-code-informed renderer.
- Current canonical remote package:
  `/home/zerlinshen/projects/reproductions/figure_packages/3d_lung_yost_author_code_figures_20260522`.
- Current Mac handoff copy:
  `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/3d_lung_yost_author_code_figures_20260522`.
- Added R renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_author_style_figures.R`.
- Output: 6 Yost-only records, each exported as PNG/SVG/PDF/TIFF:
  Yost Fig. 1-Fig. 5-style reproduction figures plus separated STK11/LKB1
  exploratory extension.
- Governance: `plotting_factory/contracts/lung_3d_author_style_profiles.json`
  now treats Yost as the only current generated author-code route; Miao Liu and
  Qiaowei Liu are record-only/source-memory routes unless explicitly reopened.
- Verification: R renderer completed; package has 6 PNG/SVG/PDF/TIFF files;
  filename check found no Qiaowei or Miao outputs; author-style profile
  validator passed with `OK profiles=3 groups=2 record_only=2`; runtime check
  found Yost `Rscript` and `python3`; `git diff --check` passed for renderer,
  contract, README, and AGENTS files; manual visual QA passed after removing
  the stray ggplot text-legend glyph and edge-label clipping.
- New journal:
  `before-every-run/journal/2026-05-22-yost-author-code-figure-package.md`.
- New human review log:
  `/home/zerlinshen/projects/reproductions/ledger/human_review/2026-05-22-yost-author-code-figure-package-human-review.md`.

---

# Latest lung 3D single-figure package update - 2026-05-22

- Completed a plotting/reporting-only redesign of the lung 3D-genome
  reproduction handoff into standalone publication-style figures.
- Current canonical project-side figure package:
  `/home/zerlinshen/projects/reproductions/figure_packages/3d_lung_stk11_single_figures_20260522`.
- Current Mac handoff copy:
  `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/3d_lung_stk11_single_figures_20260522`.
- Added plotting renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/python/report/render_lung_3d_single_figures.py`.
- Added plotting renderer README:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/python/report/README.md`.
- Output: 12 standalone figure records, each exported as PNG/SVG/PDF; package
  includes `FIGURE_MAP.tsv` and `README.md`.
- Scientific conclusion:
  - Yost 2025 Nat Genet: LUAD/LUSC processed MVS evidence supports loop
    landscape, CALDER2, lineage accessibility, gene-loop/RNA, CNV-aware axes;
    strict cell-type scATAC remains not reproduced.
  - Miao Liu 2025 Nat Genet: true single-cell 3D genome track; processed
    4DN/source-data supports progression; STK11 is a boundary/negative result.
  - Qiaowei Liu 2025 Nat Commun: TP63-MYC processed/proxy loop mechanism
    companion; not single-cell and not STK11-centric.
- Verification: `py_compile` passed, renderer wrote 12 figures, PNG QA passed
  for 12/12 nonblank high-resolution figures, and `git diff --check` passed for
  the plotting renderer.
- New journal:
  `before-every-run/journal/2026-05-22-lung-3d-single-figure-redesign.md`.
- New human review log:
  `/home/zerlinshen/projects/reproductions/ledger/human_review/2026-05-22-lung-3d-single-figure-human-review.md`.

---

# Latest governance update: human conclusion logs required - 2026-05-19

- Round9 LUSC closure completed on 2026-05-22 with three comparable production
  lanes under `/home/zerlinshen/projects/round9-singlecell-comparison`.
- Fresh runs:
  - baseline refresh:
    `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1744Z-9907a7d/python/lusc_ps01_round9_baseline_refresh_20260522_014447`
  - scDblFinder:
    `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1749Z-9907a7d/python/lusc_ps01_round9_scdblfinder_20260522_014947`
  - consensus OR:
    `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1756Z-9907a7d/python/lusc_ps01_round9_consensus_scdbl_or_20260522_015643`
- Outcome: all three lanes had `14` modules `ok`; `pseudobulk_de` skipped due
  missing explicit contrast. Baseline Scrublet still under-called (`21`
  doublets, `0.03%`, undercall ratio `206.9x`), while scDblFinder and
  consensus called `6737`/`6751` doublets (`9.30%`/`9.32%`).
- Decision: keep Scrublet as global default, but when the under-call diagnostic
  fires on heterogeneous tumor/tissue data, use scDblFinder or consensus OR
  with `scrublet_scdblfinder` as a conditional second-opinion lane.
- Downstream stability: consensus retained `65669` cells, produced `35`
  clusters, `0.0%` annotation unknown rate, and `10488` significant DE genes.
- Reports:
  - comparison:
    `/home/zerlinshen/projects/round9-singlecell-comparison/reports/round9_lusc_doublet_downstream_comparison.md`
  - curated vector figure bundle:
    `/home/zerlinshen/projects/round9-singlecell-comparison/reports/figures/round9_lusc_closure`
  - human conclusion:
    `/home/zerlinshen/projects/round9-singlecell-comparison/ledger/human_review/2026-05-22-round9-lusc-closure-human-conclusion.md`
- Figure QA: curated bundle passed `--require-vector` with `9` valid figures,
  `6` vector outputs, `3` rasters, and `0` warnings. Manual visual review is
  still required before publication use.
- New journal: `before-every-run/journal/2026-05-22-round9-lusc-closure.md`.

---

- Round9 LUSC baseline evidence run completed on 2026-05-22:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1718Z-9907a7d/python/lusc_ps01_round9_baseline_20260522_011822`.
- Purpose: second real data shape for full single-cell comparative optimization,
  with upstream and downstream modules plus Figure QA.
- Outcome: run exit 0; `14` modules `ok`; `pseudobulk_de` truthfully skipped
  because no explicit contrast was provided; final AnnData `72399 x 17267`.
- Key caution: Scrublet under-called doublets strongly (`0.03%`, `21` calls,
  undercall ratio `206.9x` vs expected `6%`). Do not promote a doublet default
  without a second-opinion/consensus lane on this data shape.
- Figure QA: `30` mechanically valid PNGs, `0` vector outputs; manuscript
  vector gate fails as expected with `missing_vector_output`.
- New project governance: `/home/zerlinshen/projects/round9-singlecell-comparison/project.yaml`.
- New journal:
  `before-every-run/journal/2026-05-22-round9-lusc-baseline.md`.

---

- Confirmed existing `before-every-run` rule already required remote run journaling, but did not explicitly require human-facing scientific decision logs.
- Added canonical rule: when a discussion changes scientific interpretation, final-run strategy, claim support, benchmark lane choice, source-of-truth status, or human-facing next actions, write a human-readable conclusion log.
- Required locations now documented:
  - `/home/zerlinshen/projects/<project-id>/ledger/human_review/<date>-<topic>.md`
  - `/home/zerlinshen/projects/<project-id>/runs/<run-id>/evidence/<topic>.md` when tied to a run
  - `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/LATEST.md` plus journal entry
- Updated suite and repo governance surfaces:
  - `/home/zerlinshen/Bioinformatics Research Pipeline/README.md`
  - `/home/zerlinshen/Bioinformatics Research Pipeline/governance/README.md`
  - `singlecell_factory/README.md`
  - `singlecell_factory/docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md`
  - `r_multiomics_factory/README.md`
  - `plotting_factory/README.md`
  - remote and local `before-every-run` skill files
- Verification: `git diff --check` passed for changed child-repo markdown/skill files.
- New journal: `before-every-run/journal/2026-05-19-governance-human-conclusion-log-rule.md`.

---

# Latest human conclusion log update - 2026-05-19

- Added human-readable conclusion log for the NC2024/Cell controller-validation discussion and final clustering strategy.
- Key decision: NC current CSS route is evidence-only/controller-validation; final real-data claims should use sparse-exact or validated hybrid clustering, with CSS retained only as smoke/control.
- Key Cell conclusion: low ARI is likely driven by reproduction-boundary and missing original annotation/parameter/multi-omic context; biological signal is not absent and needs a parity/gap report.
- Human conclusion artifacts:
  - `/home/zerlinshen/projects/nc-reproduction/ledger/human_review/2026-05-19-clustering-final-strategy-human-conclusion.md`
  - `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-19T1010Z-82f1964/evidence/human_conclusion_clustering_final_strategy_20260519.md`
  - `/home/zerlinshen/projects/wave5-trevino/ledger/human_review/2026-05-19-clustering-final-strategy-human-conclusion.md`
  - `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-19T1005Z-82f1964/evidence/human_conclusions/human_conclusion_clustering_final_strategy_20260519.md`
- New journal: `before-every-run/journal/2026-05-19-human-conclusion-clustering-final-strategy.md`.

---

# Latest two-dataset full simulation and DE/checkpoint patch update - 2026-05-19

- Cell/Trevino full public RNA controller-validation completed: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-19T1005Z-82f1964`.
- NC2024 full tumor run reached main DE evidence, then was intentionally stopped at status 143 to avoid a redundant old-code DE retry: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-19T1010Z-82f1964`.
- RAPIDS DE cause confirmed: installed `rapids-singlecell 0.13.4` does not expose `rapids_singlecell.tl.rank_genes_groups`; GPU clustering remains valid and separate.
- Main NC DE output existed before stop: `marker_genes_all.csv` 1901 lines, `marker_genes.csv` 1893 lines, `marker_top5_by_cluster.csv` 96 lines.
- Patch applied in `workflow/modular/context.py`, `workflow/modular/modules/differential_expression.py`, tests, and `README.md`: RAPIDS DE API/version gate, massive sparse CPU DE routing, correction-aware fallback, no MemoryGuard retry/false-ok, oversized substate DE skip, Dataset2D checkpoint sidecar fallback, and operator docs for sidecar-only checkpoints.
- Verification: `py_compile` passed; `git diff --check` passed; targeted pytest `9 passed`; real environment smoke confirmed `de_backend=cpu_sparse`, `de_test_actually_used=sparse_welch_fallback`, `de_correction_actually_used=bonferroni`, `de_rapids_singlecell_version=0.13.4`.
- New journal: `before-every-run/journal/2026-05-19-two-dataset-full-simulation-and-de-patch.md`.

---
# Latest Cell/Trevino module-smoke update - 2026-05-19

- Remote Cell/Trevino module smoke completed successfully from local Codex over `ssh ubuntu-tail`.
- Project root: `/home/zerlinshen/projects/wave5-trevino`.
- Evidence-only module-smoke run: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-19T0956Z-82f1964`.
- Derived smoke input: `/home/zerlinshen/projects/wave5-trevino/inputs/module_smoke_20260519_cell_rna_3000/prepared_input.h5ad`, a 3,000-cell subset of the public Trevino RNA prepared matrix with all genes retained.
- Configs used: `/home/zerlinshen/projects/wave5-trevino/configs/module_smoke_20260519/cell_brain_development_markers.json` and `cell_brain_development_signatures.json`.
- Modules all `ok`: `cellranger`, `qc`, `doublet_detection`, `clustering`, `annotation`, `cell_cycle`, `differential_expression`, `gene_signature_scoring`, `metacell`, `trajectory`, `composition`, `pathway_analysis`, `pseudobulk_de`, `cell_fate`.
- Key metrics: raw `3000 x 33355`; after QC `3000 x 22556`; after doublet removal `2997` cells; `19` clusters; annotation unknown pct `0.0`; `2279` DE genes; `15` signatures scored; `50` metacells; `8` sample groups in composition; exploratory pseudobulk completed; `5` terminal states in cell-fate fallback.
- Expected fallbacks observed: CPU scrublet fallback after RAPIDS dtype rejection, MiniBatchKMeans fallback for missing `SEACells`, composition fallback for missing `pertpy`, pathway fallback for missing `gseapy/decoupler`, pseudobulk Mann-Whitney fallback for missing `pydeseq2`, manual cell-fate fallback for missing `CellRank`.
- Governance validation after the run: overall `pass`; warnings are existing advisory/legacy Cell project fields and legacy run layout warnings, with no findings on the new module-smoke run.
- Retention: do not delete existing Cell/Trevino runs in this task. Treat the new run and derived input as `evidence-only` module-interface smoke artifacts, not the Cell/Trevino human-facing source of truth.
- Current Cell/Trevino scientific source of truth remains `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88`; current linked pipeline run remains `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88`.
- New journal: `before-every-run/journal/2026-05-19-cell-trevino-module-smoke.md`.

---

# Latest raw FASTQ/FASTA-reference local-to-remote smoke update - 2026-05-19

- Local Mac session successfully launched remote `singlecell_factory` over `ssh ubuntu-tail`.
- Project root: `/home/zerlinshen/projects/raw-fasta-smoke`.
- Canonical latest/best smoke run: `/home/zerlinshen/projects/raw-fasta-smoke/runs/2026-05-19T0931Z-82f1964`.
- Input boundary: real 10x tiny FASTQ plus Cell Ranger GRCh38 reference containing `fasta/genome.fa`; operational raw/reference smoke only, not article-level scientific truth.
- Clean run modules: `cellranger`, `qc`, `doublet_detection`, `clustering` all `ok`; raw cells `1142`, raw genes `38606`, after QC `1142 x 377`, clusters `9`.
- First attempt `2026-05-19T0928Z-82f1964` proved raw/QC/clustering but optional annotation failed due insufficient marker-gene coverage; parameter/log/manifest evidence archived under `/home/zerlinshen/projects/raw-fasta-smoke/ledger/run_records/2026-05-19T0928Z-82f1964`.
- Retention action applied: deleted superseded bulky first-run output while preserving per-run logs/params; protected raw FASTQ source and reference/FASTA.
- New journal: `before-every-run/journal/2026-05-19-raw-fasta-smoke-local-remote.md`.

---

# Latest two real-dataset final factory validation update — 2026-05-18

- Final bridge/governance validation completed for two real datasets across `singlecell_factory -> r_multiomics_factory -> plotting_factory`.
- Verdict: `PASS_TWO_REAL_DATASET_FACTORY_BRIDGE_VALIDATED`.
- Governance report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-two-real-dataset-final-validation/REPORT.md`.
- Machine report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-two-real-dataset-final-validation/two_real_dataset_final_validation.json`.
- NC2024 evidence: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`, final AnnData `5281 x 19504`, 8 modules `ok`, `batch_correction` expected `skipped`, project-root bundle and R outputs present. Scope remains P15_T1 structure/module/bridge validation, not full article-scale cohort truth.
- Cell/Trevino bridge evidence: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88`, final AnnData `55653 x 25519`, project-root bundle and R outputs present. Human-facing paper evidence remains `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88` with `conditional` public-resource quality gate.
- Fixes captured: current-project NC validator, optional `--output` when `--project-root` is provided, dataset-aware R marker selection, project-root R manifest lookup, and ggplot2 `linewidth` border styling.
- New journal: `before-every-run/journal/2026-05-18-two-real-dataset-final-validation.md`.

---

# Latest NC2024 legacy artifact cleanup update — 2026-05-18

- Off-mainline cleanup completed after the module-flow validation was recorded as `evidence-only`.
- Deleted evidence-only run: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88`.
- Deleted factory legacy/smoke residue: `/home/zerlinshen/singlecell_factory/results/test_20260516_*`, `/home/zerlinshen/singlecell_factory/results/ops`, and `/home/zerlinshen/singlecell_factory/output`.
- Preserved canonical NC2024 source of truth: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.
- Preserved module-flow validation report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-module-flow-validation/REPORT.md`.
- Preserved lightweight historical validation evidence: `/home/zerlinshen/singlecell_factory/results/small_real_validate_20260426`.
- Reclaimed bytes: `1614159948` (~1.50 GiB).
- Cleanup record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-legacy-artifact-cleanup/`.
- New journal: `before-every-run/journal/2026-05-18-nc2024-legacy-artifact-cleanup.md`.

---

# Latest NC2024 module flow validation update — 2026-05-18

- Controlled project-root validation run completed under `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88`.
- Verdict: `MODULE_FLOW_VALIDATION_PASS_WITH_TEST_CONTRACT_DRIFT`.
- Module status: 8 `ok`, 1 expected `skipped` (`batch_correction`, one `sample` batch only). Final AnnData `5281 x 19504`.
- Flow details: marker DB loaded 5 DBs / 198 entries; GPU clustering completed; context-aware annotation had 0 validation mismatches; DE used expected CPU fallback because current `rapids_singlecell.tl` lacks `rank_genes_groups`; doublet detection used fallback all-singlets after RAPIDS scrublet dtype rejection.
- Project governance validator passed for the new run: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-module-flow-validation/`.
- Module/process tests: `verify_module_references` OK (`45/45 modules verified; 2 project-local`), DAG/project-root/governance tests `29 passed`, minimal modular pipeline test `1 passed`.
- Known test-contract drift: broad targeted pytest including all `tests/test_modular.py` produced `78 passed`, `5 failed`; failures are stale DAG-size/GPU-fallback/copy-object expectations versus the current expanded DAG and M2 GPU failure policy.
- Artifact classification: new run is `evidence-only`; current NC2024 source of truth remains `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.
- New journal: `before-every-run/journal/2026-05-18-nc2024-module-flow-validation.md`.

---

# Latest paper reproduction ladder policy update — 2026-05-18

- Reproduction workflow is now faithful-first: clone/stage upstream repo, pin commit/data/license, run raw-data reproduction when possible, reproduce both data objects and figures, classify parity/resource gaps, then perform module-gap analysis and context optimization.
- Governance record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-paper-reproduction-ladder-policy/`.
- Project policies for Cell/Trevino and NC2024 now include `upstream_repository`, `raw_data_reproduction`, `data_object_reproduction`, `figure_reproduction`, `module_gap_decisions`, and `context_optimization_decisions`.
- Skills were updated both locally and on `ubuntu-tail` so future agents use this ladder before adding/updating reusable modules.
- No pipeline run was launched and no data was deleted in this policy round.

---

# Latest NC2024 project retention and Cell status update — 2026-05-18

- NC2024/cancer project retention rule adopted: keep only the latest validated project run result; older bulky project run outputs can be deleted once inventory/provenance and replacement source of truth are recorded.
- Applied to `/home/zerlinshen/projects/nc-reproduction`: kept `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/`; deleted 3 older run directories totaling `142437966218` bytes.
- Cleanup record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-project-latest-run-retention/`.
- Post-cleanup project governance validation: `pass`, severity `{}`, 1 run inspected.
- Cell/Trevino reproduction status checked: quality gate decision `conditional`; 17 evidence JSON files, 31 direct evidence panels, no unexpected false checks, no missing referenced paths, no hash mismatches.
- Cell status record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-cell-reproduction-status/REPORT.md`.

---

# Latest NC2024 real run after structure-change update — 2026-05-18

- Legacy factory NC2024 `results/` cleanup completed after user approval: 48 entries, `30667532081` bytes removed; `results/` is now `27M` and has zero NC2024/nc2024 leftovers.
- True NC2024 pipeline run completed under project root: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/`.
- Run verdict: `REAL_RUN_PASS_PROJECT_ARCHITECTURE_CONFIRMED`.
- Report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-real-run-validation/REPORT.md`.
- Module status: 8 `ok`, 1 expected `skipped` (`batch_correction`, one batch only). Final AnnData `5281 x 19504`.
- Recovery note: old `launch_wave3_nc2024.sh` generated invalid non-hex run ids under the new architecture; launcher was patched to use the current factory git short SHA.
- Rule going forward: NC2024 scientific outputs live under `/home/zerlinshen/projects/nc-reproduction/runs/`; factory keeps only code/contracts/validators/docs/governance records.

---

# Latest NC2024 architecture validation update — 2026-05-18

- New architecture validation using `/home/zerlinshen/projects/nc-reproduction` completed.
- Verdict: `PROJECT_ARCHITECTURE_PASS_WITH_FACTORY_LEGACY_DEBT`.
- Report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-architecture-validation/REPORT.md`.
- Project validator: `pass`, severity `{'info': 2, 'warning': 1}`, 3 runs inspected.
- No governance files were written under the project root; validation outputs are factory control-plane only.
- Legacy debt: 48 NC2024 result paths remain under `singlecell_factory/results/`, total `30667532081` bytes. Do not create new NC2024 scientific outputs there; perform a separate retention cleanup before deleting legacy results.

---

# Latest before/after-run note — Wave-5 Trevino public Cell reproduction

## 2026-05-18 05:56 CST — Figure 2C/2D public author-resource evidence
- Mainline: Cell paper all-figure reproduction from public earliest-computable inputs; still not FASTQ-level.
- Completed this checkpoint: Figure 2C proxy heatmap and Figure 2D GA/RNA correlation scatter from author/public resources.
- Status 2C: `GENERATED_AUTHOR_RESOURCE_GENE_ACTIVITY_RNA_HEATMAP_NOT_CRE_ACCESSIBILITY_PAIR_PARITY`.
- Status 2D: `REPRODUCED_AUTHOR_RESOURCE_GA_RNA_CORRELATION_WITH_LINK_COUNT_DISCREPANCY`.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2cd_public_author_cre_gene_resources.json`.
- Outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2C_2D_cre_gene/`.
- Metrics: significant links 64,030 vs paper expected 64,878 (delta -848); heatmap genes 400; pseudobulks 1267; GA/RNA table rows 19,173; unique genes 2,739.
- Boundary: Figure 2D uses author GA_RNA_Correlations.RDS joined to public author PeakGeneLinks_Significant.RDS. Figure 2C is a proxy heatmap using author glial KNN pseudobulk ATAC gene activity and RNA matrices for genes with significant linked CREs. Exact Figure 2C CRE-accessibility pair rows are not claimed because the locally available accessibility pseudobulk RDS lacks peak names and therefore cannot be safely joined to peak.name without additional author row-order metadata.
- Next immediate action: Figure 2E GO/gene-set enrichment, then 2F-I multiome panels.

# Latest before/after-run note — Wave-5 Trevino public Cell reproduction

## 2026-05-18 05:44 CST — Figure 2B author CCA matching generated/validated
- Mainline: Cell paper all-figure reproduction from public earliest-computable inputs; still not FASTQ-level because public FASTQ/SRA raw reads are not available in the current evidence set.
- Completed this checkpoint: Figure 2B matched RNA/ATAC cluster UMAPs from author public `CCA_Matching.RDS`.
- Status: `REPRODUCED_AUTHOR_CCA_PUBLIC_MATRIX_FIGURE2B`; all validation checks pass.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2b_public_author_cca_matching.json`.
- Outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2B_cca_matching/`.
- Metrics: CCA rows 89,172; RNA cells 55,653; ATAC cells 12,000; RNA/ATAC match-cluster coverage 1.000/1.000; RNA cell-type proxy agreement 0.808.
- Truth boundary: Author public CCA_Matching.RDS was used for RNA↔ATAC nearest-neighbor links. RNA UMAP is the current singlecell_factory public RNA pipeline embedding; ATAC UMAP is the bounded public gene-activity embedding produced for Figure 3B. This is not FASTQ/fragments-level reprocessing and not exact author peak-LSI UMAP parity.
- Verification: `python -m py_compile ...` passed; `git diff --check` passed; `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> 8 passed in 0.79s.
- Next immediate action: Figure 2C/2D CRE-gene linkage heatmap/correlation evidence using public/author peak-gene resources.

# Latest Wave-5 Figure 1E quality gate update (2026-05-18)

- Quality gate PASS after Figure 1E addition: Python syntax checks, `git diff --check`, and Trevino loader regression tests passed.
- Test evidence: `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> `8 passed in 0.79s`.
- Next immediate action: Figure 2B public RNA/ATAC integration evidence.

---

# Latest Wave-5 Figure 1E multimodal overlay update (2026-05-18)

- Figure 1E public-matrix multimodal marker overlays are available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_1E_multimodal_marker_overlays`.
- Validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure1e_public_matrix_multimodal_marker_overlays.json`; all checks pass `True`.
- Metrics: RNA cells `55653`, ATAC overlay cells `12000`, markers `SOX9, EOMES, NEUROD2, DLX2`.
- Honest caveat: RNA expression, ATAC gene activity, and chromVAR motif activity were rendered from public processed matrices. ATAC panels use the 12k public ATAC subset embedding. DLX2 uses a DLX6 motif row as a DLX-family proxy because no exact DLX2 chromVAR row was found.
- Next immediate action: Figure 2B public RNA/ATAC integration evidence.

---

# Latest Wave-5 Figure 1 quality gate update (2026-05-18)

- Quality gate PASS after Figure 1D/1F/1G/1H addition: Python syntax checks, `git diff --check`, and Trevino loader regression tests passed.
- Test evidence: `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> `8 passed in 0.77s`.
- Next immediate action: Figure 1E multimodal marker overlays or Figure 2B public integration evidence.

---

# Latest Wave-5 Figure 1 public-matrix update (2026-05-18)

- Figure 1D/1F/1G/1H public-matrix outputs are available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_1_public_panels`.
- Validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure1_public_matrix_1d_1f_1g_1h.json`; all checks pass `True`.
- Metrics: RNA `55653` cells, RNA markers `28/28`; ATAC UMAP subset `12000` cells; ATAC marker dotplot `31304` cells, ATAC markers `28/28`.
- Honest caveat: Public processed RNA and ATAC gene-activity matrices were used. ATAC UMAP panels use the existing 12k Figure 3B subset embedding; ATAC marker dotplot uses all public 31,304 ATAC metadata-aligned cells. This is not FASTQ/fragments-level or exact author embedding parity.
- Next immediate action: Figure 1E multimodal overlays or Figure 2B public integration evidence.

---

# Latest Wave-5 Figure 3 quality gate update (2026-05-18)

- Quality gate PASS after Figure 3F/3H/gap-audit additions: Python syntax checks, R parse, `git diff --check`, and Trevino loader regression tests all passed.
- Test evidence: `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> `8 passed in 0.78s`.
- Current Figure 3 status: 3A/3C/3D/3E/3F/3H generated/validated; 3B generated with weak PAX6 concordance caveat; 3G/3I honestly marked motif-synergy method-resource gaps.
- Next immediate action: continue remaining target-matrix panels outside Figure 3, prioritizing main Figure 1/2 READY_FOR_LOADER panels.

---

# Latest Wave-5 Figure 3G/3H/3I motif-correlation/synergy update (2026-05-18)

- Figure 3H motif/gene correlation heatmap outputs are available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3H_tf_motif_gene_correlation`.
- Figure 3H validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3h_public_resource_tf_motif_gene_correlation.json`; metrics: `158/682` observed gene-motif correlations, `136` strong |rho| pairs.
- Figure 3G/3I synergy audit: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3g_3i_public_resource_motif_synergy_gap.json`; statuses `{'3G': 'METHOD_RESOURCE_GAP_SYNERGY_NOT_REPRODUCED', '3I': 'METHOD_RESOURCE_GAP_MEAN_SYNERGY_NOT_REPRODUCED'}`.
- Honest caveat: 3H is author-resource motif-level correlation evidence, not exact 24 motif-cluster parity; 3G/3I are not reproduced because chromVAR/motifmatchr and the motif-cluster synergy mapping/workflow are absent locally.
- Next immediate action: run quality checks for new Figure 3 scripts, then proceed to remaining main-figure target-matrix panels.

---

# Latest Wave-5 Figure 3F TF/motif heatmap update (2026-05-18)

- Figure 3F public-matrix TF expression / motif activity heatmap outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3F_tf_motif_heatmap`.
- Validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3f_public_matrix_tf_motif_heatmap.json`; all checks pass `True`.
- Metrics: candidate pairs after filters `1668`, selected pairs `31`, unique genes `31`, unique motifs `22`; RNA cells binned `6000`, ATAC cells binned `12000`.
- Honest caveat: status `GENERATED_PUBLIC_MATRIX_TF_MOTIF_HEATMAP_DERIVED_PAIRS`; TF/motif heatmap from public RNA/chromVAR matrices binned by Figure 3A/3B pseudotime and TF_MotifExpressionCorrelation-derived pairs. Not exact author 363-pseudobulk, motif-cluster, or original 31/24 object parity.
- Next immediate action: Figure 3G/H/I resource audit or remaining Figure 3 panels, then code-quality checks for newly added scripts.

---

# Latest Wave-5 Figure 3E motif enrichment update (2026-05-18)

- Figure 3E public-matrix derived-k5 motif enrichment outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3E_motif_enrichment`.
- Validation PASS for generated outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3e_public_matrix_derived_k5_motif_enrichment.json`.
- Motif matrix readiness: `657930 x 452`; cluster count `5`; all link rows mapped to peak IDs.
- Enrichment: `1021` positive rows in `all_consensus_peaks`; Bonferroni<=0.05 positive rows `160`.
- Honest caveat: status `GENERATED_PUBLIC_MATRIX_DERIVED_K5_MOTIF_ENRICHMENT`; not exact LOLA/JASPAR/topGO environment parity and inherits Figure 3C derived-k5 caveat.
- Next immediate action: Figure 3F TF expression/motif activity heatmap resource audit/build.

---

# Latest Wave-5 Figure 3D gene-set enrichment update (2026-05-18)

- Figure 3D public-matrix local gene-set enrichment outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3D_geneset_enrichment`.
- Validation PASS for generated outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3d_public_matrix_local_geneset_enrichment.json`.
- Loaded `16` local author gene-set files/terms; tested `35` enrichment rows; FDR<=0.05 rows `3`.
- Honest caveat: status `GENERATED_PUBLIC_MATRIX_LOCAL_GENESET_ENRICHMENT_NOT_TOPGO_PARITY`; this is not exact topGO v2.36.0 / GO database-version parity.
- Next immediate action: Figure 3E motif enrichment rework on Figure 3C derived k5 clusters.

---

# Latest Wave-5 Figure 3C CRE-gene heatmap update (2026-05-18)

- Figure 3C public-matrix derived-k5 CRE-gene heatmap outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3C_cre_gene_heatmap`.
- Validation PASS for generated outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3c_public_matrix_cre_gene_heatmap.json`; PNG/SVG/PDF, link-row TSV, matrix NPZ, summary JSON, and log are present/nonempty.
- Coverage: `12927/13989` Table S3B links usable (`0.924`); peak-coordinate mapping complete; finite ATAC accessibility and RNA expression matrices.
- Honest caveat: the paper declares k=5, but local Table S3B has `10` cluster labels; this run derives k=5 from public RNA expression profiles and does not claim exact author 363-pseudobulk/CCA/fragment-level parity.
- Next immediate action: Figure 3D/3E enrichment panels using `figure3c_link_row_summary.tsv`, author motif resources, and local gene-set files.

---

# Latest Wave-5 Figure 3B ATAC pseudotime transfer update (2026-05-18)

- Figure 3B public-matrix ATAC transferred pseudotime outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3B_atac_pseudotime`.
- Generated-output validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3b_public_matrix_atac_pseudotime_transfer.json`; all expected PNG/SVG/PDF/H5AD/JSON/log files exist and are nonempty.
- Honest biological status: `GENERATED_PUBLIC_MATRIX_SUBSET_WEAK_MARKER_CONCORDANCE`, because PAX6 direction failed (`rho=+0.053877` where pre-registered expectation was negative). SOX2 and NEUROD1/NEUROD6/SLC17A7/SATB2 directions passed.
- Computed ATAC subset: `12000 / 31304` cells; transferred from Figure 3A RNA velocity subset with `6000` finite pseudotime reference cells; `25` shared transfer genes; finite transferred pseudotime cells `12000`.
- Boundary: public processed ATAC gene activity + RNA velocity transfer; not fragments/FASTQ-level scATAC and not author integration-model parity.
- Next immediate action: start Figure 3C linkage/co-accessibility evidence using Figure 3A/3B plus Table S3B/FCM/link resources.

---

# Latest Wave-5 Figure 3A velocity reproduction update (2026-05-18)

- Figure 3A public-matrix RNA velocity/pseudotime outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3A_velocity`.
- Validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3a_public_matrix_velocity_reproduction.json`; status `REPRODUCED_PUBLIC_MATRIX_SUBSET`; all output PNG/SVG/PDF files exist and nonempty.
- Computed subset `[6000, 2432]` from `55653` overlapping QC-passed public RNA cells; velocity graph `[6000, 6000]`; velocity pseudotime finite cells `6000`.
- Marker direction check passed: SOX2/PAX6 decrease with pseudotime, NEUROD1/NEUROD6/SLC17A7/SATB2 increase with pseudotime.
- Boundary: this is public processed-matrix, 6000-cell stratified-subset reproduction, not FASTQ-level or full-cell velocity graph.

---

# Latest Wave-5 RNA velocity input update (2026-05-18)

- Full velocity-ready public RNA input: `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_public_rna_velocity_prepared/prepared_input.h5ad`; shape `57868 x 32648`; layers `['ambiguous', 'spliced', 'unspliced']`.
- Layer audit: counts has `33355` genes, velocity layers have `32648` genes, intersection `32648`; `707` counts-only genes are intentionally dropped only for the velocity object.
- scVelo smoke PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/rna_velocity_smoke/scvelo_smoke_retry.json`; 3k-cell subset produced `velocity` layer and `velocity_graph` `[3000, 3000]` with scvelo `0.3.4`.
- Keep the counts-only pipeline object as the clustering/annotation baseline; use this layer-aligned object for Figure 3A/RNA velocity reproduction.

---

# Latest Wave-5 public RNA viability update (2026-05-18)

- Current public-matrix pipeline run: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/`.
- Prepared public RNA input: `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_public_rna_prepared/prepared_input.h5ad`; shape before QC `57868 x 33355`; gene-symbol var_names mapped for marker annotation (`{'mapped_to_symbol_var_names': 32366, 'fallback_to_gene_id': 989, 'unmapped_gene_ids': 976, 'duplicate_symbols_fell_back': 13}`).
- Pipeline output: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458` with `final_adata.h5ad`, `run_manifest.json`, and `module_status.csv` present.
- Module status: `cellranger -> qc -> doublet_detection -> clustering -> annotation -> differential_expression -> composition` all `ok`.
- Final object after QC/doublet filtering: `55653 x 25519`, `24` Leiden clusters; top annotations: Excitatory neuron=36503, Radial glia / neural progenitor=6973, Inhibitory interneuron=4969, Cycling progenitor=3020, Intermediate progenitor=1990.
- Truth boundary: no public FASTQ/SRA raw-read source found yet; this lane is GEO public processed scRNA matrix-level reproduction, with velocity layers deferred to a separate loader audit.

---


# Latest Wave-5.1 Final Completion Update (2026-05-17)

- Wave-5.1 Trevino Cell reproduction evidence package is complete with honest MV3 deferral.
- Final audit: `.omx/reports/wave5-v51-final-completion-audit-20260517.json`.
- v5.1 ledger: `ops/run_ledger/wave5_trevino_20260517T1005Z-7b539a5.v5.1.json` validates against `wave5_v5_1.schema.json`; v4.2/v5.0 schemas reject it.
- v5.1 comparison PDF: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1005Z-7b539a5/python/figures/v5/comparison_3fa006c876fae372_v5.1.pdf`, SHA `3fa006c876fae37289661bc5875766b76bff25634bd1346e0276c0158b239d6d`.
- MV3 F-3/F-6/F-7 status: `HONEST-DEFERRED` at `python/figures/v5/mv3/mv3_disposition.json`; no new dependencies or external GWAS/motif/genome resources were installed or fabricated.
- Final targeted gate: `75 passed, 11 deselected in 1.23s`; contract parity OK; MV2 ledger status is `PASS-CONTRACT-ONLY` (not full analytical R parity).

---

# Latest Wave-5.1 HV1 Update (2026-05-17)

- HV1 real Layer-3 metrics completed for run `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1005Z-7b539a5/`.
- Outputs: `python/figures/v5/quantitative_metrics/*.json` (`15` JSON files) and `python/figures/v5/reproduction_rate.json`.
- Real-metric figure count: `10` main/RNA-only figures with `n_scored_metrics > 0`; AC-V51-HV1-1 PASS (threshold `>=8`).
- Reproduction rate: v5.1 `0.4031615505603515` vs v5.0.4 baseline `0.338`, delta `+0.0651615505603515`; AC-V51-HV1-3 PASS.
- GEO SHA verification: all local GSE162170 files OK; manifest at `data/external/trevino_2021/geo/MANIFEST.json`.
- ATAC peak Spearman remains honest `NOT-COMPUTED`: the public ATAC counts file lacks cell-barcode headers in this local copy, so F-5 is not used to inflate the headline rate.
- Post-HV1 audit: `.omx/reports/wave5-v51-cell-mainline-gap-audit-20260517-post-hv1.json`; remaining mainline gaps are MV3 disposition, comparison_v5.1 PDF, and v5.1 ledger.

---
# Latest Wave-5.1 Trevino Cell Reproduction Update (2026-05-17)

- Current v5.1 evidence run: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1005Z-7b539a5/`.
- HV2 WNN pre-registered sweep completed: `288/288` OK rows, selected `resolution=0.5`, `n_neighbors=20`, `rna_weighted_2x`, `mutual`; ARI `0.5202309143816428`; selected config at `python/figures/v5/rna_only/selected_config.json`.
- MV1a processed-count replay completed: `python/final_adata.h5ad`, `python/run_manifest.json`, and `python/module_status.csv` exist under the project root; parity vs v4.2 PASS at `python/parity/parity_vs_v4_2.json`.
- R bundle exported under the project run at `python/bundle/` with `singlecell_r_bundle_v2.2` and `r_factory_sha_at_export=cba6162`.
- Operational caution: this adopted run id is legacy/no-hyphen (`20260517T1005Z-7b539a5`) and does not satisfy newer `workflow.modular.project_paths` run-id regex; direct project-run paths were used for bundle export.
- New journal: `ops/before_every_run/journal/2026-05-17-wave5-v51-hv2-mv1a.md`.

---

# BeforeEveryRun Latest Summary

## Current Canonical Execution State

- Canonical prepared input:
  `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- Retained fresh stage-1 evidence run:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
- Current clean full-cohort rerun source of truth:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652`
- Current clean rerun launch log:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO.launch.log`
- Current clean rerun local package:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase5_clean_full_cohort_rerun_20260424`
- Current human-facing workspace entrance:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/00_HUMAN_START_HERE.md`
- Current agent-facing run ledger:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/agent_runs/README.md`

## Current Methodology / Governance State

- As of the `2026-04-25` methodology audit, the next paper-aligned optimization
  lane is sparse-exact rather than CSS approximation:
  `/home/zerlinshen/singlecell_factory/ops/nc2024_methodology_audit/AUDIT_2026-04-25.md`
- That audit records paper-faithful Scrublet, Harmony 15-PC clustering, Leiden
  resolution `1.0`, paper-aligned DE parameters, and tumor-vs-background /
  healthy cohort splitting through `--cohort-subset`.
- Paper-aligned launcher:
  `/home/zerlinshen/singlecell_factory/scripts/run_nc2024_paper_aligned_20260425.sh`
- Sparse-exact exploratory launcher:
  `/home/zerlinshen/singlecell_factory/scripts/run_nc2024_full_cohort_sparse_exact_20260425.sh`
- Observed `2026-04-25` `NC2024_NSCLC_FULL_COHORT_SPARSE_EXACT_REAL_AUTO_*`
  result directories are not canonical successful runs unless a later audit
  finds `final_adata.h5ad`, `run_manifest.json`, and `module_status.csv`.
  Current read-only inspection saw only early mandatory outputs/checkpoints and
  a zero-byte sparse-exact launch log.
- Both Mac-led SSH orchestration and direct remote operation are valid. The Mac
  remains the review/organization/report-packaging surface; the remote repo
  remains the compute, remote-R, and run-truth surface.
- As of the `2026-04-30` architecture contract round, module hierarchy metadata
  lives in `/home/zerlinshen/singlecell_factory/workflow/modular/module_catalog.py`.
  `workflow/modular/pipeline.py` derives the legacy dependency DAG from that
  catalog, and CLI optional-module help derives from the same source.
- No full-cohort rerun was launched in the architecture contract round.
  Controller-validation smoke now lives at:
  `/home/zerlinshen/singlecell_factory/scripts/validate_nc2024_architecture_contract.py`
  It checks current NC2024 v2 run dirs, bundle manifests, run-ledger entries,
  and `bridges/local_r_pipeline_macbook/{R,R_bundle}` symlinks without opening
  the large H5ADs.

## Latest Run Verdict

- Date: `2026-04-24`
- Execution mode: `debug_massive`
- Scale mode: `massive`
- Checkpoint policy: `--scale-mode massive` (sets `mandatory_only` automatically)
- Prepare rerun: `no`
- Input: canonical prepared Zarr above.
- Result: clean rerun wrote `final_adata.h5ad`, `run_manifest.json`,
  `module_status.csv`, and `module_reconciliation.tsv`.
- Final object: `810218 x 30374`, `33G`.
- Module status: all requested modules `ok`.
- `pseudobulk_de`: `ok/completed` in the pipeline-native rerun.
- Reconciliation:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652/module_reconciliation.tsv`
- Remote R reporting:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652/r_plots/phase5_readable_20260424`

## Cleanup Truth

- User approved deleting old `results/` artifacts when they are reproduction/test
  outputs and each run is recorded.
- Pre-extended-run cleanup archive:
  `/home/zerlinshen/singlecell_factory/ops/cleanup_records/nc2024_pre_extended_real_run_20260424T051809Z`
- Deleted superseded result directories:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_LARGE_AUTO_20260423_035243`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222`
- These deleted directories must not be cited as currently existing canonical
  artifacts. Use the retained fresh stage-1 run and current clean rerun instead.
- Raw/reference data and the canonical prepared Zarr were not deleted.

## Most Important Carry-Forward Lessons

- Treat `large` as a capacity probe on this cohort, not the main completion lane.
- Use direct `massive` for actual full-cohort survivability, debugging, and
  recovery work.
- Use controller `large -> massive` only for orchestration validation.
- Do not rerun prepare unless the canonical prepared Zarr is invalid.
- Use `--scale-mode massive` for full-cohort runs; it automatically sets
  `checkpoint_policy=mandatory_only` to skip early-module AnnData checkpoints
  unless there is a specific need to preserve full per-module AnnData checkpoints.
- For reproduce-stage NC2024 work, interpret "all modules" as "all eligible
  modules", not "every imaginable module regardless of modality or contrast
  contract".
- Remote R plotting/reporting is remote-side by default on `ubuntu-tail` with:
  `/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript`
  and
  `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh`.
- The previous `/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript` runtime is
  retained as rollback, but it lacks R `arrow` and is not the default for v2
  parquet bundle plotting.
- The bridge folder name still says `local_r_pipeline_macbook`, but current
  operation is remote-first; do not describe the Mac-local R pipeline as active.
- Python still produces module-native QC/diagnostic plots. R is preferred only
  where it can render materially better publication-style figures from the same
  validated artifacts.
- Treat compact R bundle values as plotting/reporting handoff data, not DE-ready
  counts or new quantitative-expression evidence without full-object validation.
- `rna_velocity` without spliced/unspliced layers or loom/BAM+GTF should be
  treated as `skipped_missing_splicing_modality`, not as a surprising failure.
- `pseudobulk_de` should be treated as confirmatory only when an explicit
  contrast contract is provided; exploratory `group_vs_rest` output must be
  explicitly enabled.
- The winning implementation pattern for full-cohort pseudobulk is:
  keep lazy inputs through the general pipeline and materialize the counts layer
  only inside pseudobulk aggregation.
- Do not enable `cnv_inference`/`evolution` on the full cohort until a scale-safe
  CNV lane exists.
- Do not silently upgrade recovered modules: preserve original status files and
  cite recovery outputs separately.
- Use the local `reproduce-run-retention` skill when deciding whether a prior
  reproduce/test run can be deleted after evidence capture.
- Remote Codex / Claude Code agents should have the same Shenxin workflow skill
  coverage. Codex project skills live in `.codex/skills/` using standard Codex
  project skill management; Claude project skills live in `.claude/skills/`.
  `codex_skills/` remains only as a legacy compatibility mirror for historical
  project-local skills.
- For compact bundle handoff, current v2 parquet manifests validate against the
  new singlecell-to-multiomics contract. The `r_multiomics_arrow` environment can
  source `R_bundle/io_bundle.R`, validate NC2024 v2 manifests, read both current
  v2 bundles through `read_bundle()`, and plot them through
  `scripts/plot_remote_bundle_large.R`.
- `bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh` is now
  v2-aware: it reuses `bundle_manifest.json` parquet bundles when source
  `final_adata.h5ad` size/mtime and requested markers/obs/obsm are compatible,
  then delegates plotting to
  `/home/zerlinshen/r_multiomics_factory/scripts/plot_remote_bundle_large.R`.
- Remote report/figure outputs intended for Mac-side review must be pulled back
  as lightweight review copies under
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/human_review/figures/YYYY-MM-DD-<slug>/`
  with `README.md`, `TRANSFER_MANIFEST.tsv`, `checksums.sha256`, and provenance
  logs. Do not pull bulky primary compute objects by default.

## Current Scientific/Reporting State

- Phase-1 fresh massive reproduction package exists and remains the baseline
  package for the eight-module stage-1 reproduction.
- Phase-3 CSS fidelity audit remains the current controlled estimate of
  clean-vs-CSS precision:
  - clean-vs-CSS ARI `0.4190`, NMI `0.6540`
  - CSS-vs-full-massive Leiden ARI `0.5479`, NMI `0.7220`
  - exact cluster identity is weaker than broad biological interpretation.
- The clean rerun now reproduces the same broad module set but with
  `pseudobulk_de` pipeline-native success rather than a separate recovery-only
  artifact.
- Phase-5 package now exists for the clean rerun.
- The Phase-4 package now has a readable deck-style PDF plus a preserved
  `first_pass_unreadable` backup inside the same package.
- The Phase-4 package also now carries a second-pass readable R rerender with
  larger UMAP exports and split marker-dot panels.
- ELF3 should currently be framed primarily as tumor epithelial / tumor-enriched.
  In recovered pseudobulk DE:
  - `Tumor epithelial_vs_rest`: log2FC `6.24528665788722`,
    padj `3.1724278254870602e-24`
  - `Myeloid/Macro_vs_rest`: log2FC `-1.1154218654114758`,
    padj `8.227643967172654e-09`

## Read Before Next Remote Run

- `journal/2026-04-23-nc2024-stage1-recovery.md`
- `journal/2026-04-24-nc2024-fresh-massive-rerun.md`
- `journal/2026-04-24-nc2024-h5ad-schema-inspection.md`
- `journal/2026-04-24-nc2024-phase1-pdf-package.md`
- `journal/2026-04-24-nc2024-checkpoint-hardening.md`
- `journal/2026-04-24-nc2024-css-fidelity-audit.md`
- `journal/2026-04-24-nc2024-r-bridge-performance-hardening.md`
- `journal/2026-04-24-nc2024-remote-r-governance-correction.md`
- `journal/2026-04-24-nc2024-r-environment-hardening.md`
- `journal/2026-04-24-nc2024-bridge-code-review-hardening.md`
- `journal/2026-04-24-nc2024-artifact-inventory-cleanup.md`
- `journal/2026-04-24-nc2024-extended-full-cohort-real-run.md`
- `journal/2026-04-24-nc2024-phase4-report-pseudobulk-contract-hardening.md`
- `journal/2026-04-24-nc2024-clean-rerun-three-pass-resolution.md`
- `journal/2026-04-30-nc2024-architecture-contract-optimization.md`
- `journal/2026-04-30-nc2024-r-arrow-environment-validation.md`
- `journal/2026-04-30-nc2024-remote-to-mac-handoff-protocol.md`

## Current Open Operational Improvement

- Promote post-run pseudobulk recovery into a first-class pipeline recovery
  surface so future manifests can represent both original failure and recovered
  evidence cleanly.
- Add scale-safe CNV/evolution implementation before attempting full-cohort
  CNV/evolution.
- Automate the new remote-to-Mac handoff manifest/checksum generation if this
  pattern repeats across multiple figure/report packages.
- Add direct subprocess test coverage for the remote R wrapper.
- Add a first-class recovery manifest so pipeline failure and recovery can be
  linked more formally than a separate recovery directory.
- Add direct subprocess test coverage for the `r_multiomics_arrow` remote R
  wrapper path.

## Latest Hook-Confirmed Success

- Timestamp: `2026-04-24T20:15:00+08:00` final clean rerun artifacts observed.
- Success evidence:
  - current clean rerun final H5AD exists and is `33G`
  - current clean rerun `run_manifest.json` exists
  - current clean rerun `module_status.csv` exists
  - current clean rerun `module_reconciliation.tsv` exists
  - current clean rerun readable R report directory exists
  - current local Phase-5 package exists


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 06:05 CST

- Mainline remains honest reproduction of the Trevino Cell paper figures from the earliest public/computable inputs available locally; not the NC2024/NSCLC line.
- Latest completed panel: Figure 2E GPC enrichment, status `GENERATED_LOCAL_GPC_GENESET_ENRICHMENT_NOT_TOPGO_PARITY`.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2e_public_local_gpc_enrichment.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2E_gpc_enrichment`.
- Metrics: count-aligned GPCs 185 / paper 185; top TF enrichment `GO_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY` BH q=2.1e-10.
- Fresh verification: py_compile passed for Figure 2B-2E scripts; `git diff --check` passed; Trevino loader tests `8 passed in 0.78s`.
- Boundary: Figure 2E is local Fisher/BH enrichment from public resources, not exact topGO/Table-S2 parity. Continue with Figure 2F-2I multiome resource audit/panels.


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 06:20 CST

- Figure 2F/2G/2H/2I public multiome evidence generated and verified.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2fghi_public_multiome_evidence.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2F_2I_multiome`.
- Statuses: `{'2F': 'REPRODUCED_PUBLIC_AUTHOR_MULTIOME_QC_SUMMARY', '2G': 'GENERATED_PUBLIC_AUTHOR_RNA_AND_BOUNDED_ATAC_MULTIOME_PROJECTION', '2H': 'METHOD_RESOURCE_GAP_MULTIOME_LINKAGE_VENN_NOT_REPRODUCED', '2I': 'GENERATED_PUBLIC_MULTIOME_GPC_CORRELATION_CORRESPONDENCE'}`.
- Fresh verification: py_compile passed for Figure 2B-2I scripts; `git diff --check` passed; Trevino loader tests `8 passed in 0.78s`.
- Boundary: Figure 2H is a resource-gap audit, not a reproduced Venn; Figure 2G ATAC is bounded marker-NN projection, not exact author peak-LSI/uwot parity.


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 06:50 CST

- Figure 4/5 public author FCM panels generated/audited and verified.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure4_5_public_author_fcm_evidence.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_4_5_fcm_public_panels`.
- Statuses: `{'4A': 'GENERATED_AUTHOR_FCM_PUBLIC_PSEUDOBULK_OVERVIEW', '4B': 'REPRODUCED_AUTHOR_FCM_MODULE_HEATMAP_FROM_CENTERS', '4C': 'REPRODUCED_AUTHOR_FCM_SELECTED_GENE_HEATMAP', '4D': 'GENERATED_AUTHOR_FCM_MODULE_SCORE_UMAPS', '4E': 'REPRODUCED_AUTHOR_FCM_CENTROID_JACCARD_NETWORK', '4F': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4G': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4H': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4I': 'IMAGE_PANEL_NOT_COMPUTATIONAL_REPRODUCTION_LOCAL_PDF_ONLY', '4J': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4K': 'IMAGE_PANEL_NOT_COMPUTATIONAL_REPRODUCTION_LOCAL_PDF_ONLY', '5A': 'REPRODUCED_AUTHOR_FCM_ASTRO_GENE_MEMBERSHIP_EXPRESSION', '5B': 'GENERATED_AUTHOR_FCM_LINKED_PEAK_MOTIF_ENRICHMENT_M13_VS_M14', '5C': 'GENERATED_AUTHOR_FCM_AQP4_POSITIVE_RECLUSTERING_LOCAL_RULE', '5D': 'GENERATED_PUBLIC_FCM_PSEUDOBULK_DE_NOT_DESEQ2_PARITY', '5E': 'METHOD_RESOURCE_GAP_BHADURI_PRIMARY10X_NOT_STAGED', '5F': 'METHOD_RESOURCE_GAP_BHADURI_PRIMARY10X_NOT_STAGED'}`.
- Fresh verification: FCM render command passed; py_compile passed; `git diff --check` passed; 14 PNGs verified; Trevino loader tests `8 passed in 0.77s`.
- Boundary: 4I/4K are image-only IHC panels, not recomputed matrix outputs; 5D is a local pseudo-bulk Welch-test proxy, not DESeq2/Table S5 parity; 5E/5F need Bhaduri/UCSC `primary10X` staging before reproduction can be claimed.
- Next: Figure 6 chromatin-state/GPC branch resource audit/panels.


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 07:05 CST

- Figure 6 public FCM/GPC branch panels generated/audited and verified.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure6_public_fcm_gpc_branch_evidence.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_6_fcm_gpc_branch`.
- Statuses: `{'6A': 'REPRODUCED_PUBLIC_FCM_CELL_CYCLE_MODULE_CORRELATION', '6B': 'GENERATED_PUBLIC_FCM_ATAC_PROJECTION_SCHEMATIC', '6C': 'REPRODUCED_AUTHOR_FCM_ATAC_BRANCH_PROJECTION_FROM_GA_ANCHORS', '6D': 'GENERATED_PUBLIC_FCM_BRANCH_GENE_ACTIVITY_HEATMAP_LOCAL_SPECIFICITY', '6E': 'GENERATED_PUBLIC_FCM_GPC_TF_MOTIF_GENE_BRANCH_HEATMAP', '6F': 'REPRODUCED_AUTHOR_FCM_GPC_ANCHOR_REPROJECTION_BRANCH_ARROWS', '6G': 'GENERATED_BOUNDED_MULTIOME_RNA_TO_FCM_NN_PROJECTION_COLORED_BY_ATAC_CLUSTER'}`.
- Fresh verification: Figure 6 render command passed; py_compile passed; `git diff --check` passed; 7 PNGs verified; Trevino loader tests `8 passed in 0.85s`.
- Boundary: 6D/6E use local public-branch specificity/means; 6G is bounded NN multiome RNA→FCM projection, not exact author uwot parity.
- Next: Figure 7 disease/BPNet/ASD resource audit.


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 07:06 CST

- Figure 7 public S6/BPNet-ASD panels generated/audited and verified. Main Trevino Cell Figures 1-7 now have public-resource evidence/outputs or explicit resource-gap audits.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure7_public_asd_bpnet_evidence.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_7_asd_bpnet_public_panels`.
- Statuses: `{'7A': 'GENERATED_FROM_METHODS_AND_PUBLIC_TABLE_S6_COUNTS', '7B': 'GENERATED_S6D_CLUSTER_HIGH_EFFECT_SUMMARY_ON_PUBLIC_ATAC_UMAP_NOT_EXACT_AUTHOR_FISHER_DENOMINATORS', '7C': 'GENERATED_S6_COUNT_AUDIT_EXACT_OR_1_909_DENOMINATORS_NOT_PUBLIC_IN_S6', '7D': 'GENERATED_SFARI_AUDIT_CAPTION_BENCHMARK_PLUS_RECOMPUTED_COUNTS_NOT_EXACT_24_17_FROM_LOCAL_TABLES', '7E': 'REPRODUCED_FROM_PUBLIC_TABLE_S6E_MOTIF_OVERLAP_EXCESS', '7F': 'GENERATED_NFIA_S6_LOCUS_SUMMARY_NOT_EXACT_BPNET_LOGO_TRACK_PARITY', '7G': 'GENERATED_NPY_S6_LOCUS_SUMMARY_NOT_EXACT_BPNET_LOGO_TRACK_PARITY'}`.
- Fresh verification: Figure 7 render command passed; `python -m py_compile scripts/dev/wave5_render_public_figure7_asd_bpnet.py` passed; `git diff --check` passed; 7 PNGs verified; Trevino loader tests `8 passed in 0.78s`.
- Boundary: S6D duplicated `Cluster GluN6` yields 261/232 caption-like high-effect case/control vs paper 262/232; exact Figure 7B/7C Fisher denominators and 7F/7G BPNet logo/ref-alt tracks are not deposited in public S6/local Brain_ASD checkout.
- Next: final quality gate/code-review and consolidated completion handoff; no additional long compute phase is currently queued.


## Current Wave5 Cell/Trevino all-figure reproduction completion — 2026-05-18 07:13 CST

- Completed conditional public-resource all-main-figure reproduction/audit.
- Completion report: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/TREVINO_ALL_FIGURES_COMPLETION_REPORT.md`.
- Final quality gate: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/all_figures_final_quality_gate_review.json`; code review `APPROVE/CLEAR`; scientific decision `conditional` due documented public-resource gaps.
- Evidence audit: 17 JSON files, 53 panel statuses, 0 missing referenced paths, 0 hash mismatches, 0 unexpected false checks.
- Pipeline viability: public RNA pipeline run has all modules `ok` and final AnnData under `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/final_adata.h5ad`.
- Boundary: not FASTQ/raw-fragment parity; Figure 7 exact BPNet denominators/model tracks and image-only panels remain honest resource gaps.

## Current pipeline validation G002/G003 LUSC update — 2026-05-27 23:25 CST

- Started next repo-native ultragoal round because the hidden Codex goal slot is still occupied by the previous completed objective.
- G001 claim ledger created under `/home/zerlinshen/projects/pipeline-validation-20260527/`:
  - `claim_ledger.tsv`
  - `evidence_index.json`
  - `G001_claim_ledger_report.md`
- G002/G003/G004 LUSC real-data validation completed on the 92,430-cell / 9-dataset LUSC Q22 h5ad.
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/marker_retention.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/embedding_label_purity.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/cell_type_mixing.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/rare_population_preservation.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g002_g003_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/G002_G003_lusc_marker_mixing_report.md`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/G004_lusc_rare_population_report.md`
- Result:
  - Marker retention LUSC pass: 23/23 canonical/edge marker panels passed, marker coverage via `var.feature_name`, minimum AUROC 0.743973.
  - Cell-type mixing LUSC partial pass: Harmony improved same-batch neighbor fraction for 18/24 labels.
  - Purity flags for G004: `Alveolar cell type 1`, `DC mature`, and `cDC2` lost >0.10 same-label neighbor purity after Harmony.
  - Rare/edge preservation LUSC partial pass: 13/19 pass; review labels are `Alveolar cell type 1`, `DC mature`, `T cell regulatory`, `cDC1`, `cDC2`, `other`.
- Next: repeat G002/G003/G004 on Trevino before declaring global pass.

## Current pipeline validation G002/G003/G004 Trevino update — 2026-05-27 23:43 CST

- Ralph sequential validation continued with Trevino public RNA real data:
  - `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/final_adata.h5ad`
  - 55,653 cells, 25,519 genes, 8 samples, 9 `cell_type` labels.
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_marker_retention.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_embedding_label_purity.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_cell_type_mixing.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_rare_population_preservation.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_g002_g003_g004_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/G002_G003_G004_trevino_marker_mixing_report.md`
- Result:
  - Trevino marker retention pass: 9/9 canonical brain marker panels passed.
  - Trevino cell-type mixing partial pass: Harmony improved same-sample neighbor fraction by >0.05 for 8/9 labels.
  - Trevino rare/edge preservation partial pass: 5/7 labels passed.
  - Review labels: `Inhibitory interneuron` and `Intermediate progenitor` lost >0.10 same-label neighbor purity after Harmony.
- Validation: script `py_compile` passed; JSON parsed; TSV column integrity passed.
- Current combined status: G002 marker retention passed across LUSC and Trevino; G003/G004 are validation-complete but retain real biological purity/rare-population flags.
- Next: G005 downstream annotation and sample-level DE sanity on LUSC.

## Current pipeline validation G005 LUSC annotation/DE update — 2026-05-27 23:47 CST

- Ralph sequential validation completed downstream annotation/DE sanity on:
  - `/home/zerlinshen/projects/lusc-integration-gate-20260527/prepared/lusc_squamous_dataset_axis.h5ad`
  - 92,430 cells, 17,764 genes, 87 samples, 41 donors, 9 datasets.
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_annotation_summary.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_pseudobulk_sample_metadata.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_de_origin_tumor_primary_vs_normal_adjacent.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_de_tumor_stage_advanced_vs_early.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_sentinel_gene_checks.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_dataset_confounding.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/G005_lusc_annotation_de_sanity_report.md`
- Result:
  - Annotation marker/sample/dataset support: 23/24 labels pass; `other` remains review-only.
  - Annotation confidence column is unavailable in this LUSC object, so confidence cannot be claimed.
  - Origin sentinel direction sanity: 10/11 expected tumor/normal marker directions passed; `AQP5` mismatched.
  - Sample-level DE is computable but not final-claim safe: origin and tumor-stage contrasts both show top-DE dataset domination.
  - `origin:tumor_primary_vs_normal_adjacent`: top50 DE genes with dataset eta2 >0.50 = 0.82.
  - `tumor_stage:advanced_vs_early`: dataset/group Cramer's V = 0.893208 and top50 DE dataset eta2 fraction = 1.0.
- Validation: script `py_compile` passed; JSON parsed; all six TSV outputs passed column-count integrity checks.
- Current status: G005 is `conditional_pass_annotation_de_sanity_with_dataset_confounding_flags`.
- Next: G006 real 3D genome / Hi-C factory bridge validation.

## Current pipeline validation G006 Hi-C bridge update — 2026-05-27 23:55 CST

- Ralph sequential validation completed a real 3D contact factory bridge using:
  - `/home/zerlinshen/projects/ng2025-3d-genome/data/open/hic_lung/TCGA_HiChIP_hic/LUSC_H3K27ac.allValidPairs.hic`
  - hicstraw extraction: chr21 at 100 kb.
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/LUSC_H3K27ac_chr21_100000bp.contacts.tsv.gz`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/factory_run/hic_ingest/hic_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/factory_run/hic_tad/tad_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/g006_real_hic_factory_bridge.h5ad`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/singlecell_r_bundle_v22_hic_real/bundle_manifest.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/g006_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/G006_hic_factory_bridge_report.md`
- Result:
  - hicstraw records/total count: 49,616 / 1,766,659.
  - `hic_ingest` and `hic_tad` completed with status `ok`.
  - Bundle schema is `singlecell_r_bundle_v2.2`, Hi-C extension status `active`.
  - Reference checks passed: bundle contacts equal matrix nnz, and ingest total contacts equal hicstraw total contacts.
  - Existing Hi-C/bundle regression checks: `13 passed in 0.54s` with `--no-cov`.
- R-side load:
  - Passed using `/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript`.
  - stdout: `STATUS=active`, `CONTACTS_NROW=98860`, `COMPARTMENT_STATUS=low_information`.
- Scientific boundary:
  - Source is H3K27ac HiChIP, not unbiased Hi-C; TAD/compartment calls are technical module-contract evidence, not final 3D-genome biology claims.
  - Compartment status is `low_information` for chr21.
- Current status: G006 is `passed_with_hichip_scientific_boundary`.
- Next: G007 cross-factory handoff and documentation sync.

## Current pipeline validation G007 doc-sync update — 2026-05-28 00:02 CST

- Cross-factory documentation sync completed for:
  - `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory`
  - `/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory`
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/G007_cross_factory_doc_sync_report.md`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g007_doc_sync_summary.json`
- Updated singlecell surfaces: `README.md`, `AGENTS.md`, `AI_AGENT_PROTOCOL.md`, `CODEX.md`, `CLAUDE.md`, `PROTOCOL.md`, and this run memory.
- Updated R-side surfaces: `README.md`, `AGENTS.md`, `AI_AGENT_PROTOCOL.md`, `CODEX_PROFILE.md`, `CLAUDE.md`.
- Verification:
  - `repo-doc-sync --strict` for `singlecell_factory`: NO DRIFT.
  - `repo-doc-sync --strict` for `r_multiomics_factory`: NO DRIFT.
  - `git diff --check` passed for touched docs in both repos.
  - `py_compile` passed for all three validation scripts.
- Current status: G007 passed.
- Next: G008 final independent code and science review.

## Current pipeline validation G008 final review update — 2026-05-28 00:40 CST

- Final artifacts:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/G008_final_code_science_review_report.md`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g008_final_review_summary.json`
- Independent review:
  - Code-review lane: `APPROVE`.
  - Architect/science lane: `WATCH`.
- Registered verdict: `PASS_SUPPORTED_NOT_FINAL`.
- Meaning:
  - Supported: tested marker retention, bounded integration behavior with flags,
    raw-count-compatible DE computation with confounding flags, and real
    H3K27ac HiChIP bundle handoff.
  - Not final-ready: clean rare/edge preservation across all labels, final LUSC
    DE biology, unbiased Hi-C TAD/compartment biology, and global generalization.
