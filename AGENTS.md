# AGENTS.md - singlecell_factory

This file applies to `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory` and all subdirectories.
It overrides higher-level guidance where rules conflict.

## Bioinformatics Research Pipeline suite

Physical suite layout: this repository is stored at `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/`; do not recreate the retired `/home/zerlinshen/singlecell_factory` compatibility path.

This repo is the Python/global-control-plane member of **Bioinformatics Research
Pipeline**, the umbrella suite for `singlecell_factory`, `r_multiomics_factory`,
and `plotting_factory`. Treat the suite name as public/project identity; do not
merge repositories, rewrite git history, or move project data unless explicitly
planned. Governance follows **owner-by-primary-output**: this repo owns Python
upstream and global validation, `r_multiomics_factory` owns R-heavy/spatial primary scientific truth when R creates primary objects/statistics/interpretation,
and `plotting_factory` is a presentation-only plotting surface that must not own
biological conclusions.

## Authority / Read This First

Start with `AI_AGENT_PROTOCOL.md` for onboarding, read order, and task routing.
Then follow this file as the Codex/root project contract. Use `PROTOCOL.md` for
deep operational steps.

The canonical suite startup is [../QUICKSTART.md](../QUICKSTART.md). For a
complete cross-factory failure report from this checkout, use:

```bash
GATE_REPORT_ALL=1 bash ../scripts/run_all_gates.sh
```

## Mission
- Protect biological correctness, statistical validity, and reproducibility.
- Prefer small, reversible edits with verifiable evidence.
- Preserve pipeline contracts and output compatibility unless a breaking change is explicitly requested.

## Architecture (2026-05+): Factory-Project Separation

As of 2026-05, the working tree follows a three-way split. See full plan at
`/home/zerlinshen/.omc/plans/factory-project-separation.md`.

**Factories are tools, projects are data.** This repo (`singlecell_factory`) is
an immutable compute tool. All scientific artifacts must land under a project
directory resolved from `--project-root`, never inside this repo's tree.

- **Factory tools**: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/` (this repo),
  `/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory/`, and
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/` — hold only code, fixtures, and contracts.
- **Projects directory**: `/home/zerlinshen/projects/<project-id>/` (default
  `PROJECTS_ROOT`). Each project is a self-contained directory with `project.yaml`,
  `inputs/`, `runs/`, `configs/`, and `notebooks/`.
- **Bootstrap tool**: `/home/zerlinshen/projects-bootstrap/omc-new-project` —
  creates new projects outside all three factories.

### --project-root contract

All pipeline entry points accept `--project-root <path>` and `--run-id <id>`.
Outputs write to `<project-root>/runs/<run-id>/python/{,bundle/}`. R outputs
write to `<project-root>/runs/<run-id>/r/`. Manifests land at
`<project-root>/runs/<run-id>/manifest.json`.

Run-id format: `<UTC-timestamp>-<short-py-sha>`, regex
`^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{4}Z-[0-9a-f]{7}$`. Auto-generated if absent.

Legacy invocations without `--project-root` remain behind
`contracts/legacy_factory_output_deprecation.yaml`. A real compatibility launch
must emit a warning and append external JSONL access telemetry. Do not remove the
path until its 14-day monitoring window is silent, project-root migration is
complete, and owner approval is recorded; the current contract keeps both latter
conditions false.

Dirty-tree gate: the pipeline refuses to start when the factory git tree is dirty
unless `--allow-dirty` is passed. With `--allow-dirty`, `diff_sha256` is recorded
in `manifest.json` for forensics.

### Contracts vendoring

The Python-R bundle schema lives canonically at `contracts/bundle_schema.yaml`
in this repo. The vendored copy in `r_multiomics_factory/contracts/` must be
byte-identical. To update: edit the canonical file here, then run
`tools/sync-contracts.sh` — **never edit the vendored copy directly**.

### Project conventions (user-decided)

- Project ID format: lowercase-with-dashes, regex `^[a-z0-9][a-z0-9-]*$`
- Raw inputs: symlinks to `/data/raw/` — do not duplicate large files
- Migration trigger: STAGED — CLI plumbing landed (2026-05), actual data move
  deferred to a future explicit user-approved step after dry-run review

### Environment Switches (Round-1a)

Plan: `/home/zerlinshen/.omc/plans/factories-optimization-round1.md`

| Variable | Unset (default) | `=1` |
|---|---|---|
| `SC_REQUIRE_PROJECT_ROOT` | Contracted warning + external JSONL access event; legacy `results/` fallback | Hard error `sys.exit(2)`. Pass `--project-root` or unset. |

Warning and hard-error are mutually exclusive. Any later default flip must
satisfy the committed deprecation contract; the gate remains mirrored in
`scripts/export_singlecell_r_bundle.py` and `scripts/pack_run_for_mac.sh`.

### R-factory SHA fields (Round-1a)

Plan: `/home/zerlinshen/.omc/plans/factories-optimization-round1.md`

Two provenance fields track `r_multiomics_factory` git SHA at two points:

- `r_factory_sha_at_manifest_write` — written to `run_manifest.json` at pipeline manifest-write time.
- `r_factory_sha_at_export` — written to `<run-id>/python/bundle/provenance.json` at bundle export time.

The R bundle loader (`r_multiomics_factory/R_bundle/io_bundle.R`) logs a `WARNING` (not error) when these differ, surfacing both SHAs.

### Agent rule (MANDATORY)

**Never write outputs inside the factory tree. Always resolve outputs from
`--project-root`.** Any module that constructs a path under `singlecell_factory/`
for scientific artifacts is a bug.

**SC_REQUIRE_PROJECT_ROOT=1 fail-fast**: When invoking the pipeline, prefer
setting `SC_REQUIRE_PROJECT_ROOT=1` to fail-fast on a missing `--project-root`.
Never rely on the legacy `results/` fallback. It is governed by
`contracts/legacy_factory_output_deprecation.yaml`, emits a warning, and records
an access event outside the checkout. Do not remove it until migration, owner
approval, and the full 14-day silent window are all proven.

## Current State (2026-08-02)

The August pipeline-remediation round preserves the scientific claim posture
below while adding the canonical suite quickstart, 33-gate report-all route,
measured legacy-output deprecation, render-only evidence contracts, and
ground-truth promotion governance. Use `../QUICKSTART.md` and the committed
suite plan at `../docs/superpowers/plans/2026-08-02-pipeline-remediation.md` for
the executable entrypoints and acceptance record.

Integration-selection gate + batch correction hardened (comprehensive review + ralplan campaign) —
canonical detail in `docs/SCIENTIFIC_AUDIT_2026-05-15.md`:

- scVI now consumes the SAME ambient-corrected signal as clustering/DE (was STALE pre-ambient
  counts), rounded to integers for the NB likelihood, with a fail-loud integrality guard
  (contract `scvi-counts-v3-ambient-corrected-rounded-int`).
- `integration_select` requires baseline/Harmony/scVI/shuffle candidates; a degraded candidate
  set or non-firing shuffle falsifiability control fails loud. The cache key folds
  batch-label state plus all gate/scVI-affecting params, and the cache-HIT path
  re-asserts engine version pins (fail loud on drift).
- `batch_correction` CPU post-processing failure fails loud by default (audited
  `SC_ALLOW_BATCH_BACKEND_SKIP` opt-in → a distinct non-"completed" status), never silently keeping
  pre-correction clustering.

Integration-biology and multiomics validation artifacts now live under
`/home/zerlinshen/projects/pipeline-validation-20260527/` with the durable
ledger at `.omx/ultragoal/ledger-integration-biology-multiomics-validation-20260527.jsonl`.
Current claim posture:

- G002 marker retention passed on real LUSC and Trevino data with explicit
  marker label-shuffle negative controls.
- G003/G004 mixing and rare-population preservation are validation-complete with
  real purity/rare-population review flags and embedding-shuffle negative controls;
  do not collapse them into a clean pass.
- G005 annotation/DE sanity is conditional because LUSC sample-level origin and
  tumor-stage DE top genes are dataset-dominated; DE now uses only 71/87
  raw-count-compatible samples after excluding 16 fractional-count samples.
- G006 passed the real `.hic` -> `hic_ingest`/`hic_tad` -> v2.2 bundle -> R
  `load_hic_extension` bridge with all 7 technical gates, but the source is
  H3K27ac HiChIP and chr21 compartment status is `low_information`, so it is not
  unbiased Hi-C biology proof.
- 2026-05-28 follow-up evidence lives under
  `/home/zerlinshen/projects/pipeline-validation-20260528/`: targeted
  sensitivity passes for LUSC AT1/cDC2 and both Trevino flagged populations,
  while LUSC DC mature remains review; LUSC origin DE has a matched
  dataset-aware sanity pass on true-count-compatible samples, while tumor-stage
  DE is conditional single-dataset evidence; GM12878 unbiased Hi-C chr19 A/B
  compartments pass published-subcompartment concordance after sparse-tail
  handling fixes, but the current factory TAD boundary caller still fails
  B35T1NC Micro-C author TAD concordance and must remain exploratory.
- Current follow-up verdict is `PASS_SUPPORTED_NOT_FINAL_WITH_REVIEW_FLAGS`:
  evidence supports bounded tested claims, while LUSC DC mature, LUSC tumor-stage
  DE, and factory TAD boundary biology remain review/conditional. Downstream
  reports must keep the non-final claim boundaries.
- 2026-05-29 Round 3 suite evidence at
  `/home/zerlinshen/projects/pipeline-validation-20260528/` updates the
  remaining data asks: Liu BEP2D TP63 CUT&Tag count-FDR is now real-data
  supported within the GSE272822/SRP521453 claim boundary; external
  differential A/B flip scoring remains blocked by matched-contact data
  availability/storage; LUSC late-stage UICC DE remains blocked by
  multi-dataset early/late balance failure. Use
  `SUITE_LEVEL_PIPELINE_OPTIMIZATION_REPORT_R3.md`,
  `ULTRAGOAL_FINAL_QUALITY_GATE_R3.json`, and ledger key
  `round_3_2026_05_29_codex_followup`.

## Current State (2026-05-25)

The current suite-level validation posture is intentionally honest:

- The one-command suite gate is a 33-step structure/contract/science-assertion
  harness, including render-only, raw-data, documentation, real-data parity,
  and real-data science-assertion checks. Figure parity remains conditional
  where declared; the gate is not universal biological validation.
- Current bounded real-data handoff proof:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-24T2347Z-0321773`.
  It uses a retained Round9 LUSC consensus AnnData and validates current
  `singlecell_factory` compact-bundle export into `r_multiomics_factory` and
  `plotting_factory` rendering.
- Failed exploratory raw-input rerun:
  `/home/zerlinshen/projects/hgmm-smoke/runs/2026-05-24T1930Z-0321773`.
  It hung after scaffold creation and is not validation evidence.
- Both `r_multiomics` and `r_multiomics_arrow` can read/plot v2/v2.1/v2.2
  parquet bundles as of 2026-05-29. `r_multiomics` now has `r-arrow=24.0.0`;
  `r_multiomics_arrow` remains the renv-pinned reference env.

This does not prove NG2025 end-to-end reproduction or full raw-input pipeline
correctness. Treat restored real-data and curated figure-parity gates as open
release-readiness work.

## Historical State (2026-05-16)

**Wave-5 Trevino PCW21 (biology-aware validation pivot, plan v4.2)** is the most recent canonical work:

- Canonical run: `/home/zerlinshen/projects/wave5-trevino/runs/20260516T0931Z-d192836f1bb0/`
- Ledger: `ops/run_ledger/wave5_trevino_20260516T0931Z-d192836f1bb0.v4.2.json` (plan_revision=v4.2, validation_posture=biology-aware)
- Schema: `ops/run_ledger/schema/wave5_v4_2.schema.json` (Draft 2020-12, allOf+if-then on plan_revision)
- Plan v4.2 (gitignored, agent-state): `.omc/plans/wave5-completion-consensus-2026-05-16-v4.2.md`
- Outcomes: AC-VAL-3a panel_recall=0.611 CLOSED, AC-VAL-3b top50_hit_rate=0.08 CLOSED-PARTIAL per §3.4, AC-CI-1 ARI(RNA-only res=0.3)=0.674 CLOSED, AC-VAL-PLOT-1/2/3 + AC-LEDGER-1 + AC-VAL-3c CLOSED. CI gate exit 1 (CLOSED-PARTIAL overall).
- Biology recovery confirmed (OLIG2, SOX10, NKX2-2, PAX6, EOMES, VIM, HES1, GAD2, DLX1, DLX2, OLIG1 in top-1000); Trevino S2F top-K contamination diagnosed (MS4A12/FCRLA/SFTPC).

**NC2024 NSCLC full-cohort (~884k cells)** carry-over from 2026-04-26 is also canonical and completed both stages with paper-aligned parameters:

- **v1 outputs (annotation bug, do not cite)**: `results/nc2024_tumor_20260426/`, `results/nc2024_bh_20260426/`
- **v2 outputs (canonical)**: `results/nc2024_tumor_20260426_v2/`, `results/nc2024_bh_20260426_v2/`
- **Ledger entries**: `ops/run_ledger/nc2024_*_v2_*.json`
- **Methodology audit**: `ops/nc2024_methodology_audit/{AUDIT_2026-04-26_v2.md, ALIGNMENT_REPORT_v2_2026-04-26.md, SEGFAULT_TRACE_2026-04-26.md, SMALL_REAL_VALIDATION_2026-04-26.md}`
- **Publication template**: `docs/PUBLICATION_READY.md`
- **Mac transmission recipe**: `ops/MAC_PULL_RECIPE_2026-04-26.md`

Cell-type composition aligns with Sanchez-Mejias *Nat Commun* 2024 (doi:10.1038/s41467-024-48700-8). Biology-validation notebooks (4 paper findings on real v2 data) are pending — owner intends to do them on the day after this docs sync.

## Key Pipeline Surfaces (current)

- **Resource strategy**: `--scale-mode standard|large|massive` may resolve only
  `--lazy-read` and `--checkpoint-policy`. It must preserve optional modules,
  HVGs, PCs, neighbors, Leiden resolution, DE limits, doublet strategy, and
  clustering engine. Run metadata records `resource_strategy` separately.
- **Scientific profile**: canonical defaults are
  `clustering,differential_expression,annotation`. Any
  `--scientific-profile legacy-large|legacy-massive|paper-15pc` change requires
  `--acknowledge-scientific-non-equivalence`; the exact resolved scientific
  parameter diff is recorded.
- **One source of scientific truth**: `cli._CANONICAL_SCIENTIFIC_PARAMETERS` is
  authoritative. The argparse defaults read from it, and the `config.py`
  dataclass defaults are pinned to it by
  `tests/test_paper_param_alignment.py::test_dataclass_defaults_match_canonical_profile`.
  Never change a scientific default on only one of these surfaces — that is the
  exact bug fixed on 2026-07-28, where `PipelineConfig()` silently ran
  `n_pcs=15`/`resolution=1.0` while every CLI run used `40`/`0.8`. Paper-specific
  values belong in a named `_SCIENTIFIC_PROFILE_OVERRIDES` entry, not in a dataclass default.
- **Annotation strategy**: `--annotation-strategy {cluster_voting,cell_argmax}` (default `cluster_voting`; `cell_argmax` retained as fallback). The cell_argmax path drifts on >100k cohorts — do NOT use as default.
- **Memory enforcement**: `workflow/modular/_mem_guard.py` provides `MemoryEnforcer` (pre-flight budget + watchdog Event + cooperative abort via `SkipModule`). `SC_MEM_GUARD=on SC_MEM_WATCHDOG=on` enables runtime enforcement.
- **Sparse engines** (env flags): `SC_DE_ENGINE=sparse SC_CNV_ENGINE=chunked SC_CELLCOMM_ENGINE=sparse`.
- **Densify policy**: `workflow/modular/_densify_policy.py` returns `{GO, CHUNK, ABORT}` per planned-bytes budget. Any `.toarray()/.todense()` in `workflow/modular/modules/` requires a `# densify-allowed: <reason>` comment or a call through `plan_densify`. CI grep-ban via `tests/test_densify_audit.py`.
- **Shutdown cleanup**: `workflow/modular/_shutdown.py` runs cupy/torch/zarr explicit cleanup at interpreter exit (prevents the post-completion segfault that previously affected exit codes).

## Staircase Test Fixtures

Generated by `scripts/build_staircase_fixtures.py` from the prepared NC2024 zarr:

| Tier | Cells | Source | Marker |
|---|---|---|---|
| `nano` | 5k | synthetic CSR (conftest) | `pytest -m nano` |
| `small_real` | 100k | NC2024 subset | `pytest -m small_real` |
| `medium_real` | 500k | NC2024 subset | `pytest -m medium_real` |
| `full_real` | 884k | full NC2024 zarr (path file only) | `pytest -m full_real` |

Discipline: nano on every commit; small_real before any module behavior change reaches main; medium_real before phase-completion gates; full_real before publication-grade runs.

## Default Workflow (Codex)
1. Use `$ralplan` for non-trivial work.
2. Execute with focused implementation (`executor` by default).
3. Run `$review-and-validate-quality` before completion.
4. For high-risk changes (pipeline contracts, output schema, mutating modules), include `$code-review` and `$security-review`.

## Operating Modes

Both modes are valid:
- Direct remote operation inside `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory` for
  pipeline execution, run triage, remote R/reporting, and source-of-truth
  artifact inspection.
- Mac-led SSH orchestration from the local reproduction workspace when the task
  is coordination, review, handoff, or report packaging.

The remote repo remains the canonical run-truth surface in both modes. Do not
promote Mac-local summaries over `run_manifest.json`, `module_status.csv`,
`ops/before_every_run/LATEST.md`, or `ops/run_ledger/` evidence.

## Skill Routing
- Pipeline run/recovery: `$execute-and-recover-pipeline`
- Module feature/integration/refactor: `$develop-and-integrate-module`
- Performance work: `$optimize-and-guard-performance`
- Large cross-domain work: `$orchestrate-large-task`
- Final quality gate: `$review-and-validate-quality`
- Before-run memory / run journaling: `$before-every-run`
- Remote/local handoff and bridge work: `$singlecell-remote-workflow`
- Report package / handoff work: `$remote-run-report-bridge`
- Governance and README drift: `$reproduction-workspace-governance` then `$readme-sync-enforcer`
- Cleanup / retention after failed or superseded reproduce runs:
  `$post-run-failure-cleanup`, `$reproduce-run-retention`, `$test-run-cleanup-policy`

Prefer these project skills over ad-hoc prompting when they match.

Codex project skills live in `.codex/skills/` using standard Codex project
skill management. `codex_skills/` remains only as a legacy compatibility mirror
for historical project-local skills.

## Contract Rules (Must Keep)
- Mandatory module chain remains:
  `cellranger -> qc -> ambient_correction -> doublet_detection`.
- Dependency graph updates must be reflected in `workflow/modular/pipeline.py`.
- Each module writes only within its own module directory via context helpers.
- `module_status.csv` and `run_manifest.json` are source-of-truth artifacts.
- Keep output filenames stable unless a documented breaking change is requested.
- `trajectory` and `pseudo_velocity` are explicit opt-ins, not universal
  defaults. Claim-capable DPT requires an existing explicit Leiden root cluster
  plus a non-empty biological justification; otherwise its provenance and
  artifacts must remain non-claimable.
- `pseudo_velocity` is always `exploratory_proxy`, never canonical RNA velocity.
- The public GSE162170 Figure 3A handoff is a render-only boundary:
  `scripts/produce_public_rna_velocity_figure3a_artifacts.py` owns all scVelo
  analysis and emits immutable strict named CSVs plus a hash-linked manifest.
  Consumers must fail on schema/hash drift and must not recompute analysis; the
  lane remains `exploratory`.
- Explicit confirmatory pseudobulk requires a named biological sample column,
  one condition per sample, valid labels, raw counts, and sufficient biological
  replicates. Invalid contracts record non-claimable inference and fail loud.
  Only an all-pydeseq2 confirmatory result is claimable; rank fallback is
  exploratory/non-claimable. Keep `pseudobulk_counts.csv` numeric and place
  sample/group/condition plus inference fields in aligned
  `pseudobulk_metadata.csv`.
- Composition counts and proportions are grouped by biological sample.
  Confirmatory scCODA additionally requires an explicit, separate condition
  column with one condition per sample and sufficient replicates; without that
  contract the module is descriptive/non-claimable. Never use sample IDs,
  Leiden clusters, or a heuristic batch column as the model condition.
- Scale/resource strategies must never change scientific parameters or module
  selection. Use an acknowledged scientific profile for non-equivalent changes.
- `2026-04-25` sparse-exact probe directories are not canonical successful
  NC2024 evidence unless they contain `final_adata.h5ad`, `run_manifest.json`,
  and `module_status.csv`.

## Required Touchpoints For New/Changed Modules
- `workflow/modular/modules/<module>.py`
- `workflow/modular/pipeline.py` (`MODULE_DEPENDENCIES` / registry wiring)
- CLI/config wiring when new flags are added
- Tests in `tests/`
- User-facing docs when behavior changes (`README.md`, `PROTOCOL.md`)

## Verification Gate
Run the lightest set that proves correctness:
- Focused: `pytest -q tests/test_modular.py tests/test_modular_optimizations.py`
- Full regression (when needed): `pytest -q`
- Pipeline behavior checks should include relevant run artifacts and status outputs.

For long or fragile runs:
- Prefer `--checkpoint`
- Use `--resume-from <module>` for recovery

## Reporting Requirements
Final report must include:
- Changed files
- Commands run for verification
- Artifact evidence (`module_status.csv`, `run_manifest.json` when applicable)
- Remaining risks or untested paths

## Safety
- No new dependencies without explicit request.
- Do not silently alter biological interpretation logic.
- Do not weaken validation or skip reproducibility metadata.
