# AI_AGENT_PROTOCOL.md

Canonical onboarding index for AI agents entering
`/home/zerlinshen/singlecell_factory`.

## Authority / Read This First

1. Read this file first to understand where to look next.
2. Then obey the runtime-specific file for your agent:
   - Codex: `AGENTS.md`, then `CODEX.md` when Codex-specific review or
     optimization guidance is relevant
   - Claude Code: `CLAUDE.md`
3. Use `PROTOCOL.md` as the deep operational guide. It is not deprecated.
4. For any remote run, rerun, recovery, or result-truth question, read
   `ops/before_every_run/LATEST.md` before acting.

This file is an onboarding index. It does not override `AGENTS.md`,
`CLAUDE.md`, or the deep runbook in `PROTOCOL.md`.

## Architecture (2026-05+): Factory-Project Separation

As of 2026-05, the working tree follows a three-way split. Full plan:
`/home/zerlinshen/.omc/plans/factory-project-separation.md`.

**Factories are tools, projects are data.**

- **`singlecell_factory`** (this repo): pure compute tool. No scientific outputs
  land inside this tree. All artifacts write to `<project-root>/runs/<run-id>/python/`.
- **`r_multiomics_factory`** (`/home/zerlinshen/r_multiomics_factory/`): R-side
  compute tool. Writes to `<project-root>/runs/<run-id>/r/`.
- **`projects/`** (`/home/zerlinshen/projects/<project-id>/`): self-contained
  project directories. Create with `omc-new-project` from
  `/home/zerlinshen/projects-bootstrap/`.

**--project-root contract**: pass `--project-root <path>` and optional
`--run-id <id>` to any pipeline entry point. Run-id format:
`<UTC-timestamp>-<short-py-sha>`, regex
`^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{4}Z-[0-9a-f]{7}$`.

**Dirty-tree gate**: pipeline refuses if factory tree is dirty unless
`--allow-dirty` is passed; in that case `diff_sha256` is recorded in
`manifest.json`.

**Contracts**: canonical bundle schema at `contracts/bundle_schema.yaml` here;
vendored byte-identically into `r_multiomics_factory/contracts/`. Update via
`tools/sync-contracts.sh` only — never edit the vendored copy directly.

**Agent rule (MANDATORY)**: Never write outputs inside the factory tree. Always
resolve outputs from `--project-root`.

**SC_REQUIRE_PROJECT_ROOT=1 fail-fast**: Prefer setting this env var when
invoking the pipeline. Unset → `DeprecationWarning` + legacy fallback. `=1` →
hard error `sys.exit(2)`. Round-2 ADR will flip the default to required. Gate
mirrored in `scripts/export_singlecell_r_bundle.py` and
`scripts/pack_run_for_mac.sh`. Plan: `/home/zerlinshen/.omc/plans/factories-optimization-round1.md`.

### Environment Switches (Round-1a)

Plan: `/home/zerlinshen/.omc/plans/factories-optimization-round1.md`

| Variable | Unset | `=1` |
|---|---|---|
| `SC_REQUIRE_PROJECT_ROOT` | `DeprecationWarning`; legacy `output/` fallback | Hard error exit 2 |

### R-factory SHA fields (Round-1a)

Plan: `/home/zerlinshen/.omc/plans/factories-optimization-round1.md`

- `r_factory_sha_at_manifest_write` in `run_manifest.json` — resolved at pipeline end.
- `r_factory_sha_at_export` in `bundle/provenance.json` — resolved at bundle export.
- R bundle loader (`r_multiomics_factory/R_bundle/io_bundle.R`) warns (not errors) on mismatch.

### Remote Governance Control Plane

For remote factory/project governance, read
`docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md` and
`contracts/project_run_contract.yaml`. Governance validation reports written
inside `singlecell_factory` are control-plane records only; they are not
scientific outputs and must not replace project run evidence.

Before governance edits, persist branch, commit, `git status --short`, tracked
diff hash, and untracked inventory/hash in a durable
`ops/governance_records/` record. Project validation must be read-only against
`/home/zerlinshen/projects/<project-id>` unless a later plan explicitly scopes
project-root overlays.

### Multi-Cohort Annotation Policy

For NC2024-scale, multi-patient, or multi-cohort work, agents must separate
per-sample filtering from cohort-level annotation:

- per-sample/per-patient runs are valid for QC, smoke tests, and module-contract
  checks only;
- article-level annotation should be performed on a merged clean cohort,
  balanced sketch, or metacell atlas with explicit batch-aware integration;
- use `sample`, `patient`, `donor`, chemistry, or lane metadata as the batch key
  when combining cohorts;
- transfer atlas labels back to all cells before patient-level composition,
  pseudobulk, DE, or wet-lab/manuscript conclusions;
- do not present one-patient annotation as final multi-cohort evidence.

## Architecture In 60 Seconds

`singlecell_factory` is the upstream execution and run-truth surface for the
single-cell pipeline. It owns heavy compute, pipeline module execution,
checkpoint/resume, run manifests, module status, source AnnData outputs, and
remote R/reporting work when the object is too large for local handling.

The adjacent `r_multiomics_factory` is the canonical R source for broader
plotting/reporting helpers:

- `r_multiomics_factory/R/` owns R analysis and plotting modules.
- `r_multiomics_factory/R_bundle/` owns bundle-path helpers.
- `bridges/local_r_pipeline_macbook/` in this repo must remain symlinks into
  `r_multiomics_factory`; do not place real R files under the bridge.
- `workflow/modular/module_catalog.py` owns the single-cell module hierarchy:
  dependencies, architectural layers, modality tags, ownership, and bridge-ready
  flags. Pipeline compatibility constants are derived from this catalog.
- `scripts/validate_nc2024_architecture_contract.py` is the no-rerun
  controller-validation smoke for current NC2024 v2 artifacts and the
  singlecell-to-multiomics bridge contract.

## Cross-Repo Bridge (Bundle v2.1)

The `singlecell_factory -> r_multiomics_factory` bundle bridge is now at schema
`singlecell_r_bundle_v2.1` (additive, fully back-compatible with v2). The
default emitter writes v2.1; force legacy with
`scripts/export_singlecell_r_bundle.py --schema-version v2`. The R reader
accepts both via `ACCEPTED_V2_SCHEMAS`.

| Extension | Producer (Python) | Reader (R) | Contract test |
|---|---|---|---|
| `protein` | `scripts/export_singlecell_r_bundle.py::maybe_export_protein` | `r_multiomics_factory/R/protein_module.R::load_protein_extension` | `tests/test_r_bundle_contract.py` |
| `spatial` | `scripts/export_singlecell_r_bundle.py::maybe_export_spatial` | `r_multiomics_factory/R/spatial_module.R::load_spatial_extension` | `tests/test_r_bundle_contract.py` |
| `multimodal_obsm` (EXPERIMENTAL) | `scripts/export_singlecell_r_bundle.py::maybe_export_multimodal_obsm` | `r_multiomics_factory/R/integration_module.R::load_multimodal_extension` | `tests/test_r_bundle_contract.py` |

WARNING: The `multimodal_obsm` extension is EXPERIMENTAL — the R loader emits
`[multimodal_obsm extension] EXPERIMENTAL ...` at load. Embeddings published
under this extension are visualization aids only and must not be cited as
quantitative evidence for new claims.

Cross-language parity sentinels live in `tests/test_python_r_parity.py`
(top-marker Jaccard >= 0.6, PCA cosine >= 0.95). Python-to-R execution is
controlled by the `RSCRIPT_BIN` env var (default conda env
`r_multiomics_arrow`).

Both operation modes are valid:

- direct remote work inside this repo for compute, triage, remote R/reporting,
  and source-of-truth artifact inspection
- Mac-led SSH orchestration for coordination, review, handoff, and report
  packaging

In both modes, this remote repo remains the run-truth surface.

## Current State (2026-05-16)

**Wave-5 Trevino PCW21 biology-aware validation pivot (plan v4.2)** is the most recent canonical work. Run dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260516T0931Z-d192836f1bb0/`. Binding ledger: `ops/run_ledger/wave5_trevino_20260516T0931Z-d192836f1bb0.v4.2.json` (plan_revision=v4.2, validation_posture=biology-aware; schema `ops/run_ledger/schema/wave5_v4_2.schema.json`). Plan + spec live under `.omc/plans/wave5-completion-consensus-2026-05-16-v4.2.md` and `.omc/specs/deep-interview-wave5-completion.md` (both gitignored agent state). CI gate `scripts/ci/wave5_v4_2_gate.sh` exits 1 (CLOSED-PARTIAL overall: AC-VAL-3a CLOSED, AC-VAL-3b PARTIAL per §3.4, AC-CI-1 CLOSED, AC-VAL-PLOT-1/2/3 + AC-LEDGER-1 + AC-VAL-3c CLOSED). Methodology: same v3 peak-gene linkage data, comparison reference shifted from Trevino S2F string tuples (contaminated by sparse-detection artifacts MS4A12/FCRLA/SFTPC) to a SHA-pinned literature-curated PCW21 cortical marker panel at `ops/run_ledger/panels/wave5_cortical_panel_v1.json`. The session journal at `ops/before_every_run/journal/2026-05-16_wave5_trevino_pcw21_completion.md` documents the full execution arc (v3 → v4.2).

The carry-over canonical NC2024 NSCLC outputs from 2026-04-26 remain valid. The **v2** runs are at
`results/nc2024_tumor_20260426_v2/` and `results/nc2024_bh_20260426_v2/`.
Latest ledger entries live in `ops/run_ledger/nc2024_*_v2_*.json`. The matching
v1 directories (without the `_v2` suffix) shipped with a known annotation
labeling bug and must not be cited as scientific evidence.

For methodology, parameter rationale, paper alignment, and segfault fix
context, read the audit suite under `ops/nc2024_methodology_audit/`
(`AUDIT_2026-04-26_v2.md`, `ALIGNMENT_REPORT_v2_2026-04-26.md`,
`SMALL_REAL_VALIDATION_2026-04-26.md`, `SEGFAULT_TRACE_2026-04-26.md`) and the
publication template at `docs/PUBLICATION_READY.md`.

For Mac transmission of v2 bundles, see `ops/MAC_PULL_RECIPE_2026-04-26.md`.

## Source Of Truth

Prefer these over summaries or memory:

- `run_manifest.json`
- `module_status.csv`
- `ops/before_every_run/LATEST.md`
- `ops/run_ledger/`
- module output tables/figures in the relevant run directory
- launch logs and checkpoint directories for failed or recovered runs

For every project, use the project-owned latest-run retention and final-backup rule:

- keep only the latest validated project run result as the active scientific output unless the project-specific retention rule says otherwise;
- older project run directories may be deleted after a governance cleanup record captures inventory, manifest/status availability, and replacement source of truth;
- each project must record its analysis type, method family, selected/best parameters, module list, batch/integration keys, filtering thresholds, random seed, and final backup contents before old runs are removed;
- final backup should include root `manifest.json`, producer-native manifests, `module_status.csv`, launch command/log, parameter file, environment pins, final labeled object or compact atlas, figures/reports, and conclusion summary;
- never delete raw data, prepared canonical inputs, launch scripts, source code, environment definitions, governance records, or currently cited report assets;
- do not recreate or cite `singlecell_factory/results/<project>` as a canonical scientific output location. Scientific outputs belong under `/home/zerlinshen/projects/<project-id>/runs/`.

For NC2024/cancer work, the current structure-validation source of truth is `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/`; cleanup and run evidence are recorded under `ops/governance_records/` and `ops/before_every_run/`.


## Paper Reproduction Ladder

For paper reproduction, do not jump directly into the local preferred pipeline when an upstream route exists.

1. Clone or stage the upstream paper repository/method scripts and pin URL, commit/tag, DOI, data accessions, license, and environment.
2. Run raw-data reproduction first when public raw data exists and capacity allows it.
3. If raw data is absent or infeasible, start from the earliest public computable input and state that boundary.
4. Reproduce both data objects and figure panels, then write claim-level evidence.
5. Compare each claim as exact, approximate, proxy, unsupported, or resource gap.
6. Perform module-gap analysis: existing module, parameter change, new reusable module, or paper-specific script.
7. Use `$develop-and-integrate-module` / `$singlecell-factory-module-delivery` only when a paper method should become reusable in `singlecell_factory`.
8. After faithful reproduction, optimize for our context and record selected best parameters in `ledger/project_retention_policy.yaml`.

For reproduction projects, project policy must include `upstream_repository`, `raw_data_reproduction`, `data_object_reproduction`, `figure_reproduction`, `module_gap_decisions`, and `context_optimization_decisions`.

**See Also (2026-05-18, C2 reproduction methodology):** For the operational procedure (fork upstream → record provenance → run raw → identify parity gaps → route fixes per three-factory targets → record context-tuning), follow `docs/PAPER_REPRODUCTION_SOP.md`. The companion skill `paper-reproduction-from-upstream` (at `~/.claude/skills/paper-reproduction-from-upstream/SKILL.md`) provides the agent walk-through. The schema for the six policy fields above lives in `contracts/project_run_contract.yaml`, and `scripts/validate_project_governance.py` emits advisory warnings when fields are missing or malformed.

## Skill And Workflow Routing

- Before remote execution: `$before-every-run`
- Pipeline run, failure, recovery: `$execute-and-recover-pipeline`
- Module creation or module contract change: `$develop-and-integrate-module`
- Singlecell module delivery: `$singlecell-factory-module-delivery`
- Remote/local bridge or handoff: `$singlecell-remote-workflow`
- Report packaging from remote results: `$remote-run-report-bridge`
- Reproduction claim/status work: `$nsclc-baseline-reproduction`
- Cleanup after failed or superseded runs: `$post-run-failure-cleanup`,
  `$reproduce-run-retention`, `$test-run-cleanup-policy`
- Documentation or governance drift: `$reproduction-workspace-governance` then
  `$readme-sync-enforcer`
- Final quality gate: `$review-and-validate-quality`

Codex project skills live in `.codex/skills/`. `codex_skills/` is only a
legacy compatibility mirror for historical project-local skills.

## Task Routing

- If the task touches pipeline contracts, modules, AnnData structure, CLI flags,
  run recovery, or canonical outputs, start here.
- If the task is only downstream R plotting/reporting or bundle visualization,
  inspect upstream run truth here first, then work in `r_multiomics_factory`.
- If the task spans both repos, establish source-of-truth artifacts here before
  editing downstream R consumers.
- If the task changes user-facing behavior, update the relevant README/protocol
  in the same workstream.

## Verification Contract

Before claiming completion, collect evidence appropriate to the task:

- docs-only: root entrypoint grep/read checks and a concise diff summary
- module/code change: focused tests, then broader tests when risk warrants
- run/recovery: `module_status.csv`, `run_manifest.json`, logs, and recovered
  artifact paths
- bridge/R source change: verify bridge symlinks with
  `bash scripts/ci/check_bridge_symlink.sh`
- module hierarchy or bridge contract change: run
  `python3 scripts/validate_nc2024_architecture_contract.py` and focused catalog
  / bundle tests before considering the change stable

For pipeline changes, a statement without artifact paths is not evidence.

## Common Mistakes To Avoid

- Do not skip `ops/before_every_run/LATEST.md` before a remote run.
- Do not treat Mac-local summaries as more authoritative than remote artifacts.
- Do not write real R files into `bridges/local_r_pipeline_macbook/`.
- Do not demote `PROTOCOL.md`; it is still the deep operational guide.
- Do not route canonical upstream fixes into `r_multiomics_factory` first.

## Fast Links

- `AGENTS.md` - Codex/root project rules.
- `CLAUDE.md` - Claude Code project rules.
- `PROTOCOL.md` - deep pipeline and bridge runbook.
- `README.md` - long-form project reference.
- `ops/before_every_run/LATEST.md` - current run memory.
- `workflow/modular/` - modular pipeline code.
- `scripts/` - launchers, audits, verification helpers.
- `bridges/local_r_pipeline_macbook/` - symlink bridge into `r_multiomics_factory`.
