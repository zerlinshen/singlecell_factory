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

## Architecture In 60 Seconds

`singlecell_factory` is the upstream execution and run-truth surface for the
single-cell pipeline. It owns heavy compute, pipeline module execution,
checkpoint/resume, run manifests, module status, source AnnData outputs, and
remote R/reporting work when the object is too large for local handling.

The adjacent `multiomics_r_factory` is the canonical R source for broader
plotting/reporting helpers:

- `multiomics_r_factory/R/` owns R analysis and plotting modules.
- `multiomics_r_factory/R_bundle/` owns bundle-path helpers.
- `bridges/local_r_pipeline_macbook/` in this repo must remain symlinks into
  `multiomics_r_factory`; do not place real R files under the bridge.
- `workflow/modular/module_catalog.py` owns the single-cell module hierarchy:
  dependencies, architectural layers, modality tags, ownership, and bridge-ready
  flags. Pipeline compatibility constants are derived from this catalog.
- `scripts/validate_nc2024_architecture_contract.py` is the no-rerun
  controller-validation smoke for current NC2024 v2 artifacts and the
  singlecell-to-multiomics bridge contract.

## Cross-Repo Bridge (Bundle v2.1)

The `singlecell_factory -> multiomics_r_factory` bundle bridge is now at schema
`singlecell_r_bundle_v2.1` (additive, fully back-compatible with v2). The
default emitter writes v2.1; force legacy with
`scripts/export_singlecell_r_bundle.py --schema-version v2`. The R reader
accepts both via `ACCEPTED_V2_SCHEMAS`.

| Extension | Producer (Python) | Reader (R) | Contract test |
|---|---|---|---|
| `protein` | `scripts/export_singlecell_r_bundle.py::maybe_export_protein` | `multiomics_r_factory/R/protein_module.R::load_protein_extension` | `tests/test_r_bundle_contract.py` |
| `spatial` | `scripts/export_singlecell_r_bundle.py::maybe_export_spatial` | `multiomics_r_factory/R/spatial_module.R::load_spatial_extension` | `tests/test_r_bundle_contract.py` |
| `multimodal_obsm` (EXPERIMENTAL) | `scripts/export_singlecell_r_bundle.py::maybe_export_multimodal_obsm` | `multiomics_r_factory/R/integration_module.R::load_multimodal_extension` | `tests/test_r_bundle_contract.py` |

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

## Current State (2026-04-26)

The canonical NC2024 NSCLC outputs are the **v2** runs at
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

For NC2024 full-cohort work, the canonical (post-fix) results are the **v2** runs:
`results/nc2024_tumor_20260426_v2/` and `results/nc2024_bh_20260426_v2/`,
with ledger entries `ops/run_ledger/nc2024_*_v2_*.json` and the methodology
audit suite under `ops/nc2024_methodology_audit/` (`AUDIT_2026-04-26_v2.md`,
`ALIGNMENT_REPORT_v2_2026-04-26.md`, `SEGFAULT_TRACE_2026-04-26.md`,
`SMALL_REAL_VALIDATION_2026-04-26.md`). Publication-ready citation patterns
live in `docs/PUBLICATION_READY.md`. The v1 dirs (`nc2024_*_20260426/` without
the `_v2` suffix) shipped with a known annotation labeling bug and must not be
cited as scientific evidence.

Do not treat the `2026-04-25` `NC2024_NSCLC_FULL_COHORT_SPARSE_EXACT_REAL_AUTO_*`
probe directories as canonical success unless `final_adata.h5ad`,
`run_manifest.json`, and `module_status.csv` are present.

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
  inspect upstream run truth here first, then work in `multiomics_r_factory`.
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
- Do not route canonical upstream fixes into `multiomics_r_factory` first.

## Fast Links

- `AGENTS.md` - Codex/root project rules.
- `CLAUDE.md` - Claude Code project rules.
- `PROTOCOL.md` - deep pipeline and bridge runbook.
- `README.md` - long-form project reference.
- `ops/before_every_run/LATEST.md` - current run memory.
- `workflow/modular/` - modular pipeline code.
- `scripts/` - launchers, audits, verification helpers.
- `bridges/local_r_pipeline_macbook/` - symlink bridge into `multiomics_r_factory`.
