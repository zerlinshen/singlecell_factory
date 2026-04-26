# singlecell_factory: Claude Project Instructions

## Authority / Read This First

Start with `AI_AGENT_PROTOCOL.md` for onboarding, read order, and task routing.
Then follow this file as the Claude Code runtime contract. Use `PROTOCOL.md`
for deep operational steps.

## Project Precedence (Overrides Global OMC Defaults)
- This project's rules take precedence over global `~/.claude/CLAUDE.md` orchestration defaults.
- If global OMC behavior conflicts with project safety/reproducibility constraints, follow this file.
- Default to single-agent, minimal, targeted edits for optimization, bugfixes, and parameter tuning.
- Multi-agent orchestration is allowed for new module creation and cross-module major refactors when it improves delivery quality.
- For multi-agent work, require an explicit plan and clear ownership boundaries before parallel execution.
- For high-impact changes (pipeline contracts, `adata` structure, output schema), require verification evidence before claiming completion.
- Always prioritize reproducibility, deterministic outputs, and stable pipeline contracts over speed.

## Scope
- Use this repository's modular workflow as the default entry point.
- Prefer minimal, targeted edits that preserve scientific reproducibility and output compatibility.
- Default research workflow: ingest user-provided papers first, then adapt open-source reference implementations, then integrate into existing/new modules with tests.

## Operating Modes
- Direct remote operation inside `/home/zerlinshen/singlecell_factory` is valid
  for heavy compute, pipeline execution, run triage, and remote R/reporting.
- Mac-led SSH orchestration is also valid when the task is coordination,
  review, handoff, or local report packaging.
- In both modes, canonical run truth comes from remote artifacts:
  `run_manifest.json`, `module_status.csv`, `ops/before_every_run/LATEST.md`,
  and `ops/run_ledger/`.
- Do not treat the observed `2026-04-25` sparse-exact probe directories as
  successful canonical runs unless `final_adata.h5ad`, `run_manifest.json`, and
  `module_status.csv` are present.

## Core Commands
- Run pipeline:
  `python -m workflow.modular.cli --project <name> --sample-root <path> --optional-modules <modules>`
- Run with recovery:
  `python -m workflow.modular.cli ... --checkpoint`
- Resume:
  `python -m workflow.modular.cli ... --checkpoint --resume-from <module>`
- Run tests:
  `pytest -q`
- Run focused modular tests:
  `pytest -q tests/test_modular.py tests/test_modular_optimizations.py`

## Pipeline Contract
- Mandatory stages are `cellranger -> qc -> doublet_detection`.
- Optional stages are dependency-resolved from `workflow/modular/pipeline.py`.
- Each module writes only to its own output subdirectory via `ctx.set_module_dir(...)`.
- Module status and manifest are source-of-truth outputs:
  `module_status.csv`, `run_manifest.json`.

## Change Rules
- **Bridge symlink rule**: `bridges/local_r_pipeline_macbook/R` and `bridges/local_r_pipeline_macbook/R_bundle` must always be symlinks pointing to `multiomics_r_factory/R` and `multiomics_r_factory/R_bundle` respectively. Never place real R files under these bridge paths. Edit R sources in `multiomics_r_factory/` only. Verify with `bash scripts/ci/check_bridge_symlink.sh`.
- When adding a module, update all required integration points:
  - `workflow/modular/modules/<module>.py`
  - `MODULE_DEPENDENCIES` and `_build_registry()` in `workflow/modular/pipeline.py`
  - CLI/config wiring if new parameters are introduced
  - tests and README/PROTOCOL when user-facing behavior changes
- If a module mutates `adata` structure (embeddings/layers/main matrix), evaluate whether it must be in `MUTATING_MODULES`.
- Keep output filenames stable unless explicitly performing a documented breaking change.
- For paper-driven integrations, always record source paper paths, reference repo URLs, commit/tag, and license in the final summary.
- Do not paste large external code blocks directly; re-implement/adapt to this codebase style and contracts.

## Skills In This Project
- `/orchestrate-large-task`
- `/execute-and-recover-pipeline`
- `/develop-and-integrate-module`
- `/review-and-validate-quality`
- `/optimize-and-guard-performance`
- `/download-large-file`
- `/before-every-run`
- `/singlecell-remote-workflow`
- `/remote-run-report-bridge`
- `/reproduction-workspace-governance`
- `/readme-sync-enforcer`
- `/post-run-failure-cleanup`
- `/reproduce-run-retention`
- `/singlecell-factory-module-delivery`

Use these skills as a composition set, not as isolated tools.

Default combo routing:
- Large cross-domain work: `/orchestrate-large-task` then delegated specialized skills.
- Run/triage/recovery work: `/execute-and-recover-pipeline` -> `/review-and-validate-quality`.
- Feature/module work: `/develop-and-integrate-module` -> `/review-and-validate-quality`.
- Performance-sensitive changes: add `/optimize-and-guard-performance` before final quality gate.
- Remote execution or resume: `/before-every-run` -> `/singlecell-remote-workflow` -> `/execute-and-recover-pipeline`.
- Remote result to local package: `/remote-run-report-bridge` -> `/reproduction-workspace-governance` -> `/review-and-validate-quality`.
- Governance/doc drift: `/reproduction-workspace-governance` -> `/readme-sync-enforcer`.

Consolidation map (legacy -> active):
- `run-modular-pipeline`, `diagnose-modular-run`, `triage-failure` -> `/execute-and-recover-pipeline`
- `add-pipeline-module`, `integrate-paper-algorithm`, `safe-refactor`, `dependency-upgrade` -> `/develop-and-integrate-module`
- `code-review-risk`, `validate-analysis-outputs`, `docs-sync` -> `/review-and-validate-quality`
- `optimize-modular-performance`, `perf-regression-check` -> `/optimize-and-guard-performance`
- `download-large-file` -> `/download-large-file`

Skill locations:
- Global Claude reusable skills: `~/.claude/skills/`
- Project Claude skills: `.claude/skills/`
- Project Codex skills: `.codex/skills/`
- Legacy Codex compatibility mirror: `codex_skills/`
