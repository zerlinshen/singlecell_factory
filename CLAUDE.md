# singlecell_factory: Claude Project Instructions

## Scope
- Use this repository's modular workflow as the default entry point.
- Prefer minimal, targeted edits that preserve scientific reproducibility and output compatibility.
- Default research workflow: ingest user-provided papers first, then adapt open-source reference implementations, then integrate into existing/new modules with tests.

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
- `/run-modular-pipeline`
- `/diagnose-modular-run`
- `/add-pipeline-module`
- `/validate-analysis-outputs`
- `/optimize-modular-performance`
- `/integrate-paper-algorithm`
- `/triage-failure`
- `/safe-refactor`
- `/code-review-risk`
- `/dependency-upgrade`
- `/perf-regression-check`
- `/docs-sync`

Use these skills for repeatable workflows and faster context loading.
