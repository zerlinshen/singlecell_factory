# AGENTS.md - singlecell_factory

This file applies to `/home/zerlinshen/singlecell_factory` and all subdirectories.
It overrides higher-level guidance where rules conflict.

## Mission
- Protect biological correctness, statistical validity, and reproducibility.
- Prefer small, reversible edits with verifiable evidence.
- Preserve pipeline contracts and output compatibility unless a breaking change is explicitly requested.

## Default Workflow (Codex)
1. Use `$ralplan` for non-trivial work.
2. Execute with focused implementation (`executor` by default).
3. Run `$review-and-validate-quality` before completion.
4. For high-risk changes (pipeline contracts, output schema, mutating modules), include `$code-review` and `$security-review`.

## Skill Routing
- Pipeline run/recovery: `$execute-and-recover-pipeline`
- Module feature/integration/refactor: `$develop-and-integrate-module`
- Performance work: `$optimize-and-guard-performance`
- Large cross-domain work: `$orchestrate-large-task`
- Final quality gate: `$review-and-validate-quality`

Prefer these project skills over ad-hoc prompting when they match.

## Contract Rules (Must Keep)
- Mandatory module chain remains: `cellranger -> qc -> doublet_detection`.
- Dependency graph updates must be reflected in `workflow/modular/pipeline.py`.
- Each module writes only within its own module directory via context helpers.
- `module_status.csv` and `run_manifest.json` are source-of-truth artifacts.
- Keep output filenames stable unless a documented breaking change is requested.

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
