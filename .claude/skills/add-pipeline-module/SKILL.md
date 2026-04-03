---
name: add-pipeline-module
description: Add or refactor a modular analysis stage in this codebase. Use when creating a new module, changing dependencies, wiring registry/CLI support, or updating tests for module-level behavior.
argument-hint: [module-name] [depends-on]
---

Implement module changes with full pipeline integration, not just a new file.

1. Implement the module class.
- Create/update `workflow/modular/modules/<module>.py`.
- Define `<ModuleName>Module` with `name = "<module>"` and `run(self, ctx)`.
- Validate required upstream state and raise precise `ValueError` messages.
- Write outputs only under `ctx.figure_dir` / `ctx.table_dir`.
- Update `ctx.metadata` with compact, serializable summary stats.

2. Wire orchestration.
- Add dependency edges in `MODULE_DEPENDENCIES` inside `workflow/modular/pipeline.py`.
- Register the module in `_build_registry()`.
- If the module mutates embeddings/layers/structure, add it to `MUTATING_MODULES`.

3. Wire configuration and CLI (if needed).
- Add dataclass config in `workflow/modular/config.py` when new params are required.
- Add CLI args in `workflow/modular/cli.py` and pass into `PipelineConfig`.
- Update help text listing available optional modules.

4. Add or update tests.
- Dependency resolution test in `tests/test_modular.py`.
- Registry inclusion test in `tests/test_modular.py`.
- Numerical/parity or behavior test in a focused test file when algorithmic.

5. Update docs where behavior is user-facing.
- Update module list, dependency DAG, and output structure in `README.md`.
- Update `PROTOCOL.md` only if workflow instructions or recommended commands changed.

6. Verify before finishing.
- Run targeted tests first, then broader `pytest` when feasible.
- Confirm no stale references to old module names remain.
