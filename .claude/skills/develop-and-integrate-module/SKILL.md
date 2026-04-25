---
name: develop-and-integrate-module
description: Implement or evolve modular analysis capabilities, including new modules, paper-method integration, safe refactors, and dependency upgrades, while preserving pipeline contracts and reproducibility. Use for feature development in workflow/modular and related config/CLI/tests wiring.
---

Deliver module-level changes with full integration safety.

1. Define scope and contract before edits.
- State module/API behavior and expected inputs/outputs.
- Mark backward-compatibility constraints and acceptable breakage.

2. Implement with full orchestration wiring.
- Update module code under `workflow/modular/modules/`.
- Update registry and dependency graph in `workflow/modular/pipeline.py`.
- Update config/CLI plumbing when introducing parameters.
- Evaluate whether mutating behavior requires `MUTATING_MODULES` updates.

3. Integrate external methods safely (paper sub-workflow).
- Ingest paper: extract objective, key algorithm, assumptions, parameters, failure modes.
- Find reference implementation: locate repo URL, identify relevant source files, note commit/tag and license.
- Evaluate fit: compare reference API/data structures with pipeline contracts (AnnData obs/obsm/uns, module I/O pattern).
- Adapt, don't paste: re-implement core logic to match codebase style, use existing helpers (`_scanpy_compat`, `_gpu_utils`, `score_gene_sets`).
- Wire integration: add module to `MODULE_DEPENDENCIES`, `_build_registry()`, CLI args if needed.
- Record provenance: source paper path, reference repo URL, commit/tag, license, and adaptation notes in module docstring and `__references__`.
- Validate: run `module-contract-validator` agent, then `/scaffold-module-test` for isolated test.

4. Handle refactors and dependency upgrades with risk control.
- Keep edits minimal and scoped; avoid unrelated churn.
- Validate compatibility impact of dependency version changes.
- Preserve output schema and filenames unless explicitly changing contract.

5. Ship with verification.
- Add/update targeted tests first, then broader tests when feasible.
- Summarize behavior changes, migration notes, and residual risks.
- Route to `/review-and-validate-quality` for final quality gate.
