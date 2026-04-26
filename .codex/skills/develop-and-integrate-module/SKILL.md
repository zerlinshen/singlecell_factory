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

3. Integrate external methods safely.
- For paper-driven changes, extract objective, assumptions, parameters, failure modes.
- Re-implement/adapt to repository style; avoid raw code dumps.
- Record source paper path, reference repo URL, commit/tag, and license.

4. Handle refactors and dependency upgrades with risk control.
- Keep edits minimal and scoped; avoid unrelated churn.
- Validate compatibility impact of dependency version changes.
- Preserve output schema and filenames unless explicitly changing contract.

5. Ship with verification.
- Add/update targeted tests first, then broader tests when feasible.
- Summarize behavior changes, migration notes, and residual risks.
- Route to `/review-and-validate-quality` for final quality gate.
