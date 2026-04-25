---
name: integrate-paper-algorithm
description: Integrate a new algorithm from papers (PDFs) and open-source implementations into this modular scRNA-seq pipeline. Use when the user provides article paths and wants Claude to fuse the method into existing modules or add a new module.
argument-hint: [paper-path-or-folder] [target-module-or-new-module]
---

Integrate literature methods with reproducible engineering steps.

1. Ingest source materials first.
- Read all user-provided PDFs or article notes from the given path.
- Build a concise method brief: objective, inputs, outputs, assumptions, core equations/steps, hyperparameters, expected failure modes.
- Record uncertainty explicitly; do not invent missing details.

2. Harvest implementation references from open source.
- Find 1-3 credible implementations (GitHub or maintained packages).
- Capture repo URL, license, commit/tag used, and relevant files/functions.
- Prefer implementations with tests, active maintenance, and permissive license compatibility.
- Never copy large code blocks verbatim; adapt concepts into project-native implementation.

3. Decide integration strategy before coding.
- Path A (fuse into existing module): if input/output contracts already match current stage semantics.
- Path B (new module): if algorithm introduces a new stage, artifacts, or dependency edge.
- Map dependencies against `workflow/modular/pipeline.py` DAG and decide whether the stage is mutating (`MUTATING_MODULES`) or appending.

4. Implement with repository contracts.
- Keep module interface: `name` + `run(self, ctx)` and use `ctx.figure_dir` / `ctx.table_dir`.
- Update orchestration (`MODULE_DEPENDENCIES`, `_build_registry`, optional CLI/config wiring).
- Add compact metadata in `ctx.metadata` for traceability and QC.
- Preserve output naming conventions and run-manifest consistency.

5. Verify scientific and software correctness.
- Add tests for: dependency resolution, registry wiring, and algorithm-specific behavior/parity.
- Validate on a small run first (`--parallel-workers 1`, limited optional modules), then larger integration run.
- Check `module_status.csv` and `run_manifest.json` for completion and metadata plausibility.

6. Produce an audit-ready integration report.
- Paper summary and what was implemented.
- Open-source references used (URL, commit/tag, license).
- What was fused vs newly added.
- Validation evidence (tests + run outputs).
- Known limitations and next-step improvements.

7. Ask for explicit path inputs when missing.
- If article location is not provided, request the exact folder or file paths before proceeding.
