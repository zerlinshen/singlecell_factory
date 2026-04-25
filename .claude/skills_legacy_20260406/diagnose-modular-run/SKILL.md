---
name: diagnose-modular-run
description: Diagnose failed or suspicious modular pipeline runs. Use when a run fails, module_status contains errors, expected artifacts are missing, or resume strategy is needed.
argument-hint: [run-dir-or-project]
---

Diagnose modular run failures with evidence-first triage.

1. Locate the target run directory.
- If the user gives a run dir, use it directly.
- Otherwise inspect `results/` and choose the newest timestamped run.

2. Read canonical run records first.
- `<run_dir>/module_status.csv`
- `<run_dir>/run_manifest.json`
- `<run_dir>/.checkpoints/` (if present)

3. Determine failure class.
- Dependency/setup failure: missing input matrix, missing optional package, API/network dependency.
- Data validity failure: required columns absent in `adata.obs`/`adata.obsm`, empty marker sets, no significant genes.
- Resource failure: out-of-memory, timeout, over-aggressive parallelism.
- Module contract failure: prior module did not write outputs needed by downstream stage.

4. Produce a concrete recovery command.
- Prefer `--checkpoint --resume-from <module>` if checkpoint exists.
- Otherwise rerun from start with a smaller module set and `--parallel-workers 1`.
- If a specific optional backend failed, suggest fallback-friendly rerun settings (for example omit the failing module temporarily).

5. Report in this format.
- Failure point: module + exact message.
- Root-cause hypothesis: one primary cause, one backup cause.
- Evidence: files/fields inspected.
- Recovery command: exact CLI command.
- Prevention: one or two repository-level hardening actions (tests, guards, docs).
