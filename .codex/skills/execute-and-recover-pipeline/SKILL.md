---
name: execute-and-recover-pipeline
description: Execute modular pipeline runs and handle failures end-to-end with checkpoint/resume, evidence-based diagnosis, and concrete recovery actions. Use for run setup, module selection, run-time triage, failed module recovery, and post-run sanity checks.
---

Run and recover modular workflows with reproducible commands.

1. Validate run prerequisites.
- Confirm sample root and expected matrix paths.
- Choose optional modules intentionally.
- Enable checkpointing by default for non-trivial runs.

2. Execute with reproducibility.
- Build a single explicit CLI command.
- Start with conservative parallelism when debugging.
- Persist command and environment assumptions in the summary.

3. Diagnose failures with evidence-first triage.
- Inspect `module_status.csv`, `run_manifest.json`, and module artifacts.
- Identify first failing module and nearest causal error.
- Distinguish root cause from downstream cascade errors.

4. Recover instead of restarting blindly.
- Prefer `--resume-from <failed-module>` using same project/module set.
- Provide exact recovery command and expected next checks.
- If checkpoint state is invalid, justify full rerun.

5. Close run with minimal quality gate.
- Confirm all modules reach expected status.
- Verify key output artifacts exist for requested analysis scope.
- Escalate to `/review-and-validate-quality` for full scientific/reporting checks.
