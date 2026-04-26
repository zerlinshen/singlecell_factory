---
name: skill-usage-playbook
description: Use when deciding which installed skill should be invoked for a bioinformatics workflow step. Helps route tasks to the most effective skill with low friction.
---

Use this skill when choosing which existing skill should drive the next step.

1. Prefer `before-every-run` before any remote execution or resume on `ubuntu-tail`.
2. Prefer `bioinformatics-autosteer` when the task spans multiple skills and needs the Shenxin routing chain selected.
3. Prefer `singlecell-remote-workflow` for remote analysis, local artifact organization, and bridge work.
4. Prefer `execute-and-recover-pipeline` for failed runs, reruns, checkpoint recovery, and status diagnosis.
5. Prefer `remote-run-report-bridge` when a remote result should become a local organized report pack.
6. Prefer `review-and-validate-quality` when a workflow or result needs a quality gate before being treated as stable.
7. Prefer `test-run-cleanup-policy` after official 10x test runs, repeated benchmarks, or superseded pipeline experiments.
8. Prefer `nsclc-baseline-reproduction` when paper-level support status needs to be updated.
9. Prefer `readme-sync-enforcer` whenever code or workflow behavior changes.
10. Prefer `develop-and-integrate-module` for pipeline code changes, module additions, CLI/config changes, and performance-oriented refactors.
