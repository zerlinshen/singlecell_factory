---
name: project-onboarding
description: Project-specific Kimi onboarding for the singlecell_factory upstream execution workspace.
---

# singlecell_factory Kimi Onboarding

Use this project skill when Kimi starts inside
`/home/zerlinshen/singlecell_factory` or handles upstream single-cell workflow
tasks.

## Required First Reads

1. `AI_AGENT_PROTOCOL.md`
2. `README.md`
3. `ops/before_every_run/LATEST.md`
4. Relevant run directory `run_manifest.json` and `module_status.csv`

If the `bioinformatics-workflow` Kimi plugin is installed, use its read-only
tools to read protocol and run-memory files before considering direct shell
inspection.

## Workspace Contract

- This project is the upstream source of truth for single-cell execution,
  recovery, and run evidence.
- Do not treat downstream R/report output as replacing upstream run truth.
- Before launching, resuming, or rerunning a pipeline, apply the
  `before-every-run` memory protocol.
- Verify with concrete files: `run_manifest.json`, `module_status.csv`,
  expected `*.h5ad` outputs, logs, and module artifacts.
- Keep bulky failed or superseded outputs governed by the retention and cleanup
  skills; do not delete raw data, references, or current best-run evidence.

## Operation Modes

- Mac-led mode: use the Mac as the command center and operate this project over
  SSH on `ubuntu-tail`.
- Direct remote mode: operate directly from this project root after logging
  into `ubuntu-tail`.

Both modes are valid. In either mode, preserve the same run-memory and
verification requirements.
