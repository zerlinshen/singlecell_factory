---
name: reproduce-run-retention
description: Use when a bioinformatics project is still in reproduction or test mode and older bulky run outputs may be deleted after required provenance, validation evidence, source-of-truth replacement, and cleanup records are captured.
---

# Reproduce Run Retention

Use this skill when a workspace is still in reproduction, validation, benchmark,
or test mode rather than a final long-term production archive.

The goal is to keep reproduce-stage storage honest and legible:
- keep the current source of truth
- keep the evidence needed to explain what was done
- remove superseded bulky outputs only after retention requirements are met

## When To Use

- A newer reproduce/test run has replaced an older bulky run.
- The older run is reproducible enough that keeping every full result directory
  is no longer worth the storage cost.
- The user has approved cleanup of old reproduce/test outputs.
- You need a consistent checklist before deleting prior `results/` directories.

## Required Retention Before Deletion

Do not delete an old reproduce/test run until all of the following exist in a
reachable location:

- run command or wrapper path
- execution mode
- input artifact path(s)
- run directory path
- `module_status.csv`
- `run_manifest.json`
- launch log
- final high-value outputs if they are still being cited
- file inventory or equivalent directory snapshot
- cleanup record explaining what was deleted and why
- explicit replacement source of truth

For NC2024-style work, prefer storing this in:
- `before-every-run/`
- `agent_runs/YYYY-MM-DD-<slug>/RUN.md`
- remote cleanup archive under `ops/cleanup_records/`

## Hard Never-Delete List

Never delete these through reproduce-stage retention cleanup:

- raw data
- prepared canonical input
- source code
- environment definitions
- the current retained source-of-truth final object
- currently cited local report/package assets
- the latest run memory and cleanup records

## Default Workflow

1. Confirm the project is still in reproduce/test mode.
2. Identify the current source of truth and the run(s) being superseded.
3. Capture missing provenance before deleting anything.
4. Write or update:
   - run ledger entry
   - cleanup record
   - latest run memory
5. Delete only the superseded bulky outputs.
6. Record reclaimed disk space.
7. Update local and remote README/run-memory surfaces if the source of truth changed.

## Output Contract

Every cleanup performed under this skill should report:

- what was preserved
- what was deleted
- why deletion was safe
- where the retained evidence lives
- what the new source of truth is
- how much disk space was reclaimed

## Pairing

Use with:
- `before-every-run`
- `post-run-failure-cleanup`
- `test-run-cleanup-policy`
- `reproduction-workspace-governance`

Use `reproduce-run-retention` when the cleanup question is specifically about
reproduction-stage storage governance, not just failed-run debris.
