---
name: post-run-failure-cleanup
description: Clean up useless artifacts after a failed remote pipeline round on `ubuntu-tail`. Use after meaningful failed runs in `singlecell_factory` to classify failed outputs, preserve canonical evidence, delete clearly superseded partial artifacts, and keep the workspace from filling with misleading stale runs.
---

# Post Run Failure Cleanup

Use this skill immediately after a meaningful **failed** remote run when the run has already produced stale, partial, or clearly superseded artifacts.

This skill is project-specific and should be treated as mandatory after failed rounds that generated non-trivial outputs.

## Read First

Before cleanup, read:

- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run/LATEST.md`
- the newest relevant run journal entry
- `/Users/zerlinshen/Downloads/1. Codex/SKILL_PROTOCOL.md`

## Preserve First

Always preserve:

- raw/reference data
- canonical prepared input
- canonical direct successful `massive` run
- canonical controller-fallback successful run
- one evidence-rich failed run if it explains the winning fix

## Cleanup Workflow

1. Identify the failed round precisely.
- exact run directory
- exact first failing module
- whether it still contains unique debugging evidence

2. Classify the round.
- `failed exploratory`
- `evidence-only`
- `superseded`

3. Decide what is safe to remove.
- delete only clearly superseded temp/partial artifacts
- keep anything still cited by current handoff, run memory, or report files

4. Clean both local and remote surfaces when relevant.
- remote:
  - stale failed run directories
  - tiny redundant launch outputs when they add no value
- local:
  - partial handoff bundles
  - temporary render artifacts
  - stale copies of superseded outputs

5. Write back the result.
- update the journal
- update `LATEST.md` if the cleanup changes what should be treated as current
- summarize:
  - what was preserved
  - what was deleted
  - reclaimed size when practical

## Required Rule

If a failed round produced heavy partial outputs and a better later run already exists, do not leave the failed round sitting around indefinitely.

Clean it or explicitly mark why it is still being preserved.
