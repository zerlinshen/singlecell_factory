---
name: before-every-run
description: Persistent preflight and post-run journaling for remote work on `ubuntu-tail`. Use before every remote pipeline execution, rerun, resume, recovery, benchmark, or large transfer related to `singlecell_factory` or the NC2024 reproduction. Read the latest prior run notes first, preferably the full journal when manageable, then after the run update the shared local/remote log with problems, actions, resolutions, cautions, and improvement ideas.
---

# Before Every Run

Use this skill whenever a task is about to launch or resume a **remote** job on `ubuntu-tail`.

This skill exists to stop the same failure from happening twice. The rule is simple:

- before every remote run, read the previous notes
- after every remote run, write back what happened

## Canonical Locations

### Remote

- root:
  - `/home/zerlinshen/singlecell_factory/ops/before_every_run`
- latest summary:
  - `/home/zerlinshen/singlecell_factory/ops/before_every_run/LATEST.md`
- journal entries:
  - `/home/zerlinshen/singlecell_factory/ops/before_every_run/journal/`

### Local mirror

- root:
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run`
- latest summary:
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run/LATEST.md`
- journal entries:
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run/journal/`

## Required Workflow

## 1. Before launching a remote task

Always read:

- remote `LATEST.md`
- at least the newest journal entry

Prefer reading the **full journal** when:

- the task is part of the same workflow family
- the task touches the same data or modules
- there have been multiple recent failures
- you are about to do another full-cohort run

Minimum acceptable behavior:

- read the latest summary
- read the latest entry
- state the relevant carry-forward lessons before running anything

## 2. While the run is active

Track:

- run name or run directory
- launch command or controlling script
- first real blocker
- recovery actions taken
- any non-obvious behavior worth preserving

Do not wait until much later to reconstruct this from memory.

## 3. After the run finishes or fails

Update the remote journal first, then refresh the local mirror.

Every new run note must record:

- task or run name
- date/time
- objective
- what was attempted
- what succeeded
- what failed
- root cause or best current diagnosis
- what changed
- what remains risky
- what the next operator should remember
- which artifacts now look:
  - `canonical`
  - `evidence-only`
  - `superseded`
  - `failed exploratory`

Also identify obvious cleanup candidates conservatively:

- suggest deletion only for clearly superseded temp/partial artifacts
- never suggest deleting raw/reference data
- never suggest deleting canonical successful runs
- keep at least one evidence-rich failed run if it explains the winning fix

Then update `LATEST.md` so the next run can start with the newest truth quickly.

## 4. If the run repeated an old problem

Say so explicitly.

Use wording like:

- `Regression of previous issue`
- `Same failure family as prior lazy-backed matrix problem`
- `Previously solved; reappeared after new path`

That makes pattern detection much easier across sessions.

## Journal Quality Rules

- Prefer concrete paths over vague descriptions.
- Prefer real run IDs over “the latest run”.
- Separate:
  - observed facts
  - inferred diagnosis
  - fixes applied
  - unresolved risks
- If a failure was superseded later, say which later run proved the fix.
- If a fix is only a workaround, say so directly.

## Keep/Ignore Guidance

The journal should mention when something is now canonical and when something is outdated.

Mark artifacts clearly as one of:

- `canonical`
- `evidence-only`
- `superseded`
- `failed exploratory`

## Promotion Rule

Use the journal for short-term operational memory.

If a lesson remains true across **2 or more runs**, promote it into the project wiki under:

- `/Users/zerlinshen/Downloads/1. Codex/.omx/wiki/`

Good candidates for promotion:

- canonical execution lane rules
- repeated failure families
- stable controller/orchestration lessons
- cleanup and preservation rules that are no longer one-off

This helps later cleanup and prevents new agents from anchoring on stale runs.

## Entry Format

Use the template in:

- [references/journal-template.md](references/journal-template.md)

Use the helper scripts when helpful:

- [scripts/sync_before_every_run.sh](scripts/sync_before_every_run.sh)
  - refresh the local mirror from the remote canonical folder
- [scripts/new_run_note.py](scripts/new_run_note.py)
  - scaffold a new journal entry from the standard template

For the latest summary, keep it shorter:

- current canonical run(s)
- newest major lesson
- current top unresolved risk
- immediate next recommended action

The canonical Shenxin routing contract is documented in:

- [references/shenxin-routing-contract.md](references/shenxin-routing-contract.md)

## When This Skill Pairs Well With Others

- Start with `before-every-run`, then follow the Shenxin routing contract.
- Use with `singlecell-remote-workflow` for remote analysis plus local handoff.
- Use with `execute-and-recover-pipeline` when debugging or resuming runs.
- Use with `remote-run-report-bridge` after a successful remote run that needs local reporting.
- Use with `review-and-validate-quality` after a workflow/result change that should be treated as stable.
- Use with `test-run-cleanup-policy` after repeated benchmarks, superseded runs, or result-harvested experiments.
- Use with `post-run-failure-cleanup` after failed rounds that left behind stale partial artifacts.

## Non-Negotiable Rule

Do not launch a new remote run blind.

At minimum, read the previous note first.

## Execution Mode Defaults For NC2024-Style Runs

When the task is about full-cohort NC2024-style single-cell work, state the execution mode explicitly before running:

- `execution_mode = debug_massive`
- `execution_mode = controller_validation`
- `execution_mode = benchmark`

Default choices:

- `debug_massive`
  - direct `massive`
  - stop at first failing module
  - patch only that module
  - rerun direct `massive`
- `controller_validation`
  - use official orchestration path
  - expect `large` to behave as a capacity probe
  - verify truthful fallback to `massive`
- `benchmark`
  - use only when intentionally comparing fidelity/performance lanes

Carry-forward defaults for this project:

- prefer sparse/lazy preservation
- do not rerun prepare if the canonical prepared input is still valid
- treat `large` as a capacity probe on the full cohort unless there is an explicit reason not to

## Draft Outputs

The post-run hook can create automatic first-pass drafts under:

- local:
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run/journal/drafts/`
- remote:
  - `/home/zerlinshen/singlecell_factory/ops/before_every_run/journal/drafts/`

These are scaffolds, not final journal entries. Promote useful content into the real journal and `LATEST.md`.

When the post-run hook detects explicit success evidence, it may also update the bounded `Latest Hook-Confirmed Success` section in `LATEST.md`. Keep the manual canonical sections above it as the higher-authority operator summary.
