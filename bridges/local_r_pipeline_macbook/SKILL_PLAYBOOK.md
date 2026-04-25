# SKILL_PLAYBOOK.md

This is the practical guide for when to call which skill in this project.

## Core rule

If the work is about remote single-cell analysis plus remote R plotting, start with:
- `singlecell-remote-workflow`

R plotting/reporting is remote-side. Treat "local R plotting" language in old
notes as historical unless a user explicitly revives a local Mac plotting lane.

## When to call each skill

### `singlecell-remote-workflow`
Use when:
- running or checking remote `singlecell_factory`
- bridging remote outputs into remote R plots
- validating compact R bundles and remote plot outputs

### `execute-and-recover-pipeline`
Use when:
- a run failed
- a checkpoint exists and you want recovery instead of blind rerun
- you need to diagnose what completed and what is still usable

### `develop-and-integrate-module`
Use when:
- changing pipeline code
- adding or adjusting modules
- changing CLI/config behavior
- optimizing performance or memory usage

### `review-and-validate-quality`
Use when:
- deciding whether a workflow change is stable
- checking whether outputs are really complete and interpretable
- doing a quality gate before treating a run as final

### `remote-run-report-bridge`
Use when:
- a remote run has usable outputs
- you want remote R plots and a summary PDF/report package
- you need a bridge bundle rather than a giant raw object transfer

### `readme-sync-enforcer`
Use when:
- pipeline behavior changed
- fallback logic changed
- new modes, flags, or caveats were introduced
- you updated workflow logic and need the README to match

### `test-run-cleanup-policy`
Use when:
- you finished official 10x test runs
- you repeated benchmark runs
- you want to clean stale result folders while preserving raw data/reference assets

### `bioinformatics-autosteer`
Use when:
- the next step is ambiguous
- multiple skills seem relevant
- you want a low-friction nudge toward the right workflow pattern

## Recommended patterns

### Pattern A: remote run -> remote plots -> report
1. `singlecell-remote-workflow`
2. `remote-run-report-bridge`
3. `review-and-validate-quality`

### Pattern B: failed run -> diagnosis -> rerun
1. `execute-and-recover-pipeline`
2. `develop-and-integrate-module` if code changes are needed
3. `readme-sync-enforcer` if behavior changed

### Pattern C: big test dataset iteration
1. `singlecell-remote-workflow`
2. `develop-and-integrate-module`
3. `test-run-cleanup-policy`

## Practical note

The best skill is usually the one that reduces repeated manual decisions, not the one with the fanciest name.
