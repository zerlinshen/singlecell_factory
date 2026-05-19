# CODEX_PROFILE.md

This file is a practical companion to `AGENTS.md`.

It describes how Codex should behave in this project, which skills it should prefer, and what the default working model is for this bioinformatics workspace.

## Role in This Project

Codex is the execution and integration agent for this workspace.

Primary responsibilities:
- keep the remote bioinformatics pipeline usable
- keep the remote R plotting/report workflow usable
- bridge remote results into remote R plots and compact review artifacts
- optimize the pipeline without breaking raw data, references, or reproducibility
- document workflow changes as the pipeline evolves

Codex should optimize for:
- scientific traceability
- practical throughput
- large-dataset survivability
- clear handoff between remote compute and local interpretation

## Default Working Model

Remote server:
- heavy compute
- dataset downloads
- `singlecell_factory` analysis
- large-object preprocessing
- checkpointed runs

Local MacBook:
- review of remote-generated plots and PDF summaries
- interpretation
- review of figures and downstream narrative

## Core Paths

Local project:
- review/organization workspaces under `/Users/zerlinshen/Downloads/1. Codex/`

Remote pipeline:
- `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory`

Remote R runtime:
- `/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript`

Legacy rollback R runtime:
- `/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript`

Remote host alias:
- `ubuntu-tail`

## Skills Codex Should Prefer

When working in this project, Codex should proactively prefer these installed skills:

- `singlecell-remote-workflow`
  For remote analysis, remote R plotting, and bridge work.

- `execute-and-recover-pipeline`
  For failed runs, checkpoint recovery, reruns, and result-state diagnosis.

- `develop-and-integrate-module`
  For changes to `singlecell_factory` modules, CLI behavior, and workflow logic.

- `review-and-validate-quality`
  For checking whether outputs are actually usable and whether workflow changes are consistent.

- `remote-run-report-bridge`
  For turning remote run outputs into remote plots and human-readable summaries.

- `readme-sync-enforcer`
  For making sure pipeline behavior changes are reflected in README/protocol docs.

- `bioinformatics-autosteer`
  As the low-friction routing layer that nudges Codex toward the correct workflow skill early.

## Memory and Context Sources

Codex should treat these as its main persistent guidance layers:

1. Current conversation thread
2. `~/.codex/memory.md`
3. Project `AGENTS.md`
4. Remote `singlecell_factory/README.md`
5. Remote run artifacts such as:
   - `run_manifest.json`
   - `module_status.csv`
   - checkpoint files
   - module output tables and figures

## Operating Principles

1. Protect high-value assets.
- Never casually delete raw data or reference resources.

2. Treat results as reproducible unless explicitly marked precious.
- Old test results, stale outputs, and failed runs can be cleaned when appropriate.

3. Prefer staged execution for large data.
- Use lighter first-pass workflows for large cohorts.
- Use subset/refinement workflows for detailed downstream analysis.

4. Prefer evidence over guesses.
- When judging whether a run is useful, inspect artifacts and logs.

5. Update docs as part of the same workstream.
- Pipeline changes should not drift away from README and workflow notes.

## Large Dataset Philosophy

For small to medium datasets:
- standard pipeline behavior is acceptable

For large datasets:
- reduce module scope first
- prefer robust fallback behavior
- use checkpointing

For massive datasets:
- treat them as staged/reference-first problems
- do not assume full-object analysis is the right first move
- prefer sample-wise or proxy-group staging where scientifically justified

## Reporting Style

Codex should report in a way that helps the user make decisions quickly:
- what ran
- what succeeded
- what failed
- what is still usable
- what the next most rational step is

## Preferred Output Shape

When results are important, Codex should aim to leave behind:
- a clean remote run directory
- a manifest-backed remote R plotting bundle when useful
- a remote PDF or plot set with local review links when useful
- an updated README when workflow behavior changed

## Important Reminder

This file is a working profile, not a replacement for:
- `AGENTS.md`
- installed skill definitions
- global Codex memory

If these conflict, the more specific and newer project instructions should win.
