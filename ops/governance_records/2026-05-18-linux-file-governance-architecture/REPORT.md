# Linux File Governance Architecture Report

Control-plane governance record only. This is not scientific evidence.

## Objective

Clarify the Linux file architecture for `/home/zerlinshen` without modifying
the bioinformatics pipeline, launching runs, deleting artifacts, or moving
canonical data.

## Baseline Observed

- `/home/zerlinshen/projects` contains two governed project roots:
  `nc-reproduction` and `wave5-trevino`.
- Both projects already have `project.yaml`, `README.md`, and
  `ledger/project_retention_policy.yaml`.
- `/home/zerlinshen/singlecell_factory/results` is about `27M` and contains
  legacy/smoke remnants, including `small_real_validate_20260426` and small
  `test_20260516_*` directories.
- `/home/zerlinshen/projects` is about `28G`, with `wave5-trevino` at about
  `27G` and `nc-reproduction` at about `1.6G`.
- The remote factory working tree was already dirty before this round, with
  many tracked and untracked changes unrelated to file governance.

## Strategy Chosen

Use safe documentation and index updates, not artifact migration.

```text
factory repos = tools and control plane
project roots = scientific data and run evidence
run directories = immutable evidence bundles
human-facing = conclusions, reports, figures, caveats
agent-facing = manifests, logs, module status, validation JSON
legacy/smoke = old factory results unless explicitly promoted
```

## Files Added Or Updated

- `/home/zerlinshen/projects/README.md`
- `/home/zerlinshen/projects/FILE_GOVERNANCE.md`
- `/home/zerlinshen/projects/nc-reproduction/runs/README.md`
- `/home/zerlinshen/projects/wave5-trevino/runs/README.md`
- `/home/zerlinshen/singlecell_factory/docs/LINUX_FILE_GOVERNANCE.md`
- `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-linux-file-governance-architecture/REPORT.md`

## Source-Of-Truth Pointers Preserved

- NC2024:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`
- Trevino:
  `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88`

## Non-Actions

- No `workflow/modular/` files changed.
- No Python or R pipeline source files changed.
- No project artifacts moved.
- No data deleted.
- No pipeline run launched.
- No retention cleanup performed.

## Verification Plan

- Confirm every touched file exists and is nonempty.
- Confirm required governance terms appear in the new docs.
- Confirm canonical source-of-truth run directories still exist.
- Confirm `git status --short` for pipeline implementation paths did not change
  due to this round.
- Confirm the final tree is understandable from both
  `/home/zerlinshen/projects/README.md` and
  `/home/zerlinshen/singlecell_factory/docs/LINUX_FILE_GOVERNANCE.md`.

## Verification Evidence

timestamp_local: 2026-05-18T17:55+08:00

Passed:

- `test -s` for all six touched governance files.
- Required terms found across project and factory governance docs:
  `singlecell_factory`, `multiomics_r_factory`,
  `/home/zerlinshen/projects`, `human-facing`, `agent-facing`,
  `legacy/smoke`, and `do not modify pipeline`.
- Source-of-truth directories exist:
  - `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`
  - `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88`
- `git -C /home/zerlinshen/singlecell_factory diff --check -- docs/LINUX_FILE_GOVERNANCE.md ops/governance_records/2026-05-18-linux-file-governance-architecture/REPORT.md`
  returned cleanly.
- `git status --short -- docs/LINUX_FILE_GOVERNANCE.md ops/governance_records/2026-05-18-linux-file-governance-architecture/REPORT.md`
  shows only the two expected untracked factory governance files.
- Deployment used `rsync` from a local staging directory containing only:
  - `projects/README.md`
  - `projects/FILE_GOVERNANCE.md`
  - `projects/nc-reproduction/runs/README.md`
  - `projects/wave5-trevino/runs/README.md`
  - `singlecell_factory/docs/LINUX_FILE_GOVERNANCE.md`
  - this report

Residual risk:

- The remote factory working tree was dirty before this round. This report does
  not claim that pre-existing pipeline or module changes are clean; it only
  records that this governance round did not add pipeline changes.

## Review Cycle 1 Closure

timestamp_local: 2026-05-18T18:00+08:00

Code-review returned `COMMENT` and architecture review returned `WATCH` because
the project root READMEs did not link the shared project governance policy.

Doc-only correction applied:

- `/home/zerlinshen/projects/nc-reproduction/README.md` now links
  `/home/zerlinshen/projects/FILE_GOVERNANCE.md`.
- `/home/zerlinshen/projects/wave5-trevino/README.md` now links
  `/home/zerlinshen/projects/FILE_GOVERNANCE.md`.

No pipeline files, run outputs, or retention state were changed for this review
cycle.
