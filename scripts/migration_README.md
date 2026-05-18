# Migration: Factory Output → Project Layout

Parent plan: `/home/zerlinshen/.omc/plans/factory-project-separation.md` (§4, §4a)

This document describes the two-phase workflow for migrating existing scientific
data from `singlecell_factory/output/` and `r_multiomics_factory/output/` into
the new `${PROJECTS_ROOT}/<project-id>/runs/<run-id>/` layout.

**No data is moved by any script in this directory without explicit human approval
and removal of the `SAFETY_DISABLE` code constant.**

---

## Two-Phase Workflow

### Phase 1 — Inventory (read-only)

```bash
python scripts/migration_inventory.py
```

Walks both `output/` trees. Emits:

- `/home/zerlinshen/.omc/state/migration_manifest.csv` — one row per file with:
  `src_path, size_bytes, mtime, project_id_guess, run_timestamp_guess, target_relative_path`
- Summary printed to stdout: total files, total GiB, count per project.

**Review the CSV.** Correct any wrong `project_id_guess` values before proceeding.
Rows with unknown project land in `_unsorted/` for triage.

### Phase 2 — Dry-run Apply

```bash
python scripts/migration_apply.py
```

Reads the CSV and prints what *would* be moved — no data is touched.

Review the output. When satisfied:

```bash
MIGRATION_APPROVED=1 python scripts/migration_apply.py \
    --execute \
    --i-have-reviewed-the-dry-run
```

This runs all pre-flight checks (see below). If pre-flights pass, it prints a
SAFETY_DISABLE notice and exits 0. **Actual movement is still disabled** until a
human removes `SAFETY_DISABLE = True` from `migration_apply.py` and re-runs.

---

## Pre-Flight Gates

All five gates run before any data movement attempt. A single failure aborts.

| Gate | What it checks | Why |
|------|---------------|-----|
| **Device check** | `df --output=source` on src and dst parent | Determines `mv` vs rsync two-phase strategy |
| **Free-space gate** | `df --output=avail` >= 1.5 × total src size | Prevents out-of-space mid-migration |
| **Projects-root collision** | `${PROJECTS_ROOT}` empty or absent | Prevents overwriting an existing project tree |
| **Dirty-tree gate** | `git status --porcelain` on both factory repos | Ensures factory SHAs recorded in manifests are unambiguous |
| **Contracts parity** | `sha256sum` of both `contracts/bundle_schema.yaml` | Ensures Python and R factories agree on schema before migration |

Override `${PROJECTS_ROOT}` (default `/home/zerlinshen/projects`) via env var:

```bash
PROJECTS_ROOT=/mnt/data/projects MIGRATION_APPROVED=1 python scripts/migration_apply.py \
    --execute --i-have-reviewed-the-dry-run
```

---

## Safety Layers

Three independent safety layers prevent accidental execution:

1. **`SAFETY_DISABLE = True`** in `migration_apply.py` — code-level gate; requires
   a human to edit the file and understand the consequences.
2. **`--execute --i-have-reviewed-the-dry-run`** — two CLI flags that must both be
   present; neither alone is sufficient.
3. **`MIGRATION_APPROVED=1`** env var — must be set in the shell; cannot be passed
   as a CLI argument.

All three must be satisfied simultaneously. This is intentional: automated tooling
(Claude, scripts, CI) can never accidentally trigger execution.

---

## After Execution

Once data is moved (after SAFETY_DISABLE is removed):

1. For each migrated run, `manifest.json` is written by `manifest_writer.py`
   recording `original_path`, `migrated_at`, factory SHAs, and `bundle_sha256`.
2. Legacy `output/<thing>` paths are replaced by wrapper scripts that log access
   events to `<project-root>/runs/<run-id>/logs/legacy_access.log`.
3. Shims are retired after a 14-day silent window (zero logged accesses).

See plan §4b for the full step-by-step sequence.
