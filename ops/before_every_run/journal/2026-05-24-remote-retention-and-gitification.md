# Remote Project Retention Cleanup And Suite Gitification - 2026-05-24

## Scope

- Task: execute the approved Ralph plan combining project retention cleanup and
  pipeline/suite gitification.
- Execution mode: controller_validation for governance/retention only.
- Remote host: `/home/zerlinshen` on `ubuntu-tail`.

## What Changed

- Created project cleanup manifest:
  `/home/zerlinshen/projects/ledger/cleanup_candidates/2026-05-24-project-retention-gitification.tsv`.
- Created missing retention policies for:
  - `/home/zerlinshen/projects/round9-singlecell-comparison/ledger/project_retention_policy.yaml`
  - `/home/zerlinshen/projects/hgmm-smoke/ledger/project_retention_policy.yaml`
- Preserved metadata/evidence records under per-project
  `ledger/cleanup_records/2026-05-24-project-retention-gitification/`.
- Deleted only explicit superseded candidates with preservation records:
  - `wave5-trevino/runs/20260516T0931Z-d192836f1bb0`
  - `wave5-trevino/runs/20260517T1005Z-7b539a5`
  - `hgmm-smoke/runs/2026-05-21T1625Z-9907a7d`
  - `hgmm-smoke/runs/2026-05-21T1626Z-9907a7d`
  - `hgmm-smoke/runs/2026-05-21T1630Z-9907a7d`
- Reclaimed approximately 6.07 GB; `/home/zerlinshen/projects` dropped from
  about 104G to about 98G.
- Updated suite root docs/ignore for a lightweight governance git repository.
- Remediated retired paths in `/home/zerlinshen/projects-bootstrap` and smoke
  tested `omc-new-project`.

## Preserved / Deferred

- Preserved Wave5 source-of-truth run `20260517T1436Z-13c2c88` and linked
  pipeline run `2026-05-17T2004Z-13c2c88`.
- Preserved Round9 comparison lanes and initial evidence-only baseline.
- Preserved NG2025 data/open/staging and source-of-truth run family.
- Deferred large Wave5 2026-05-19 runs pending a dedicated Cell cleanup
  inventory.

## Verification

- Retention checks passed with manifest coverage, source-of-truth existence,
  preservation records, and non-git project/data root checks.
- Suite gate and git commit/push status are recorded in the active Codex
  session final report.
