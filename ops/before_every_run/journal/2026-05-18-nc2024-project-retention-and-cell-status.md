# NC2024 latest-run retention and Cell status — 2026-05-18

Date/time UTC: `2026-05-18T09:20:48Z`
Execution mode: `storage_governance_cleanup_and_status_review`

## Objective

Apply user-approved cancer/NC2024 latest-run retention and check current Cell/Trevino reproduction status.

## What succeeded

- Updated README, AI agent protocol, and remote governance docs with reproduction run retention policy.
- Replaced stale NC2024 canonical references in `AI_AGENT_PROTOCOL.md` that pointed to deleted factory `results/nc2024_*_v2` outputs.
- Deleted older NC2024 project run directories:
  - `2026-05-14T0542Z-d192836`
  - `2026-05-15T1815Z-wave4-200k-preflight`
  - `2026-05-15T1833Z-wave4-800k-clustering`
- Kept current source of truth:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/`
- Reclaimed/deleted target bytes: `142437966218`.
- Post-cleanup governance validation: `pass`, severity `{}`, one run inspected.
- Cell/Trevino status reviewed from `all_figures_final_quality_gate_review.json`: decision `conditional`, expected resource gaps only, no unexpected false checks, missing paths, or hash mismatches.

## Artifact classification

canonical:
- `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/`
- `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/all_figures_final_quality_gate_review.json`

evidence-only:
- `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-project-latest-run-retention/`
- `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-cell-reproduction-status/`

superseded/deleted:
- older NC2024 project run directories listed in the cleanup manifest

failed exploratory:
- none in this cleanup/status pass

## Next operator should remember

For NC2024/cancer reproduction, keep the latest validated project run result only. Do not recreate old factory `results/NC2024*` outputs or preserve superseded bulky project runs after provenance and cleanup records exist.
