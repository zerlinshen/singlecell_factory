# NC2024 architecture validation — 2026-05-18

Task: validate the new factory/project governance architecture using the NC2024 project dataset.
Date/time UTC: `2026-05-18T08:57:43Z`
Execution mode: `architecture_validation` (read-only project validation, not a new pipeline run)

## Objective

Verify that NC2024 can be governed under the new architecture:

- factory owns contracts, validators, docs, and governance records
- project root owns run/evidence/scientific artifacts
- governance reports remain control-plane only
- no new scientific outputs are written into the factory or project run directories during validation

## What was attempted

- Read prior before-every-run memory and latest real journal entry.
- Validated `/home/zerlinshen/projects/nc-reproduction` with `scripts/validate_project_governance.py`.
- Created factory-side control-plane report under:
  `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-architecture-validation/`
- Audited factory legacy NC2024 outputs under `singlecell_factory/results/`.

## What succeeded

- Project governance validator status: `pass`.
- Severity counts: `{'info': 2, 'warning': 1}`.
- Three NC project runs were inspected.
- Two runs have root `manifest.json`; the 200k preflight run is legacy/partial and missing root manifest.
- The 800k clustering run has Python `run_manifest.json` and `module_status.csv`.
- No governance report was written under `/home/zerlinshen/projects/nc-reproduction`.
- `ops/run_records` is absent after the prior cleanup; new governance records live under `ops/governance_records`.

## What remains risky

- The factory is not fully pure yet: 48 legacy NC2024 paths remain under `singlecell_factory/results/`, totaling `30667532081` bytes.
- These are pre-project-architecture scientific/test outputs. They were not deleted in this validation pass.
- Latest older NC2024 references in `ops/before_every_run/LATEST.md` may still mention factory `results/` as canonical and should be updated after retention cleanup.

## Artifact classification

canonical:
- `/home/zerlinshen/projects/nc-reproduction/project.yaml`
- `/home/zerlinshen/projects/nc-reproduction/runs/`
- `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-architecture-validation/REPORT.md`
- `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-architecture-validation/architecture_validation.json`

evidence-only:
- `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-architecture-validation/factory_legacy_nc2024_results.txt`

superseded / cleanup candidates:
- `singlecell_factory/results/NC2024*`
- `singlecell_factory/results/nc2024_*`

failed exploratory:
- none in this validation pass

## Next operator should remember

Treat new NC2024 work as project-owned. Do not create new NC2024 scientific outputs under `singlecell_factory/results/`. Run a separate retention cleanup for the 48 legacy factory result paths only after confirming no local/report references still require them.
