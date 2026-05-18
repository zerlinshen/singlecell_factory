# Two Project-Specific Retention And Final Backup Policies

Created: `2026-05-18T09:30:11Z`

This control-plane record captures the first concrete project-level policies under the remote factory/project governance contract.

## Projects

- Cell/Trevino reproduction: `/home/zerlinshen/projects/wave5-trevino/ledger/project_retention_policy.yaml`
- NC2024 NSCLC reproduction: `/home/zerlinshen/projects/nc-reproduction/ledger/project_retention_policy.yaml`

## Shared Rule

Each project owns its own retention rule. The default is latest validated run only, but the project policy decides when an older run is still protected, what best parameters define the source of truth, and what must be included in final backup.

## Cell/Trevino Rule

Current source of truth: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88`.

Keep the latest conditional public-resource all-figure evidence run. Do not delete older Cell runs until a Cell-specific cleanup inventory proves that all cited figure evidence, resource manifests, quality gates, and reports are preserved in the latest run or final backup. The scientific boundary remains public processed matrices/author resources, not FASTQ/fragments/BPNet exact parity.

## NC2024 Rule

Current source of truth: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.

Keep only the latest validated project run result. Older bulky project runs were already deleted after inventory. Future article-scale work must filter per sample/patient first, then merge clean smaller objects for batch-aware cohort-level annotation and label transfer.

## Verification Plan

- Parse both policy YAML files.
- Run project governance validator for both project roots.
- Run factory governance validator unit tests.
- Run `git diff --check` on governance files.
