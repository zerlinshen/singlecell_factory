# 2026-04-26 Skill And README Governance

## Objective

Align local/remote README governance and remote Codex / Claude Code skill
availability after the `2026-04-25` NC2024 sparse-exact and paper-aligned
workflow changes.

## What Changed

- Clarified that both direct remote operation and Mac-led SSH orchestration are
  valid.
- Kept the Mac role as coordination, review, handoff, and report packaging.
- Kept `/home/zerlinshen/singlecell_factory` as the heavy-compute,
  source-of-truth execution, and remote-R/reporting surface.
- Documented `2026-04-25` sparse-exact / paper-aligned methodology changes and
  pointed operators to:
  `/home/zerlinshen/singlecell_factory/ops/nc2024_methodology_audit/AUDIT_2026-04-25.md`
- Marked observed `2026-04-25`
  `NC2024_NSCLC_FULL_COHORT_SPARSE_EXACT_REAL_AUTO_*` directories as
  non-canonical probes unless final run artifacts are present.
- Synchronized selected workflow/governance skills to remote Codex and Claude
  Code global/project surfaces.

## Artifact Classification

- canonical:
  - `/home/zerlinshen/singlecell_factory/ops/before_every_run/LATEST.md`
  - `/home/zerlinshen/singlecell_factory/README.md`
  - `/home/zerlinshen/singlecell_factory/AGENTS.md`
  - `/home/zerlinshen/singlecell_factory/CLAUDE.md`
  - `/home/zerlinshen/singlecell_factory/BEST_PRACTICES.md`
- evidence-only:
  - `/home/zerlinshen/singlecell_factory/ops/nc2024_methodology_audit/AUDIT_2026-04-25.md`
- failed exploratory / non-canonical until proven otherwise:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_SPARSE_EXACT_REAL_AUTO_20260425_211822`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_SPARSE_EXACT_REAL_AUTO_20260425_212725`

## Residual Risk

- No pipeline run was launched in this governance round.
- Remote Codex project-scope skill management now uses `.codex/skills/`.
  `codex_skills/` remains only as a legacy compatibility mirror for historical
  project-local skills.
