---
name: reproduction-workspace-governance
description: Use when a bioinformatics reproduction workspace is becoming cluttered and needs a stable split between human-facing deliverables, figure/comparison/ROI evidence, and agent-facing run memory, handoff, and per-run provenance. Use before or after report packaging, remote-run handoff, README cleanup, or skill/playbook updates that change how future agents should navigate the workspace.
---

# Reproduction Workspace Governance

## Purpose

Keep reproduction workspaces legible for both humans and agents without moving or deleting canonical artifacts by default.

The core split is:
- human-facing evidence: final reports, figures, figure comparisons, fidelity/gap summaries, and ROI topic pages
- agent-facing memory: run parameters, remote paths, commands, outcomes, issues, cleanup notes, and next-handoff instructions

## Default Workflow

1. Identify the workspace root and nearest scoped `AGENTS.md`, `README.md`, `CODEX_PROFILE.md`, and `SKILL_PLAYBOOK.md`.
2. Preserve canonical artifacts in place unless the user explicitly approves a migration.
3. Add or refresh these entry surfaces:
   - `00_HUMAN_START_HERE.md`
   - `human_review/README.md`
   - `human_review/figures/README.md`
   - `human_review/comparisons/README.md`
   - `human_review/roi/README.md`
   - `agent_runs/README.md`
   - `agent_runs/TEMPLATE_RUN_README.md`
4. Create or update one per-run record under `agent_runs/YYYY-MM-DD-<slug>/RUN.md` for the current workstream.
5. Update the root `README.md`, scoped `AGENTS.md`, `CODEX_PROFILE.md`, and `SKILL_PLAYBOOK.md` so the new split is discoverable.
6. If remote execution or remote-side workflow behavior changed, also update `before-every-run` memory and the remote README.
7. Run the bundled checker:
   - `python3 /Users/zerlinshen/.codex/skills/reproduction-workspace-governance/scripts/check_reproduction_workspace.py <workspace-root>`

## Governance Rules

- Do not make the root directory a dumping ground for every PDF, plot, CSV, log, and scratch memo.
- Do not force humans to read agent handoff logs before seeing the current result.
- Do not force agents to infer run provenance from human-facing report prose.
- Prefer indexes and links over disruptive moves during an active reproduction.
- Keep old artifacts reachable; label them as first-pass, historical, superseded, or archive material rather than silently hiding them.
- For ROI topics such as `ELF3`, create one short topic page that links evidence, open questions, and next analyses.
- For each meaningful remote round, record data source, execution mode, command or wrapper, output path, outcome, validation evidence, and residual risk.

## When To Read More

Read `references/nc-reproduction-layout.md` when designing a new NC2024-style layout or explaining why the human/agent split exists.
