---
name: architecture-governance-sync
description: Use when a project or reproduction workspace structure, workflow boundary, artifact layout, architecture, or operational contract changes and Codex must refresh architecture governance docs, update the nearest README in the same workstream, and regenerate a workspace structure diagram PNG.
---

# Architecture Governance Sync

## Purpose

Keep a workspace understandable after structural or workflow changes. Use this
skill after creating folders, moving artifacts, adding run lanes, changing
remote/local boundaries, introducing new scripts/modules, or updating handoff
contracts.

This skill complements narrower skills:
- Use `readme-sync-enforcer` for ordinary README drift after behavior changes.
- Use `reproduction-workspace-governance` for NC2024-style human/agent layout.
- Use this skill when the architecture map itself must be refreshed.

## Default Workflow

1. Identify the workspace root and nearest scoped `README.md`.
2. Inspect the current top-level structure and important governance files:
   `AGENTS.md`, `CODEX_PROFILE.md`, `SKILL_PLAYBOOK.md`, `before-every-run/`,
   `agent_runs/`, `human_review/`, `docs/`, `scripts/`, `workflow/`,
   `bridges/`, `results/`, and report/figure folders.
3. Run the bundled updater:

```bash
python3 /Users/zerlinshen/.codex/skills/architecture-governance-sync/scripts/sync_architecture_governance.py <workspace-root>
```

4. Review the generated README section and diagram before claiming completion.
5. If this is a remote workflow, update the relevant remote README or run-memory
   surface as well.
6. Report the changed README path, generated `structure.png`, and verification
   evidence.

## Output Contract

The script updates or creates:
- `README.md` with a marked `Architecture Governance` block.
- `docs/architecture/ARCHITECTURE_GOVERNANCE.md` with the same snapshot.
- `docs/architecture/structure.svg` for editable structure diagrams.
- `docs/architecture/structure.png` for human-facing visual review.

The README block is bounded by:

```markdown
<!-- ARCHITECTURE-GOVERNANCE:START -->
...
<!-- ARCHITECTURE-GOVERNANCE:END -->
```

Future runs replace only this block by default.

## Governance Rules

- Preserve canonical artifacts unless the user explicitly asks to move or
  delete them.
- Prefer indexes, README blocks, and diagrams over disruptive directory moves.
- Keep human-facing result surfaces separate from agent-facing run memory.
- Show remote/local boundaries explicitly when a workflow spans machines.
- Treat `structure.png` as a required artifact whenever architecture governance
  changes.
- Do not hide uncertainty: if a folder's role is inferred from its name, say it
  is inferred in the snapshot.

## Verification

After running the updater, verify:

```bash
test -s <workspace-root>/README.md
test -s <workspace-root>/docs/architecture/structure.png
test -s <workspace-root>/docs/architecture/ARCHITECTURE_GOVERNANCE.md
python3 /Users/zerlinshen/.codex/skills/architecture-governance-sync/scripts/sync_architecture_governance.py <workspace-root> --check
```

If `--check` reports drift, rerun without `--check`, inspect the diff or changed
files, then rerun verification.
