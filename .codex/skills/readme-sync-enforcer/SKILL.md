---
name: readme-sync-enforcer
description: Use whenever code, workflow behavior, execution modes, environment assumptions, or operational guidance change. Ensure the relevant README is updated in the same workstream; if none exists, create one.
---

Use this skill whenever a pipeline, module, workflow, or project behavior changes.

1. Assume docs drift unless checked.
- If code or workflow behavior changes, inspect the nearest relevant README.
- If no README exists for the changed surface, create one.

2. Update docs in the same workstream.
- Do not postpone README updates to a vague later step.
- Record new flags, fallback behavior, environment caveats, and changed outputs.

3. Keep the README operational.
- Prefer concrete execution guidance, caveats, and examples.
- Mention important failure modes and workarounds.

4. For bioinformatics projects, explicitly document:
- dataset-scale modes
- GPU vs CPU behavior
- fallback rules
- staging strategies for large cohorts
- remote-to-local bridge/report workflow when relevant
- run-memory and handoff surfaces when operational defaults changed
- `SKILL_PROTOCOL.md` and `SKILL_PLAYBOOK.md` when skill-routing behavior changed

5. Final check.
- Before declaring work complete, ask: did behavior change enough that a user would be surprised if the README stayed old?
- If yes, update the README before finishing.
