---
name: safe-refactor
description: Refactor code while preserving behavior through contract checks and staged validation. Use when restructuring modules, renaming internals, or extracting shared logic.
argument-hint: [target-path-or-module]
---

Refactor with behavior safety guarantees.

1. Define contracts before editing.
- Public APIs, expected inputs/outputs, side effects, file formats, CLI flags.
- Existing tests that represent those contracts.

2. Plan minimal-risk sequence.
- Small commits/patches by concern (extract, rename, move, cleanup).
- Keep compatibility shims when needed during transition.

3. Edit in stages with checks.
- After each stage, run the narrowest relevant tests/lints.
- Stop and fix immediately on first regression.

4. Preserve external behavior.
- Keep user-facing interfaces stable unless change is explicitly requested.
- If behavior must change, update tests/docs in the same patch.

5. Output final assurance summary.
- What changed structurally.
- What contracts were verified.
- Remaining risks and follow-up checks.
