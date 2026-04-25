---
name: docs-sync
description: Keep documentation aligned with code changes by updating affected README/protocol/CLI sections in the same change. Use after implementing features, flags, module changes, or behavior updates.
argument-hint: [changed-files-or-feature]
---

Synchronize docs with shipped behavior.

1. Identify doc-impacting changes.
- CLI arguments and defaults.
- Module names/dependencies.
- Output files/paths.
- Setup/run instructions and examples.

2. Update authoritative docs only.
- Edit existing docs instead of creating redundant new files.
- Keep examples executable and consistent with current code.

3. Validate documentation correctness.
- Verify command snippets against actual CLI/options.
- Remove stale references to deleted behavior.

4. Summarize doc updates.
- What user-facing behavior changed.
- Which docs were updated.
- Any known documentation gaps to follow up.
