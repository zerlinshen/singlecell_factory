---
name: dependency-upgrade
description: Upgrade dependencies safely with compatibility checks, changelog review, and validation runs. Use when bumping Python packages, CLI tools, or pinned versions.
argument-hint: [dependency-name-or-file]
---

Upgrade dependencies with controlled risk.

1. Scope the upgrade.
- Identify target package(s), current version(s), desired version(s), and motivation.
- Prefer incremental upgrades over multi-major jumps.

2. Check upstream change impact.
- Review release notes/changelog for breaking changes.
- Identify deprecated APIs used in this codebase.

3. Apply upgrade minimally.
- Update lock/spec files only as required.
- Avoid unrelated dependency churn.

4. Validate compatibility.
- Run targeted tests first, then broader suite if needed.
- Run a representative runtime command for high-risk dependencies.

5. Record upgrade evidence.
- Versions before/after.
- Breaking changes considered.
- Test/runtime results.
- Rollback command or version pin if issues appear.
