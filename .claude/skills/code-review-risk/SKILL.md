---
name: code-review-risk
description: Perform risk-first code review focused on bugs, regressions, and missing tests. Use when reviewing PRs, diffs, or local changes before merge.
argument-hint: [diff-or-path]
---

Review changes by severity, not by style preference.

1. Focus on correctness and regressions first.
- Functional bugs.
- Contract breaks and backward-compatibility risks.
- Concurrency/state/resource issues.
- Security or data integrity impact.

2. Validate change coverage.
- Missing tests for new behavior.
- Inadequate edge-case handling.
- Weak failure-path assertions.

3. Report findings in ranked order.
- Critical: can break production or corrupt results.
- High: likely incorrect behavior under normal scenarios.
- Medium: correctness risk under edge conditions.
- Low: maintainability concerns with limited immediate impact.

4. Each finding must include.
- File/path reference.
- Why it is risky.
- Concrete fix direction.

5. If no findings exist.
- State "No correctness findings."
- List residual risks or untested areas.
