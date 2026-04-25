---
name: review-and-validate-quality
description: Perform the quality gate across code risk review, scientific output validation, and documentation synchronization. Use before merge/reporting to catch correctness regressions, missing artifacts, weak tests, and stale docs.
---

Apply one integrated quality gate before declaring success.

1. Run risk-first code review.
- Prioritize correctness, regression, data integrity, and contract breaks.
- Rank findings by severity with path references and concrete fixes.
- If no correctness issues are found, state that explicitly.

2. Validate analysis outputs.
- Check run structure, module completeness, and key artifacts.
- Verify metadata sanity and scientific plausibility signals.
- Flag missing/invalid outputs with explicit impact.

3. Check test adequacy.
- Ensure new behavior has focused tests.
- Confirm failure paths and edge cases are covered where risk is high.

4. Synchronize docs with shipped behavior.
- Update README/PROTOCOL/CLI docs impacted by code changes.
- Keep module lists, dependencies, and outputs aligned with current behavior.

5. Produce gate decision.
- `pass`: no blocking correctness/scientific issues.
- `conditional`: non-blocking issues with explicit follow-ups.
- `fail`: blocking risks that must be fixed before merge/report.
