---
name: review-and-validate-quality
description: Perform the quality gate across code risk review, scientific output validation, and documentation synchronization. Use before merge/reporting to catch correctness regressions, missing artifacts, weak tests, and stale docs.
---

Apply one integrated quality gate before declaring success.

1. Run risk-first code review.
- Prioritize correctness, regression, data integrity, and contract breaks.
- Rank findings by severity with path references and concrete fixes.
- If no correctness issues are found, state that explicitly.

2. Validate analysis outputs (including biological plausibility).
- Check run structure, module completeness, and key artifacts.
- Verify metadata sanity and scientific plausibility signals.
- Flag missing/invalid outputs with explicit impact.
- Biological plausibility checklist:
  - Marker genes: are assigned cell types consistent with known tissue markers?
  - Cluster count: reasonable for dataset size (typically 5-30 for 3K cells)?
  - DE results: number of significant genes per cluster in expected range (50-500)?
  - QC pass rate: was cell filtering neither too aggressive (<30% remaining) nor too lax (>95%)?
  - Trajectory: does pseudotime direction align with expected biological progression?
  - Batch correction: is biological variation preserved post-correction?
- When biological concerns are found, suggest launching `bio-reviewer` agent for detailed domain review.

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
