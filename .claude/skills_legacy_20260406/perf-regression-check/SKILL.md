---
name: perf-regression-check
description: Detect and explain runtime or memory regressions using before/after benchmarks and acceptance thresholds. Use when performance may have changed after code edits.
argument-hint: [benchmark-command-or-target]
---

Run performance checks with explicit pass/fail criteria.

1. Establish baseline and candidate runs.
- Use same dataset/input and similar environment.
- Run baseline (before) and candidate (after) using reproducible commands.

2. Capture core metrics.
- Wall time.
- CPU time/utilization (if available).
- Peak memory.
- Accuracy/parity metrics when algorithmic output must match.

3. Compare against thresholds.
- Define acceptable deltas before analysis.
- Mark pass/fail for each metric, not just aggregate opinion.

4. Investigate regressions.
- Identify likely bottleneck stage/function.
- Suggest targeted optimization, not broad rewrites.

5. Report in concise table format.
- Metric | before | after | delta | threshold | pass/fail.
