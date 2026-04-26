---
name: test-run-cleanup-policy
description: Use when working with official 10x test datasets, smoke tests, benchmark reruns, or repeated pipeline experiments. After useful outputs are harvested, remove redundant or superseded run results while preserving raw data, references, and the current best run.
---

Use this skill for repeated test and benchmark workflows, especially with official 10x data.

1. Preserve high-value assets.
- Always keep raw downloaded data and reference resources.
- Do not delete pipeline source code or environment definitions.

2. Treat repeated benchmark/test outputs as reproducible.
- Old reruns, failed runs, stale result folders, and superseded experiment outputs can be cleaned after a better or newer run exists.
- Keep at least one best or currently referenced run unless the user explicitly says everything may be removed.
- For NC2024-style work, preserve:
  - canonical prepared input
  - canonical direct successful `massive` run
  - canonical controller-fallback successful run
  - one evidence-rich fallback log

3. Before cleanup, identify:
- which runs are test-only or official-dataset benchmarks
- which run is the current best reference
- which files are still needed for local reporting or bridge work

4. Safe cleanup targets usually include:
- `results/`
- `output/`
- transient bundles
- caches created only for test runs

5. Report cleanup clearly.
- State what was preserved
- State what was deleted
- State reclaimed disk space
