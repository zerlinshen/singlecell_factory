---
name: optimize-modular-performance
description: Improve runtime and memory efficiency of the modular pipeline while preserving analytical correctness. Use when profiling bottlenecks, vectorizing module code, tuning parallel workers, or verifying performance regressions.
argument-hint: [target-module-or-stage]
---

Optimize performance with correctness gates.

1. Baseline before changes.
- Run existing perf checks (`scripts/benchmark_modular_performance.py` when dataset is available).
- Capture wall time, peak RSS, and key output parity metrics.

2. Prefer safe optimization patterns already used in this repo.
- Vectorize loops in NumPy/Pandas/Scipy.
- Avoid unnecessary `adata` copies and dense conversions.
- Reuse precomputed matrices for repeated scoring operations.
- Keep plotting I/O asynchronous where supported by `PipelineContext`.

3. Preserve module contracts.
- Do not change required `adata` fields without migration updates.
- Keep output filenames stable unless docs/tests are updated in the same change.
- For parallelization, verify branch-merge behavior for `obs`, `obsm`, and `uns`.

4. Validate correctness after optimization.
- Run focused tests (for example `tests/test_modular_optimizations.py`).
- Run dependency and integration tests touching orchestration.
- Require numerical parity or bounded tolerance for algorithm changes.

5. Report optimization outcome.
- Before/after runtime and memory metrics.
- Accuracy/parity result.
- Risk notes (dataset-specific gains, backend variability, optional dependency effects).
