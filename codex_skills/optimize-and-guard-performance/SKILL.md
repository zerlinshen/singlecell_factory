---
name: optimize-and-guard-performance
description: Improve runtime and memory efficiency while enforcing performance regression guards and analytical parity. Use when profiling bottlenecks, tuning parallelism, reducing memory pressure, or validating before/after benchmark changes.
---

Optimize performance with explicit guardrails.

1. Baseline first.
- Capture wall-time, memory, and relevant throughput metrics.
- Use reproducible commands and comparable inputs/environment.

2. Optimize with bounded scope.
- Target highest-impact bottlenecks first.
- Prefer algorithmic/vectorization/dataflow improvements over cosmetic micro-tuning.
- Keep changes small enough for clear attribution.

3. Guard correctness and parity.
- Validate key outputs against baseline expectations.
- Define acceptable tolerance for numerical drift when needed.
- Treat unexplained result divergence as a blocker.

4. Run regression check.
- Compare before/after with explicit acceptance thresholds.
- Report pass/fail with command-level reproducibility.

5. Integrate with quality gate.
- If optimization changes behavior exposure, route to `/review-and-validate-quality`.
- Document tradeoffs and residual performance risks.
