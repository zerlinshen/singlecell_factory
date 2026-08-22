---
name: optimize-and-guard-performance
description: Improve runtime, memory, GPU utilization, or throughput while guarding analytical parity, fallback behavior, and reproducibility. Use for profiling, scaling, CPU/GPU comparisons, parallelism, compilation/cache effects, or before/after performance validation; not for unsupported publication-speed claims from a single synthetic smoke.
---

# Optimize and Guard Performance

Optimize only after defining the result that must remain true.

## Establish a comparable baseline

- Record code/input identity, command, environment/lock, host resources,
  thread/process settings, seed, warm-up policy, cache state, and data shape.
- Capture wall time, peak RSS, peak VRAM where relevant, throughput, output
  size, and the domain metric that defines correctness.
- Separate environment/launch failures from code, algorithm, and input failures.
  Retry one safe transient once; investigate repeated failures.
- For JIT, CUDA, distributed, or cached paths, separate first-run compilation
  from steady state. Use repeated trials or explain why one trial is sufficient.

## Optimize with bounded attribution

- Profile first and target the dominant bottleneck. Prefer algorithm, sparsity,
  dataflow, batching, vectorization, and avoided materialization before
  micro-tuning.
- Change one causal surface at a time when practical. If a dependency or whole
  environment changes, treat it as a compatibility migration and use
  `upstream-fit-assessment`.
- Keep incompatible analysis, plotting, or reporting consumers in separate
  lanes. Faster GPU execution does not justify a suite-wide environment switch.
- For large data, stage scale deliberately and treat swap as a safety buffer,
  not additional working memory.

## Prove execution and parity

- Show that the optimized backend actually ran: device residency, native-kernel
  evidence, backend metadata, or another direct signal. GPU availability and a
  zero exit do not exclude CPU fallback.
- Define acceptance thresholds before comparison. Use exact identity only when
  the method is deterministic; otherwise use justified tolerances and seeds.
- Compare scientific outputs at the level users consume: assignments,
  neighborhoods, markers/rankings, statistics, or other domain invariants, not
  only array shape and finiteness.
- Use synthetic fixtures for wiring and falsifiability, then a bounded
  representative dataset before production integration. Do not generalize an
  easy planted fixture into real-cohort or publication-performance evidence.
- Treat unexplained divergence, missing fallback provenance, OOM/segfault, or a
  failed claim-critical negative control as a blocker.

## Regression and promotion

1. Run a fast functional smoke and a meaningful before/after benchmark on the
   same inputs and environment boundary.
2. Report raw measurements, repeated-summary statistics when available,
   correctness/parity results, thresholds, and pass/fail. Keep failed attempts
   that explain the final launch contract.
3. Check relevant broader gates for side effects. Classify unrelated governance
   or visual-reference failures separately instead of weakening baselines.
4. Promote only the exercised lane and retain a rollback path. Record exact
   environment locks and locally built artifact hashes when applicable.
5. Route behavior or claim changes through `review-and-validate-quality`; use
   independent adversarial review when the project contract or material claim
   requires it.

## Output

Return the bottleneck, baseline and candidate identities, benchmark protocol,
raw and summarized metrics, analytical parity, execution/fallback evidence,
decision, promoted scope, rollback, and residual risks.
