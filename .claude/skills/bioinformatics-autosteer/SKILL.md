---
name: bioinformatics-autosteer
description: Route single-cell, spatial, WGS, Hi-C/3C, proteomics, large-data execution, recovery, and artifact handoff through the canonical project and the smallest relevant bioinformatics workflow. Use when ownership, execution lane, environment, or evidence requirements must be established before work begins.
---

# Bioinformatics Autosteer

Establish ownership, lane, and proof requirements before selecting a helper.
Do not turn routing into a second orchestrator.

## Start from current evidence

1. Resolve the canonical factory or project from the portfolio and nearest
   `AGENTS.md`; a similarly named checkout is not ownership evidence.
2. For the factory suite, keep code in the owning factory and scientific
   outputs under `/home/zerlinshen/projects/<project-id>/runs/<run-id>/` with a
   manifest and explicit claim boundary.
3. Read the current README, environment declaration, lock, and latest run
   journal. Do not hard-code a historical GPU environment or provider version
   in this skill; resolve the documented current lane at execution time.
4. Classify the request as module development, unchanged execution, recovery,
   performance work, upstream selection, remote handoff, or final review.

## Route to the smallest useful workflow

- Remote execution, resume, benchmark, or large transfer on `ubuntu-tail`:
  `before-every-run` first, then the registered entrypoint. Read
  [the routing contract](../before-every-run/references/shenxin-routing-contract.md)
  when remote lane ordering matters.
- Unchanged mature pipeline run: `execute-and-recover-pipeline`; use
  `singlecell-remote-workflow` when remote-to-local handoff is part of the ask.
- New or evolved `singlecell_factory` module/contract:
  `singlecell-factory-module-delivery`. Use `develop-and-integrate-module` for
  broader modular or cross-factory implementation.
- New package, upstream method, SDK, or native/GPU stack:
  `upstream-fit-assessment` before installation or custom code.
- Runtime, memory, CPU/GPU, or throughput work:
  `optimize-and-guard-performance` with analytical parity and fallback checks.
- Final substantive handoff: `review-and-validate-quality`; use the standalone
  Grok evidence review only when requested or required by the project contract.

Do not stack helpers that repeat the same planning or review function. Dev OS
remains an explicit risk-based development choice, not an automatic consequence
of working on bioinformatics.

## Execution and environment boundaries

- Record the actual host, interpreter, environment prefix, package versions,
  driver/runtime/compute capability for GPU work, command, input identity,
  seed, and output location.
- Prove device residency or native-kernel execution where CPU fallback is
  possible. An import, `--version`, available GPU, or zero exit is insufficient.
- Keep incompatible consumers in separate lanes instead of weakening a valid
  solve. A GPU analysis lane does not automatically become the plotting,
  report, or suite-wide default.
- Treat swap as a crash buffer, not RAM. Stage large runs, preserve checkpoints,
  and select the smallest module set that answers the current question.
- Preserve raw data, verified references/indexes, manifests, best runs, and
  claim-critical evidence. Caches and superseded outputs are removable only
  under the applicable retention policy.

## Evidence boundary

- Use synthetic fixtures for wiring, fallback, and bounded parity checks.
- Require representative real data before production integration or biological,
  real-cohort, or publication-performance claims.
- Separate technical success, analytical parity, scientific validity, and
  suite-wide governance status. Report each independently, including blocked or
  untested lanes.
