# 2026-08-22 - RAPIDS 26.08 GPU stack promotion

## Objective

- Restore reliable GPU execution for clustering, marker DE, doublet detection,
  scVI, and CellRank without mutating `sc_gpu` or `sc_gpu_stable`.

## Starting Context

- Canonical factory: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory`.
- Prior evidence showed basic CUDA worked but RAPIDS PCA/DE could fail with
  `CUSOLVER_STATUS_INTERNAL_ERROR` or `CUBLAS_STATUS_NOT_INITIALIZED`.
- The existing `sc_gpu` imported CellRank 2.2.0 but failed because `pygpcca`
  was absent.

## What Was Run

- Host: `zerlinshen-MS-7E49` (`ubuntu-tail`).
- Working directory: canonical `singlecell_factory` plus governed project output.
- Controlling script:
  `/home/zerlinshen/projects/gpu-stack-validation-20260822/runs/2026-08-21T1924Z-1137750/python/gpu_stack_smoke.py`.
- Run directory:
  `/home/zerlinshen/projects/gpu-stack-validation-20260822/runs/2026-08-21T1924Z-1137750/`.
- Execution mode: `benchmark` (bounded technical comparison, not a cohort run).

## Outcome

- `success` for the bounded core GPU lane.
- Promoted environment: `sc_gpu_rapids2608`.
- `sc_gpu` and `sc_gpu_stable` retained unchanged as rollback/full-optional lanes.

## Evidence

- Driver 595.84 loaded on kernel 7.0.0-29; RTX 5090 D v2; CUDA capability 13.2.
- CuPy cuSOLVER and PyTorch CUDA smokes passed.
- RAPIDS normalize/HVG/PCA/neighbors/UMAP/Leiden, GPU DE, and GPU Scrublet passed.
- Historical 50k x 20k RAPIDS PCA reproducer passed without cuSOLVER/cuBLAS failure.
- scVI recorded both trainer and module on `cuda:0` for two epochs; CellRank
  produced a row-stochastic transition matrix and completed GPCCA Schur/
  macrostate computation.
- CPU/GPU clustering ARI 1.0; DE top-30 marker Jaccard 1.0 for both groups.
- Dedicated 50k doublet fixture after installing Annoy: 99.936% CPU/GPU call
  agreement; CPU/GPU calls 2,468/2,500; truth F1 0.9936/1.0.
- Focused factory tests: 44 passed, 43 deselected.
- The full 34-gate suite under the established launcher environment passed
  32/34. The only failures were the license gate and its governance self-test,
  both caused by pre-existing unregistered external clones `DiffHiChIP` and
  `JASPAR2024`; all single-cell, plotting, R, stress, and promotion gates passed.
- Running the suite with `sc_gpu_rapids2608` on `PATH` passed 31/34: the same
  two governance failures plus ten RMS plotting pixel regressions in the
  environment containing Matplotlib 3.11. Matplotlib is the most plausible
  additional cause, but the environments were not a controlled one-variable
  A/B. This confirms that the promoted environment is a GPU lane, not a
  replacement for the suite's reference-render environment.
- Human report: run `evidence/VALIDATION_REPORT.md`; machine results: run `ops/*.json`.

## Problems Encountered

- Tailscale control traffic initially failed through the previous proxy route;
  resolved separately with a Tailscale-only Clash rule before this run.
- `rapids-singlecell 0.16.1` is source-only on PyPI and initially failed because
  CMake could not find nvcc. Adding the CUDA 13.2 `cuda-nvcc` entry and setting
  the active toolkit path resolved the build.
- CellRank failure in `sc_gpu` was the missing `pygpcca`, not a GPU failure.
- A compact doublet fixture yielded misleading low score correlation because
  CPU Scrublet called zero doublets. It was superseded for the decision by the
  dedicated 50k injected-doublet harness.
- LIANA 1.9 requires `pandas<3`, while RAPIDS 26.08 requires `pandas>=3`.
  A dry-run resolver confirmed they cannot share this environment cleanly.
- The suite license gate is presently blocked independently of this change by
  unregistered external clones under `/home/zerlinshen/downloads/external_refs`.
  They were neither deleted nor registered because that requires a separate
  provenance/license decision.
- An isolated focused pytest initially failed three CPU-DE tests because
  `tests/conftest.py` defaults `NUMBA_DISABLE_JIT=1`. A fresh NumPy/Numba cache
  reproduced the same failure in the established launcher environment. The
  repository's documented subprocess contract, `NUMBA_DISABLE_JIT=0`, restored
  the expected 44-pass result without changing the production dependency
  stack. Failed diagnostic logs are retained beside the passing log.
- The historical 50k PCA script also requires an activated environment (or an
  explicit CUDA root at `targets/x86_64-linux`) so CuPy can find toolkit
  headers. Both failed launch attempts and the corrected passing transcript are
  retained.

## What Was Changed

- Added `environment_gpu_rapids2608.yml` and exact Conda/pip locks.
- Updated README, PROTOCOL, AGENTS, and AI_AGENT_PROTOCOL GPU guidance.
- Created the governed validation project and promoted clone
  `/home/zerlinshen/conda/envs/sc_gpu_rapids2608`.
- Removed the temporary `sc_gpu_rapids2608_py312_candidate` diagnostic
  environment after proving that the focused-test failures came from the JIT
  launch contract rather than Python 3.13 or Numba 0.64. The original
  `sc_gpu_rapids2608_candidate`, promoted lane, and both rollback lanes remain.

## Resolution

- Use `sc_gpu_rapids2608` for clustering, marker DE, Scrublet, scVI, and CellRank.
- Keep LIANA and other pandas<3 consumers in their existing compatible lane.
- Retain the pre-import CUDA warm-up because it remains cheap and defensive,
  even though the updated stack passed the old-order stress reproducer.

## Initial cautions (historical; suite-status caution superseded below)

- This is bounded synthetic/technical validation, not real-cohort biological certification.
- Do not force-install LIANA into the RAPIDS 26.08 environment.
- Do not use the promoted environment to update plotting reference images.
- At this point in the run, suite-level 34/34 was not yet claimable because
  `DiffHiChIP` and `JASPAR2024` still needed external-reference disposition.
  The final governance resolution below supersedes this status.
- Use `PYTHONNOUSERSITE=1` for validation to prevent user-site package leakage.
- Set `NUMBA_DISABLE_JIT=0` for isolated factory pytest validation.
- The candidate environment is retained until a later authorized cleanup; do
  not delete it while this run is the promotion evidence.

## Improvement Ideas

- Add the bounded GPU smoke as an opt-in hardware gate, not to data-free CI.
- Add an environment-lane resolver so LIANA routing is explicit rather than operator memory.
- Repeat clustering/DE/doublet parity on a bounded real dataset before any
  publication-grade performance claim.

## Classification

- Run and report: `canonical` technical environment evidence.
- Compact `ops/parity.json` doublet-score interpretation: `evidence-only`,
  superseded for the doublet decision by `ops/doublet-parity-50k.json`.
- No scientific output or biological claim was produced.

## Final governance resolution (supersedes the earlier suite cautions)

- A later same-day GPU-first rerun passed the complete canonical suite: 34/34.
- The plotting gate now invokes the validated reference-render interpreter
  explicitly, so the GPU lane's Matplotlib does not redefine visual truth.
- `DiffHiChIP`, active JASPAR 2024, and the JASPAR 2026 candidate are registered
  with license, provenance, scope, and exact local identity. The corrected
  JASPAR 2024 consumer sidecar matches the independently recomputed digest.
- Gate 18 proves retained produced-artifact integrity for 8/8 figures. The
  separate original-paper reference audit has no active reference set, reports
  `not_run`, and contributes no pass to 34/34.
- Therefore, the earlier 32/34 and 31/34 results and the caution against a
  34/34 claim remain diagnostic history only. They are superseded by
  `evidence/VALIDATION_REPORT.md` and
  `logs/full-suite-gpu-path-artifact-integrity-final-20260822.log` in the
  governed validation run.
