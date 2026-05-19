# Component E — GPU DE Parity Report

**Status:** ⚠️ **G-E2 PARTIAL FAIL — GPU DE INADMISSIBLE at current pin**
**Component:** E (GPU DE environment validation)
**Plan:** `~/.omc/plans/nc-cell-clustering-final-strategy-plan.md` (APPROVED v4)
**Date:** 2026-05-19/20 (Tail-host session)
**Author:** zerlinshen + Claude (Opus 4.7)

## Headline result

**GPU DE via `rapids_singlecell.tl.rank_genes_groups` is INADMISSIBLE for final-claim runs at this version pinning.** Component D's methodology MUST use CPU Wilcoxon DE per F-3 contract. CPU Wilcoxon is the science-grade default; GPU DE was a performance optimization, not a scientific requirement.

## Gate outcomes

| Gate | Result | Evidence |
|---|---|---|
| **G-E1 (isolation invariant)** | ✅ **PASS** | `sc_gpu_stable` SHA256 unchanged across all 7 build attempts (v1 clone, v2-v7 fresh env). `ops/env_validation/G-E1-sc_gpu_stable-{before,after}-*.sha256` hashes match. Final post-attempt check: `sc_gpu_stable` still has `rapids-singlecell 0.13.4` + `cuml 25.10` intact. |
| **G-E2 (API presence)** | ⚠️ **PARTIAL FAIL** | rapids-singlecell **0.14.1** installed in `sc_gpu_de_test` env (Python 3.12, RAPIDS 26.04 base). Package fails to import at `__init__.py` due to missing runtime dependency `cuvs` (`ModuleNotFoundError: No module named 'cuvs'`). Installing `cuvs` requires upgrading the entire RAPIDS stack from CUDA 12.9 → CUDA 13.2 (3 GB download, 37 package upgrades, 27 changes), which would destabilize the env and create CUDA version mismatch with the host driver (currently CUDA 13.2 driver). Without `cuvs`, `rapids_singlecell.tl.rank_genes_groups` cannot be reached. |
| **G-E3 (parity)** | ❌ **N/A** | Cannot run parity tests if `rapids_singlecell` cannot be imported. Recorded as N/A pending G-E2 resolution. |

## Final env state

- `sc_gpu_de_test` exists at `/home/zerlinshen/conda/envs/sc_gpu_de_test`
- Python 3.12.13, RAPIDS 26.04 (cuML / cudf / cugraph / rmm), cuda-version 12.9
- `rapids-singlecell-0.14.1` installed via pip (with cuda-nvcc + cuda-cudart-dev for build)
- `cuvs` NOT installed (the unmet runtime dependency)
- `sc_gpu_stable` UNCHANGED from start of Component E work

## Build attempts (chronological)

1. **v1 (clone of sc_gpu_stable, pip upgrade to rapids-singlecell>=0.14):** FAIL — Python 3.11 incompatible (rapids-singlecell ≥0.14 requires Python ≥3.12).
2. **v2 (background nohup with set -e; mamba create fresh):** orphaned + raced with v3.
3. **v3 (fresh mamba env, Python 3.12 + RAPIDS 25.10 stack):** mamba create succeeded but pip install rapids-singlecell failed — wheel build needed CUDA toolkit (nvcc), absent from pip build sandbox.
4. **v4 (mamba install rapids-singlecell from rapidsai/conda-forge channels):** FAIL — rapids-singlecell is not on conda channels (verified via `mamba search`). Falls back to pip + cuda-nvcc install.
5. **v5 (install cuda-nvcc + cuda-cudart-dev into env, then pip with CUDA env vars):** SUCCEEDED for pip install of rapids-singlecell-0.14.1, but G-E2 smoke FAILED with `No module named 'cuvs'`.
6. **v6 (mamba install cuvs + pylibcuvs):** FAIL — `pylibcuvs` package does not exist on any channel.
7. **v7 (mamba install cuvs alone):** mid-transaction proposed CUDA 12→13 upgrade across 37+27 packages with 3 GB download; killed before commit to preserve env stability.

## Why this is acceptable per the canonical plan

The parent plan's Component E (`~/.omc/plans/nc-cell-clustering-final-strategy-plan.md` Component E + the E executor brief at `singlecell_factory/ops/briefs/E-gpu-de-env-executor-brief.md`) anticipates this outcome:

> "On G-E3 fail at any scale tier: GPU DE excluded from final-claim paths until a documented fix; Component D uses CPU DE only."

Generalizing to G-E2 partial fail: same effect — Component D uses CPU Wilcoxon DE only (which is already the F-3-enforced default per Principle 6: "No silent algorithmic substitution"). The factory's existing CPU Wilcoxon path on the 877k NC tumor cohort is the legitimate, peer-reviewable, defensible DE method.

The `sc_gpu_de_test` env at its current state is a useful artifact for future investigation (when rapids-singlecell's `cuvs` runtime dependency relaxes, or when a CUDA 13 RAPIDS upgrade is taken on as a separate plan).

## Follow-up suggestions (out of scope for this Component E completion)

1. **Wait for rapids-singlecell ≥0.15-stable** (currently 0.15.0rc7 is latest pre-release). Check whether 0.15-stable softens the `cuvs` import to lazy / optional.
2. **CUDA 12 → 13 RAPIDS upgrade plan**: a separate strategic plan that audits ALL existing factory consumers (sc10x, sc_rna_velocity_pseudotime, sc_gpu_stable, sc_gpu, r_multiomics envs) for CUDA compatibility. Not within Component E's scope.
3. **Alternative GPU DE path**: investigate whether `cudf-pandas` or `cupy` numpy backend can accelerate scanpy's CPU Wilcoxon without rapids-singlecell. Out of scope.

## Recommended D methodology binding

D's methodology-binding subsection (gated on G-G2 + ambient module integration + G-C3) should commit to: **DE engine = CPU Wilcoxon via `scanpy.tl.rank_genes_groups(method="wilcoxon")` on the F-3-hardened path.** No `--allow-welch-fallback`. No GPU DE.

This decision is recorded here as Component E's final output for the current planning round.
