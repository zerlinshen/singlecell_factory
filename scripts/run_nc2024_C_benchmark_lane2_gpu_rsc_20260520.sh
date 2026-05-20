#!/usr/bin/env bash
# C benchmark Lane-2: GPU rapids-singlecell drop-in (Principle 8 compatible scanpy semantics).
#
# Inputs (G-C0 A1): prepared_input.zarr (884k barcodes; --cohort-subset condition=tumor → 877,023 cells)
# Implementation: clustering.py:_run_gpu uses rsc.pp.neighbors + rsc.tl.umap + rsc.tl.leiden (drop-ins)
# Determinism: Lane-2 GPU determinism probe PASSED STRICT (ARI=1.0, 2026-05-20)
# Contract: G-C0 v2 APPROVED 2026-05-20 (approvals ledger G-C0-2026-05-20.md)
#
# Estimated wall clock: 1.5-3 hours on RTX 5090.
# Output: /home/zerlinshen/projects/nc-reproduction/runs/<timestamp>/

set -euo pipefail

REPO_ROOT="/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory"
PROJECT_ROOT="/home/zerlinshen/projects/nc-reproduction"
RAW_ROOT="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526"
SAMPLE_ROOT="${RAW_ROOT}/full_cohort"

PROJECT_NAME="NC2024_C_LANE2_GPU_RSC_20260520"
PROBE_LOG="${PROJECT_ROOT}/runs/${PROJECT_NAME}.launch.log"

# Lane-2 uses GPU; sc_gpu_stable has rapids-singlecell 0.13.4 + scanpy_external
CONDA_RUN="/home/zerlinshen/conda/bin/conda run -n sc_gpu_stable"

# C benchmark: clustering + batch_correction ONLY
MODULES="batch_correction,clustering"

mkdir -p "$(dirname "${PROBE_LOG}")"
: > "${PROBE_LOG}"

cd "${REPO_ROOT}"

set +e
# IMPORTANT: SC_CLUSTERING_ENGINE NOT set (defaults to "auto" → routes to GPU when --gpu-mode=force)
# IMPORTANT: --gpu-mode force — Lane-2 must run on GPU; silent CPU fallback would invalidate the benchmark
SC_MEM_GUARD=on \
SC_MEM_WATCHDOG=on \
SC_DENSIFY_HARD_CAP_GB=24 \
SC_DENSIFY_SOFT_CAP_GB=12 \
${CONDA_RUN} python -m workflow.modular.cli \
  --project-root "${PROJECT_ROOT}" \
  --project "${PROJECT_NAME}" \
  --sample-root "${SAMPLE_ROOT}" \
  --optional-modules "${MODULES}" \
  --batch-key sample \
  --batch-method harmony \
  --cohort-subset condition=tumor \
  --gpu-mode force \
  --scale-mode large \
  --n-pcs 40 \
  --n-neighbors 15 \
  --leiden-resolution 0.8 \
  --random-state 42 \
  --checkpoint \
  --parallel-workers 1 \
  --no-run-cellranger-if-missing \
  --allow-dirty > "${PROBE_LOG}" 2>&1
EC=$?
set -e

echo "EC=${EC}" | tee -a "${PROBE_LOG}"
echo "LANE2_LOG=${PROBE_LOG}"
exit ${EC}
