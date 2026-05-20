#!/usr/bin/env bash
# C benchmark Lane-1: CPU sparse_exact arpack PCA + scanpy CPU neighbors/UMAP/Leiden.
#
# Inputs (G-C0 A1): prepared_input.zarr (884k barcodes; --cohort-subset condition=tumor → 877,023 cells)
# Implementation: clustering.py:474-563 `_run_sparse_exact` (F-2 done, commit e34e837)
# Contract: G-C0 v2 APPROVED 2026-05-20 (approvals ledger G-C0-2026-05-20.md)
#
# Estimated wall clock: 3-6 hours on 96GB RAM workstation (CPU-bound).
# Output: /home/zerlinshen/projects/nc-reproduction/runs/<timestamp>/

set -euo pipefail

REPO_ROOT="/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory"
PROJECT_ROOT="/home/zerlinshen/projects/nc-reproduction"
RAW_ROOT="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526"
SAMPLE_ROOT="${RAW_ROOT}/full_cohort"

PROJECT_NAME="NC2024_C_LANE1_SPARSE_EXACT_20260520"
PROBE_LOG="${PROJECT_ROOT}/runs/${PROJECT_NAME}.launch.log"

# Lane-1 uses CPU; sc10x env has scanpy 1.11.5 + scanpy_external + harmonypy
CONDA_RUN="/home/zerlinshen/conda/bin/conda run -n sc10x"

# C benchmark: clustering + batch_correction ONLY (downstream modules are Component D's work, not C's)
MODULES="batch_correction,clustering"

mkdir -p "$(dirname "${PROBE_LOG}")"
: > "${PROBE_LOG}"

cd "${REPO_ROOT}"

set +e
SC_CLUSTERING_ENGINE=sparse_exact \
SC_DE_ENGINE=sparse \
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
  --gpu-mode off \
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
echo "LANE1_LOG=${PROBE_LOG}"
exit ${EC}
