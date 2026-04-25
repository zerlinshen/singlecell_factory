#!/usr/bin/env bash
# NC2024 paper-faithful two-stage run: tumor cohort + B/H cohort
# Phase 7B — paper-aligned params, capability flags explicit (no --scale-mode)
# Usage: bash scripts/run_nc2024_paper_aligned_20260425.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${NC2024_DATA_ROOT:-${REPO_ROOT}/data/raw/emtab13526}"
RESULTS_ROOT="${NC2024_RESULTS_ROOT:-${REPO_ROOT}/results}"
TIMESTAMP="20260425"

# --- Paper-aligned capability flags (explicit, no --scale-mode) ---
CAPABILITY_FLAGS=(
    --doublet-strategy grouped
    --clustering-engine sparse_exact
    --checkpoint-policy mandatory_only
    --lazy-read auto
)

# --- Paper-aligned DE params (override defaults) ---
PAPER_DE_FLAGS=(
    --de-correction bonferroni
    --de-min-pct 0.30
    --de-method wilcoxon
    --de-logfc-threshold 0.0
)

# --- Env: sparse engines + memory guards ---
export SC_DE_ENGINE=sparse
export SC_CNV_ENGINE=chunked
export SC_CELLCOMM_ENGINE=sparse
export SC_MEM_GUARD=on
export SC_MEM_WATCHDOG=on

# ============================================================
# Stage 1: Tumor cohort (lung adenocarcinoma + squamous cell)
# ============================================================
TUMOR_RUN_DIR="${RESULTS_ROOT}/nc2024_tumor_${TIMESTAMP}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage 1: Tumor cohort -> ${TUMOR_RUN_DIR}"

python -m workflow.modular.cli \
    --project nc2024_tumor \
    --sample-root "${DATA_ROOT}" \
    --output-dir "${TUMOR_RUN_DIR}" \
    --optional-modules "clustering,batch_correction,differential_expression,annotation,pseudobulk_de" \
    --cohort-subset "disease=lung_adenocarcinoma,lung_squamous_cell_carcinoma" \
    --n-pcs 15 \
    --leiden-resolution 1.0 \
    --expected-doublet-rate 0.06 \
    --batch-method harmony \
    --batch-key sample \
    --checkpoint \
    "${CAPABILITY_FLAGS[@]}" \
    "${PAPER_DE_FLAGS[@]}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage 1 complete."

# ============================================================
# Stage 2: Background + Healthy cohort (normal tissue)
# ============================================================
BH_RUN_DIR="${RESULTS_ROOT}/nc2024_bh_${TIMESTAMP}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage 2: B+H cohort -> ${BH_RUN_DIR}"

python -m workflow.modular.cli \
    --project nc2024_bh \
    --sample-root "${DATA_ROOT}" \
    --output-dir "${BH_RUN_DIR}" \
    --optional-modules "clustering,batch_correction,differential_expression,annotation,pseudobulk_de" \
    --cohort-subset "disease=normal" \
    --n-pcs 15 \
    --leiden-resolution 1.0 \
    --expected-doublet-rate 0.06 \
    --batch-method harmony \
    --batch-key sample \
    --checkpoint \
    "${CAPABILITY_FLAGS[@]}" \
    "${PAPER_DE_FLAGS[@]}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage 2 complete."
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Both stages finished. Results in:"
echo "  Tumor:       ${TUMOR_RUN_DIR}"
echo "  Background:  ${BH_RUN_DIR}"
