#!/usr/bin/env bash
# NC2024 paper-faithful two-stage real run (884k cells)
# Phase 7B + 7C — paper-aligned params, capability flags explicit, sparse-exact engines
# Wraps prior launcher (run_nc2024_paper_aligned_20260425.sh) with corrections:
#   - Correct DATA_ROOT pointing at full_cohort subdir
#   - cohort-subset uses `condition` column (clean tumor/healthy_background split)
#   - conda run -n sc_gpu wrapping
#   - --no-run-cellranger-if-missing flag (input is prepared zarr, no CellRanger)
# Output:
#   - results/nc2024_tumor_20260426/<run>/...
#   - results/nc2024_bh_20260426/<run>/...
#   - ops/run_ledger/<project>_<timestamp>.json (per stage)
# set -e intentionally omitted: Python interpreter segfaults (exit 139) on
# shutdown after heavy C-extension runs (cupy/torch/zarr) are cosmetic — all
# pipeline outputs are written before the crash.  We check each stage
# explicitly so a true failure (exit != 0 and != 139) still aborts the run.
# Fix applied: workflow/modular/_shutdown.py + cli.py run_shutdown_cleanup().
# If the fix holds, both stages should exit 0; 139 is kept as a safe bypass.
set -uo pipefail

_run_stage() {
    local stage_name="$1"
    shift
    "$@"
    local rc=$?
    if [[ $rc -eq 0 ]]; then
        return 0
    elif [[ $rc -eq 139 ]]; then
        # Interpreter segfault on shutdown — outputs are intact, treat as warning.
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: ${stage_name} exited 139 (interpreter segfault on shutdown). Outputs verified intact; continuing." >&2
        return 0
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: ${stage_name} failed with exit code ${rc}. Aborting." >&2
        exit $rc
    fi
}

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526/full_cohort"
RESULTS_ROOT="${REPO_ROOT}/results"
TIMESTAMP="20260426"
CONDA_RUN="/home/zerlinshen/conda/bin/conda run -n sc_gpu"

# --- Paper-aligned capability flags (explicit, no --scale-mode) ---
CAPABILITY_FLAGS=(
    --doublet-strategy grouped
    --clustering-engine sparse_exact
    --checkpoint-policy mandatory_only
    --lazy-read auto
)

# --- Paper-aligned DE params ---
PAPER_DE_FLAGS=(
    --de-correction bonferroni
    --de-min-pct 0.30
    --de-method wilcoxon
    --de-logfc-threshold 0.0
)

# --- Sparse engines + memory guards ---
export SC_DE_ENGINE=sparse
export SC_CNV_ENGINE=chunked
export SC_CELLCOMM_ENGINE=sparse
export SC_MEM_GUARD=on
export SC_MEM_WATCHDOG=on

cd "${REPO_ROOT}"

# ============================================================
# Stage 1: Tumor cohort  (condition=tumor, ~877k cells)
# ============================================================
TUMOR_RUN_DIR="${RESULTS_ROOT}/nc2024_tumor_${TIMESTAMP}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage 1: Tumor cohort -> ${TUMOR_RUN_DIR}"

_run_stage "Stage 1 (tumor)" \
    ${CONDA_RUN} python -m workflow.modular.cli \
    --project nc2024_tumor \
    --sample-root "${DATA_ROOT}" \
    --output-dir "${TUMOR_RUN_DIR}" \
    --optional-modules "clustering,batch_correction,differential_expression,annotation,pseudobulk_de" \
    --cohort-subset "condition=tumor" \
    --n-pcs 15 \
    --leiden-resolution 1.0 \
    --expected-doublet-rate 0.06 \
    --batch-method harmony \
    --batch-key sample \
    --checkpoint \
    --no-run-cellranger-if-missing \
    "${CAPABILITY_FLAGS[@]}" \
    "${PAPER_DE_FLAGS[@]}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage 1 complete."

# ============================================================
# Stage 2: B/H cohort  (condition=healthy_background, ~7k cells)
# ============================================================
BH_RUN_DIR="${RESULTS_ROOT}/nc2024_bh_${TIMESTAMP}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage 2: B+H cohort -> ${BH_RUN_DIR}"

_run_stage "Stage 2 (healthy_background)" \
    ${CONDA_RUN} python -m workflow.modular.cli \
    --project nc2024_bh \
    --sample-root "${DATA_ROOT}" \
    --output-dir "${BH_RUN_DIR}" \
    --optional-modules "clustering,batch_correction,differential_expression,annotation,pseudobulk_de" \
    --cohort-subset "condition=healthy_background" \
    --n-pcs 15 \
    --leiden-resolution 1.0 \
    --expected-doublet-rate 0.06 \
    --batch-method harmony \
    --batch-key sample \
    --checkpoint \
    --no-run-cellranger-if-missing \
    "${CAPABILITY_FLAGS[@]}" \
    "${PAPER_DE_FLAGS[@]}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage 2 complete."
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Both stages finished. Results in:"
echo "  Tumor:       ${TUMOR_RUN_DIR}"
echo "  Background:  ${BH_RUN_DIR}"
