#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/home/zerlinshen/singlecell_factory"
CONDA_RUN="/home/zerlinshen/conda/bin/conda run -n sc_gpu"
RESULTS_DIR="${REPO_ROOT}/results"

LUAD_RAW_DEFAULT="${REPO_ROOT}/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/luad/cell_communication/cell_communication_liana.csv"
LUSC_RAW_DEFAULT="${REPO_ROOT}/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/lusc/cell_communication/cell_communication_liana.csv"
LUAD_TOP20_DEFAULT="${REPO_ROOT}/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_adenocarcinoma_checkpoint_top20.csv"
LUSC_TOP20_DEFAULT="${REPO_ROOT}/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_squamous_cell_carcinoma_checkpoint_top20.csv"

LUAD_RAW="${LUAD_RAW_DEFAULT}"
LUSC_RAW="${LUSC_RAW_DEFAULT}"
LUAD_TOP20="${LUAD_TOP20_DEFAULT}"
LUSC_TOP20="${LUSC_TOP20_DEFAULT}"
OUTPUT_DIR="${RESULTS_DIR}"
PROJECT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --luad-raw-csv) LUAD_RAW="$2"; shift 2 ;;
    --lusc-raw-csv) LUSC_RAW="$2"; shift 2 ;;
    --luad-top20-csv) LUAD_TOP20="$2"; shift 2 ;;
    --lusc-top20-csv) LUSC_TOP20="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --project) PROJECT="$2"; shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "${PROJECT}" ]]; then
  PROJECT="NC2024_NSCLC_SUBTYPE_CHECKPOINT_AUDIT_VERIFIED_AUTO_$(date -u +%Y%m%d_%H%M%S)"
fi

cd "${REPO_ROOT}"
${CONDA_RUN} python -m pytest -q -o addopts='' tests/test_nc2024_subtype_checkpoint_audit_script.py
AUDIT_DIR="$(${CONDA_RUN} python scripts/audit_nc2024_subtype_checkpoint_pairs.py \
  --luad-raw-csv "${LUAD_RAW}" \
  --lusc-raw-csv "${LUSC_RAW}" \
  --luad-top20-csv "${LUAD_TOP20}" \
  --lusc-top20-csv "${LUSC_TOP20}" \
  --output-dir "${OUTPUT_DIR}" \
  --project "${PROJECT}")"
${CONDA_RUN} python scripts/summarize_nc2024_checkpoint_hierarchy.py --audit-dir "${AUDIT_DIR}"
echo "RUN_DIR=${AUDIT_DIR}"
