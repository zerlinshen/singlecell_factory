#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_ROOT="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526"
FULL_ROOT="${RAW_ROOT}/full_cohort"
RESULTS="${REPO_ROOT}/results"
PREP_SCRIPT="${REPO_ROOT}/scripts/prepare_emtab13526_full_cohort_zarr.py"
STAGE1_LAUNCHER="${REPO_ROOT}/scripts/run_emtab13526_full_cohort_stage1.sh"
PREPARED_ZARR="${FULL_ROOT}/prepared_input.zarr"
READY_SENTINEL="${FULL_ROOT}/prepared_input.ready"
PREP_LOG="${RESULTS}/NC2024_NSCLC_FULL_COHORT_PREPARE.launch.log"
CTRL_LOG="${RESULTS}/NC2024_NSCLC_FULL_COHORT_STAGE1_WITH_FALLBACK.launch.log"
LARGE_PROJECT="NC2024_NSCLC_FULL_COHORT_STAGE1_LARGE_AUTO"
MASSIVE_PROJECT="NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO"
CONDA_RUN="/home/zerlinshen/conda/bin/conda run -n sc_gpu"

timestamp_utc() {
  date -u +"%Y%m%d_%H%M%S"
}

log() {
  printf '[%s] %s\n' "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" | tee -a "${CTRL_LOG}"
}

capture_preflight_inventory() {
  local ts inventory
  ts="$(timestamp_utc)"
  inventory="${RESULTS}/NC2024_NSCLC_FULL_COHORT_PREFLIGHT_${ts}.inventory.txt"
  {
    echo "## RAW_ROOT"
    find "${RAW_ROOT}" -maxdepth 2 | sort
    echo
    echo "## FULL_COHORT"
    find "${FULL_ROOT}" -maxdepth 3 2>/dev/null | sort || true
    echo
    echo "## RESULTS"
    find "${RESULTS}" -maxdepth 1 -mindepth 1 | sort
    echo
    echo "## TMP"
    find /tmp -maxdepth 1 \( -name 'prepare_emtab13526*' -o -name 'find_bad_emtab*' -o -name 'test_lazy_to_memory*' \) | sort
    echo
    echo "## SCRIPTS"
    find "${REPO_ROOT}/scripts" -maxdepth 1 \( -name '*emtab13526*' -o -name '*nc2024*' \) | sort
  } > "${inventory}"
  log "WROTE_PREFLIGHT_INVENTORY=${inventory}"
}

cleanup_derived_artifacts() {
  log "START_CLEANUP"
  rm -rf \
    "${RAW_ROOT}/_prepared_parts" \
    "${RAW_ROOT}/_prepared_parts_clean" \
    "${RAW_ROOT}/tumor/prepared_input.zarr" \
    "${RAW_ROOT}/tumor/prepared_input.ready" \
    "${FULL_ROOT}/_prepared_parts_zarr" \
    "${FULL_ROOT}/prepared_input.zarr" \
    "${FULL_ROOT}/prepared_input.ready"

  find "${RESULTS}" -maxdepth 1 -type d -name 'NC2024_NSCLC_TUMOR_*' -exec rm -rf {} +
  find "${RESULTS}" -maxdepth 1 -type f \( \
    -name 'NC2024_NSCLC_TUMOR*.launch.log' -o \
    -name 'NC2024_NSCLC_TUMOR_PIPELINE.launch.log' \
  \) -delete
  find /tmp -maxdepth 1 \( \
    -name 'prepare_emtab13526*' -o \
    -name 'find_bad_emtab*' -o \
    -name 'test_lazy_to_memory*' \
  \) -exec rm -rf {} +
  log "END_CLEANUP"
}

validate_prepare() {
  local tmp
  tmp="$(mktemp)"
  set +e
  ${CONDA_RUN} python "${PREP_SCRIPT}" --validate-only >"${tmp}" 2>&1
  local ec=$?
  set -e
  cat "${tmp}" >> "${PREP_LOG}"
  rm -f "${tmp}"
  return "${ec}"
}

run_prepare() {
  log "RUN_PREPARE"
  : > "${PREP_LOG}"
  ${CONDA_RUN} python "${PREP_SCRIPT}" >> "${PREP_LOG}" 2>&1
  if [[ ! -f "${READY_SENTINEL}" ]]; then
    log "FINAL_STATUS=STOP_NO_FALLBACK"
    log "PREPARE_FAILED_MISSING_READY"
    exit 1
  fi
  if ! validate_prepare; then
    log "FINAL_STATUS=STOP_NO_FALLBACK"
    log "PREPARE_VALIDATE_ONLY_FAILED"
    exit 1
  fi
}

launch_mode() {
  local mode="$1"
  local project="$2"
  local tmp run_dir ec

  tmp="$(mktemp)"
  set +e
  "${STAGE1_LAUNCHER}" "${mode}" "${project}" >"${tmp}" 2>&1
  ec=$?
  set -e
  cat "${tmp}" >> "${CTRL_LOG}"
  run_dir="$(grep '^RUN_DIR=' "${tmp}" | tail -n 1 | cut -d= -f2- || true)"
  rm -f "${tmp}"
  printf '%s\n' "${run_dir}"
  return "${ec}"
}

run_successful() {
  local run_dir="$1"
  if [[ -z "${run_dir}" || ! -d "${run_dir}" ]]; then
    return 1
  fi

  RUN_DIR_TO_VALIDATE="${run_dir}" ${CONDA_RUN} python - <<'PY'
from pathlib import Path
import pandas as pd
import sys
import os

run_dir = Path(os.environ["RUN_DIR_TO_VALIDATE"])
required_modules = [
    "cellranger",
    "qc",
    "doublet_detection",
    "clustering",
    "annotation",
    "composition",
    "immune_phenotyping",
    "tumor_microenvironment",
]
required_dirs = [
    "clustering",
    "annotation",
    "composition",
    "immune_phenotyping",
    "tumor_microenvironment",
]
required_tables = [
    ("annotation", "cell_type_annotation.csv"),
    ("composition", "composition_proportions.csv"),
    ("immune_phenotyping", "immune_phenotyping.csv"),
    ("tumor_microenvironment", "tme_scores_per_cell.csv"),
]

manifest = run_dir / "run_manifest.json"
status_csv = run_dir / "module_status.csv"
if not manifest.exists() or not status_csv.exists():
    sys.exit(1)

df = pd.read_csv(status_csv)
status_map = dict(zip(df["module"].astype(str), df["status"].astype(str)))
for module in required_modules:
    if status_map.get(module) != "ok":
        sys.exit(1)

for module in required_dirs:
    if not (run_dir / module).is_dir():
        sys.exit(1)

for module, filename in required_tables:
    if not (run_dir / module / filename).exists():
        sys.exit(1)
PY
}

is_capacity_failure() {
  local run_dir="$1"
  local launch_log="$2"

  if [[ -n "${launch_log}" && -f "${launch_log}" ]] && grep -Eq 'Killed|MemoryError|memory exhaustion|out of memory|Cannot allocate memory|allocator' "${launch_log}"; then
    return 0
  fi

  if [[ -n "${run_dir}" && -d "${run_dir}" ]]; then
    if grep -RIEq 'Killed|MemoryError|memory exhaustion|out of memory|Cannot allocate memory|allocator' "${run_dir}" 2>/dev/null; then
      return 0
    fi
  fi

  return 1
}

main() {
  local large_run_dir massive_run_dir large_ec massive_ec
  mkdir -p "${RESULTS}"
  : > "${CTRL_LOG}"

  capture_preflight_inventory
  cleanup_derived_artifacts

  if [[ ! -d "${PREPARED_ZARR}" ]] || [[ ! -f "${READY_SENTINEL}" ]]; then
    run_prepare
  elif ! validate_prepare; then
    run_prepare
  else
    log "REUSING_EXISTING_PREPARED_INPUT=${PREPARED_ZARR}"
  fi

  if [[ ! -d "${PREPARED_ZARR}" ]]; then
    log "FINAL_STATUS=STOP_NO_FALLBACK"
    log "PREPARED_INPUT_MISSING_AFTER_PREPARE"
    exit 1
  fi
  if [[ ! -f "${READY_SENTINEL}" ]]; then
    log "FINAL_STATUS=STOP_NO_FALLBACK"
    log "READY_SENTINEL_MISSING_AFTER_PREPARE"
    exit 1
  fi

  set +e
  large_run_dir="$(launch_mode large "${LARGE_PROJECT}")"
  large_ec=$?
  set -e

  if [[ ${large_ec} -eq 0 ]] && run_successful "${large_run_dir}"; then
    log "FINAL_STATUS=SUCCESS_LARGE"
    exit 0
  fi

  if [[ ${large_ec} -eq 137 || ${large_ec} -eq 143 ]] || is_capacity_failure "${large_run_dir}" "${RESULTS}/${LARGE_PROJECT}.launch.log"; then
    log "LARGE_FAILED_CAPACITY_PROMOTING_TO_MASSIVE"
    set +e
    massive_run_dir="$(launch_mode massive "${MASSIVE_PROJECT}")"
    massive_ec=$?
    set -e
    if [[ ${massive_ec} -eq 0 ]] && run_successful "${massive_run_dir}"; then
      log "FINAL_STATUS=SUCCESS_MASSIVE"
      exit 0
    fi
    log "FINAL_STATUS=STOP_AFTER_MASSIVE_FAILURE"
    exit "${massive_ec:-1}"
  fi

  log "FINAL_STATUS=STOP_NO_FALLBACK"
  exit "${large_ec:-1}"
}

main "$@"
