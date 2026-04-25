#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 <scale_mode> <project_name>" >&2
  exit 2
fi

SCALE_MODE="$1"
PROJECT_NAME="$2"
RAW_ROOT="/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526"
SAMPLE_ROOT="${RAW_ROOT}/full_cohort"
RESULTS="/home/zerlinshen/singlecell_factory/results"
REPO_ROOT="/home/zerlinshen/singlecell_factory"
CONDA_RUN="/home/zerlinshen/conda/bin/conda run -n sc_gpu"
PROJECT_LOG="${RESULTS}/${PROJECT_NAME}.launch.log"
BEFORE_FILE="$(mktemp)"
AFTER_FILE="$(mktemp)"

case "${SCALE_MODE}" in
  large|massive)
    ;;
  *)
    echo "invalid scale_mode: ${SCALE_MODE}" >&2
    exit 2
    ;;
esac

mkdir -p "${RESULTS}"
: > "${PROJECT_LOG}"
find "${RESULTS}" -maxdepth 1 -type d -name "${PROJECT_NAME}_*" | sort > "${BEFORE_FILE}"

set +e
cd "${REPO_ROOT}"
${CONDA_RUN} python -m workflow.modular.cli \
  --project "${PROJECT_NAME}" \
  --sample-root "${SAMPLE_ROOT}" \
  --output-dir "${RESULTS}" \
  --optional-modules clustering,annotation,composition,immune_phenotyping,tumor_microenvironment \
  --batch-key sample \
  --checkpoint \
  --gpu-mode auto \
  --scale-mode "${SCALE_MODE}" > "${PROJECT_LOG}" 2>&1
EC=$?
set -e

find "${RESULTS}" -maxdepth 1 -type d -name "${PROJECT_NAME}_*" | sort > "${AFTER_FILE}"
RUN_DIR="$(comm -13 "${BEFORE_FILE}" "${AFTER_FILE}" | tail -n 1 || true)"
if [[ -z "${RUN_DIR}" && ${EC} -eq 0 ]]; then
  RUN_DIR="$(tail -n 1 "${AFTER_FILE}" || true)"
fi
rm -f "${BEFORE_FILE}" "${AFTER_FILE}"
echo "RUN_DIR=${RUN_DIR}"
exit "${EC}"
