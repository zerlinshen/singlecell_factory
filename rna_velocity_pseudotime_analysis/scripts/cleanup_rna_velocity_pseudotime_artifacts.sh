#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUNTIME_DIR="${ROOT_DIR}/rna_velocity_pseudotime_analysis/runtime"
KEEP_DAYS=7
PURGE_ALL=false

while [ $# -gt 0 ]; do
  case "$1" in
    --keep-days) KEEP_DAYS="$2"; shift ;;
    --keep-days=*) KEEP_DAYS="${1#*=}" ;;
    --all) PURGE_ALL=true ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done

if [ ! -d "${RUNTIME_DIR}" ]; then
  echo "No runtime directory found: ${RUNTIME_DIR}"
  exit 0
fi

if [ "${PURGE_ALL}" = true ]; then
  find "${RUNTIME_DIR}" -mindepth 1 -maxdepth 1 ! -name ".gitkeep" -exec rm -rf {} +
  echo "All RNA Velocity and Pseudotime artifacts removed."
  exit 0
fi

find "${RUNTIME_DIR}" -type f -mtime +"${KEEP_DAYS}" -delete
find "${RUNTIME_DIR}" -type d -empty -delete
echo "Artifacts older than ${KEEP_DAYS} days were cleaned."
