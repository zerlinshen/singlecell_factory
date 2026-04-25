#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <remote_result_dir> <remote_out_dir_or_tag> [group_by] [cluster_by] [max_cells]" >&2
  echo "This compatibility wrapper now runs R plotting on the remote server; it no longer pulls data to a local Mac R pipeline." >&2
  exit 1
fi

REMOTE_RESULT_DIR="${1%/}"
OUT_ARG="$2"
GROUP_BY="${3:-cell_type}"
CLUSTER_BY="${4:-leiden}"
MAX_CELLS="${5:-200000}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

if [[ "$OUT_ARG" = /* ]]; then
  REMOTE_OUT_DIR="$OUT_ARG"
else
  REMOTE_OUT_DIR="${REMOTE_RESULT_DIR}/r_plots/${OUT_ARG}"
fi

echo "Compatibility note: plotting is remote-only now. Delegating to run_remote_bundle_plot.sh" >&2
exec "$SCRIPT_DIR/run_remote_bundle_plot.sh" \
  "$REMOTE_RESULT_DIR" \
  "$REMOTE_OUT_DIR" \
  "$GROUP_BY" \
  "$CLUSTER_BY" \
  "$MAX_CELLS"
