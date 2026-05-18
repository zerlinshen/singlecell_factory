#!/usr/bin/env bash
# US-W4-11 — Wave 4 freeze sentinel.
#
# Writes .omc/state/wave4_800k_dispatch.lock when 800k dispatch enters; removes
# on completion. Consumed by scripts/ci/check_freeze.sh in pre-commit / CI to
# refuse changes touching workflow/, pyproject.toml, environments/ during an
# in-flight 800k run.
#
# Usage:
#   scripts/ci/wave4_dispatch_freeze.sh start <run_id>   # write lock
#   scripts/ci/wave4_dispatch_freeze.sh end <run_id>     # remove lock
set -euo pipefail

CMD="${1:-status}"
RUN_ID="${2:-}"
LOCK_FILE=".omc/state/wave4_800k_dispatch.lock"

mkdir -p .omc/state

case "$CMD" in
  start)
    if [[ -z "$RUN_ID" ]]; then echo "usage: $0 start <run_id>" >&2; exit 2; fi
    now=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    exp=$(date -u -d '+6 hours' +%Y-%m-%dT%H:%M:%SZ)
    printf '{"run_id":"%s","started_at":"%s","expected_done_by":"%s"}\n' "$RUN_ID" "$now" "$exp" > "$LOCK_FILE"
    echo "Wave 4 freeze STARTED: $RUN_ID (expires $exp)"
    ;;
  end)
    if [[ -f "$LOCK_FILE" ]]; then
      rm -f "$LOCK_FILE"
      echo "Wave 4 freeze CLEARED"
    else
      echo "Wave 4 freeze already clear"
    fi
    ;;
  status)
    if [[ -f "$LOCK_FILE" ]]; then
      cat "$LOCK_FILE"
    else
      echo "no active 800k dispatch lock"
    fi
    ;;
  *)
    echo "usage: $0 {start <run_id> | end <run_id> | status}" >&2; exit 2;;
esac
