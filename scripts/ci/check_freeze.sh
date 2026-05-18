#!/usr/bin/env bash
# US-W4-11 — Wave 4 freeze enforcement.
#
# Refuses commits / CI runs that touch workflow/, pyproject.toml, or
# environments/ while .omc/state/wave4_800k_dispatch.lock exists. Docs / tests
# / bridges changes remain allowed.
#
# Usage (pre-commit hook):
#   scripts/ci/check_freeze.sh                       # checks staged files
#   scripts/ci/check_freeze.sh path1 path2 ...       # checks given paths
set -euo pipefail

LOCK_FILE=".omc/state/wave4_800k_dispatch.lock"
if [[ ! -f "$LOCK_FILE" ]]; then
  exit 0  # no freeze active
fi

if [[ $# -eq 0 ]]; then
  mapfile -t paths < <(git diff --cached --name-only --diff-filter=ACMR 2>/dev/null || true)
else
  paths=("$@")
fi

violators=()
for p in "${paths[@]}"; do
  case "$p" in
    workflow/*|pyproject.toml|environments/*)
      violators+=("$p")
      ;;
  esac
done

if (( ${#violators[@]} > 0 )); then
  echo "ERROR: Wave 4 freeze active. The following files cannot be modified during 800k dispatch:" >&2
  for v in "${violators[@]}"; do echo "  - $v" >&2; done
  echo "" >&2
  echo "Lock contents:" >&2
  cat "$LOCK_FILE" >&2
  echo "" >&2
  echo "Wait for dispatch end (scripts/ci/wave4_dispatch_freeze.sh end <run_id>)" >&2
  echo "or commit to docs/tests/bridges only." >&2
  exit 1
fi
exit 0
