#!/usr/bin/env bash
set -euo pipefail

ALLOW_FALLBACK=0
for arg in "$@"; do
    case "$arg" in
        --allow-fallback) ALLOW_FALLBACK=1 ;;
        *) ;;
    esac
done

if [[ "$ALLOW_FALLBACK" -eq 0 ]]; then
    echo "Lockfile hook hard-stop default. Use --allow-fallback to engage pip-freeze fallback (requires journal entry)."
    exit 1
fi

ALLOWLIST="environments/wave5_test_deps_allowlist.json"
LOCKFILE="environments/sc_gpu_stable.lock.txt"
FREEZE_OUT="/tmp/wave5_v4_2_freeze_candidate.txt"
JOURNAL_DIR="ops/before_every_run/journal"

if [[ ! -f "$ALLOWLIST" ]]; then
    echo "FAIL allowlist not found: $ALLOWLIST" >&2
    exit 2
fi
if [[ ! -f "$LOCKFILE" ]]; then
    echo "FAIL lockfile not found: $LOCKFILE" >&2
    exit 2
fi

pip freeze --exclude-editable > "$FREEZE_OUT"

python - "$ALLOWLIST" "$LOCKFILE" "$FREEZE_OUT" "$JOURNAL_DIR" <<'PY'
import json
import os
import re
import sys
from pathlib import Path

allowlist_path, lockfile_path, freeze_path, journal_dir = sys.argv[1:5]

with open(allowlist_path) as f:
    allowlist = json.load(f)

direct = set(allowlist.get("direct_deps", []) or [])
transitive = set(allowlist.get("transitive_deps", []) or [])
allowed = {name.lower() for name in (direct | transitive)}

def parse_pkg_name(line: str) -> str:
    line = line.strip()
    if not line or line.startswith("#"):
        return ""
    # Strip env markers, extras.
    line = line.split(";", 1)[0].strip()
    # Split on common version operators.
    for op in ("==", ">=", "<=", "~=", "!=", ">", "<", "==="):
        if op in line:
            return line.split(op, 1)[0].strip().lower()
    return line.lower()

def pkg_set(path: str) -> set:
    out = set()
    with open(path) as f:
        for raw in f:
            name = parse_pkg_name(raw)
            if name:
                out.add(name)
    return out

lock_pkgs = pkg_set(lockfile_path)
freeze_pkgs = pkg_set(freeze_path)
delta = freeze_pkgs - lock_pkgs

unauthorized = sorted(p for p in delta if p not in allowed)
if unauthorized:
    first = unauthorized[0]
    print(f"FAIL package {first} not in wave5_test_deps_allowlist.json; open follow-up ticket per Plan v4.2 Section 8 Scenario C.", file=sys.stderr)
    if len(unauthorized) > 1:
        print(f"FAIL additional unauthorized deltas: {unauthorized[1:]}", file=sys.stderr)
    sys.exit(2)

# Check journal rationale.
jdir = Path(journal_dir)
rationale_found = False
if jdir.is_dir():
    entries = sorted(
        (p for p in jdir.iterdir() if p.is_file()),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )[:5]
    for entry in entries:
        try:
            text = entry.read_text(errors="ignore")
        except OSError:
            continue
        if "--allow-fallback" in text:
            rationale_found = True
            break

if not rationale_found:
    print("FAIL --allow-fallback used without journal rationale entry.", file=sys.stderr)
    sys.exit(2)

print("PASS lockfile fallback within allowlist; rationale documented.")
sys.exit(0)
PY
