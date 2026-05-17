#!/usr/bin/env bash
# Wave-5 v5.0 ledger schema gate (plan v5.0.4 §3.3 + Stage 0.5).
#
# Usage: wave5_v5_0_schema_gate.sh <ledger.json>
#
# Mirrors the v4.2 gate. Dispatches via scripts/ci/select_ledger_schema.py
# so a v4.2-revision record routed here will EXIT NON-ZERO with a clear
# mutual-rejection message (and vice versa for the v4.2 gate).
#
# Exit codes:
#   0  ledger validates against the dispatched schema
#   1  schema validation failure (mismatch fields, missing required, etc.)
#   2  dispatcher / I/O error (cannot resolve schema)
#   3  mutual rejection — record is v4.2 revision but this gate is v5.0-only

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <ledger.json>" >&2
  exit 2
fi

LEDGER="$1"
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DISPATCHER="$REPO_ROOT/scripts/ci/select_ledger_schema.py"
V50_SCHEMA="$REPO_ROOT/ops/run_ledger/schema/wave5_v5_0.schema.json"

if [[ ! -f "$LEDGER" ]]; then
  echo "ERROR: ledger not found: $LEDGER" >&2
  exit 2
fi

DISPATCHED_SCHEMA="$(python3 "$DISPATCHER" "$LEDGER" --print-schema-path-only 2>&1)" || {
  echo "ERROR: dispatcher failed: $DISPATCHED_SCHEMA" >&2
  exit 2
}

if [[ "$DISPATCHED_SCHEMA" != "$V50_SCHEMA" ]]; then
  echo "MUTUAL REJECTION: ledger dispatched to schema=$DISPATCHED_SCHEMA but this gate enforces v5.0 schema." >&2
  echo "Use scripts/ci/wave5_v4_2_gate.sh for v4.x records." >&2
  exit 3
fi

python3 -c "
import json, sys
import jsonschema
with open('$LEDGER') as f: rec = json.load(f)
with open('$V50_SCHEMA') as f: schema = json.load(f)
validator = jsonschema.Draft202012Validator(schema)
errors = sorted(validator.iter_errors(rec), key=lambda e: e.path)
if not errors:
    print('PASS: ledger validates against wave5_v5_0.schema.json')
    sys.exit(0)
for err in errors:
    path = '/'.join(map(str, err.path))
    print(f'FAIL: {path}: {err.message}', file=sys.stderr)
sys.exit(1)
"
