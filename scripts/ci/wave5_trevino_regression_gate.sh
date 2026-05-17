#!/usr/bin/env bash
# Wave-5 Trevino regression gate (US-W5-9).
#
# Reads a metrics.json produced by the Wave-5 pipeline run and enforces three
# thresholds. Exits 0 on PASS, 1 on threshold breach, 2 on input error.
#
# Usage:
#   bash scripts/ci/wave5_trevino_regression_gate.sh [path/to/metrics.json]
#
# Defaults to .omc/research/wave5/metrics.json. Companion allowlist for the
# R1 high-IF reference grep check lives at scripts/ci/wave5_new_refs.txt.
#
# Wave-5 spec: .omc/specs/deep-interview-wave5-completion.md
# Wave-5 plan: .omc/plans/wave5-completion-consensus-2026-05-16.md
set -euo pipefail

M="${1:-.omc/research/wave5/metrics.json}"
if [[ ! -f "$M" ]]; then
    echo "FAIL: metrics file not found: $M" >&2
    exit 2
fi

python - "$M" <<'PY'
import json
import sys

path = sys.argv[1]
m = json.load(open(path))

checks = [
    ("cell_type_ari",              0.70, "US-W5-5 cell-type ARI (joint Leiden vs Trevino)"),
    ("peak_gene_overlap_top1000",  0.70, "US-W5-7 peak-gene top-1000 overlap vs Trevino S2F"),
    ("pseudotime_spearman",        0.70, "US-W5-8 pseudotime Spearman vs Trevino Fig 4D"),
]

failed = []
missing = []
for key, threshold, label in checks:
    value = m.get(key)
    if value is None:
        missing.append((key, label))
        continue
    if value < threshold:
        failed.append((key, value, threshold, label))

if missing:
    for key, label in missing:
        print(f"MISSING {key}  ({label})")
    sys.exit(2)

if failed:
    for key, value, threshold, label in failed:
        print(f"FAIL    {key}={value:.3f} < {threshold:.2f}  ({label})")
    sys.exit(1)

reported = {k: round(m[k], 3) for k, _, _ in checks}
print(f"PASS    {reported}")
sys.exit(0)
PY
