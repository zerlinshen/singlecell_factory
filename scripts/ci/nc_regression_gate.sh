#!/usr/bin/env bash
#
# nc_regression_gate.sh — Wave 1 P-CC.S26 / Wave 2+ US-009
#
# Compares two NC pipeline runs (baseline vs candidate) on:
#   gate-1  cluster_label_ari        (sklearn.metrics.adjusted_rand_score)
#   gate-2  top_de_markers_jaccard   (set-Jaccard over top-K markers per cluster)
#   gate-3  structural               (h5ad shape/dtype/obs schema parity)
#   gate-4  bundle_schema            (manifest schema_version accepted set)
#
# Thresholds: scripts/ci/nc_regression_tolerance_spec.yaml
# Engine:     scripts/ci/nc_regression_compare.py
#
# Usage:
#   nc_regression_gate.sh --baseline-run-dir PATH --candidate-run-dir PATH
#
# Exit codes (from spec):
#   0   all four sub-gates pass
#   10  ARI sub-gate failed
#   11  Jaccard sub-gate failed
#   12  structural sub-gate failed
#   13  bundle schema sub-gate failed
#   2   missing inputs

set -euo pipefail

usage() {
    cat <<'EOF'
nc_regression_gate.sh — compare baseline vs candidate NC pipeline runs.

Required:
  --baseline-run-dir PATH    Run directory containing the baseline outputs.
  --candidate-run-dir PATH   Run directory containing the candidate outputs.

Optional:
  --conda-env NAME           Python env (default: sc_gpu).
  --help                     Show this help.
EOF
}

BASELINE=""
CANDIDATE=""
CONDA_ENV="sc_gpu"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --baseline-run-dir)   shift; BASELINE="${1:?--baseline-run-dir requires PATH}" ;;
        --candidate-run-dir)  shift; CANDIDATE="${1:?--candidate-run-dir requires PATH}" ;;
        --conda-env)          shift; CONDA_ENV="${1:?--conda-env requires NAME}" ;;
        --help|-h)            usage; exit 0 ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage; exit 2
            ;;
    esac
    shift
done

if [[ -z "$BASELINE" || -z "$CANDIDATE" ]]; then
    echo "ERROR: --baseline-run-dir and --candidate-run-dir are required" >&2
    usage
    exit 2
fi

if [[ ! -d "$BASELINE" ]]; then
    echo "ERROR: baseline run dir does not exist: $BASELINE" >&2
    exit 2
fi
if [[ ! -d "$CANDIDATE" ]]; then
    echo "ERROR: candidate run dir does not exist: $CANDIDATE" >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPARE_PY="$SCRIPT_DIR/nc_regression_compare.py"

if [[ ! -f "$COMPARE_PY" ]]; then
    echo "ERROR: helper missing: $COMPARE_PY" >&2
    exit 2
fi

echo "=== nc_regression_gate ==="
echo "baseline:  $BASELINE"
echo "candidate: $CANDIDATE"
echo "spec:      $SCRIPT_DIR/nc_regression_tolerance_spec.yaml"
echo "engine:    $COMPARE_PY"
echo

# Invoke the python comparator. It owns the spec parsing + per-gate exit codes.
exec conda run -n "$CONDA_ENV" --live-stream python "$COMPARE_PY" \
    --baseline-run-dir "$BASELINE" \
    --candidate-run-dir "$CANDIDATE"
