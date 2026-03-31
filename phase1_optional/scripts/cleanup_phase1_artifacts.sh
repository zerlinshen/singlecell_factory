#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec bash "${ROOT_DIR}/rna_velocity_pseudotime_analysis/scripts/cleanup_rna_velocity_pseudotime_artifacts.sh" "$@"
