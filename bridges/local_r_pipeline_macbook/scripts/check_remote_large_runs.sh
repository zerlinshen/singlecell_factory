#!/usr/bin/env bash
set -euo pipefail
REMOTE_SINGLECELL_ROOT="${REMOTE_SINGLECELL_ROOT:-/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory}"
printf -v REMOTE_STATUS_SCRIPT '%q' "${REMOTE_SINGLECELL_ROOT}/scripts/print_large_data_status.sh"
exec ssh ubuntu-tail "bash ${REMOTE_STATUS_SCRIPT}"
