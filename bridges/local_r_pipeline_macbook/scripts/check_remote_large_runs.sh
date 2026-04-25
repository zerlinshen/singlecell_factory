#!/usr/bin/env bash
set -euo pipefail
exec ssh ubuntu-tail 'bash /home/zerlinshen/singlecell_factory/scripts/print_large_data_status.sh'
