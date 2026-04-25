#!/usr/bin/env bash
set -euo pipefail
exec ssh -tt ubuntu-tail "$@"
