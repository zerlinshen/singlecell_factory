#!/usr/bin/env bash
set -euo pipefail

REMOTE_ROOT="/home/zerlinshen/singlecell_factory/ops/before_every_run"
LOCAL_ROOT="/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run"

mkdir -p "$LOCAL_ROOT"
rsync -av --delete "ubuntu-tail:${REMOTE_ROOT}/" "${LOCAL_ROOT}/"
echo "Synced ${REMOTE_ROOT} -> ${LOCAL_ROOT}"
