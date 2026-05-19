#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="${REPO_ROOT}/data/raw/10x_nsclc_900k_flex"
TAR="${ROOT}/16plex_900k_32_NSCLC_multiplex_count_filtered_feature_bc_matrix.tar.gz"
OUTDIR="${ROOT}/outs/filtered_feature_bc_matrix"
mkdir -p "$OUTDIR"
if [[ ! -f "$TAR" ]]; then
  echo "Missing tarball: $TAR" >&2
  exit 1
fi
tar -tzf "$TAR" >/dev/null
if [[ -f "$OUTDIR/matrix.mtx.gz" ]]; then
  rm -rf "$OUTDIR"
  mkdir -p "$OUTDIR"
fi
tar -xzf "$TAR" -C "$OUTDIR" --strip-components=1
ls -lh "$OUTDIR"
