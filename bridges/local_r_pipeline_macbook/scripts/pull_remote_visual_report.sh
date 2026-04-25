#!/usr/bin/env bash
# pull_remote_visual_report.sh
#
# Without --bundle (original behavior):
#   pull individual visual report assets via scp.
#   Usage: $0 <remote_result_dir> <local_tag>
#
# With --bundle mode:
#   pull the 3 mac_assets tarballs from <run_dir>/mac_assets/ and verify SHA256.
#   Usage: $0 --bundle --remote-run <remote_run_dir> [--remote-host <host>]
#            [--local-dir <local_dir>]
set -euo pipefail

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
BUNDLE_MODE=0
REMOTE_HOST="ubuntu-tail"
REMOTE_RUN=""
LOCAL_DIR_OVERRIDE=""

# Detect bundle vs legacy positional mode
if [[ $# -ge 1 ]] && [[ "$1" == "--bundle" ]]; then
  BUNDLE_MODE=1
  shift
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --remote-host) REMOTE_HOST="$2"; shift 2 ;;
      --remote-run)  REMOTE_RUN="$2";  shift 2 ;;
      --local-dir)   LOCAL_DIR_OVERRIDE="$2"; shift 2 ;;
      *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
  done
fi

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

# ---------------------------------------------------------------------------
# Bundle mode
# ---------------------------------------------------------------------------
if [[ "$BUNDLE_MODE" -eq 1 ]]; then
  if [[ -z "$REMOTE_RUN" ]]; then
    echo "Usage: $0 --bundle --remote-run <remote_run_dir> [--remote-host <host>] [--local-dir <dir>]" >&2
    exit 1
  fi

  REMOTE_ASSETS="${REMOTE_RUN}/mac_assets"
  RUN_BASENAME="$(basename "$REMOTE_RUN")"

  if [[ -n "$LOCAL_DIR_OVERRIDE" ]]; then
    LOCAL_DIR="$LOCAL_DIR_OVERRIDE"
  else
    LOCAL_DIR="${ROOT}/output/${RUN_BASENAME}_bundle_assets"
  fi
  mkdir -p "$LOCAL_DIR"

  echo "[pull_remote] Bundle mode: pulling from ${REMOTE_HOST}:${REMOTE_ASSETS}"

  # Pull each tarball and its SHA256 sidecar
  TARBALLS=(bundle.tgz figures.tgz report.tgz)
  PULLED=()

  for tgz in "${TARBALLS[@]}"; do
    remote_path="${REMOTE_ASSETS}/${tgz}"
    sidecar_path="${REMOTE_ASSETS}/${tgz}.sha256"
    local_path="${LOCAL_DIR}/${tgz}"
    local_sidecar="${LOCAL_DIR}/${tgz}.sha256"

    # Pull tarball (skip silently if absent on remote)
    if scp "${REMOTE_HOST}:${remote_path}" "$local_path" 2>/dev/null; then
      PULLED+=("$tgz")
      # Pull sidecar
      scp "${REMOTE_HOST}:${sidecar_path}" "$local_sidecar" 2>/dev/null || true
    else
      echo "[pull_remote] WARNING: ${tgz} not found on remote, skipping." >&2
    fi
  done

  if [[ ${#PULLED[@]} -eq 0 ]]; then
    echo "[pull_remote] ERROR: no tarballs pulled from ${REMOTE_HOST}:${REMOTE_ASSETS}" >&2
    exit 1
  fi

  # Verify SHA256 of each pulled tarball
  echo "[pull_remote] Verifying SHA256 ..."
  FAIL=0
  for tgz in "${PULLED[@]}"; do
    local_path="${LOCAL_DIR}/${tgz}"
    local_sidecar="${LOCAL_DIR}/${tgz}.sha256"
    if [[ -f "$local_sidecar" ]]; then
      expected="$(awk '{print $1}' "$local_sidecar")"
      actual="$(sha256sum "$local_path" | awk '{print $1}')"
      if [[ "$expected" == "$actual" ]]; then
        echo "  OK    $tgz  ($actual)"
      else
        echo "  FAIL  $tgz  expected=$expected  actual=$actual" >&2
        FAIL=1
      fi
    else
      echo "  SKIP  $tgz  (no sidecar)"
    fi
  done

  if [[ "$FAIL" -ne 0 ]]; then
    echo "[pull_remote] SHA256 verification FAILED for one or more tarballs." >&2
    exit 1
  fi

  echo "[pull_remote] Done. Assets in: $LOCAL_DIR"
  echo "$LOCAL_DIR"
  exit 0
fi

# ---------------------------------------------------------------------------
# Legacy positional mode (original behavior — unchanged)
# ---------------------------------------------------------------------------
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <remote_result_dir> <local_tag>" >&2
  exit 1
fi
REMOTE_RESULT_DIR="$1"
LOCAL_TAG="$2"
LOCAL_DIR="$ROOT/output/${LOCAL_TAG}_remote_report_assets"
mkdir -p "$LOCAL_DIR"
for rel in \
  run_manifest.json \
  module_status.csv \
  clustering/umap_leiden.png \
  annotation/umap_cell_type.png \
  annotation/umap_annotation_confidence.png \
  annotation/cell_type_composition.png \
  composition/composition_barplot.png \
  composition/composition_boxplot.png \
  immune_phenotyping/umap_immune_subtype.png \
  immune_phenotyping/immune_signature_heatmap.png \
  tumor_microenvironment/tme_cyt_umap.png \
  tumor_microenvironment/tme_tis_umap.png \
  tumor_microenvironment/tme_signature_heatmap.png \
  qc/qc_violin_post_filter.png \
  qc/qc_scatter_post_filter.png
  do
  mkdir -p "$LOCAL_DIR/$(dirname "$rel")"
  scp "ubuntu-tail:${REMOTE_RESULT_DIR}/${rel}" "$LOCAL_DIR/${rel}" 2>/dev/null || true
 done
 echo "$LOCAL_DIR"
