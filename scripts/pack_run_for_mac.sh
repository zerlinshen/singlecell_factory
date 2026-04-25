#!/usr/bin/env bash
# pack_run_for_mac.sh <run_dir> [--schema-version v1|v2]
#
# Exports a v2 R bundle (default) from the run's final_adata.h5ad, then
# produces three tarballs + SHA256 sidecars under <run_dir>/mac_assets/:
#   bundle.tgz   — the exported R bundle directory
#   figures.tgz  — all PNG and PDF files found under the run dir
#   report.tgz   — module_status.csv, run_manifest.json, ledger JSON(s),
#                  parity_summary.json (if present)
#
# Also writes <run_dir>/mac_assets/mac_pull_command.sh — run this on the Mac
# to pull the tarballs via scp.
set -euo pipefail

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <run_dir> [--schema-version v1|v2]" >&2
  exit 1
fi

RUN_DIR="$(realpath "$1")"
SCHEMA_VERSION="v2"

shift
while [[ $# -gt 0 ]]; do
  case "$1" in
    --schema-version)
      SCHEMA_VERSION="$2"; shift 2 ;;
    *)
      echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

if [[ ! -d "$RUN_DIR" ]]; then
  echo "Run directory not found: $RUN_DIR" >&2
  exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"

ASSETS_DIR="${RUN_DIR}/mac_assets"
BUNDLE_DIR="${RUN_DIR}/r_bundle"
H5AD="${RUN_DIR}/final_adata.h5ad"

mkdir -p "$ASSETS_DIR"

# ---------------------------------------------------------------------------
# Step 1 — export R bundle
# ---------------------------------------------------------------------------
if [[ ! -f "$H5AD" ]]; then
  echo "WARNING: final_adata.h5ad not found at ${H5AD}. Skipping bundle export." >&2
else
  echo "[pack_run_for_mac] Exporting R bundle (schema=${SCHEMA_VERSION}) ..."
  python "${REPO_ROOT}/scripts/export_singlecell_r_bundle.py" \
    --input    "$H5AD" \
    --output   "$BUNDLE_DIR" \
    --schema-version "$SCHEMA_VERSION" \
    --format   auto
  echo "[pack_run_for_mac] Bundle exported to: $BUNDLE_DIR"
fi

# ---------------------------------------------------------------------------
# Step 2 — tarball helper
# ---------------------------------------------------------------------------
make_tarball() {
  local name="$1"       # e.g. bundle
  local tarball="${ASSETS_DIR}/${name}.tgz"
  shift
  # Remaining args are passed directly to tar
  tar -czf "$tarball" "$@"
  sha256sum "$tarball" > "${tarball}.sha256"
  echo "[pack_run_for_mac] ${name}.tgz  $(du -sh "$tarball" | cut -f1)  $(cat "${tarball}.sha256" | awk '{print $1}')"
}

# ---------------------------------------------------------------------------
# Step 3 — bundle.tgz
# ---------------------------------------------------------------------------
if [[ -d "$BUNDLE_DIR" ]]; then
  make_tarball bundle -C "$RUN_DIR" r_bundle
else
  echo "[pack_run_for_mac] WARNING: bundle dir not found, skipping bundle.tgz" >&2
fi

# ---------------------------------------------------------------------------
# Step 4 — figures.tgz
# ---------------------------------------------------------------------------
FIGURE_LIST_FILE="$(mktemp)"
find "$RUN_DIR" -not -path "${ASSETS_DIR}/*" \( -name "*.png" -o -name "*.pdf" \) \
  > "$FIGURE_LIST_FILE" 2>/dev/null || true

if [[ -s "$FIGURE_LIST_FILE" ]]; then
  tar -czf "${ASSETS_DIR}/figures.tgz" -T "$FIGURE_LIST_FILE" 2>/dev/null || \
    tar -czf "${ASSETS_DIR}/figures.tgz" --files-from="$FIGURE_LIST_FILE" 2>/dev/null || true
  sha256sum "${ASSETS_DIR}/figures.tgz" > "${ASSETS_DIR}/figures.tgz.sha256"
  echo "[pack_run_for_mac] figures.tgz  $(du -sh "${ASSETS_DIR}/figures.tgz" | cut -f1)  $(cat "${ASSETS_DIR}/figures.tgz.sha256" | awk '{print $1}')"
else
  echo "[pack_run_for_mac] WARNING: no figures found under $RUN_DIR" >&2
fi
rm -f "$FIGURE_LIST_FILE"

# ---------------------------------------------------------------------------
# Step 5 — report.tgz
# ---------------------------------------------------------------------------
REPORT_STAGING="$(mktemp -d)"
for f in module_status.csv run_manifest.json parity_summary.json; do
  [[ -f "${RUN_DIR}/${f}" ]] && cp "${RUN_DIR}/${f}" "${REPORT_STAGING}/" || true
done

# ops/run_ledger/<project>_*.json
LEDGER_DIR="${REPO_ROOT}/ops/run_ledger"
if [[ -d "$LEDGER_DIR" ]]; then
  # Infer project name from run dir basename (strip trailing timestamp if present)
  PROJECT_NAME="$(basename "$RUN_DIR" | sed 's/_[0-9]\{8\}T[0-9]\{6\}.*$//' | sed 's/_[0-9]\{8\}$//')"
  find "$LEDGER_DIR" -name "${PROJECT_NAME}_*.json" -exec cp {} "${REPORT_STAGING}/" \; 2>/dev/null || true
fi

if ls "${REPORT_STAGING}"/* &>/dev/null; then
  tar -czf "${ASSETS_DIR}/report.tgz" -C "$REPORT_STAGING" .
  sha256sum "${ASSETS_DIR}/report.tgz" > "${ASSETS_DIR}/report.tgz.sha256"
  echo "[pack_run_for_mac] report.tgz   $(du -sh "${ASSETS_DIR}/report.tgz" | cut -f1)  $(cat "${ASSETS_DIR}/report.tgz.sha256" | awk '{print $1}')"
else
  echo "[pack_run_for_mac] WARNING: no report files found, skipping report.tgz" >&2
fi
rm -rf "$REPORT_STAGING"

# ---------------------------------------------------------------------------
# Step 6 — generate mac_pull_command.sh
# ---------------------------------------------------------------------------
HOSTNAME_FULL="$(hostname -f 2>/dev/null || hostname)"
REMOTE_USER="${USER:-ubuntu}"

# Build file list of what actually exists
PULL_FILES=()
for f in bundle.tgz figures.tgz report.tgz; do
  [[ -f "${ASSETS_DIR}/${f}" ]] && PULL_FILES+=("$f" "${f}.sha256")
done

PULL_CMD_PATH="${ASSETS_DIR}/mac_pull_command.sh"
cat > "$PULL_CMD_PATH" <<PULLEOF
#!/usr/bin/env bash
# Run this script ON THE MAC to pull the packed assets from the Linux host.
# Generated by pack_run_for_mac.sh on $(date -u +"%Y-%m-%dT%H:%M:%SZ")
set -euo pipefail

REMOTE_HOST="\${1:-${REMOTE_USER}@${HOSTNAME_FULL}}"
REMOTE_ASSETS="${ASSETS_DIR}"
LOCAL_DIR="\${2:-\${HOME}/scrna_remote_pulls/$(basename "$RUN_DIR")}"

mkdir -p "\$LOCAL_DIR"

FILES=(${PULL_FILES[*]:-})
for f in "\${FILES[@]}"; do
  echo "[mac_pull] Pulling \$f ..."
  scp "\${REMOTE_HOST}:\${REMOTE_ASSETS}/\${f}" "\${LOCAL_DIR}/\${f}"
done

echo "[mac_pull] Verifying SHA256 ..."
for tgz in bundle.tgz figures.tgz report.tgz; do
  sidecar="\${LOCAL_DIR}/\${tgz}.sha256"
  tarball="\${LOCAL_DIR}/\${tgz}"
  if [[ -f "\$sidecar" ]] && [[ -f "\$tarball" ]]; then
    expected="\$(awk '{print \$1}' "\$sidecar")"
    actual="\$(sha256sum "\$tarball" | awk '{print \$1}')"
    if [[ "\$expected" == "\$actual" ]]; then
      echo "  OK  \$tgz"
    else
      echo "  FAIL \$tgz: expected \$expected, got \$actual" >&2
      exit 1
    fi
  fi
done

echo "[mac_pull] All files pulled to: \$LOCAL_DIR"
PULLEOF

chmod +x "$PULL_CMD_PATH"
echo "[pack_run_for_mac] Pull command written to: $PULL_CMD_PATH"
echo ""
echo "To pull on Mac, run:"
echo "  bash ${PULL_CMD_PATH} [remote_host] [local_dir]"
echo ""
echo "[pack_run_for_mac] Done. Assets in: $ASSETS_DIR"
