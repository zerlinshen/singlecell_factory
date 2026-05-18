#!/usr/bin/env bash
# sync-contracts.sh — Propagate canonical contract schemas from
# singlecell_factory/contracts/ to r_multiomics_factory/contracts/ and
# plotting_factory/contracts/, and update all .expected_sha256 sentinel files.
#
# Schemas synced:
#   bundle_schema.yaml        → r_multiomics_factory/contracts/ (byte-identical)
#   figure_bundle_schema.yaml → r_multiomics_factory/contracts/ AND
#                               plotting_factory/contracts/ (byte-identical)
#
# Usage:
#   bash tools/sync-contracts.sh
#
# Pre-commit hook installation:
#   Add the following to singlecell_factory/.git/hooks/pre-commit:
#
#     #!/usr/bin/env bash
#     set -euo pipefail
#     bash "$(git rev-parse --show-toplevel)/tools/sync-contracts.sh" || exit 1
#
# The script exits 0 on success and non-zero on any error, including:
#   - missing sibling repos
#   - sha256sum mismatch after copy (guards against fs issues)
#
# CI parity check (no sibling required):
#   Each repo's CI job can run the local parity check independently:
#     actual=$(sha256sum contracts/bundle_schema.yaml | awk '{print $1}')
#     expected=$(grep '^bundle_schema ' contracts/.expected_sha256 | awk '{print $2}')
#     [ "$actual" = "$expected" ] || { echo "contracts parity FAIL: bundle_schema"; exit 1; }

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ---------------------------------------------------------------------------
# Locate sibling repos
# ---------------------------------------------------------------------------
SIBLING_R="$(cd "${REPO_ROOT}/../r_multiomics_factory" 2>/dev/null && pwd)" || true
if [[ -z "${SIBLING_R}" || ! -d "${SIBLING_R}" ]]; then
  echo "ERROR: sibling repo not found at ${REPO_ROOT}/../r_multiomics_factory" >&2
  echo "       Clone r_multiomics_factory alongside singlecell_factory and retry." >&2
  exit 1
fi

SIBLING_PF="$(cd "${REPO_ROOT}/../plotting_factory" 2>/dev/null && pwd)" || true
if [[ -z "${SIBLING_PF}" || ! -d "${SIBLING_PF}" ]]; then
  echo "ERROR: sibling repo not found at ${REPO_ROOT}/../plotting_factory" >&2
  echo "       Clone plotting_factory alongside singlecell_factory and retry." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Helper: copy_and_verify <src> <dst> <label>
# Copies src -> dst, verifies sha256 equality, returns hash.
# ---------------------------------------------------------------------------
copy_and_verify() {
  local src="$1"
  local dst="$2"
  local label="$3"

  if [[ ! -f "${src}" ]]; then
    echo "ERROR: canonical schema not found: ${src}" >&2
    exit 1
  fi

  mkdir -p "$(dirname "${dst}")"
  cp "${src}" "${dst}"

  local src_hash dst_hash
  src_hash=$(sha256sum "${src}" | awk '{print $1}')
  dst_hash=$(sha256sum "${dst}" | awk '{print $1}')

  if [[ "${src_hash}" != "${dst_hash}" ]]; then
    echo "ERROR: sha256 mismatch after copy of ${label} — filesystem issue?" >&2
    echo "  src: ${src_hash}" >&2
    echo "  dst: ${dst_hash}" >&2
    exit 1
  fi

  echo "${src_hash}"
}

# ---------------------------------------------------------------------------
# Sync bundle_schema.yaml → r_multiomics_factory only
# ---------------------------------------------------------------------------
BUNDLE_CANONICAL="${REPO_ROOT}/contracts/bundle_schema.yaml"
BUNDLE_VENDORED_R="${SIBLING_R}/contracts/bundle_schema.yaml"

BUNDLE_HASH=$(copy_and_verify "${BUNDLE_CANONICAL}" "${BUNDLE_VENDORED_R}" "bundle_schema")
echo "bundle_schema sync OK: ${BUNDLE_HASH}"
echo "  canonical: ${BUNDLE_CANONICAL}"
echo "  vendored:  ${BUNDLE_VENDORED_R}"

# ---------------------------------------------------------------------------
# Sync figure_bundle_schema.yaml → r_multiomics_factory AND plotting_factory
# ---------------------------------------------------------------------------
FIG_CANONICAL="${REPO_ROOT}/contracts/figure_bundle_schema.yaml"
FIG_VENDORED_R="${SIBLING_R}/contracts/figure_bundle_schema.yaml"
FIG_VENDORED_PF="${SIBLING_PF}/contracts/figure_bundle_schema.yaml"

FIG_HASH_R=$(copy_and_verify "${FIG_CANONICAL}" "${FIG_VENDORED_R}" "figure_bundle_schema→r_multiomics_factory")
FIG_HASH_PF=$(copy_and_verify "${FIG_CANONICAL}" "${FIG_VENDORED_PF}" "figure_bundle_schema→plotting_factory")

# Paranoia: all three copies must have the same hash.
FIG_HASH_CANONICAL=$(sha256sum "${FIG_CANONICAL}" | awk '{print $1}')
if [[ "${FIG_HASH_CANONICAL}" != "${FIG_HASH_R}" || "${FIG_HASH_CANONICAL}" != "${FIG_HASH_PF}" ]]; then
  echo "ERROR: figure_bundle_schema three-way parity FAIL" >&2
  echo "  canonical:          ${FIG_HASH_CANONICAL}" >&2
  echo "  r_multiomics_factory: ${FIG_HASH_R}" >&2
  echo "  plotting_factory:     ${FIG_HASH_PF}" >&2
  exit 1
fi

echo "figure_bundle_schema sync OK (3-way parity): ${FIG_HASH_CANONICAL}"
echo "  canonical: ${FIG_CANONICAL}"
echo "  vendored:  ${FIG_VENDORED_R}"
echo "  vendored:  ${FIG_VENDORED_PF}"

# ---------------------------------------------------------------------------
# Update .expected_sha256 sentinels in all repos
#
# Format: one line per tracked schema — "<stem> <sha256>"
# This multi-schema format supersedes the single-hash format used in Phase 2.
# Readers checking only bundle_schema may still grep for 'bundle_schema'.
# ---------------------------------------------------------------------------
write_sha_file() {
  local path="$1"
  local bundle_hash="$2"
  local fig_hash="$3"
  printf 'bundle_schema %s\nfigure_bundle_schema %s\n' "${bundle_hash}" "${fig_hash}" > "${path}"
}

write_sha_file "${REPO_ROOT}/contracts/.expected_sha256"       "${BUNDLE_HASH}" "${FIG_HASH_CANONICAL}"
write_sha_file "${SIBLING_R}/contracts/.expected_sha256"       "${BUNDLE_HASH}" "${FIG_HASH_R}"
write_sha_file "${SIBLING_PF}/contracts/.expected_sha256"      "${BUNDLE_HASH}" "${FIG_HASH_PF}"

echo "contracts sync complete — all .expected_sha256 sentinels updated."
