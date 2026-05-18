#!/usr/bin/env bash
# sync-contracts.sh — Propagate the canonical bundle_schema.yaml from
# singlecell_factory/contracts/ to multiomics_r_factory/contracts/ and
# update both .expected_sha256 sentinel files.
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
#   - missing sibling multiomics_r_factory repo
#   - sha256sum mismatch after copy (should never happen; guards against fs issues)
#
# CI parity check (no sibling required):
#   Each repo's CI job can run the local parity check independently:
#     actual=$(sha256sum contracts/bundle_schema.yaml | awk '{print $1}')
#     expected=$(cat contracts/.expected_sha256)
#     [ "$actual" = "$expected" ] || { echo "contracts parity FAIL"; exit 1; }

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CANONICAL="${REPO_ROOT}/contracts/bundle_schema.yaml"
CANONICAL_SHA_FILE="${REPO_ROOT}/contracts/.expected_sha256"

# Locate the sibling multiomics_r_factory repo.
SIBLING_ROOT="$(cd "${REPO_ROOT}/../multiomics_r_factory" 2>/dev/null && pwd)" || true
if [[ -z "${SIBLING_ROOT}" || ! -d "${SIBLING_ROOT}" ]]; then
  echo "ERROR: sibling repo not found at ${REPO_ROOT}/../multiomics_r_factory" >&2
  echo "       Clone multiomics_r_factory alongside singlecell_factory and retry." >&2
  exit 1
fi

VENDORED="${SIBLING_ROOT}/contracts/bundle_schema.yaml"
VENDORED_SHA_FILE="${SIBLING_ROOT}/contracts/.expected_sha256"

if [[ ! -f "${CANONICAL}" ]]; then
  echo "ERROR: canonical schema not found: ${CANONICAL}" >&2
  exit 1
fi

# Ensure destination contracts/ directory exists.
mkdir -p "$(dirname "${VENDORED}")"

# Copy canonical -> vendored (byte-for-byte).
cp "${CANONICAL}" "${VENDORED}"

# Compute sha256 and verify the copy is identical.
CANONICAL_HASH=$(sha256sum "${CANONICAL}" | awk '{print $1}')
VENDORED_HASH=$(sha256sum "${VENDORED}"  | awk '{print $1}')

if [[ "${CANONICAL_HASH}" != "${VENDORED_HASH}" ]]; then
  echo "ERROR: sha256 mismatch after copy — filesystem issue?" >&2
  echo "  canonical: ${CANONICAL_HASH}" >&2
  echo "  vendored:  ${VENDORED_HASH}" >&2
  exit 1
fi

# Update both .expected_sha256 sentinels.
printf '%s\n' "${CANONICAL_HASH}" > "${CANONICAL_SHA_FILE}"
printf '%s\n' "${VENDORED_HASH}"  > "${VENDORED_SHA_FILE}"

echo "contracts sync OK: ${CANONICAL_HASH}"
echo "  canonical: ${CANONICAL}"
echo "  vendored:  ${VENDORED}"
