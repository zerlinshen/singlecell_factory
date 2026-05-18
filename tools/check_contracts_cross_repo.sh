#!/usr/bin/env bash
# check_contracts_cross_repo.sh — Cross-repo contracts parity gate.
#
# Exit codes:
#   0  success or sibling-missing warning (non-fatal)
#   1  parity failure (sha mismatch)
#   2  internal error (missing file, sha256sum unavailable)
#
# Canonical home: singlecell_factory/tools/
# Vendored copy:  multiomics_r_factory/tools/
# Both copies MUST remain byte-identical.

set -euo pipefail

THIS_SCRIPT="$(realpath "${BASH_SOURCE[0]}")"
THIS_TOOLS_DIR="$(dirname "$THIS_SCRIPT")"
THIS_REPO="$(dirname "$THIS_TOOLS_DIR")"

# ---------------------------------------------------------------------------
# Locate sibling repo
# ---------------------------------------------------------------------------
if [[ -n "${SIBLING_REPO:-}" ]]; then
    SIBLING="$(realpath "$SIBLING_REPO")"
else
    # Determine which repo we are in by inspecting the basename, then derive
    # the sibling from the absolute path (no symlink-bridge resolution per scope).
    THIS_REPO_NAME="$(basename "$THIS_REPO")"
    PARENT_DIR="$(dirname "$THIS_REPO")"
    if [[ "$THIS_REPO_NAME" == "singlecell_factory" ]]; then
        if [[ -d "$PARENT_DIR/r_multiomics_factory" ]]; then
            SIBLING="$PARENT_DIR/r_multiomics_factory"
        elif [[ -d "$PARENT_DIR/multiomics_r_factory" ]]; then
            SIBLING="$PARENT_DIR/multiomics_r_factory"  # legacy name, accepted during transition
        fi
    elif [[ "$THIS_REPO_NAME" == "r_multiomics_factory" ]] || [[ "$THIS_REPO_NAME" == "multiomics_r_factory" ]]; then
        SIBLING="$PARENT_DIR/singlecell_factory"
    else
        echo "WARNING: check_contracts_cross_repo.sh: cannot determine sibling repo from repo name '$THIS_REPO_NAME'; cross-repo parity check skipped" >&2
        exit 0
    fi
fi

if [[ ! -d "$SIBLING" ]]; then
    echo "WARNING: check_contracts_cross_repo.sh: sibling repo not found at $SIBLING; cross-repo parity check skipped" >&2
    exit 0
fi

# ---------------------------------------------------------------------------
# Verify sha256sum is available
# ---------------------------------------------------------------------------
if ! command -v sha256sum >/dev/null 2>&1; then
    echo "ERROR: check_contracts_cross_repo.sh: sha256sum not found on PATH" >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
THIS_SCHEMA="$THIS_REPO/contracts/bundle_schema.yaml"
THIS_EXPECTED="$THIS_REPO/contracts/.expected_sha256"
SIBLING_SCHEMA="$SIBLING/contracts/bundle_schema.yaml"
SIBLING_EXPECTED="$SIBLING/contracts/.expected_sha256"

# ---------------------------------------------------------------------------
# Verify required files exist
# ---------------------------------------------------------------------------
for f in "$THIS_SCHEMA" "$THIS_EXPECTED" "$SIBLING_SCHEMA" "$SIBLING_EXPECTED"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: check_contracts_cross_repo.sh: required file missing: $f" >&2
        exit 2
    fi
done

# ---------------------------------------------------------------------------
# Compute sha256 hashes
# ---------------------------------------------------------------------------
THIS_HASH="$(sha256sum "$THIS_SCHEMA" | cut -d' ' -f1)"
SIBLING_HASH="$(sha256sum "$SIBLING_SCHEMA" | cut -d' ' -f1)"
THIS_EXPECTED_HASH="$(tr -d '[:space:]' < "$THIS_EXPECTED")"
SIBLING_EXPECTED_HASH="$(tr -d '[:space:]' < "$SIBLING_EXPECTED")"

# ---------------------------------------------------------------------------
# Parity checks
# ---------------------------------------------------------------------------
FAIL=0

if [[ "$THIS_HASH" != "$SIBLING_HASH" ]]; then
    echo "contracts parity FAIL: bundle_schema.yaml hash mismatch between repos" >&2
    echo "  $THIS_REPO:    $THIS_HASH" >&2
    echo "  $SIBLING: $SIBLING_HASH" >&2
    FAIL=1
fi

if [[ "$THIS_HASH" != "$THIS_EXPECTED_HASH" ]]; then
    echo "contracts parity FAIL: $THIS_REPO hash does not match .expected_sha256" >&2
    echo "  computed:  $THIS_HASH" >&2
    echo "  expected:  $THIS_EXPECTED_HASH" >&2
    FAIL=1
fi

if [[ "$SIBLING_HASH" != "$SIBLING_EXPECTED_HASH" ]]; then
    echo "contracts parity FAIL: $SIBLING hash does not match .expected_sha256" >&2
    echo "  computed:  $SIBLING_HASH" >&2
    echo "  expected:  $SIBLING_EXPECTED_HASH" >&2
    FAIL=1
fi

if [[ "$FAIL" -ne 0 ]]; then
    exit 1
fi

echo "contracts parity OK: both repos sha256=$THIS_HASH"
exit 0
