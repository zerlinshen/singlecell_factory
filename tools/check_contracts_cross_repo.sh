#!/usr/bin/env bash
# check_contracts_cross_repo.sh — Cross-repo contracts parity gate.
#
# Exit codes:
#   0  success or sibling-missing warning (non-fatal)
#   1  parity failure (sha mismatch)
#   2  internal error (missing file, sha256sum unavailable)
#
# Canonical home: singlecell_factory/tools/
# Vendored copy:  r_multiomics_factory/tools/
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
        else
            echo "WARNING: check_contracts_cross_repo.sh: expected sibling r_multiomics_factory not found under $PARENT_DIR" >&2
            exit 0
        fi
    elif [[ "$THIS_REPO_NAME" == "r_multiomics_factory" ]]; then
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
# Helpers
# ---------------------------------------------------------------------------
expected_hash_for_key() {
    local expected_file="$1"
    local key="$2"
    local value

    value="$(awk -v key="$key" '$1 == key { print $2 }' "$expected_file")"
    if [[ -n "$value" ]]; then
        printf '%s\n' "$value"
        return 0
    fi

    # Legacy fallback: the Phase 2 sentinel was a single bare sha for
    # bundle_schema only. Figure schema governance requires named sentinels.
    if [[ "$key" == "bundle_schema" ]]; then
        value="$(tr -d '[:space:]' < "$expected_file")"
        if [[ "$value" =~ ^[0-9a-fA-F]{64}$ ]]; then
            printf '%s\n' "$value"
            return 0
        fi
    fi

    echo "ERROR: check_contracts_cross_repo.sh: missing '$key' entry in $expected_file" >&2
    return 1
}

check_schema_key() {
    local key="$1"
    local filename="$2"
    local this_schema="$THIS_REPO/contracts/$filename"
    local sibling_schema="$SIBLING/contracts/$filename"
    local this_expected="$THIS_REPO/contracts/.expected_sha256"
    local sibling_expected="$SIBLING/contracts/.expected_sha256"

    for f in "$this_schema" "$sibling_schema" "$this_expected" "$sibling_expected"; do
        if [[ ! -f "$f" ]]; then
            echo "ERROR: check_contracts_cross_repo.sh: required file missing for $key: $f" >&2
            exit 2
        fi
    done

    local this_hash sibling_hash this_expected_hash sibling_expected_hash
    this_hash="$(sha256sum "$this_schema" | cut -d' ' -f1)"
    sibling_hash="$(sha256sum "$sibling_schema" | cut -d' ' -f1)"
    if ! this_expected_hash="$(expected_hash_for_key "$this_expected" "$key")"; then
        exit 2
    fi
    if ! sibling_expected_hash="$(expected_hash_for_key "$sibling_expected" "$key")"; then
        exit 2
    fi

    if [[ "$this_hash" != "$sibling_hash" ]]; then
        echo "contracts parity FAIL: $key ($filename) hash mismatch between repos" >&2
        echo "  $THIS_REPO: $this_hash" >&2
        echo "  $SIBLING: $sibling_hash" >&2
        FAIL=1
    fi

    if [[ "$this_hash" != "$this_expected_hash" ]]; then
        echo "contracts parity FAIL: $THIS_REPO $key hash does not match .expected_sha256" >&2
        echo "  computed: $this_hash" >&2
        echo "  expected: $this_expected_hash" >&2
        FAIL=1
    fi

    if [[ "$sibling_hash" != "$sibling_expected_hash" ]]; then
        echo "contracts parity FAIL: $SIBLING $key hash does not match .expected_sha256" >&2
        echo "  computed: $sibling_hash" >&2
        echo "  expected: $sibling_expected_hash" >&2
        FAIL=1
    fi

    if [[ "$this_hash" == "$sibling_hash" && "$this_hash" == "$this_expected_hash" && "$sibling_hash" == "$sibling_expected_hash" ]]; then
        echo "contracts parity OK: $key sha256=$this_hash"
    fi
}

# ---------------------------------------------------------------------------
# Parity checks
# ---------------------------------------------------------------------------
FAIL=0

check_schema_key "bundle_schema" "bundle_schema.yaml"
check_schema_key "figure_bundle_schema" "figure_bundle_schema.yaml"

if [[ "$FAIL" -ne 0 ]]; then
    exit 1
fi

exit 0
