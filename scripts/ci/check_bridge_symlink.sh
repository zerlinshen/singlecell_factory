#!/usr/bin/env bash
# CI guard: assert that bridges/local_r_pipeline_macbook/R and R_bundle are symlinks.
# Exits non-zero with a clear error message if either path is a real directory/file.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BRIDGE_DIR="$REPO_ROOT/bridges/local_r_pipeline_macbook"

R_LINK="$BRIDGE_DIR/R"
R_BUNDLE_LINK="$BRIDGE_DIR/R_bundle"

PASS=true

check_symlink() {
    local path="$1"
    local label="$2"
    if [ ! -e "$path" ] && [ ! -L "$path" ]; then
        echo "ERROR: $label does not exist: $path" >&2
        PASS=false
    elif [ ! -L "$path" ]; then
        echo "ERROR: $label exists but is NOT a symlink (found real file/dir): $path" >&2
        PASS=false
    else
        local target
        target="$(readlink "$path")"
        echo "OK: $label -> $target"
    fi
}

check_symlink "$R_LINK"        "bridges/local_r_pipeline_macbook/R"
check_symlink "$R_BUNDLE_LINK" "bridges/local_r_pipeline_macbook/R_bundle"

if [ "$PASS" = "false" ]; then
    echo ""
    echo "FAIL: Bridge symlink check failed." >&2
    echo "      Run Phase 1 setup to create symlinks from bridges/local_r_pipeline_macbook/" >&2
    echo "      to multiomics_r_factory/R and multiomics_r_factory/R_bundle." >&2
    exit 1
fi

echo ""
echo "PASS: Bridge symlink check succeeded."
