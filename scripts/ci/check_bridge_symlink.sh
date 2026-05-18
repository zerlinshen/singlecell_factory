#!/usr/bin/env bash
# CI guard: assert that all bridge symlinks are valid, non-dangling symlinks
# resolving into expected sibling repos.
#
# Exits non-zero with a clear error message on any failure.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BRIDGE_BASE="$REPO_ROOT/bridges"
PARENT_DIR="$(dirname "$REPO_ROOT")"

PASS=true

check_bridge_symlink() {
    local path="$1"
    local label="$2"
    local expected_prefix="$3"

    # Must exist as a path entry
    if [ ! -e "$path" ] && [ ! -L "$path" ]; then
        echo "ERROR: $label does not exist: $path" >&2
        PASS=false
        return
    fi

    # Must be a symlink, not a real file or directory
    if [ ! -L "$path" ]; then
        echo "ERROR: $label exists but is NOT a symlink (found real file/dir): $path" >&2
        PASS=false
        return
    fi

    local raw_target
    raw_target="$(readlink "$path")"

    # Must not be dangling (target must resolve to an existing path)
    if [ ! -e "$path" ]; then
        echo "ERROR: $label is a dangling symlink (target does not exist): $raw_target" >&2
        PASS=false
        return
    fi

    # Resolve to absolute path
    local resolved
    resolved="$(readlink -f "$path")"

    # Must not point back into the same repo (self-bridge)
    if [[ "$resolved" == "$REPO_ROOT"* ]]; then
        echo "ERROR: $label is a self-bridge (resolves inside singlecell_factory): $resolved" >&2
        PASS=false
        return
    fi

    # Must resolve to a directory
    if [ ! -d "$resolved" ]; then
        echo "ERROR: $label resolves to a non-directory: $resolved" >&2
        PASS=false
        return
    fi

    # Must resolve inside the expected sibling repo prefix
    local expected_abs="$PARENT_DIR/$expected_prefix"
    if [[ "$resolved" != "$expected_abs"* ]]; then
        echo "ERROR: $label resolves to $resolved but expected prefix $expected_abs" >&2
        PASS=false
        return
    fi

    echo "OK: $label -> $raw_target (resolves: $resolved)"
}

# local_r_pipeline_macbook bridge
check_bridge_symlink \
    "$BRIDGE_BASE/local_r_pipeline_macbook/R" \
    "bridges/local_r_pipeline_macbook/R" \
    "r_multiomics_factory/R"

check_bridge_symlink \
    "$BRIDGE_BASE/local_r_pipeline_macbook/R_bundle" \
    "bridges/local_r_pipeline_macbook/R_bundle" \
    "r_multiomics_factory/R_bundle"

# local_plot_pipeline bridge
check_bridge_symlink \
    "$BRIDGE_BASE/local_plot_pipeline/r" \
    "bridges/local_plot_pipeline/r" \
    "plotting_factory/r"

check_bridge_symlink \
    "$BRIDGE_BASE/local_plot_pipeline/python" \
    "bridges/local_plot_pipeline/python" \
    "plotting_factory/python"

if [ "$PASS" = "false" ]; then
    echo ""
    echo "FAIL: Bridge symlink check failed." >&2
    exit 1
fi

echo ""
echo "PASS: Bridge symlink check succeeded."
