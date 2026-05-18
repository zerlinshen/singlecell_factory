#!/usr/bin/env bash
# cross_factory_contract_gate.sh — Cross-factory contract integrity gate.
#
# Three-factory context (post-2026-05-18):
# Verifies that the three-factory layout (singlecell_factory, r_multiomics_factory,
# plotting_factory) is internally consistent: shared schemas are byte-identical,
# bridge symlinks resolve correctly, and no legacy name leakage exists.
#
# Exit codes:
#   0  all checks pass
#   1  parity / bridge / grep failure
#   2  missing required file or tool
#
# Advisory in Phase 2; promote to blocking in Phase 3 §3.4.

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PARENT_DIR="$(dirname "$REPO_ROOT")"
SC_REPO="$REPO_ROOT"
R_REPO="$PARENT_DIR/r_multiomics_factory"

FAIL=0

# ---------------------------------------------------------------------------
# 1. bundle_schema.yaml byte-identical between SC and r_multiomics_factory
# ---------------------------------------------------------------------------
echo "[cross-factory] Checking bundle_schema.yaml parity..."

SC_SCHEMA="$SC_REPO/contracts/bundle_schema.yaml"
R_SCHEMA="$R_REPO/contracts/bundle_schema.yaml"

if [[ ! -f "$SC_SCHEMA" ]]; then
    echo "ERROR: missing $SC_SCHEMA" >&2
    exit 2
fi
if [[ ! -f "$R_SCHEMA" ]]; then
    echo "ERROR: missing $R_SCHEMA" >&2
    exit 2
fi

if ! command -v sha256sum >/dev/null 2>&1; then
    echo "ERROR: sha256sum not found on PATH" >&2
    exit 2
fi

SC_HASH="$(sha256sum "$SC_SCHEMA" | cut -d' ' -f1)"
R_HASH="$(sha256sum "$R_SCHEMA" | cut -d' ' -f1)"

if [[ "$SC_HASH" != "$R_HASH" ]]; then
    echo "FAIL: bundle_schema.yaml hash mismatch" >&2
    echo "  SC:  $SC_HASH" >&2
    echo "  R:   $R_HASH" >&2
    FAIL=1
else
    echo "OK: bundle_schema.yaml sha256=$SC_HASH"
fi

# ---------------------------------------------------------------------------
# 2. Phase 3 placeholder: figure_bundle_schema.yaml parity
# ---------------------------------------------------------------------------
# Phase 3 placeholder: will verify figure_bundle_schema.yaml byte-identical
# across all three repos (SC, r_multiomics_factory, plotting_factory).

# ---------------------------------------------------------------------------
# 3. Bridge symlinks resolve via readlink -f
# ---------------------------------------------------------------------------
echo "[cross-factory] Checking bridge symlinks..."

check_symlink() {
    local path="$1"
    local label="$2"
    local expected_abs="$3"

    if [[ ! -L "$path" ]]; then
        echo "FAIL: $label is not a symlink: $path" >&2
        FAIL=1
        return
    fi
    local resolved
    resolved="$(readlink -f "$path")"
    if [[ ! -e "$resolved" ]]; then
        echo "FAIL: $label is a dangling symlink -> $resolved" >&2
        FAIL=1
        return
    fi
    if [[ "$resolved" != "$expected_abs"* ]]; then
        echo "FAIL: $label resolves to $resolved (expected prefix $expected_abs)" >&2
        FAIL=1
        return
    fi
    echo "OK: $label -> $resolved"
}

check_symlink \
    "$SC_REPO/bridges/local_r_pipeline_macbook/R" \
    "bridges/local_r_pipeline_macbook/R" \
    "$R_REPO/R"

check_symlink \
    "$SC_REPO/bridges/local_r_pipeline_macbook/R_bundle" \
    "bridges/local_r_pipeline_macbook/R_bundle" \
    "$R_REPO/R_bundle"

check_symlink \
    "$SC_REPO/bridges/local_plot_pipeline/r" \
    "bridges/local_plot_pipeline/r" \
    "$PARENT_DIR/plotting_factory/r"

check_symlink \
    "$SC_REPO/bridges/local_plot_pipeline/python" \
    "bridges/local_plot_pipeline/python" \
    "$PARENT_DIR/plotting_factory/python"

# ---------------------------------------------------------------------------
# 4. No real (non-symlink) file at the named bridge mount points themselves
# ---------------------------------------------------------------------------
echo "[cross-factory] Checking bridge mount points are symlinks (not real files)..."

for bridge_path in \
    "$SC_REPO/bridges/local_r_pipeline_macbook/R" \
    "$SC_REPO/bridges/local_r_pipeline_macbook/R_bundle" \
    "$SC_REPO/bridges/local_plot_pipeline/r" \
    "$SC_REPO/bridges/local_plot_pipeline/python"; do
    if [[ -e "$bridge_path" ]] && [[ ! -L "$bridge_path" ]]; then
        echo "FAIL: bridge mount point is a real file/dir (not a symlink): $bridge_path" >&2
        FAIL=1
    fi
done

if [[ "$FAIL" -eq 0 ]]; then
    echo "OK: all bridge mount points are symlinks"
fi

# ---------------------------------------------------------------------------
# 5. No non-allowlisted multiomics_r_factory references in active source
# ---------------------------------------------------------------------------
echo "[cross-factory] Checking for legacy multiomics_r_factory references..."

LEGACY_HITS="$(grep -rn "multiomics_r_factory" "$SC_REPO/" \
    --include="*.md" --include="*.py" --include="*.sh" \
    --include="*.yaml" --include="*.json" \
    | grep -vE "RENAME_NOTE|ops/governance_records|ops/before_every_run/journal|ops/cleanup_records|previous_name|legacy name accepted during transition" \
    | grep -vE "^$SC_REPO/\.omc/|^$SC_REPO/\.omx/" \
    | grep -vE "cross_factory_contract_gate\.sh" \
    | grep -vE "renamed from|renamed 2026|R factory, renamed" \
    | wc -l)"

if [[ "$LEGACY_HITS" -gt 0 ]]; then
    echo "FAIL: $LEGACY_HITS non-allowlisted multiomics_r_factory reference(s) found:" >&2
    grep -rn "multiomics_r_factory" "$SC_REPO/" \
        --include="*.md" --include="*.py" --include="*.sh" \
        --include="*.yaml" --include="*.json" \
        | grep -vE "RENAME_NOTE|ops/governance_records|ops/before_every_run/journal|ops/cleanup_records|previous_name|legacy name accepted during transition" \
        | grep -vE "^$SC_REPO/\.omc/|^$SC_REPO/\.omx/" \
        | grep -vE "cross_factory_contract_gate\.sh" \
        | grep -vE "renamed from|renamed 2026|R factory, renamed" >&2
    FAIL=1
else
    echo "OK: 0 non-allowlisted multiomics_r_factory references"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
if [[ "$FAIL" -ne 0 ]]; then
    echo "FAIL: cross_factory_contract_gate — one or more checks failed." >&2
    exit 1
fi

echo "PASS: cross_factory_contract_gate — all checks passed."
exit 0
