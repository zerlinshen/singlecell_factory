#!/usr/bin/env bash
# CI guard: assert that both bundle_schema.yaml copies are byte-identical.
#
# Canonical:    singlecell_factory/contracts/bundle_schema.yaml
# Vendored copy: multiomics_r_factory/contracts/bundle_schema.yaml
#
# Run on every PR that touches either schema copy. Wire into the same CI lane
# as check_bridge_symlink.sh.
#
# Exit 0 = files are identical.
# Exit 1 = files differ (or one is missing).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CANONICAL="$REPO_ROOT/contracts/bundle_schema.yaml"
VENDORED="$REPO_ROOT/../multiomics_r_factory/contracts/bundle_schema.yaml"

if [ ! -f "$CANONICAL" ]; then
    echo "ERROR: canonical schema not found: $CANONICAL" >&2
    exit 1
fi

if [ ! -f "$VENDORED" ]; then
    echo "ERROR: vendored schema not found: $VENDORED" >&2
    echo "       Run tools/sync-contracts.sh from singlecell_factory/ to sync." >&2
    exit 1
fi

if cmp -s "$CANONICAL" "$VENDORED"; then
    DIGEST=$(sha256sum "$CANONICAL" | awk '{print $1}')
    echo "PASS: bundle_schema.yaml copies are byte-identical (sha256: $DIGEST)"
    exit 0
else
    echo "FAIL: bundle_schema.yaml copies differ." >&2
    echo "      Canonical: $CANONICAL" >&2
    echo "      Vendored:  $VENDORED" >&2
    echo "      Run tools/sync-contracts.sh from singlecell_factory/ to sync." >&2
    diff "$CANONICAL" "$VENDORED" >&2 || true
    exit 1
fi
