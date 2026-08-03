from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# --- r-contract fixtures ---

@pytest.fixture(scope="session")
def rscript_path():
    """Return the absolute path to Rscript, or skip if not found or missing R deps.

    Resolution order:
      1. ``RSCRIPT_BIN`` environment variable (honored verbatim if executable)
      2. ``r_multiomics_arrow`` conda env (production reference; matches
         bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh)
      3. ``r_multiomics`` conda env (parquet-capable fallback as of 2026-05-29)
      4. ``Rscript`` on PATH (last resort)

    Tests are skipped when no candidate is executable, or when the chosen
    Rscript is missing required packages (arrow, jsonlite, Matrix).
    """
    rscript: str | None = None

    env_override = os.environ.get("RSCRIPT_BIN")
    if env_override and os.path.isfile(env_override) and os.access(env_override, os.X_OK):
        rscript = env_override
    else:
        # H3: pin r_multiomics_arrow first so the test suite tracks the
        # renv-pinned production runner (run_remote_bundle_plot.sh) by default.
        candidates = [
            "/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript",
            "/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript",
        ]
        for candidate in candidates:
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                rscript = candidate
                break
        if rscript is None:
            rscript = shutil.which("Rscript")

    if rscript is None:
        pytest.skip("Rscript not found; skipping r_contract tests")

    # Check that required R packages are available.
    probe = subprocess.run(
        [rscript, "-e",
         "pkgs <- c('arrow','jsonlite','Matrix'); "
         "missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; "
         "if (length(missing)) { cat('MISSING:', paste(missing, collapse=','), '\\n'); quit(status=1) } else { cat('OK\\n') }"],
        capture_output=True, text=True, timeout=30,
    )
    if probe.returncode != 0 or "MISSING" in probe.stdout:
        missing_info = probe.stdout.strip() or probe.stderr.strip()
        pytest.skip(f"Required R packages not available ({missing_info}); skipping r_contract tests")
    return rscript
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Stabilize scanpy/numba imports in this CI/runtime environment.
#
# KNOWN BUG this setdefault can expose (found 2026-08-03, scanpy 1.12 /
# numba 0.61.2 / numpy 2.2.6 / python 3.13): with JIT disabled, an
# `@njit`-decorated function falls back to plain CPython execution, and
# scanpy's `_normalize_csr` (scanpy/preprocessing/_normalization.py:65)
# unconditionally `return`s `counts_per_cell, counts_per_cols`, but
# `counts_per_cols` is only ever assigned inside
# `if exclude_highly_expressed:`. Any `sc.pp.normalize_total(...)` call on
# a CSR matrix with the default `exclude_highly_expressed=False`, run
# anywhere JIT is off (including a subprocess that inherits this env var),
# raises `UnboundLocalError: cannot access local variable 'counts_per_cols'`.
# Invisible with JIT on (numba's SSA typing does not hit the same
# CPython-semantics error). If you spawn a subprocess from a pytest test
# that needs real JIT (e.g. a heavy scanpy stress), explicitly set
# `NUMBA_DISABLE_JIT=0` in that child's env rather than inheriting -- see
# tests/test_clustering_memory_regression.py::wave3_200k_stress_result for
# the working pattern and the (separate, JIT-vs-pytest) reason it needs a
# subprocess in the first place. Not yet filed upstream; if scanpy bumps
# past 1.12 and this stops reproducing, this comment can be pruned.
MPL_DIR = ROOT / ".mplconfig"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
