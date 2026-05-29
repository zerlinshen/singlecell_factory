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
MPL_DIR = ROOT / ".mplconfig"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
