"""Cross-language bundle SHA parity test (MISS-1).

Asserts that Python's _bundle_sha256_concat() and R's compute_bundle_sha()
produce identical hex SHA256 strings for the same bundle directory.

Parity failure means the on-disk cache key and the provenance.json manifest
field reference different bundle identities — a precision bug.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.export_singlecell_r_bundle import ExportConfig, export_bundle, _bundle_sha256_concat

RSCRIPT_BIN = "/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript"
BUNDLE_CACHE_R = ROOT.parent / "r_multiomics_factory" / "R_bundle" / "bundle_cache.R"


# ---------------------------------------------------------------------------
# Synthetic bundle fixture
# ---------------------------------------------------------------------------

def _make_synthetic_adata(n_cells: int = 80, n_genes: int = 20, seed: int = 7) -> AnnData:
    rng = np.random.default_rng(seed)
    markers = ["CD3E", "LYZ", "MS4A1", "NKG7"]
    other = [f"GENE_{i}" for i in range(n_genes - len(markers))]
    var_names = markers + other

    X = sp.csr_matrix(rng.uniform(0, 5, size=(n_cells, n_genes)).astype(np.float32))
    adata = AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = var_names
    adata.obs["leiden"] = pd.Categorical(["0"] * (n_cells // 2) + ["1"] * (n_cells - n_cells // 2))
    adata.obs["cell_type"] = pd.Categorical(["T"] * n_cells)
    adata.obsm["X_umap"] = rng.normal(size=(n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"]  = rng.normal(size=(n_cells, 10)).astype(np.float32)
    return adata


@pytest.fixture(scope="module")
def synthetic_bundle(tmp_path_factory):
    """Export a small synthetic v2 bundle and return its path."""
    tmp = tmp_path_factory.mktemp("sha_parity_bundle")
    h5ad = tmp / "test.h5ad"
    adata = _make_synthetic_adata()
    adata.write_h5ad(str(h5ad))

    bundle_dir = tmp / "bundle"
    cfg = ExportConfig(
        input_h5ad=h5ad,
        output_dir=bundle_dir,
        schema_version="v2",
        format="parquet",
        markers=("CD3E", "LYZ", "MS4A1", "NKG7"),
        compression="gzip",
    )
    export_bundle(cfg)
    return bundle_dir


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_bundle_sha_python_r_parity(synthetic_bundle):
    """Python and R must produce identical sha256 for the same bundle."""
    import os
    import shutil

    rscript = RSCRIPT_BIN
    if not Path(rscript).is_file():
        pytest.skip(f"Rscript not found at {rscript}")

    bundle_cache_r = BUNDLE_CACHE_R
    if not bundle_cache_r.is_file():
        pytest.skip(f"bundle_cache.R not found at {bundle_cache_r}")

    # Python SHA
    python_sha = _bundle_sha256_concat(synthetic_bundle)
    assert len(python_sha) == 64, f"Python SHA unexpected length: {python_sha!r}"

    # R SHA via subprocess
    r_script = textwrap.dedent(f"""\
        source({repr(str(bundle_cache_r))})
        sha <- compute_bundle_sha({repr(str(synthetic_bundle))})
        cat(sha, "\\n")
    """)

    result = subprocess.run(
        [rscript, "--vanilla", "-e", r_script],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if result.returncode != 0:
        pytest.fail(
            f"R compute_bundle_sha failed (rc={result.returncode}):\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )

    r_sha = result.stdout.strip()
    assert len(r_sha) == 64, f"R SHA unexpected length: {r_sha!r}\nstderr: {result.stderr}"

    assert python_sha == r_sha, (
        f"SHA parity FAIL: Python={python_sha!r}  R={r_sha!r}\n"
        f"Bundle: {synthetic_bundle}\n"
        f"R stderr: {result.stderr}"
    )


def test_bundle_sha_matches_provenance_json(synthetic_bundle):
    """SHA in provenance.json must match Python's _bundle_sha256_concat."""
    provenance_path = synthetic_bundle / "provenance.json"
    if not provenance_path.exists():
        pytest.skip("provenance.json not written by this export (optional field)")

    with provenance_path.open() as fh:
        prov = json.load(fh)

    recorded_sha = prov.get("bundle_sha256", "")
    if not recorded_sha:
        pytest.skip("bundle_sha256 not present in provenance.json")

    recomputed = _bundle_sha256_concat(synthetic_bundle)
    assert recorded_sha == recomputed, (
        f"provenance.json bundle_sha256 mismatch:\n"
        f"  recorded:   {recorded_sha}\n"
        f"  recomputed: {recomputed}"
    )
