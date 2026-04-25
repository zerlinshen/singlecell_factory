"""Phase 7D: Python <-> R bundle v2 round-trip subprocess contract tests.

These tests export a synthetic v2 bundle from Python then call io_bundle.R via
subprocess to read it back, asserting that:
  - R sources io_bundle.R without error
  - obs row count matches the exported cell count
  - obsm dimensions are intact
  - marker expression matrix dimensions are intact

All tests are opt-in via @pytest.mark.r_contract and are excluded from the
default CI run (see pyproject.toml addopts).  They are skipped automatically
when Rscript is not available.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ad = pytest.importorskip("anndata")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.export_singlecell_r_bundle import ExportConfig, export_bundle

IO_BUNDLE_R = "/home/zerlinshen/multiomics_r_factory/R_bundle/io_bundle.R"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_tiny_h5ad(path: Path, n_cells: int = 6, n_genes: int = 5) -> None:
    rng = np.random.default_rng(42)
    x = sparse.csr_matrix(
        rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    )
    obs = pd.DataFrame(
        {
            "cell_type": ["T", "Myeloid", "Tumor", "NK", "B", "Epi"][:n_cells],
            "leiden": [str(i % 3) for i in range(n_cells)],
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=["CD3E", "LYZ", "MS4A1", "NKG7", "ELF3"][:n_genes])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.random((n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((n_cells, 3)).astype(np.float32)
    adata.write_h5ad(path)


def _run_r(rscript: str, script: str, timeout: int = 60) -> subprocess.CompletedProcess:
    return subprocess.run(
        [rscript, "-e", script],
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _parse_kv_stdout(stdout: str) -> dict[str, str]:
    """Parse KEY=VALUE lines emitted by cat() in R scripts."""
    result: dict[str, str] = {}
    for line in stdout.splitlines():
        if "=" in line:
            key, _, val = line.partition("=")
            result[key.strip()] = val.strip()
    return result


def _assert_r_ok(proc: subprocess.CompletedProcess, label: str) -> None:
    assert proc.returncode == 0, (
        f"R script [{label}] exited with code {proc.returncode}.\n"
        f"--- stdout ---\n{proc.stdout}\n"
        f"--- stderr ---\n{proc.stderr}"
    )
    # Any R-side stop() message ends up in stderr; treat non-empty stderr as a
    # warning but allow informational messages (arrow/jsonlite load messages).
    # Hard-fail on "Error" in stderr.
    if "Error" in proc.stderr:
        pytest.fail(
            f"R script [{label}] has Error in stderr:\n{proc.stderr}\nstdout:\n{proc.stdout}"
        )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.r_contract
def test_r_sources_io_bundle_cleanly(rscript_path):
    """io_bundle.R sources without error and prints 'ok'."""
    proc = _run_r(
        rscript_path,
        f"source('{IO_BUNDLE_R}'); cat('ok\\n')",
    )
    _assert_r_ok(proc, "source_io_bundle")
    assert "ok" in proc.stdout, f"Expected 'ok' in stdout, got: {proc.stdout!r}"


@pytest.mark.r_contract
def test_r_reads_v2_obs_nrow(rscript_path, tmp_path):
    """R reads v2 bundle obs.parquet and reports correct cell count."""
    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_v2"
    n_cells = 6
    _make_tiny_h5ad(h5ad, n_cells=n_cells)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ", "MS4A1"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2",
            format="parquet",
        )
    )

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"cat(sprintf('OBS_NROW=%d\\n', nrow(b$obs)))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "read_v2_obs_nrow")

    kv = _parse_kv_stdout(proc.stdout)
    assert "OBS_NROW" in kv, f"OBS_NROW not found in stdout: {proc.stdout!r}"
    assert int(kv["OBS_NROW"]) == n_cells, (
        f"Expected {n_cells} obs rows, R reported {kv['OBS_NROW']}"
    )


@pytest.mark.r_contract
def test_r_reads_v2_obsm_dims(rscript_path, tmp_path):
    """R reads v2 bundle obsm parquets and reports correct matrix dimensions."""
    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_v2"
    n_cells = 6
    _make_tiny_h5ad(h5ad, n_cells=n_cells)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2",
            format="parquet",
        )
    )

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"cat(sprintf('UMAP_NROW=%d\\n', nrow(b$obsm[['X_umap']]))); "
        f"cat(sprintf('UMAP_NCOL=%d\\n', ncol(b$obsm[['X_umap']]))); "
        f"cat(sprintf('PCA_NROW=%d\\n', nrow(b$obsm[['X_pca']]))); "
        f"cat(sprintf('PCA_NCOL=%d\\n', ncol(b$obsm[['X_pca']])))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "read_v2_obsm_dims")

    kv = _parse_kv_stdout(proc.stdout)
    assert int(kv["UMAP_NROW"]) == n_cells
    assert int(kv["UMAP_NCOL"]) == 2
    assert int(kv["PCA_NROW"]) == n_cells
    assert int(kv["PCA_NCOL"]) == 3


@pytest.mark.r_contract
def test_r_reads_v2_marker_expr_dims(rscript_path, tmp_path):
    """R reads v2 bundle marker_expr.parquet and reports correct expression dims."""
    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_v2"
    n_cells = 6
    markers = ("CD3E", "LYZ", "MS4A1")
    _make_tiny_h5ad(h5ad, n_cells=n_cells, n_genes=5)

    manifest = export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=markers,
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2",
            format="parquet",
        )
    )
    n_markers_present = len(manifest["bundle"]["markers_present"])

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"cat(sprintf('EXPR_NROW=%d\\n', nrow(b$expr_sparse))); "
        f"cat(sprintf('EXPR_NCOL=%d\\n', ncol(b$expr_sparse)))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "read_v2_marker_expr_dims")

    kv = _parse_kv_stdout(proc.stdout)
    assert int(kv["EXPR_NROW"]) == n_cells, (
        f"Expected {n_cells} expr rows, got {kv['EXPR_NROW']}"
    )
    assert int(kv["EXPR_NCOL"]) == n_markers_present, (
        f"Expected {n_markers_present} marker cols, got {kv['EXPR_NCOL']}"
    )
