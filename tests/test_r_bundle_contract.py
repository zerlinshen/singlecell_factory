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
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ad = pytest.importorskip("anndata")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
R_FACTORY_ROOT = ROOT.parent / "r_multiomics_factory"

from scripts.export_singlecell_r_bundle import ExportConfig, export_bundle

IO_BUNDLE_R = str(R_FACTORY_ROOT / "R_bundle" / "io_bundle.R")
ATAC_MODULE_R = str(R_FACTORY_ROOT / "R" / "atac_module.R")
HIC_MODULE_R = str(R_FACTORY_ROOT / "R" / "hic_module.R")

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


def _make_tiny_atac_h5ad(
    path: Path,
    n_cells: int = 6,
    n_genes: int = 5,
    n_peaks: int = 4,
    n_lsi: int = 2,
) -> None:
    rng = np.random.default_rng(44)
    x = sparse.csr_matrix(
        rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    )
    obs = pd.DataFrame(
        {
            "cell_type": [["T", "B", "Myeloid"][i % 3] for i in range(n_cells)],
            "leiden": [str(i % 3) for i in range(n_cells)],
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=["CD3E", "LYZ", "MS4A1", "NKG7", "ELF3"][:n_genes])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.random((n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((n_cells, 3)).astype(np.float32)
    adata.obsm["X_lsi"] = rng.normal(size=(n_cells, n_lsi)).astype(np.float32)
    adata.obsm["atac_peaks"] = sparse.csr_matrix(
        rng.integers(0, 4, size=(n_cells, n_peaks)).astype(np.float32)
    )
    starts = np.arange(1000, 1000 + n_peaks * 200, 200)
    peaks = pd.DataFrame(
        {
            "peak_id": [f"chr1:{start}-{start + 100}" for start in starts],
            "chrom": ["chr1"] * n_peaks,
            "start": starts.astype("int64"),
            "end": (starts + 100).astype("int64"),
        }
    )
    adata.uns["atac_peaks"] = peaks
    adata.write_h5ad(path)


def _make_tiny_hic_h5ad(
    path: Path,
    n_cells: int = 6,
    n_genes: int = 5,
    n_bins: int = 6,
) -> None:
    rng = np.random.default_rng(45)
    x = sparse.csr_matrix(
        rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    )
    obs = pd.DataFrame(
        {
            "cell_type": [["T", "B", "Myeloid"][i % 3] for i in range(n_cells)],
            "leiden": [str(i % 3) for i in range(n_cells)],
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=["CD3E", "LYZ", "MS4A1", "NKG7", "ELF3"][:n_genes])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.random((n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((n_cells, 3)).astype(np.float32)

    starts = np.arange(0, n_bins * 1000, 1000)
    adata.uns["hic_bins"] = pd.DataFrame(
        {
            "bin_id": np.arange(n_bins, dtype=np.int64),
            "chrom": ["chr1"] * n_bins,
            "start": starts.astype(np.int64),
            "end": (starts + 1000).astype(np.int64),
        }
    )
    mat = np.zeros((n_bins, n_bins), dtype=np.float32)
    for i in range(n_bins):
        mat[i, i] = 10 - i
        if i + 1 < n_bins:
            mat[i, i + 1] = 2
            mat[i + 1, i] = 2
    adata.uns["hic_contact_matrix"] = sparse.csr_matrix(mat)
    adata.uns["hic_tad_boundaries"] = pd.DataFrame(
        {
            "bin_id": np.arange(n_bins, dtype=np.int64),
            "chrom": ["chr1"] * n_bins,
            "start": starts.astype(np.int64),
            "end": (starts + 1000).astype(np.int64),
            "insulation": np.array([0.2, 0.1, -1.2, 0.0, -0.8, 0.1], dtype=np.float32)[:n_bins],
            "is_boundary": [False, False, True, False, True, False][:n_bins],
        }
    )
    adata.uns["hic_compartments"] = pd.DataFrame(
        {
            "bin_id": np.arange(n_bins, dtype=np.int64),
            "chrom": ["chr1"] * n_bins,
            "start": starts.astype(np.int64),
            "end": (starts + 1000).astype(np.int64),
            "eigenvector_1": np.array([0.5, 0.4, 0.2, -0.1, -0.2, -0.3], dtype=np.float32)[:n_bins],
            "compartment": ["A", "A", "A", "B", "B", "B"][:n_bins],
        }
    )
    adata.uns["hic_ingest_metadata"] = {
        "resolution_bp": 1000,
        "matrix_format": "csr_sparse",
        "input_format": "tsv_contact_pairs",
        "normalization": "raw_counts_or_unbalanced",
    }
    adata.uns["hic_tad_metadata"] = {
        "compartment_status": "confident",
        "low_information_chromosomes": [],
        "compartment_status_by_chrom": {"chr1": "confident"},
    }
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


# ---------------------------------------------------------------------------
# C2 / H2 / M1 / M2 companion tests
# ---------------------------------------------------------------------------

def test_export_is_atomic_no_partial_manifest_visible(tmp_path, monkeypatch):
    """C2: while a bundle is being exported, the final output_dir must NEVER
    contain bundle_manifest.json until the temp-dir rename completes.

    We instrument _export_bundle_v2 with a sleep+sentinel: a poller thread
    samples final_output_dir during the export and asserts it is either
    absent or fully populated (manifest present implies all sibling files).
    """
    import scripts.export_singlecell_r_bundle as exp_mod

    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_v2"
    n_cells = 6
    _make_tiny_h5ad(h5ad, n_cells=n_cells)

    real_export_v2 = exp_mod._export_bundle_v2

    def slow_export_v2(adata, cell_idx, config):
        result = real_export_v2(adata, cell_idx, config)
        # Hold the temp dir populated for a beat so the poller can observe.
        time.sleep(0.4)
        return result

    monkeypatch.setattr(exp_mod, "_export_bundle_v2", slow_export_v2)

    observations: list[tuple[bool, bool]] = []
    stop_flag = threading.Event()

    def poll():
        while not stop_flag.is_set():
            final_present = bundle_dir.exists()
            manifest_present = (bundle_dir / "bundle_manifest.json").exists()
            observations.append((final_present, manifest_present))
            time.sleep(0.02)

    poller = threading.Thread(target=poll, daemon=True)
    poller.start()
    try:
        exp_mod.export_bundle(
            exp_mod.ExportConfig(
                input_h5ad=h5ad,
                output_dir=bundle_dir,
                obs_columns=("cell_type", "leiden"),
                markers=("CD3E", "LYZ", "MS4A1"),
                obsm_keys=("X_umap", "X_pca"),
                schema_version="v2",
                format="parquet",
            )
        )
    finally:
        stop_flag.set()
        poller.join(timeout=2.0)

    # Final dir must exist and be fully populated post-export.
    assert (bundle_dir / "bundle_manifest.json").exists()

    # No observation may show final_dir-present-but-manifest-missing while
    # we were sampling; that would indicate a non-atomic publish.
    bad = [(fp, mp) for (fp, mp) in observations if fp and not mp]
    assert not bad, (
        f"Non-atomic publish detected: {len(bad)} samples saw the final "
        f"output_dir before the manifest landed. Total samples={len(observations)}"
    )


@pytest.mark.r_contract
def test_claim_guard_surfaces_in_r(rscript_path, tmp_path):
    """H2: read_bundle_v2 must surface claim_guard= via message() so callers
    that wrap the reader in a logging shell see the use-restriction string.
    """
    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_v2"
    _make_tiny_h5ad(h5ad, n_cells=6)

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
        f"cat('done\\n')"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "claim_guard_surface")
    combined = proc.stdout + proc.stderr
    assert "claim_guard=" in combined, (
        f"Expected 'claim_guard=' in R output, got:\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
    )


def test_categorical_round_trip_preserves_ordered_levels(tmp_path):
    """M2: an ordered Categorical with explicit levels must survive the parquet
    round-trip with both ordered=True and the original level order intact."""
    import anndata as ad_local

    h5ad = tmp_path / "ord.h5ad"
    bundle_dir = tmp_path / "bundle_v2_ord"
    n_cells = 4
    rng = np.random.default_rng(0)
    x = sparse.csr_matrix(rng.integers(0, 5, size=(n_cells, 3)).astype(np.float32))
    stage_levels = ["I", "II", "III", "IV"]
    obs = pd.DataFrame(
        {
            "stage": pd.Categorical(
                ["II", "I", "IV", "III"], categories=stage_levels, ordered=True
            ),
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=["CD3E", "LYZ", "MS4A1"])
    adata = ad_local.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.random((n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((n_cells, 2)).astype(np.float32)
    adata.write_h5ad(h5ad)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("stage",),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2",
            format="parquet",
        )
    )

    df = pd.read_parquet(bundle_dir / "obs.parquet")
    assert "stage" in df.columns
    stage_back = df["stage"]
    assert isinstance(stage_back.dtype, pd.CategoricalDtype), (
        f"Expected categorical dtype after round-trip, got {stage_back.dtype}"
    )
    assert stage_back.cat.ordered is True, "ordered flag was lost in round-trip"
    assert list(stage_back.cat.categories) == stage_levels, (
        f"Expected levels {stage_levels}, got {list(stage_back.cat.categories)}"
    )


def test_cell_column_present_in_all_parquets(tmp_path):
    """M1: every parquet written by the v2 exporter must expose a `cell`
    column literally (not just an unnamed index, and not __index_level_0__)."""
    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_v2_cellcol"
    _make_tiny_h5ad(h5ad, n_cells=6, n_genes=5)

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

    targets = [
        bundle_dir / "obs.parquet",
        bundle_dir / "obsm" / "X_umap.parquet",
        bundle_dir / "obsm" / "X_pca.parquet",
        bundle_dir / "marker_expr.parquet",
    ]
    for path in targets:
        assert path.exists(), f"Expected parquet at {path}"
        df = pd.read_parquet(path)
        assert "cell" in df.columns, (
            f"{path.name} is missing literal 'cell' column; got {list(df.columns)}"
        )
        assert "__index_level_0__" not in df.columns, (
            f"{path.name} still has legacy __index_level_0__ column"
        )


# ---------------------------------------------------------------------------
# Phase B: bundle v2.1 schema tests (additive, backwards-compatible)
# ---------------------------------------------------------------------------

def test_v2_1_schema_default(tmp_path):
    """Default writer emits singlecell_r_bundle_v2.1 with empty extensions={}."""
    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_default"
    _make_tiny_h5ad(h5ad, n_cells=6)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            format="parquet",
        )
    )

    manifest_path = bundle_dir / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["schema_version"] == "singlecell_r_bundle_v2.1", (
        f"Default writer should emit v2.1, got: {manifest['schema_version']!r}"
    )
    assert manifest.get("extensions") == {}, (
        f"v2.1 default should have empty extensions={{}}, got: {manifest.get('extensions')!r}"
    )


def test_v2_schema_explicit_legacy(tmp_path):
    """--schema-version v2 emits the legacy schema string and OMITS extensions field entirely."""
    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_v2_legacy"
    _make_tiny_h5ad(h5ad, n_cells=6)

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

    manifest_path = bundle_dir / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["schema_version"] == "singlecell_r_bundle_v2", (
        f"v2 emitter should produce singlecell_r_bundle_v2, got: {manifest['schema_version']!r}"
    )
    assert "extensions" not in manifest, (
        "v2 manifest must NOT carry an extensions field (omitted, not silently empty)"
    )


@pytest.mark.r_contract
def test_v2_1_reader_accepts_v2_bundle(rscript_path, tmp_path):
    """R reader (post-v2.1) still reads legacy v2 bundles cleanly."""
    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_v2_back_compat"
    n_cells = 6
    _make_tiny_h5ad(h5ad, n_cells=n_cells)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type",),
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
        f"cat(sprintf('SCHEMA=%s\\n', b$manifest$schema_version)); "
        f"cat(sprintf('OBS_NROW=%d\\n', nrow(b$obs))); "
        f"cat(sprintf('HAS_EXT=%s\\n', as.character(!is.null(b$extensions))))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "v2.1_reader_accepts_v2")
    kv = _parse_kv_stdout(proc.stdout)
    assert kv["SCHEMA"] == "singlecell_r_bundle_v2"
    assert int(kv["OBS_NROW"]) == n_cells
    # extensions key always exists in the result list, but is empty for v2 input
    assert kv["HAS_EXT"] == "TRUE"


@pytest.mark.r_contract
def test_v2_1_unknown_extension_skipped(rscript_path, tmp_path):
    """R reader emits 'skipping unknown extension' message and does not error."""
    h5ad = tmp_path / "tiny.h5ad"
    bundle_dir = tmp_path / "bundle_unknown_ext"
    _make_tiny_h5ad(h5ad, n_cells=6)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type",),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            format="parquet",
        )
    )

    # Hand-craft a manifest with a fake extension key that the R reader doesn't know.
    manifest_path = bundle_dir / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["schema_version"] == "singlecell_r_bundle_v2.1"
    manifest["extensions"]["future_modality_xyz"] = {
        "version": "1.0",
        "files": [],
        "claim_guard": "not_for_de_or_new_quantitative_claims_without_full_object_validation",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"cat('READ_OK\\n')"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "v2.1_unknown_extension_skipped")
    assert "READ_OK" in proc.stdout
    assert "skipping unknown extension 'future_modality_xyz'" in proc.stderr, (
        f"Expected skip message in stderr, got:\n{proc.stderr}"
    )


def test_add_extension_helper():
    """Python add_extension API: idempotent, validates inputs, populates fields."""
    from scripts.export_singlecell_r_bundle import add_extension, CLAIM_GUARD

    manifest: dict = {}

    # First registration creates manifest['extensions'] and stores the entry.
    entry = add_extension(
        manifest,
        "protein",
        version="1.0",
        files=["protein_expr"],
        normalization="CLR",
    )
    assert manifest["extensions"]["protein"]["version"] == "1.0"
    assert manifest["extensions"]["protein"]["files"] == ["protein_expr"]
    assert manifest["extensions"]["protein"]["claim_guard"] == CLAIM_GUARD
    assert manifest["extensions"]["protein"]["normalization"] == "CLR"
    assert entry is manifest["extensions"]["protein"]

    # Idempotent for same name: a second call replaces the prior entry.
    add_extension(
        manifest,
        "protein",
        version="1.1",
        files=["protein_expr", "isotype_ctrl"],
        normalization="DSB",
    )
    assert manifest["extensions"]["protein"]["version"] == "1.1"
    assert manifest["extensions"]["protein"]["files"] == ["protein_expr", "isotype_ctrl"]
    assert manifest["extensions"]["protein"]["normalization"] == "DSB"
    # Only one extension key remains after idempotent replace.
    assert list(manifest["extensions"].keys()) == ["protein"]

    # Input validation: rejects non-dict manifest, empty name, empty version, non-list files,
    # empty claim_guard.
    with pytest.raises(TypeError, match="manifest to be a dict"):
        add_extension([], "protein", version="1.0", files=["x"])
    with pytest.raises(TypeError, match="non-empty string `name`"):
        add_extension({}, "", version="1.0", files=["x"])
    with pytest.raises(TypeError, match="non-empty string `version`"):
        add_extension({}, "protein", version="", files=["x"])
    with pytest.raises(TypeError, match="`files` to be a list/tuple"):
        add_extension({}, "protein", version="1.0", files="x")  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="non-empty string `claim_guard`"):
        add_extension({}, "protein", version="1.0", files=["x"], claim_guard="")


# ---------------------------------------------------------------------------
# Phase B (v2.1) protein / ADT extension tests
# ---------------------------------------------------------------------------

def _make_tiny_protein_h5ad(
    path: Path,
    *,
    n_cells: int = 8,
    n_genes: int = 5,
    n_proteins: int = 10,
    isotypes: tuple[str, ...] = ("ISO_IgG1", "ISO_IgG2"),
) -> tuple[list[str], np.ndarray]:
    """Synthesize a tiny AnnData with both an RNA matrix and an ADT modality.

    Returns (protein_names, clr_matrix) so tests can assert round-trip values.
    """
    rng = np.random.default_rng(7)
    x = sparse.csr_matrix(
        rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    )
    obs = pd.DataFrame(
        {
            "cell_type": [["T", "B", "Myeloid"][i % 3] for i in range(n_cells)],
            "leiden": [str(i % 3) for i in range(n_cells)],
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=["CD3E", "LYZ", "MS4A1", "NKG7", "ELF3"][:n_genes])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.random((n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((n_cells, 3)).astype(np.float32)

    # Build a protein panel with two isotype controls at the front.
    base = [f"ADT_{i + 1}" for i in range(n_proteins - len(isotypes))]
    protein_names = list(isotypes) + base
    clr = rng.random((n_cells, n_proteins)).astype(np.float32)
    adata.obsm["protein_clr"] = clr
    adata.uns["protein_names"] = protein_names

    adata.write_h5ad(path)
    return protein_names, clr


@pytest.mark.r_contract
def test_protein_extension_round_trip(rscript_path, tmp_path):
    """Python writes protein.parquet; R reads it back with matching shape and rownames."""
    h5ad = tmp_path / "tiny_protein.h5ad"
    bundle_dir = tmp_path / "bundle_v2_protein"
    n_cells = 8
    n_proteins = 10
    isotypes = ("ISO_IgG1", "ISO_IgG2")
    protein_names, _ = _make_tiny_protein_h5ad(
        h5ad, n_cells=n_cells, n_proteins=n_proteins, isotypes=isotypes,
    )

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2.1",
            format="parquet",
            include_protein=True,
            protein_obsm_key="protein_clr",
            protein_isotype_controls=isotypes,
        )
    )

    # Manifest assertions on the Python side.
    manifest = json.loads((bundle_dir / "bundle_manifest.json").read_text())
    assert manifest["extensions"]["protein"]["normalization"] == "CLR"
    assert manifest["extensions"]["protein"]["files"] == ["protein"]
    assert manifest["extensions"]["protein"]["isotype_controls"] == list(isotypes)
    assert manifest["extensions"]["protein"]["claim_guard"] == (
        "not_for_de_or_new_quantitative_claims_without_full_object_validation"
    )
    assert (bundle_dir / "protein.parquet").exists()
    assert "protein" in manifest["files"]

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"m <- b$extensions$protein$metadata; "
        f"d <- b$extensions$protein$data; "
        f"cat(sprintf('NORM=%s\\n', m$normalization)); "
        f"cat(sprintf('GUARD=%s\\n', m$claim_guard)); "
        f"cat(sprintf('NROW=%d\\n', nrow(d))); "
        f"cat(sprintf('NCOL=%d\\n', ncol(d))); "
        f"cat(sprintf('ROWS_MATCH=%s\\n', as.character(identical(rownames(d), rownames(b$obs))))); "
        f"cat(sprintf('FIRST_COL=%s\\n', colnames(d)[1]))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "protein_round_trip")
    kv = _parse_kv_stdout(proc.stdout)
    assert kv["NORM"] == "CLR"
    assert kv["GUARD"] == (
        "not_for_de_or_new_quantitative_claims_without_full_object_validation"
    )
    assert int(kv["NROW"]) == n_cells
    assert int(kv["NCOL"]) == n_proteins
    assert kv["ROWS_MATCH"] == "TRUE"
    assert kv["FIRST_COL"] == protein_names[0]


def test_protein_module_python_unit(tmp_path):
    """ProteinADTModule.run wires obsm/obs/uns correctly and respects the DSB stub."""
    from workflow.modular.modules.protein_adt import (
        ProteinADTConfig,
        ProteinADTModule,
    )
    from workflow.modular.context import PipelineContext
    from workflow.modular.config import CellRangerConfig, PipelineConfig

    n_cells = 8
    n_proteins = 10
    isotypes = ("ISO_IgG1", "ISO_IgG2")

    rng = np.random.default_rng(11)
    x = sparse.csr_matrix(rng.integers(0, 5, size=(n_cells, 4)).astype(np.float32))
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=["A", "B", "C", "D"])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    counts = sparse.csr_matrix(
        rng.integers(0, 50, size=(n_cells, n_proteins)).astype(np.int32)
    )
    adata.obsm["protein_counts"] = counts
    protein_names = list(isotypes) + [f"ADT_{i + 1}" for i in range(n_proteins - len(isotypes))]
    adata.uns["protein_names"] = protein_names

    cfg = PipelineConfig(
        project="protein_unit",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path / "outs"),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True)
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    ctx.adata = adata
    ctx.set_module_dir("protein_adt")

    mod = ProteinADTModule(ProteinADTConfig(isotype_controls=isotypes, normalization="CLR"))
    mod.run(ctx)

    # obsm key created with the right shape and dtype.
    assert "protein_clr" in adata.obsm
    clr = adata.obsm["protein_clr"]
    assert clr.shape == (n_cells, n_proteins)
    assert clr.dtype == np.float32
    # CLR rows sum to ~0 by construction (mean-centered in log space).
    np.testing.assert_allclose(clr.sum(axis=1), np.zeros(n_cells), atol=1e-4)

    # Per-cell QC + multimodal flag.
    assert "protein_total_counts" in adata.obs.columns
    assert "protein_n_detected" in adata.obs.columns
    assert "X_wnn_prep" in adata.obs.columns
    assert adata.obs["X_wnn_prep"].dtype == bool

    # Isotype flags present and aligned with the input names.
    qc = adata.uns["protein_qc"]
    assert qc["normalization"] == "CLR"
    assert qc["panel_size"] == n_proteins
    assert qc["isotype_flags"][:2] == [True, True]
    assert sum(qc["isotype_flags"]) == len(isotypes)

    # QC plot was written.
    assert (run_dir / "protein_adt" / "protein_detection_rate_violin.png").exists()

    # DSB request must raise NotImplementedError -- never silently degrade.
    adata2 = ad.AnnData(
        X=sparse.csr_matrix(np.zeros((4, 2), dtype=np.float32)),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(4)]),
        var=pd.DataFrame(index=["A", "B"]),
    )
    adata2.obsm["protein_counts"] = sparse.csr_matrix(
        np.array([[1, 2], [3, 4], [5, 6], [7, 8]], dtype=np.int32)
    )
    ctx2 = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    ctx2.adata = adata2
    ctx2.set_module_dir("protein_adt_dsb")
    mod_dsb = ProteinADTModule(ProteinADTConfig(normalization="DSB"))
    with pytest.raises(NotImplementedError, match="DSB normalization"):
        mod_dsb.run(ctx2)


def test_protein_extension_omitted_when_flag_off(tmp_path):
    """Default export path produces no protein extension or parquet, even when adata has the modality."""
    h5ad = tmp_path / "tiny_protein.h5ad"
    bundle_dir = tmp_path / "bundle_no_protein"
    _make_tiny_protein_h5ad(h5ad, n_cells=8, n_proteins=10)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            format="parquet",
            # include_protein defaults to False -- explicit assertion.
        )
    )

    manifest = json.loads((bundle_dir / "bundle_manifest.json").read_text())
    # v2.1 default still emits an empty extensions field.
    assert manifest["schema_version"] == "singlecell_r_bundle_v2.1"
    assert manifest.get("extensions") == {}
    assert "protein" not in manifest["files"]
    assert not (bundle_dir / "protein.parquet").exists()


# ---------------------------------------------------------------------------
# Phase B (v2.1) spatial transcriptomics extension tests
# ---------------------------------------------------------------------------

def _make_tiny_spatial_h5ad(
    path: Path,
    *,
    n_cells: int = 200,
    n_genes: int = 10,
    libraries: tuple[str, ...] = ("sample_A", "sample_B"),
    image_paths: dict[str, str] | None = None,
) -> tuple[np.ndarray, list[str]]:
    """Synthesize an AnnData with a populated spatial obsm + RNA matrix.

    Returns (coords, library_assignment) so tests can assert round-trip values.
    """
    rng = np.random.default_rng(13)
    x = sparse.csr_matrix(
        rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    )
    gene_names = [f"GENE_{i + 1}" for i in range(n_genes)]
    # Make sure CD3E / LYZ exist for the v2 marker_expr requirement.
    gene_names[0] = "CD3E"
    gene_names[1] = "LYZ"
    obs = pd.DataFrame(
        {
            "cell_type": [["T", "B", "Myeloid"][i % 3] for i in range(n_cells)],
            "leiden": [str(i % 3) for i in range(n_cells)],
        },
        index=[f"spot_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=gene_names)
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.random((n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((n_cells, 3)).astype(np.float32)

    coords = rng.uniform(0, 1000, size=(n_cells, 2)).astype(np.float32)
    adata.obsm["spatial"] = coords

    library_assignment = [libraries[i % len(libraries)] for i in range(n_cells)]
    adata.obs["spatial_sample"] = pd.Categorical(library_assignment)
    adata.obs["spatial_library_id"] = pd.Categorical(library_assignment)

    adata.uns["spatial_coord_system"] = {
        "platform": "visium",
        "units": "pixel",
        "in_tissue_only": True,
    }
    if image_paths is not None:
        adata.uns["spatial"] = {"library_id": dict(image_paths)}

    adata.write_h5ad(path)
    return coords, library_assignment


@pytest.mark.r_contract
def test_spatial_extension_round_trip(rscript_path, tmp_path):
    """Python writes spatial.parquet; R reads it back with matching shape and rownames."""
    h5ad = tmp_path / "tiny_spatial.h5ad"
    bundle_dir = tmp_path / "bundle_v2_spatial"
    n_cells = 200
    image_paths = {
        "sample_A": "/data/visium/sampleA/tissue_hires_image.png",
        "sample_B": "/data/visium/sampleB/tissue_hires_image.png",
    }
    _make_tiny_spatial_h5ad(h5ad, n_cells=n_cells, image_paths=image_paths)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2.1",
            format="parquet",
            include_spatial=True,
            spatial_obsm_key="spatial",
        )
    )

    # Python-side manifest assertions.
    manifest = json.loads((bundle_dir / "bundle_manifest.json").read_text())
    spatial_ext = manifest["extensions"]["spatial"]
    assert spatial_ext["files"] == ["spatial"]
    assert spatial_ext["claim_guard"] == (
        "not_for_de_or_new_quantitative_claims_without_full_object_validation"
    )
    assert spatial_ext["coord_system"]["platform"] == "visium"
    assert spatial_ext["coord_system"]["units"] == "pixel"
    assert spatial_ext["library_image_paths"] == image_paths
    assert (bundle_dir / "spatial.parquet").exists()
    assert "spatial" in manifest["files"]

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"m <- b$extensions$spatial$metadata; "
        f"d <- b$extensions$spatial$data; "
        f"cat(sprintf('PLATFORM=%s\\n', m$coord_system$platform)); "
        f"cat(sprintf('UNITS=%s\\n', m$coord_system$units)); "
        f"cat(sprintf('GUARD=%s\\n', m$claim_guard)); "
        f"cat(sprintf('NROW=%d\\n', nrow(d$coords))); "
        f"cat(sprintf('NCOL=%d\\n', ncol(d$coords))); "
        f"cat(sprintf('ROWS_MATCH=%s\\n', as.character(identical(rownames(d$coords), rownames(b$obs))))); "
        f"cat(sprintf('HAS_LIB_PATHS=%s\\n', as.character(!is.null(d$library_image_paths)))); "
        f"cat(sprintf('LIB_A=%s\\n', d$library_image_paths$sample_A))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "spatial_round_trip")
    kv = _parse_kv_stdout(proc.stdout)
    assert kv["PLATFORM"] == "visium"
    assert kv["UNITS"] == "pixel"
    assert kv["GUARD"] == (
        "not_for_de_or_new_quantitative_claims_without_full_object_validation"
    )
    assert int(kv["NROW"]) == n_cells
    assert int(kv["NCOL"]) == 2
    assert kv["ROWS_MATCH"] == "TRUE"
    assert kv["HAS_LIB_PATHS"] == "TRUE"
    assert kv["LIB_A"] == image_paths["sample_A"]


def test_spatial_module_python_unit(tmp_path):
    """SpatialIngestModule preserves obsm and writes coord-system metadata."""
    from workflow.modular.modules.spatial_ingest import (
        SpatialIngestConfig,
        SpatialIngestModule,
    )
    from workflow.modular.context import PipelineContext
    from workflow.modular.config import CellRangerConfig, PipelineConfig

    n_cells = 12
    rng = np.random.default_rng(3)
    x = sparse.csr_matrix(rng.integers(0, 5, size=(n_cells, 4)).astype(np.float32))
    obs = pd.DataFrame(index=[f"spot_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=["A", "B", "C", "D"])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    coords = rng.uniform(0, 100, size=(n_cells, 2)).astype(np.float32)
    adata.obsm["spatial"] = coords

    cfg = PipelineConfig(
        project="spatial_unit",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path / "outs"),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True)
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    ctx.adata = adata
    ctx.set_module_dir("spatial_ingest")

    image_paths = {"libA": "/data/libA/image.png"}
    mod = SpatialIngestModule(SpatialIngestConfig(
        platform="visium", library_image_paths=image_paths,
    ))
    mod.run(ctx)

    # obsm preserved.
    assert "spatial" in adata.obsm
    np.testing.assert_array_equal(np.asarray(adata.obsm["spatial"]), coords)
    # uns metadata correct.
    assert adata.uns["spatial_coord_system"]["platform"] == "visium"
    assert adata.uns["spatial_coord_system"]["units"] == "pixel"
    assert adata.uns["spatial_coord_system"]["in_tissue_only"] is True
    assert adata.uns["spatial"]["library_id"] == {"libA": "/data/libA/image.png"}
    assert ctx.metadata["spatial_status"] == "ok"
    assert ctx.metadata["spatial_n_spots"] == n_cells
    assert ctx.metadata["spatial_platform"] == "visium"


def test_spatial_squidpy_skip_clean(tmp_path, monkeypatch):
    """SpatialNeighborhoodsModule exits cleanly when squidpy is not importable."""
    import builtins

    from workflow.modular.modules.spatial_neighborhoods import (
        SpatialNeighborhoodsConfig,
        SpatialNeighborhoodsModule,
    )
    from workflow.modular.context import PipelineContext
    from workflow.modular.config import CellRangerConfig, PipelineConfig

    n_cells = 8
    rng = np.random.default_rng(5)
    x = sparse.csr_matrix(rng.integers(0, 5, size=(n_cells, 4)).astype(np.float32))
    obs = pd.DataFrame(
        {"leiden": [str(i % 2) for i in range(n_cells)]},
        index=[f"spot_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=["A", "B", "C", "D"])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["spatial"] = rng.uniform(0, 100, size=(n_cells, 2)).astype(np.float32)

    cfg = PipelineConfig(
        project="spatial_skip",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path / "outs"),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True)
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    ctx.adata = adata
    ctx.set_module_dir("spatial_neighborhoods")

    # Force ImportError for squidpy regardless of whether it's actually installed.
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "squidpy" or name.startswith("squidpy."):
            raise ImportError("simulated: squidpy not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    mod = SpatialNeighborhoodsModule(SpatialNeighborhoodsConfig(group_by="leiden"))
    # Must not raise.
    mod.run(ctx)

    assert ctx.metadata["spatial_neighborhoods_status"] == "skipped_no_squidpy"
    # Status entry recorded as skipped.
    statuses = [s for s in ctx.module_status if s.get("module") == "spatial_neighborhoods"]
    assert statuses, "expected a module_status entry for spatial_neighborhoods"
    assert statuses[-1]["status"] == "skipped"


def test_spatial_extension_image_paths_string_only(tmp_path):
    """Tissue images must NEVER be serialized into the bundle; only string paths."""
    h5ad = tmp_path / "tiny_spatial_imgs.h5ad"
    bundle_dir = tmp_path / "bundle_spatial_imgs"
    image_paths = {
        "sample_A": "/data/visium/sampleA/tissue_hires_image.png",
        "sample_B": "/data/visium/sampleB/tissue_hires_image.tiff",
    }
    _make_tiny_spatial_h5ad(h5ad, n_cells=64, image_paths=image_paths)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2.1",
            format="parquet",
            include_spatial=True,
        )
    )

    # Walk the bundle dir and assert no image files were materialized.
    image_suffixes = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".webp", ".bmp"}
    for path in bundle_dir.rglob("*"):
        if path.is_file():
            assert path.suffix.lower() not in image_suffixes, (
                f"Bundle must not contain image files; found: {path}"
            )

    # Manifest must have STRING paths only (not bytes/None).
    manifest = json.loads((bundle_dir / "bundle_manifest.json").read_text())
    lib_paths = manifest["extensions"]["spatial"]["library_image_paths"]
    assert lib_paths == image_paths
    for value in lib_paths.values():
        assert isinstance(value, str), f"library_image_paths must be strings; got {type(value)}"


def test_spatial_extension_omitted_when_flag_off(tmp_path):
    """Default export path produces no spatial extension or parquet."""
    h5ad = tmp_path / "tiny_spatial_off.h5ad"
    bundle_dir = tmp_path / "bundle_no_spatial"
    _make_tiny_spatial_h5ad(h5ad, n_cells=32)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            format="parquet",
            # include_spatial defaults to False -- explicit assertion.
        )
    )

    manifest = json.loads((bundle_dir / "bundle_manifest.json").read_text())
    assert manifest["schema_version"] == "singlecell_r_bundle_v2.1"
    assert manifest.get("extensions") == {}
    assert "spatial" not in manifest["files"]
    assert not (bundle_dir / "spatial.parquet").exists()


# ---------------------------------------------------------------------------
# Phase B (v2.1) multimodal_obsm extension tests (EXPERIMENTAL)
# ---------------------------------------------------------------------------

def _make_tiny_multimodal_h5ad(
    path: Path,
    *,
    n_cells: int = 8,
    n_genes: int = 5,
    wnn_dims: int = 2,
) -> np.ndarray:
    """Synthesize an AnnData with a populated X_wnn obsm + RNA matrix.

    The WNN matrix is generated directly (no Seurat call); this is a unit
    fixture for the round-trip wiring, not an integration test.
    """
    rng = np.random.default_rng(17)
    x = sparse.csr_matrix(
        rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    )
    obs = pd.DataFrame(
        {
            "cell_type": [["T", "B", "Myeloid"][i % 3] for i in range(n_cells)],
            "leiden": [str(i % 3) for i in range(n_cells)],
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=["CD3E", "LYZ", "MS4A1", "NKG7", "ELF3"][:n_genes])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.random((n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((n_cells, 3)).astype(np.float32)
    wnn = rng.random((n_cells, wnn_dims)).astype(np.float32)
    adata.obsm["X_wnn"] = wnn
    adata.uns["multimodal_status"] = {"engine": "wnn", "status": "ok"}
    adata.write_h5ad(path)
    return wnn


@pytest.mark.r_contract
def test_multimodal_extension_round_trip(rscript_path, tmp_path):
    """Python writes multimodal_obsm_X_wnn.parquet; R reads it back with matching shape."""
    h5ad = tmp_path / "tiny_mm.h5ad"
    bundle_dir = tmp_path / "bundle_v2_mm"
    n_cells = 8
    _make_tiny_multimodal_h5ad(h5ad, n_cells=n_cells)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2.1",
            format="parquet",
            include_multimodal_obsm=True,
            multimodal_obsm_keys=("X_wnn", "X_mofa"),
        )
    )

    # Python-side manifest assertions.
    manifest = json.loads((bundle_dir / "bundle_manifest.json").read_text())
    mm_ext = manifest["extensions"]["multimodal_obsm"]
    assert mm_ext["files"] == ["multimodal_obsm_X_wnn"]
    assert mm_ext["engine_used"] == "wnn"
    assert mm_ext["experimental"] is True
    assert mm_ext["claim_guard"] == (
        "not_for_de_or_new_quantitative_claims_without_full_object_validation"
    )
    assert (bundle_dir / "multimodal_obsm_X_wnn.parquet").exists()
    assert "multimodal_obsm_X_wnn" in manifest["files"]
    # Mofa absent on the adata -> no parquet, no extra entry
    assert not (bundle_dir / "multimodal_obsm_X_mofa.parquet").exists()

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"m <- b$extensions$multimodal_obsm$metadata; "
        f"d <- b$extensions$multimodal_obsm$data; "
        f"cat(sprintf('ENGINE=%s\\n', m$engine_used)); "
        f"cat(sprintf('GUARD=%s\\n', m$claim_guard)); "
        f"cat(sprintf('NROW=%d\\n', nrow(d$X_wnn))); "
        f"cat(sprintf('NCOL=%d\\n', ncol(d$X_wnn))); "
        f"cat(sprintf('ROWS_MATCH=%s\\n', as.character(identical(rownames(d$X_wnn), rownames(b$obs)))))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "multimodal_round_trip")
    kv = _parse_kv_stdout(proc.stdout)
    assert kv["ENGINE"] == "wnn"
    assert kv["GUARD"] == (
        "not_for_de_or_new_quantitative_claims_without_full_object_validation"
    )
    assert int(kv["NROW"]) == n_cells
    assert int(kv["NCOL"]) == 2
    assert kv["ROWS_MATCH"] == "TRUE"


@pytest.mark.r_contract
def test_atac_extension_round_trip(rscript_path, tmp_path):
    """Python writes v2.2 ATAC LSI/peaks; R attaches the active ATAC payload."""
    h5ad = tmp_path / "tiny_atac.h5ad"
    bundle_dir = tmp_path / "bundle_v2_atac"
    n_cells = 6
    n_peaks = 4
    n_lsi = 2
    _make_tiny_atac_h5ad(h5ad, n_cells=n_cells, n_peaks=n_peaks, n_lsi=n_lsi)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2.2",
            format="parquet",
            include_atac=True,
            atac_lsi_obsm_key="X_lsi",
            atac_peaks_uns_key="atac_peaks",
        )
    )

    manifest = json.loads((bundle_dir / "bundle_manifest.json").read_text())
    assert manifest["schema_version"] == "singlecell_r_bundle_v2.2"
    atac_ext = manifest["extensions"]["atac"]
    assert atac_ext["status"] == "active"
    assert atac_ext["files"] == ["atac_lsi", "atac_peaks"]
    assert atac_ext["n_cells"] == n_cells
    assert atac_ext["n_components"] == n_lsi
    assert atac_ext["n_peaks"] == n_peaks
    assert atac_ext["table"]["lsi_index_column"] == "cell"
    assert atac_ext["table"]["peak_id_column"] == "peak_id"
    peaks_df = pd.read_parquet(bundle_dir / "extensions" / "atac" / "peaks.parquet")
    assert list(peaks_df.columns) == ["peak_id", "chrom", "start", "end"]

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"source('{ATAC_MODULE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"m <- b$extensions$atac$metadata; "
        f"d <- b$extensions$atac$data; "
        f"cat(sprintf('STATUS=%s\\n', m$status)); "
        f"cat(sprintf('LSI_NROW=%d\\n', nrow(d$lsi))); "
        f"cat(sprintf('LSI_NCOL=%d\\n', ncol(d$lsi))); "
        f"cat(sprintf('PEAKS_NROW=%d\\n', nrow(d$peaks))); "
        f"cat(sprintf('ROWS_MATCH=%s\\n', as.character(identical(rownames(d$lsi), rownames(b$obs))))); "
        f"cat(sprintf('PEAK_ID=%s\\n', d$peaks$peak_id[1])); "
        f"p <- plot_atac_lsi(b, group_by='cell_type'); "
        f"cols <- unique(ggplot2::ggplot_build(p)$data[[1]]$colour); "
        f"cat(sprintf('PLOT_COLORS=%d\\n', length(cols)))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "atac_round_trip")
    kv = _parse_kv_stdout(proc.stdout)
    assert kv["STATUS"] == "active"
    assert int(kv["LSI_NROW"]) == n_cells
    assert int(kv["LSI_NCOL"]) == n_lsi
    assert int(kv["PEAKS_NROW"]) == n_peaks
    assert kv["ROWS_MATCH"] == "TRUE"
    assert kv["PEAK_ID"].startswith("chr1:")
    assert int(kv["PLOT_COLORS"]) >= 2


@pytest.mark.r_contract
@pytest.mark.parametrize(
    ("extension", "loader_env"),
    [
        ("protein", "PROTEIN_MODULE_R_PATH"),
        ("spatial", "SPATIAL_MODULE_R_PATH"),
        ("multimodal_obsm", "INTEGRATION_MODULE_R_PATH"),
        ("atac", "ATAC_MODULE_R_PATH"),
        ("hic", "HIC_MODULE_R_PATH"),
    ],
)
def test_requested_active_extension_loader_failure_is_hard_error(
    rscript_path, tmp_path, extension, loader_env
):
    """An advertised active extension may not degrade to metadata-only success."""
    h5ad = tmp_path / f"broken_{extension}_loader.h5ad"
    bundle_dir = tmp_path / f"bundle_broken_{extension}_loader"
    export_kwargs: dict[str, object] = {
        "input_h5ad": h5ad,
        "output_dir": bundle_dir,
        "obs_columns": ("cell_type", "leiden"),
        "markers": ("CD3E", "LYZ"),
        "obsm_keys": ("X_umap", "X_pca"),
        "schema_version": "v2.1",
        "format": "parquet",
    }
    if extension == "protein":
        _make_tiny_protein_h5ad(h5ad)
        export_kwargs["include_protein"] = True
        export_kwargs["protein_obsm_key"] = "protein_clr"
    elif extension == "spatial":
        _make_tiny_spatial_h5ad(h5ad, n_cells=8)
        export_kwargs["include_spatial"] = True
        export_kwargs["spatial_obsm_key"] = "spatial"
    elif extension == "multimodal_obsm":
        _make_tiny_multimodal_h5ad(h5ad, n_cells=8)
        export_kwargs["include_multimodal_obsm"] = True
        export_kwargs["multimodal_obsm_keys"] = ("X_wnn",)
    elif extension == "atac":
        _make_tiny_atac_h5ad(h5ad)
        export_kwargs["schema_version"] = "v2.2"
        export_kwargs["include_atac"] = True
        export_kwargs["atac_lsi_obsm_key"] = "X_lsi"
        export_kwargs["atac_peaks_uns_key"] = "atac_peaks"
    else:
        _make_tiny_hic_h5ad(h5ad)
        export_kwargs["schema_version"] = "v2.2"
        export_kwargs["include_hic"] = True
        export_kwargs["hic_max_contacts"] = 100

    export_bundle(ExportConfig(**export_kwargs))
    broken_loader = tmp_path / f"broken_{extension}_module.R"
    broken_loader.write_text(
        f"stop('intentional {extension} loader failure')\n",
        encoding="utf-8",
    )
    env = dict(os.environ)
    env[loader_env] = str(broken_loader)
    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    proc = subprocess.run(
        [
            rscript_path,
            "-e",
            f"source('{IO_BUNDLE_R}'); read_bundle_v2('{bundle_dir_r}')",
        ],
        capture_output=True,
        text=True,
        timeout=180,
        env=env,
    )

    combined = proc.stdout + proc.stderr
    assert proc.returncode != 0, (
        f"requested {extension} loader failure was silently accepted:\n{combined}"
    )
    assert f"requested {extension} extension loader failed" in combined


@pytest.mark.r_contract
def test_hic_extension_round_trip(rscript_path, tmp_path):
    """Python writes v2.2 HIC tables; R attaches active contacts/bins/boundaries/compartments."""
    h5ad = tmp_path / "tiny_hic.h5ad"
    bundle_dir = tmp_path / "bundle_v2_hic"
    n_bins = 6
    _make_tiny_hic_h5ad(h5ad, n_bins=n_bins)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2.2",
            format="parquet",
            include_hic=True,
            hic_max_contacts=100,
        )
    )

    manifest = json.loads((bundle_dir / "bundle_manifest.json").read_text())
    assert manifest["schema_version"] == "singlecell_r_bundle_v2.2"
    hic_ext = manifest["extensions"]["hic"]
    assert hic_ext["status"] == "active"
    assert hic_ext["files"] == ["hic_bins", "hic_contacts", "hic_boundaries", "hic_compartments"]
    assert hic_ext["n_bins"] == n_bins
    assert hic_ext["n_boundary_bins"] == 2
    assert hic_ext["table"]["contacts_path"] == "extensions/hic/contacts.parquet"
    assert hic_ext["hic_tad_metadata"]["compartment_status"] == "confident"
    assert hic_ext["hic_tad_metadata"]["compartment_status_by_chrom"] == {"chr1": "confident"}
    boundaries_df = pd.read_parquet(bundle_dir / "extensions" / "hic" / "boundaries.parquet")
    assert "position" in boundaries_df.columns
    assert int(boundaries_df.loc[boundaries_df["bin_id"] == 2, "position"].item()) == 2500

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"source('{HIC_MODULE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"m <- b$extensions$hic$metadata; "
        f"d <- b$extensions$hic$data; "
        f"h <- load_hic_extension(b); "
        f"cat(sprintf('STATUS=%s\\n', m$status)); "
        f"cat(sprintf('COMPARTMENT_STATUS=%s\\n', m$hic_tad_metadata$compartment_status)); "
        f"cat(sprintf('LOW_INFO_N=%d\\n', length(m$hic_tad_metadata$low_information_chromosomes))); "
        f"cat(sprintf('CHR1_STATUS=%s\\n', m$hic_tad_metadata$compartment_status_by_chrom$chr1)); "
        f"cat(sprintf('LOAD_COMPARTMENT_STATUS=%s\\n', h$metadata$hic_tad_metadata$compartment_status)); "
        f"cat(sprintf('BINS_NROW=%d\\n', nrow(d$bins))); "
        f"cat(sprintf('CONTACTS_NROW=%d\\n', nrow(d$contacts))); "
        f"cat(sprintf('BOUNDARY_BINS=%d\\n', sum(d$boundaries$is_boundary))); "
        f"cat(sprintf('COMPARTMENTS=%s\\n', paste(sort(unique(d$compartments$compartment)), collapse=''))); "
        f"p <- plot_contact_heatmap(b, max_bins=10); "
        f"cat(sprintf('PLOT_CLASS=%s\\n', class(p)[1]))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "hic_round_trip")
    kv = _parse_kv_stdout(proc.stdout)
    assert kv["STATUS"] == "active"
    assert kv["COMPARTMENT_STATUS"] == "confident"
    assert int(kv["LOW_INFO_N"]) == 0
    assert kv["CHR1_STATUS"] == "confident"
    assert kv["LOAD_COMPARTMENT_STATUS"] == "confident"
    assert int(kv["BINS_NROW"]) == n_bins
    assert int(kv["CONTACTS_NROW"]) == hic_ext["n_contacts"]
    assert int(kv["BOUNDARY_BINS"]) == 2
    assert kv["COMPARTMENTS"] == "AB"
    assert "ggplot" in kv["PLOT_CLASS"]


@pytest.mark.r_contract
def test_hic_ab_plot_honors_low_information_metadata(rscript_path):
    """A/B plotting must not render low-information or unknown rows as A/B calls."""
    script = (
        f"source('{HIC_MODULE_R}'); "
        "bundle <- list(extensions=list(hic=list("
        "metadata=list(status='active', hic_tad_metadata=list("
        "compartment_status='partial_low_information', "
        "low_information_chromosomes=c('chr2'), "
        "compartment_status_by_chrom=list(chr1='confident', chr2='low_information'))), "
        "data=list(compartments=data.frame("
        "bin_id=0:3, chrom=c('chr1','chr1','chr2','chr2'), "
        "eigenvector_1=c(0.3,-0.2,0,0), "
        "compartment=c('A','B','low_information','low_information'), "
        "stringsAsFactors=FALSE))))); "
        "p <- plot_ab_compartments(bundle); "
        "cat(sprintf('AB_COLORS=%s\\n', paste(sort(unique(p$data$ab_color)), collapse=','))); "
        "cat(sprintf('LOW_INFO_N=%d\\n', sum(p$data$ab_color == 'low_information'))); "
        "unknown_bundle <- list(extensions=list(hic=list("
        "metadata=list(status='active', hic_tad_metadata=list("
        "compartment_status='unknown_or_unvalidated', "
        "low_information_chromosomes=c(), "
        "compartment_status_by_chrom=list())), "
        "data=list(compartments=data.frame("
        "bin_id=0:1, eigenvector_1=c(0.3,-0.2), compartment=c('A','B'), "
        "stringsAsFactors=FALSE))))); "
        "p2 <- plot_ab_compartments(unknown_bundle); "
        "cat(sprintf('UNKNOWN_COLORS=%s\\n', paste(sort(unique(p2$data$ab_color)), collapse=',')))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "hic_ab_plot_low_information")
    kv = _parse_kv_stdout(proc.stdout)
    assert kv["AB_COLORS"] == "A,B,low_information"
    assert int(kv["LOW_INFO_N"]) == 2
    assert kv["UNKNOWN_COLORS"] == "unknown_or_unvalidated"


@pytest.mark.r_contract
def test_exported_hic_low_information_round_trip_to_ab_plot(rscript_path, tmp_path):
    """Python-exported low-information HIC compartments stay non-A/B in R plots."""
    h5ad = tmp_path / "tiny_hic_low_information.h5ad"
    bundle_dir = tmp_path / "bundle_v2_hic_low_information"
    _make_tiny_hic_h5ad(h5ad, n_bins=6)
    adata = ad.read_h5ad(h5ad)
    compartments = adata.uns["hic_compartments"].copy()
    compartments["eigenvector_1"] = 0.0
    compartments["compartment"] = "low_information"
    adata.uns["hic_compartments"] = compartments
    adata.uns["hic_tad_metadata"] = {
        "compartment_status": "low_information",
        "low_information_chromosomes": ["chr1"],
        "compartment_status_by_chrom": {"chr1": "low_information"},
    }
    adata.write_h5ad(h5ad)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            schema_version="v2.2",
            format="parquet",
            include_hic=True,
            hic_max_contacts=100,
        )
    )

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    script = (
        f"source('{IO_BUNDLE_R}'); "
        f"source('{HIC_MODULE_R}'); "
        f"b <- read_bundle_v2('{bundle_dir_r}'); "
        f"h <- load_hic_extension(b); "
        f"p <- plot_ab_compartments(b); "
        f"cat(sprintf('LOAD_COMPARTMENT_STATUS=%s\\n', h$metadata$hic_tad_metadata$compartment_status)); "
        f"cat(sprintf('PARQUET_COMPARTMENTS=%s\\n', paste(sort(unique(h$compartments$compartment)), collapse=','))); "
        f"cat(sprintf('PLOT_COLORS=%s\\n', paste(sort(unique(p$data$ab_color)), collapse=','))); "
        f"cat(sprintf('LOW_INFO_N=%d\\n', sum(p$data$ab_color == 'low_information')))"
    )
    proc = _run_r(rscript_path, script)
    _assert_r_ok(proc, "exported_hic_low_information_plot")
    kv = _parse_kv_stdout(proc.stdout)
    assert kv["LOAD_COMPARTMENT_STATUS"] == "low_information"
    assert kv["PARQUET_COMPARTMENTS"] == "low_information"
    assert kv["PLOT_COLORS"] == "low_information"
    assert int(kv["LOW_INFO_N"]) == 6


def test_multimodal_module_off_by_default(tmp_path):
    """Default config (engine='off') is a no-op: no obsm keys added, status='off'."""
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    from workflow.modular.context import PipelineContext
    from workflow.modular.config import CellRangerConfig, PipelineConfig

    n_cells = 6
    rng = np.random.default_rng(19)
    x = sparse.csr_matrix(rng.integers(0, 5, size=(n_cells, 4)).astype(np.float32))
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=["A", "B", "C", "D"])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    # Even with both latents present, engine='off' must not act.
    adata.obsm["X_pca"] = rng.random((n_cells, 3)).astype(np.float32)
    adata.obsm["protein_clr"] = rng.random((n_cells, 4)).astype(np.float32)

    cfg = PipelineConfig(
        project="mm_off",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path / "outs"),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True)
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    ctx.adata = adata
    ctx.set_module_dir("multimodal_integration")

    mod = MultimodalIntegrationModule(MultimodalIntegrationConfig())
    # Make sure env override does not flip the default during this test.
    import os as _os
    prior_env = _os.environ.pop("SC_MULTIMODAL_ENGINE", None)
    try:
        mod.run(ctx)
    finally:
        if prior_env is not None:
            _os.environ["SC_MULTIMODAL_ENGINE"] = prior_env

    assert "X_wnn" not in adata.obsm
    assert "X_mofa" not in adata.obsm
    status = adata.uns["multimodal_status"]
    assert status["engine"] == "off"
    assert status["status"] == "skipped"
    assert ctx.metadata["multimodal_status"] == "off"


def test_multimodal_wnn_skips_when_r_missing(tmp_path, monkeypatch):
    """WNN engine skips cleanly when the R subprocess fails (return code != 0)."""
    import subprocess as _subprocess
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    from workflow.modular.context import PipelineContext
    from workflow.modular.config import CellRangerConfig, PipelineConfig

    n_cells = 6
    rng = np.random.default_rng(21)
    x = sparse.csr_matrix(rng.integers(0, 5, size=(n_cells, 4)).astype(np.float32))
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=["A", "B", "C", "D"])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_pca"] = rng.random((n_cells, 3)).astype(np.float32)
    adata.obsm["protein_clr"] = rng.random((n_cells, 4)).astype(np.float32)

    cfg = PipelineConfig(
        project="mm_wnn_skip",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path / "outs"),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True)
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    ctx.adata = adata
    ctx.set_module_dir("multimodal_integration")

    # Force the module to think Rscript exists, but make subprocess.run return code 127.
    import workflow.modular.modules.multimodal_integration as mm_mod

    monkeypatch.setattr(mm_mod.MultimodalIntegrationModule, "_resolve_rscript",
                        staticmethod(lambda cfg: "/usr/bin/Rscript"))
    # Also make sure the run_wnn.R driver "exists" so we go down the subprocess path.
    monkeypatch.setattr(mm_mod.Path, "exists", lambda self: True)

    class _Result:
        returncode = 127
        stderr = "Rscript: command not found"
        stdout = ""

    def fake_run(*args, **kwargs):
        return _Result()

    monkeypatch.setattr(mm_mod.subprocess, "run", fake_run)

    mod = MultimodalIntegrationModule(MultimodalIntegrationConfig(engine="wnn"))
    mod.run(ctx)

    assert "X_wnn" not in adata.obsm
    status = adata.uns["multimodal_status"]
    assert status["engine"] == "wnn"
    assert status["status"] == "skipped"
    assert "127" in status["reason"]
    assert ctx.metadata["multimodal_status"] == "skipped_driver_failed"


def test_multimodal_mofa_skips_when_muon_missing(tmp_path, monkeypatch):
    """MOFA engine skips cleanly when neither mofapy2 nor muon are importable."""
    import builtins as _builtins
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    from workflow.modular.context import PipelineContext
    from workflow.modular.config import CellRangerConfig, PipelineConfig

    n_cells = 6
    rng = np.random.default_rng(23)
    x = sparse.csr_matrix(rng.integers(0, 5, size=(n_cells, 4)).astype(np.float32))
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=["A", "B", "C", "D"])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_pca"] = rng.random((n_cells, 3)).astype(np.float32)
    adata.obsm["protein_clr"] = rng.random((n_cells, 4)).astype(np.float32)

    cfg = PipelineConfig(
        project="mm_mofa_skip",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path / "outs"),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True)
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    ctx.adata = adata
    ctx.set_module_dir("multimodal_integration")

    real_import = _builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "mofapy2" or name.startswith("mofapy2."):
            raise ImportError("simulated: mofapy2 not installed")
        if name == "muon" or name.startswith("muon."):
            raise ImportError("simulated: muon not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(_builtins, "__import__", fake_import)

    mod = MultimodalIntegrationModule(MultimodalIntegrationConfig(engine="mofa"))
    mod.run(ctx)

    assert "X_mofa" not in adata.obsm
    status = adata.uns["multimodal_status"]
    assert status["engine"] == "mofa"
    assert status["status"] == "skipped"
    assert "mofapy2" in status["reason"] or "muon" in status["reason"]
    assert ctx.metadata["multimodal_status"] == "skipped_no_mofa"


def test_multimodal_extension_omitted_when_flag_off(tmp_path):
    """Default export path produces no multimodal_obsm extension or parquets."""
    h5ad = tmp_path / "tiny_mm_off.h5ad"
    bundle_dir = tmp_path / "bundle_no_mm"
    _make_tiny_multimodal_h5ad(h5ad, n_cells=6)

    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("cell_type", "leiden"),
            markers=("CD3E", "LYZ"),
            obsm_keys=("X_umap", "X_pca"),
            format="parquet",
            # include_multimodal_obsm defaults to False -- explicit assertion.
        )
    )

    manifest = json.loads((bundle_dir / "bundle_manifest.json").read_text())
    assert manifest["schema_version"] == "singlecell_r_bundle_v2.1"
    assert manifest.get("extensions") == {}
    assert "multimodal_obsm" not in manifest.get("extensions", {})
    assert not (bundle_dir / "multimodal_obsm_X_wnn.parquet").exists()
    assert not (bundle_dir / "multimodal_obsm_X_mofa.parquet").exists()
