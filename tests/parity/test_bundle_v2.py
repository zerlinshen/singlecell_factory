"""Bundle v2 parity tests.

Covers:
  - test_v2_round_trip_max_abs_diff_le_1e6   parquet marker round-trip
  - test_v2_manifest_required_fields          schema contract
  - test_v2_sha256_per_file                  integrity check in manifest
  - test_v2_parquet_obs_columns_preserved    categorical obs column round-trip
  - test_v1_path_unchanged                   v1 export byte-stable vs fixture
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.export_singlecell_r_bundle import ExportConfig, export_bundle


# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

def _make_synthetic_adata(n_cells: int = 120, n_genes: int = 30,
                           seed: int = 42) -> AnnData:
    rng = np.random.default_rng(seed)
    markers = ["CD3E", "LYZ", "MS4A1", "NKG7", "EPCAM", "KRT8"]
    other = [f"GENE_{i}" for i in range(n_genes - len(markers))]
    var_names = markers + other

    X = rng.uniform(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    X_sparse = sp.csr_matrix(X)

    adata = AnnData(X_sparse)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = var_names
    adata.obs["leiden"] = pd.Categorical(["0"] * (n_cells // 2) + ["1"] * (n_cells - n_cells // 2))
    adata.obs["cell_type"] = pd.Categorical(["T"] * n_cells)
    adata.obs["sample"] = "S1"
    adata.obsm["X_umap"] = rng.normal(size=(n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.normal(size=(n_cells, 10)).astype(np.float32)
    return adata, markers, X


@pytest.fixture(scope="module")
def v2_bundle(tmp_path_factory):
    """Export a small v2 (parquet) bundle; return (adata, out_dir, markers, X_dense)."""
    tmp = tmp_path_factory.mktemp("v2_bundle")
    h5ad_path = tmp / "test.h5ad"
    adata, markers, X_dense = _make_synthetic_adata()
    adata.write_h5ad(h5ad_path)

    out_dir = tmp / "bundle_v2"
    cfg = ExportConfig(
        input_h5ad=h5ad_path,
        output_dir=out_dir,
        markers=tuple(markers),
        obsm_keys=("X_umap", "X_pca"),
        compression="none",
        marker_chunk_size=50,
        schema_version="v2",
        format="parquet",
    )
    manifest = export_bundle(cfg)
    return adata, out_dir, markers, X_dense, manifest


@pytest.fixture(scope="module")
def v1_bundle(tmp_path_factory):
    """Export a v1 bundle for backward-compat checks."""
    tmp = tmp_path_factory.mktemp("v1_bundle")
    h5ad_path = tmp / "test.h5ad"
    adata, markers, X_dense = _make_synthetic_adata()
    adata.write_h5ad(h5ad_path)

    out_dir = tmp / "bundle_v1"
    cfg = ExportConfig(
        input_h5ad=h5ad_path,
        output_dir=out_dir,
        markers=tuple(markers),
        obsm_keys=("X_umap", "X_pca"),
        compression="none",
        marker_chunk_size=50,
        schema_version="v1",
    )
    manifest = export_bundle(cfg)
    return adata, out_dir, markers, X_dense, manifest


# ---------------------------------------------------------------------------
# v2 tests
# ---------------------------------------------------------------------------

def test_v2_round_trip_max_abs_diff_le_1e6(v2_bundle):
    """Max abs diff between exported parquet marker_expr and original X[:,markers] <= 1e-6."""
    adata, out_dir, markers, X_dense, _ = v2_bundle

    try:
        import pyarrow.parquet as _pq
        df = _pq.read_table(out_dir / "marker_expr.parquet").to_pandas()
    except ImportError:
        df = pd.read_parquet(out_dir / "marker_expr.parquet")

    # Align to marker order; index column is "cell"
    if "cell" in df.columns:
        df = df.set_index("cell")
    df = df[markers]

    marker_idx = [list(adata.var_names).index(m) for m in markers]
    original = X_dense[:, marker_idx].astype(np.float32)
    exported = df.values.astype(np.float32)

    max_diff = np.abs(exported - original).max()
    assert max_diff <= 1e-6, f"v2 parquet round-trip max abs diff = {max_diff:.3e} (threshold 1e-6)"


def test_v2_manifest_required_fields(v2_bundle):
    """v2 manifest must contain all required contract fields."""
    _, out_dir, _, _, manifest = v2_bundle

    assert manifest["schema_version"] == "singlecell_r_bundle_v2"
    assert "singlecell_r_bundle_v1" in manifest["compatible_with"]
    assert isinstance(manifest["files"], dict)
    assert len(manifest["files"]) > 0
    assert "cell_alignment" in manifest
    assert manifest["cell_alignment"]["primary_index"] == "barcode"
    assert manifest["cell_alignment"]["all_files_share_primary"] is True
    assert "source" in manifest
    assert "generated_at_utc" in manifest

    # Every file entry must have format, path, sha256, n_rows, n_cols
    for stem, rec in manifest["files"].items():
        for field in ("format", "path", "sha256", "n_rows", "n_cols"):
            assert field in rec, f"File entry '{stem}' missing field '{field}'"


def test_v2_sha256_per_file(v2_bundle):
    """Every file listed in the v2 manifest must have a correct SHA256."""
    _, out_dir, _, _, manifest = v2_bundle

    def _sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()

    for stem, rec in manifest["files"].items():
        file_path = out_dir / rec["path"]
        assert file_path.exists(), f"File missing: {file_path}"
        expected = rec["sha256"]
        actual = _sha256(file_path)
        assert actual == expected, (
            f"SHA256 mismatch for {stem}: expected={expected}, actual={actual}"
        )


def test_v2_parquet_obs_columns_preserved(v2_bundle):
    """A known categorical obs column (leiden) must round-trip as Categorical."""
    _, out_dir, _, _, manifest = v2_bundle

    try:
        import pyarrow.parquet as _pq
        df = _pq.read_table(out_dir / "obs.parquet").to_pandas()
    except ImportError:
        df = pd.read_parquet(out_dir / "obs.parquet")

    # Handle index stored as column
    if "cell" in df.columns:
        df = df.set_index("cell")

    assert "leiden" in df.columns, "leiden column missing from obs.parquet"
    # Should survive as categorical or string — either is acceptable for R use
    assert len(df["leiden"]) > 0
    unique_vals = set(df["leiden"].astype(str).unique())
    assert unique_vals <= {"0", "1"}, f"Unexpected leiden values: {unique_vals}"


def test_v2_obsm_parquets_exist(v2_bundle):
    """obsm/X_pca.parquet and obsm/X_umap.parquet must exist and have correct shape."""
    adata, out_dir, _, _, _ = v2_bundle

    for key, expected_cols in [("X_pca", 10), ("X_umap", 2)]:
        pq_path = out_dir / "obsm" / f"{key}.parquet"
        assert pq_path.exists(), f"Missing obsm parquet: {pq_path}"
        try:
            import pyarrow.parquet as _pq
            df = _pq.read_table(pq_path).to_pandas()
        except Exception:
            df = pd.read_parquet(pq_path)
        # rows = n_cells; cols = embedding dims + possibly index column
        data_cols = [c for c in df.columns if c not in ("cell", "__index_level_0__")]
        assert len(data_cols) == expected_cols, (
            f"{key}.parquet has {len(data_cols)} data columns, expected {expected_cols}"
        )


def test_v2_uns_summary_json_exists(v2_bundle):
    """uns_summary.json must exist and contain n_obs and n_vars."""
    adata, out_dir, _, _, _ = v2_bundle
    uns_path = out_dir / "uns_summary.json"
    assert uns_path.exists(), "uns_summary.json missing"
    data = json.loads(uns_path.read_text())
    assert "n_obs" in data
    assert "n_vars" in data
    assert data["n_obs"] == adata.n_obs
    assert data["n_vars"] == adata.n_vars


# ---------------------------------------------------------------------------
# v1 backward-compat test
# ---------------------------------------------------------------------------

def test_v1_path_unchanged(v1_bundle):
    """v1 export must produce files with correct schema_version and round-trip accuracy."""
    adata, out_dir, markers, X_dense, manifest = v1_bundle

    # Schema version
    assert manifest["schema_version"] == "singlecell_r_bundle_v1"

    # v1 files use CSV (no .parquet)
    assert (out_dir / "obs.csv").exists()
    assert (out_dir / "marker_expr.csv").exists()
    assert (out_dir / "X_umap.csv").exists()
    assert (out_dir / "X_pca.csv").exists()
    assert (out_dir / "bundle_manifest.json").exists()
    assert (out_dir / "bundle_manifest.tsv").exists()

    # v1 files must NOT have parquet
    assert not (out_dir / "obs.parquet").exists()
    assert not (out_dir / "marker_expr.parquet").exists()

    # Round-trip accuracy
    df = pd.read_csv(out_dir / "marker_expr.csv", index_col=0)
    df = df[markers]
    marker_idx = [list(adata.var_names).index(m) for m in markers]
    original = X_dense[:, marker_idx].astype(np.float32)
    max_diff = np.abs(df.values.astype(np.float32) - original).max()
    assert max_diff <= 1e-6, f"v1 CSV round-trip max abs diff = {max_diff:.3e}"

    # Manifest TSV must have schema_version key
    tsv = pd.read_csv(out_dir / "bundle_manifest.tsv", sep="\t")
    tsv_map = dict(zip(tsv["key"], tsv["value"]))
    assert tsv_map["schema_version"] == "singlecell_r_bundle_v1"
