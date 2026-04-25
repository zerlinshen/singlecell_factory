"""Bundle round-trip parity gate.

Exports a synthetic AnnData via export_bundle (v1 schema), reads back the
marker_expr CSV, and verifies max abs diff vs original X[:,markers] <= 1e-6.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from scripts.export_singlecell_r_bundle import ExportConfig, export_bundle


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def bundle_adata_and_dir(tmp_path_factory):
    """Create a small AnnData, export it, return (adata, output_dir)."""
    rng = np.random.default_rng(42)
    n_cells, n_genes = 200, 50

    markers = ["CD3E", "LYZ", "MS4A1", "NKG7", "EPCAM"]
    other_genes = [f"GENE_{i}" for i in range(n_genes - len(markers))]
    var_names = markers + other_genes

    X = rng.uniform(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    X_sparse = sp.csr_matrix(X)

    adata = AnnData(X_sparse)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = var_names
    adata.obs["leiden"] = pd.Categorical(["0"] * n_cells)
    adata.obsm["X_umap"] = rng.normal(size=(n_cells, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.normal(size=(n_cells, 10)).astype(np.float32)

    # Write to temp h5ad
    tmp = tmp_path_factory.mktemp("bundle_rt")
    h5ad_path = tmp / "test.h5ad"
    adata.write_h5ad(h5ad_path)

    out_dir = tmp / "bundle_out"
    cfg = ExportConfig(
        input_h5ad=h5ad_path,
        output_dir=out_dir,
        markers=tuple(markers),
        obsm_keys=("X_umap", "X_pca"),
        compression="none",
        marker_chunk_size=100,
    )
    export_bundle(cfg)

    return adata, out_dir, markers, X


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_bundle_files_exist(bundle_adata_and_dir):
    _, out_dir, _, _ = bundle_adata_and_dir
    assert (out_dir / "obs.csv").exists()
    assert (out_dir / "marker_expr.csv").exists()
    assert (out_dir / "X_umap.csv").exists()
    assert (out_dir / "X_pca.csv").exists()
    assert (out_dir / "bundle_manifest.tsv").exists()


def test_bundle_roundtrip_max_abs_diff(bundle_adata_and_dir):
    """Max abs diff between exported marker_expr and original X[:,markers] <= 1e-6."""
    adata, out_dir, markers, X_dense = bundle_adata_and_dir

    exported = pd.read_csv(out_dir / "marker_expr.csv", index_col=0)
    # align columns to marker order
    exported = exported[markers]

    marker_idx = [list(adata.var_names).index(m) for m in markers]
    original = X_dense[:, marker_idx].astype(np.float32)

    exported_arr = exported.values.astype(np.float32)
    max_diff = np.abs(exported_arr - original).max()
    assert max_diff <= 1e-6, f"Bundle round-trip max abs diff = {max_diff:.3e} (threshold 1e-6)"


def test_bundle_roundtrip_cell_alignment(bundle_adata_and_dir):
    """Cell order in marker_expr must match obs order from the exported obs file."""
    adata, out_dir, markers, _ = bundle_adata_and_dir
    obs_df = pd.read_csv(out_dir / "obs.csv", index_col=0)
    expr_df = pd.read_csv(out_dir / "marker_expr.csv", index_col=0)
    assert list(obs_df.index) == list(expr_df.index), "Cell order mismatch between obs and marker_expr"


def test_bundle_manifest_schema_version(bundle_adata_and_dir):
    """Bundle manifest must record schema_version = singlecell_r_bundle_v1."""
    _, out_dir, _, _ = bundle_adata_and_dir
    manifest = pd.read_csv(out_dir / "bundle_manifest.tsv", sep="\t", index_col=0)
    version = manifest.loc["schema_version", "value"]
    assert version == "singlecell_r_bundle_v1"
