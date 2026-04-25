from __future__ import annotations

import gzip
import json
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


def read_bundle_csv(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return pd.read_csv(handle, index_col=0)


def make_tiny_h5ad(path: Path, include_pca: bool = True) -> None:
    x = sparse.csr_matrix(
        np.array(
            [
                [0, 1, 2, 0],
                [3, 0, 0, 4],
                [0, 5, 0, 0],
                [6, 0, 7, 0],
            ],
            dtype=np.float32,
        )
    )
    obs = pd.DataFrame(
        {
            "cell_type": ["T", "Myeloid", "Tumor", "Tumor"],
            "leiden": ["0", "1", "1", "2"],
            "unused": ["a", "b", "c", "d"],
        },
        index=["c1", "c2", "c3", "c4"],
    )
    var = pd.DataFrame(index=["CD3E", "LYZ", "ELF3", "KRT8"])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = np.arange(8, dtype=np.float32).reshape(4, 2)
    if include_pca:
        adata.obsm["X_pca"] = np.arange(12, dtype=np.float32).reshape(4, 3)
    adata.write_h5ad(path)


def test_export_bundle_writes_sparse_safe_compact_contract(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny.h5ad"
    out_dir = tmp_path / "bundle"
    make_tiny_h5ad(input_path)

    manifest = export_bundle(
        ExportConfig(
            input_h5ad=input_path,
            output_dir=out_dir,
            source_run_dir=tmp_path,
            obs_columns=("cell_type", "leiden", "missing_col"),
            markers=("ELF3", "CD3E", "MISSING"),
            obsm_keys=("X_umap", "X_pca"),
            marker_chunk_size=2,
        )
    )

    assert manifest["schema_version"] == "singlecell_r_bundle_v1"
    assert manifest["source"]["n_obs"] == 4
    assert manifest["source"]["n_vars"] == 4
    assert manifest["bundle"]["markers_present"] == ["ELF3", "CD3E"]
    assert manifest["bundle"]["required_files"] == ["X_pca", "X_umap", "marker_expr", "obs"]
    assert manifest["expression"]["source_slot"] == "X"
    assert manifest["expression"]["value_scale"] == "source_X_as_stored"
    assert manifest["expression"]["export_dtype"] == "float32"
    assert manifest["expression"]["intended_use"] == "plotting_and_visual_summary_only"
    assert (
        manifest["expression"]["claim_guard"]
        == "not_for_de_or_new_quantitative_claims_without_full_object_validation"
    )
    assert "full expression matrix is not converted to dense" in manifest["precision_policy"]

    obs = read_bundle_csv(out_dir / "obs.csv.gz")
    marker = read_bundle_csv(out_dir / "marker_expr.csv.gz")
    umap = read_bundle_csv(out_dir / "X_umap.csv.gz")
    pca = read_bundle_csv(out_dir / "X_pca.csv.gz")

    assert obs.columns.tolist() == ["cell_type", "leiden"]
    assert marker.loc["c2", "CD3E"] == 3
    assert marker.loc["c3", "ELF3"] == 0
    assert umap.shape == (4, 2)
    assert pca.shape == (4, 3)

    manifest_disk = json.loads((out_dir / "bundle_manifest.json").read_text())
    assert manifest_disk["files"]["marker_expr"]["sha256"]
    tsv = pd.read_csv(out_dir / "bundle_manifest.tsv", sep="\t")
    assert "precision_policy" in set(tsv["key"])
    tsv_map = dict(zip(tsv["key"], tsv["value"]))
    assert tsv_map["expression_value_scale"] == "source_X_as_stored"
    assert tsv_map["file_marker_expr_n_rows"] == "4"
    assert tsv_map["file_marker_expr_n_cols"] == "2"


def test_export_bundle_can_deterministically_subset_cells(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny.h5ad"
    out_dir = tmp_path / "bundle"
    make_tiny_h5ad(input_path)

    export_bundle(
        ExportConfig(
            input_h5ad=input_path,
            output_dir=out_dir,
            markers=("CD3E",),
            max_cells=2,
            seed=7,
        )
    )

    obs = read_bundle_csv(out_dir / "obs.csv.gz")
    marker = read_bundle_csv(out_dir / "marker_expr.csv.gz")
    assert obs.shape[0] == 2
    assert marker.shape == (2, 1)


def test_export_bundle_fails_when_required_embedding_missing(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny_missing_pca.h5ad"
    out_dir = tmp_path / "bundle"
    make_tiny_h5ad(input_path, include_pca=False)

    with pytest.raises(ValueError, match="X_pca"):
        export_bundle(ExportConfig(input_h5ad=input_path, output_dir=out_dir))


def test_export_bundle_rejects_invalid_marker_chunk_size(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny.h5ad"
    out_dir = tmp_path / "bundle"
    make_tiny_h5ad(input_path)

    with pytest.raises(ValueError, match="marker-chunk-size"):
        export_bundle(ExportConfig(input_h5ad=input_path, output_dir=out_dir, marker_chunk_size=0))
