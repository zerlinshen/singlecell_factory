"""Contract definition test for the RNA-seq pipeline AnnData output schema.

This test verifies that a synthetic adata object can be constructed and
validated against the keys + dtypes named in:
  .omc/research/rna-multiomics-contract.md

This is a CONTRACT DEFINITION test, not a pipeline integration test.
The pipeline itself is not executed; keys are populated manually so that
downstream multi-omics modules (VDJ/ATAC/Hi-C/Ribo) can rely on the schema
being stable.

US-MO1a — Wave 3 (2026-05-15)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Fixture: synthetic 100-cell adata with all contracted keys populated
# ---------------------------------------------------------------------------

def _make_rna_contract_adata():
    """Build a minimal 100-cell × 200-gene AnnData with all contracted keys."""
    import anndata as ad

    n_obs = 100
    n_vars = 200
    rng = np.random.default_rng(42)

    # Raw count matrix (sparse, integer-like)
    X = sp.random(n_obs, n_vars, density=0.3, format="csr", dtype=np.float32,
                  random_state=rng)
    X.data = np.round(X.data * 10).astype(np.float32)

    obs_names = pd.Index([f"CELL_{i:04d}-1" for i in range(n_obs)], name=None)
    var_names = pd.Index([f"GENE_{j:04d}" for j in range(n_vars)])

    # --- obs columns ---
    leiden_labels = [str(i % 5) for i in range(n_obs)]
    cell_types = [
        "T cell", "B cell", "Myeloid/Macro", "Fibroblast", "Endothelial"
    ]
    obs = pd.DataFrame(
        {
            "sample": pd.array(
                [f"sample_{i % 3}" for i in range(n_obs)], dtype="object"
            ),
            "leiden": pd.Categorical(leiden_labels),
            "cell_type": pd.array(
                [cell_types[int(l) % len(cell_types)] for l in leiden_labels],
                dtype="object",
            ),
            "context_aware_celltype": pd.array(
                [cell_types[int(l) % len(cell_types)] for l in leiden_labels],
                dtype="object",
            ),
            "context_aware_substate": pd.array(
                [None if i % 3 != 0 else "CD8+ effector" for i in range(n_obs)],
                dtype="object",
            ),
            "annotation_confidence": rng.uniform(0.1, 1.0, n_obs).astype(np.float64),
        },
        index=obs_names,
    )

    # --- var columns ---
    var = pd.DataFrame(
        {
            "highly_variable": rng.choice([True, False], n_vars),
            "feature_type": pd.array(
                ["Gene Expression"] * n_vars, dtype="object"
            ),
        },
        index=var_names,
    )

    # --- obsm ---
    obsm = {
        "X_pca": rng.standard_normal((n_obs, 50)).astype(np.float32),
        "X_pca_harmony": rng.standard_normal((n_obs, 50)).astype(np.float32),
        "X_umap": rng.standard_normal((n_obs, 2)).astype(np.float32),
    }

    # --- uns ---
    uns = {
        "pca": {
            "variance_ratio": rng.dirichlet(np.ones(50)).astype(np.float64),
            "variance": rng.uniform(1.0, 5.0, 50).astype(np.float64),
        },
        "rank_genes_groups": {
            "params": {"groupby": "leiden", "method": "wilcoxon"},
            "names": {str(i): [f"GENE_{j:04d}" for j in range(10)] for i in range(5)},
        },
        "marker_db_index": {
            "CellMarker2": pd.DataFrame(
                {
                    "cell_type": ["T cell", "B cell"],
                    "marker_gene": ["CD3D", "CD79A"],
                    "marker_type": ["positive", "positive"],
                }
            )
        },
    }

    adata = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns)

    # Set adata.raw to simulate pre-normalization counts preservation
    adata.raw = adata.copy()

    return adata


@pytest.fixture(scope="module")
def rna_adata():
    return _make_rna_contract_adata()


# ---------------------------------------------------------------------------
# Section 1: obs keys
# ---------------------------------------------------------------------------

class TestObsContract:
    def test_sample_present(self, rna_adata):
        assert "sample" in rna_adata.obs.columns, "obs['sample'] must be present"

    def test_sample_dtype(self, rna_adata):
        assert rna_adata.obs["sample"].dtype == object, (
            "obs['sample'] must be object/str dtype"
        )

    def test_leiden_present(self, rna_adata):
        assert "leiden" in rna_adata.obs.columns, "obs['leiden'] must be present"

    def test_leiden_dtype(self, rna_adata):
        assert hasattr(rna_adata.obs["leiden"], "cat"), (
            "obs['leiden'] must be Categorical"
        )

    def test_leiden_labels_are_strings(self, rna_adata):
        cats = rna_adata.obs["leiden"].cat.categories.tolist()
        assert all(isinstance(c, str) for c in cats), (
            "obs['leiden'] categories must all be str"
        )

    def test_cell_type_present(self, rna_adata):
        assert "cell_type" in rna_adata.obs.columns, "obs['cell_type'] must be present"

    def test_cell_type_dtype(self, rna_adata):
        assert rna_adata.obs["cell_type"].dtype == object, (
            "obs['cell_type'] must be object/str dtype"
        )

    def test_context_aware_celltype_present(self, rna_adata):
        assert "context_aware_celltype" in rna_adata.obs.columns, (
            "obs['context_aware_celltype'] must be present"
        )

    def test_context_aware_celltype_dtype(self, rna_adata):
        assert rna_adata.obs["context_aware_celltype"].dtype == object, (
            "obs['context_aware_celltype'] must be object/str dtype"
        )

    def test_context_aware_substate_present(self, rna_adata):
        assert "context_aware_substate" in rna_adata.obs.columns, (
            "obs['context_aware_substate'] must be present when sub-states detected"
        )

    def test_context_aware_substate_dtype(self, rna_adata):
        # Column may be object with None values for cells with no substate
        assert rna_adata.obs["context_aware_substate"].dtype == object, (
            "obs['context_aware_substate'] must be object dtype"
        )

    def test_annotation_confidence_present(self, rna_adata):
        assert "annotation_confidence" in rna_adata.obs.columns, (
            "obs['annotation_confidence'] must be present"
        )

    def test_annotation_confidence_dtype(self, rna_adata):
        assert rna_adata.obs["annotation_confidence"].dtype == np.float64, (
            "obs['annotation_confidence'] must be float64"
        )

    def test_annotation_confidence_range(self, rna_adata):
        vals = rna_adata.obs["annotation_confidence"].dropna()
        assert (vals >= 0).all(), "annotation_confidence must be non-negative"


# ---------------------------------------------------------------------------
# Section 2: var keys
# ---------------------------------------------------------------------------

class TestVarContract:
    def test_highly_variable_present(self, rna_adata):
        assert "highly_variable" in rna_adata.var.columns, (
            "var['highly_variable'] must be present"
        )

    def test_highly_variable_dtype(self, rna_adata):
        assert rna_adata.var["highly_variable"].dtype == bool, (
            "var['highly_variable'] must be bool dtype"
        )

    def test_feature_type_present(self, rna_adata):
        assert "feature_type" in rna_adata.var.columns, (
            "var['feature_type'] must be present"
        )

    def test_feature_type_dtype(self, rna_adata):
        assert rna_adata.var["feature_type"].dtype == object, (
            "var['feature_type'] must be object/str dtype"
        )


# ---------------------------------------------------------------------------
# Section 3: obsm keys
# ---------------------------------------------------------------------------

class TestObsmContract:
    def test_X_pca_present(self, rna_adata):
        assert "X_pca" in rna_adata.obsm, "obsm['X_pca'] must be present"

    def test_X_pca_shape(self, rna_adata):
        arr = rna_adata.obsm["X_pca"]
        assert arr.ndim == 2, "obsm['X_pca'] must be 2-D"
        assert arr.shape[0] == rna_adata.n_obs, (
            "obsm['X_pca'] row count must equal n_obs"
        )
        assert arr.shape[1] >= 1, "obsm['X_pca'] must have at least 1 component"

    def test_X_pca_dtype(self, rna_adata):
        assert rna_adata.obsm["X_pca"].dtype == np.float32, (
            "obsm['X_pca'] must be float32"
        )

    def test_X_pca_harmony_present(self, rna_adata):
        assert "X_pca_harmony" in rna_adata.obsm, (
            "obsm['X_pca_harmony'] must be present (post-Harmony)"
        )

    def test_X_pca_harmony_shape(self, rna_adata):
        pca = rna_adata.obsm["X_pca"]
        harmony = rna_adata.obsm["X_pca_harmony"]
        assert harmony.shape == pca.shape, (
            "obsm['X_pca_harmony'] must have the same shape as obsm['X_pca']"
        )

    def test_X_pca_harmony_dtype(self, rna_adata):
        assert rna_adata.obsm["X_pca_harmony"].dtype == np.float32, (
            "obsm['X_pca_harmony'] must be float32"
        )

    def test_X_umap_present(self, rna_adata):
        assert "X_umap" in rna_adata.obsm, "obsm['X_umap'] must be present"

    def test_X_umap_shape(self, rna_adata):
        arr = rna_adata.obsm["X_umap"]
        assert arr.ndim == 2, "obsm['X_umap'] must be 2-D"
        assert arr.shape[0] == rna_adata.n_obs, (
            "obsm['X_umap'] row count must equal n_obs"
        )
        assert arr.shape[1] == 2, "obsm['X_umap'] must have exactly 2 columns"

    def test_X_umap_dtype(self, rna_adata):
        assert rna_adata.obsm["X_umap"].dtype == np.float32, (
            "obsm['X_umap'] must be float32"
        )


# ---------------------------------------------------------------------------
# Section 4: uns keys
# ---------------------------------------------------------------------------

class TestUnsContract:
    def test_pca_present(self, rna_adata):
        assert "pca" in rna_adata.uns, "uns['pca'] must be present"

    def test_pca_has_variance_ratio(self, rna_adata):
        assert "variance_ratio" in rna_adata.uns["pca"], (
            "uns['pca']['variance_ratio'] must be present"
        )

    def test_pca_variance_ratio_shape(self, rna_adata):
        vr = np.asarray(rna_adata.uns["pca"]["variance_ratio"])
        assert vr.ndim == 1, "uns['pca']['variance_ratio'] must be 1-D"
        assert len(vr) == rna_adata.obsm["X_pca"].shape[1], (
            "uns['pca']['variance_ratio'] length must equal n_pcs"
        )

    def test_rank_genes_groups_present(self, rna_adata):
        assert "rank_genes_groups" in rna_adata.uns, (
            "uns['rank_genes_groups'] must be present"
        )

    def test_rank_genes_groups_is_dict(self, rna_adata):
        assert isinstance(rna_adata.uns["rank_genes_groups"], dict), (
            "uns['rank_genes_groups'] must be a dict"
        )

    def test_marker_db_index_present(self, rna_adata):
        assert "marker_db_index" in rna_adata.uns, (
            "uns['marker_db_index'] must be present"
        )

    def test_marker_db_index_is_dict(self, rna_adata):
        assert isinstance(rna_adata.uns["marker_db_index"], dict), (
            "uns['marker_db_index'] must be a dict"
        )

    def test_marker_db_index_dataframe_schema(self, rna_adata):
        required_cols = {"cell_type", "marker_gene", "marker_type"}
        for db_name, df in rna_adata.uns["marker_db_index"].items():
            assert isinstance(df, pd.DataFrame), (
                f"uns['marker_db_index']['{db_name}'] must be a DataFrame"
            )
            missing = required_cols - set(df.columns)
            assert not missing, (
                f"uns['marker_db_index']['{db_name}'] missing columns: {missing}"
            )


# ---------------------------------------------------------------------------
# Section 5: raw counts location
# ---------------------------------------------------------------------------

class TestRawCountsContract:
    def test_raw_is_set(self, rna_adata):
        assert rna_adata.raw is not None, (
            "adata.raw must be set (current canonical raw counts location)"
        )

    def test_raw_var_names_match(self, rna_adata):
        assert set(rna_adata.raw.var_names) == set(rna_adata.var_names), (
            "adata.raw.var_names must match adata.var_names"
        )

    def test_raw_n_obs_matches(self, rna_adata):
        assert rna_adata.raw.X.shape[0] == rna_adata.n_obs, (
            "adata.raw.X row count must equal n_obs"
        )


# ---------------------------------------------------------------------------
# Section 6: cross-omics join key
# ---------------------------------------------------------------------------

class TestJoinKeyContract:
    def test_obs_index_is_unique(self, rna_adata):
        assert rna_adata.obs.index.is_unique, (
            "obs.index (cell barcodes) must be unique — required for cross-omics joins"
        )

    def test_obs_index_non_empty(self, rna_adata):
        assert len(rna_adata.obs.index) > 0, "obs.index must not be empty"

    def test_obs_index_no_nulls(self, rna_adata):
        assert not rna_adata.obs.index.isna().any(), (
            "obs.index must not contain NaN values"
        )

    def test_cross_omics_join_simulation(self, rna_adata):
        """Simulate a VDJ/ATAC subset join to verify barcode alignment works."""
        rng = np.random.default_rng(0)
        subset_barcodes = rna_adata.obs.index[
            rng.choice(rna_adata.n_obs, size=40, replace=False)
        ]
        other_obs = pd.DataFrame(
            {"vdj_chain": ["TRA"] * 40},
            index=subset_barcodes,
        )
        joined = rna_adata.obs.join(other_obs, how="left")
        assert len(joined) == rna_adata.n_obs, (
            "Left join on obs.index must preserve all RNA cells"
        )
        matched = joined["vdj_chain"].notna().sum()
        assert matched == 40, (
            f"Expected 40 matched barcodes in join, got {matched}"
        )
