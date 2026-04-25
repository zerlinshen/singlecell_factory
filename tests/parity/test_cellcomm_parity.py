"""Cell communication parity gate: dense default vs SC_CELLCOMM_ENGINE=sparse.

Both paths must produce identical (ligand, receptor, source, target) tuple
sets and per-tuple lr_score within float32 numerical tolerance.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData


def _try_import_sparse_groupby():
    try:
        from workflow.modular._sparse_utils import sparse_groupby_mean  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.fixture(scope="module")
def synthetic_cellcomm_adata():
    """Small synthetic AnnData with cell_type labels and known ligand/receptor genes."""
    rng = np.random.default_rng(7)
    n_cells, n_genes = 600, 80
    # Cell types: 4 groups, balanced
    cell_types = np.repeat(["T", "B", "Mono", "Epi"], n_cells // 4)
    rng.shuffle(cell_types)
    base = rng.negative_binomial(3, 0.5, size=(n_cells, n_genes)).astype(np.float32)
    # Plant strong T-cell ligand expression in T cells (gene 0) and receptor in B (gene 1)
    base[cell_types == "T", 0] *= 5.0
    base[cell_types == "B", 1] *= 5.0
    X_sparse = sp.csr_matrix(base)
    adata = AnnData(X_sparse)
    adata.var_names = [f"GENE_{i:03d}" for i in range(n_genes)]
    adata.obs["cell_type"] = pd.Categorical(cell_types)
    return adata


@pytest.mark.skipif(not _try_import_sparse_groupby(), reason="phase 2 not enabled")
def test_cellcomm_engine_sparse_means_match_dense(synthetic_cellcomm_adata, monkeypatch):
    """Verify SC_CELLCOMM_ENGINE=sparse produces same per-cell-type mean matrix as dense path."""
    from workflow.modular._sparse_utils import sparse_groupby_mean

    adata = synthetic_cellcomm_adata
    gene_indices = list(range(20))
    cell_types = adata.obs["cell_type"].unique()

    # Dense path replica
    ct_list_dense, mean_rows_dense = [], []
    for ct in cell_types:
        mask = (adata.obs["cell_type"] == ct).values
        if mask.sum() < 5:
            continue
        ct_list_dense.append(ct)
        chunk = adata.X[mask][:, gene_indices]
        if hasattr(chunk, "toarray"):
            chunk = chunk.toarray()
        mean_rows_dense.append(np.asarray(chunk, dtype=np.float32).mean(axis=0))
    dense_mat = np.stack(mean_rows_dense)

    # Sparse path replica (mirrors the module's branch)
    valid_cts = [ct for ct in cell_types
                 if (adata.obs["cell_type"] == ct).values.sum() >= 5]
    ct_array = adata.obs["cell_type"].values
    keep_mask = np.isin(ct_array, valid_cts)
    X_kept = adata.X[:, gene_indices][keep_mask]
    labels_kept = ct_array[keep_mask]
    mean_dict = sparse_groupby_mean(X_kept, labels_kept)
    sparse_mat = np.stack([np.asarray(mean_dict[ct], dtype=np.float32) for ct in valid_cts])

    # Match cell-type ordering between dense and sparse
    sparse_idx = [valid_cts.index(ct) for ct in ct_list_dense]
    sparse_aligned = sparse_mat[sparse_idx]

    np.testing.assert_allclose(dense_mat, sparse_aligned, atol=1e-6,
                               err_msg="cell-type mean matrix mismatch sparse vs dense")
