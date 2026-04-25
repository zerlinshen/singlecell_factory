"""QC parity gate: sparse vs dense fallback path must produce fully equal metrics.

Threshold: np.array_equal for all five QC columns (exact equality, no tolerance).
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from anndata import AnnData

from workflow.modular.modules.qc import QCModule


@pytest.fixture(scope="module")
def qc_adata_pair():
    """Return (adata_sparse, adata_dense) with identical expression data."""
    rng = np.random.default_rng(7)
    n_cells, n_genes = 500, 300
    base = rng.negative_binomial(3, 0.4, size=(n_cells, n_genes)).astype(np.float64)

    # Add mt/ribo/hb genes at known positions
    mt_idx = list(range(10))
    ribo_idx = list(range(10, 30))
    hb_idx = list(range(30, 35))
    gene_names = (
        [f"MT-G{i}" for i in mt_idx]
        + [f"RPL{i}" for i in range(len(ribo_idx))]
        + [f"HBA{i}" for i in range(len(hb_idx))]
        + [f"GENE_{i}" for i in range(n_genes - len(mt_idx) - len(ribo_idx) - len(hb_idx))]
    )

    def _make_adata(X):
        adata = AnnData(X.copy())
        adata.var_names = gene_names
        adata.var["mt"] = [n.startswith("MT-") for n in gene_names]
        adata.var["ribo"] = [n.startswith(("RPS", "RPL")) for n in gene_names]
        adata.var["hb"] = [n.startswith(("HBA", "HBB")) for n in gene_names]
        return adata

    adata_sparse = _make_adata(sp.csr_matrix(base))
    adata_dense = _make_adata(base)
    return adata_sparse, adata_dense


QC_COLS = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]


def test_qc_fallback_sparse_vs_dense_exact(qc_adata_pair):
    """sparse and dense fallback QC metrics must be exactly equal (np.array_equal)."""
    adata_sparse, adata_dense = qc_adata_pair

    QCModule._calculate_qc_metrics_fallback(adata_sparse)
    QCModule._calculate_qc_metrics_fallback(adata_dense)

    for col in QC_COLS:
        sparse_vals = adata_sparse.obs[col].values
        dense_vals = adata_dense.obs[col].values
        assert np.array_equal(sparse_vals, dense_vals), (
            f"QC column '{col}' differs between sparse and dense paths.\n"
            f"  max abs diff = {np.abs(sparse_vals - dense_vals).max():.3e}"
        )


def test_qc_fallback_total_counts_correct(qc_adata_pair):
    """total_counts must equal row sums of the raw matrix."""
    adata_sparse, _ = qc_adata_pair
    QCModule._calculate_qc_metrics_fallback(adata_sparse)
    X = adata_sparse.X
    expected = np.asarray(X.sum(axis=1)).ravel()
    np.testing.assert_array_equal(adata_sparse.obs["total_counts"].values, expected)


def test_qc_fallback_pct_mt_in_range(qc_adata_pair):
    """pct_counts_mt must be in [0, 100]."""
    adata_sparse, _ = qc_adata_pair
    QCModule._calculate_qc_metrics_fallback(adata_sparse)
    pct = adata_sparse.obs["pct_counts_mt"].values
    assert (pct >= 0).all() and (pct <= 100).all()


def test_qc_fallback_n_genes_correct(qc_adata_pair):
    """n_genes_by_counts must equal count of non-zero columns per row."""
    adata_sparse, _ = qc_adata_pair
    QCModule._calculate_qc_metrics_fallback(adata_sparse)
    X = adata_sparse.X
    expected = np.asarray((X > 0).sum(axis=1)).ravel()
    np.testing.assert_array_equal(adata_sparse.obs["n_genes_by_counts"].values, expected)
