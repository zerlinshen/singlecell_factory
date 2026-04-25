import numpy as np
import pytest
import scipy.sparse as sp
import scipy.stats as stats

from workflow.modular._sparse_utils import (
    MemoryGuardError,
    chunked_row_densify,
    estimate_dense_bytes,
    safe_densify,
    sparse_groupby_mean,
    sparse_welch_t,
)


def _make_sparse(rng, n_cells, n_genes, density=0.3):
    data = rng.random((n_cells, n_genes)) * (rng.random((n_cells, n_genes)) < density)
    return sp.csr_matrix(data)


def test_estimate_dense_bytes():
    assert estimate_dense_bytes((100, 200), np.float32) == 100 * 200 * 4
    assert estimate_dense_bytes((50, 50), np.float64) == 50 * 50 * 8


def test_safe_densify_dense_input():
    arr = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
    result = safe_densify(arr, reason="test", dtype=np.float64)
    assert result.dtype == np.float64
    np.testing.assert_allclose(result, arr)


def test_safe_densify_sparse_small():
    rng = np.random.default_rng(0)
    X = _make_sparse(rng, 20, 10)
    result = safe_densify(X, reason="test", dtype=np.float64)
    assert result.dtype == np.float64
    np.testing.assert_allclose(result, X.toarray())


def test_safe_densify_raises_on_abort():
    from unittest.mock import patch
    from workflow.modular._densify_policy import DensifyDecision
    with patch("workflow.modular._sparse_utils.plan_densify", return_value=DensifyDecision.ABORT):
        X = sp.csr_matrix(np.ones((100, 100)))
        with pytest.raises(MemoryGuardError):
            safe_densify(X, reason="oom_test", dtype=np.float64)


def test_safe_densify_chunk_no_allow():
    from unittest.mock import patch
    from workflow.modular._densify_policy import DensifyDecision
    with patch("workflow.modular._sparse_utils.plan_densify", return_value=DensifyDecision.CHUNK):
        X = sp.csr_matrix(np.ones((50, 50)))
        with pytest.raises(MemoryGuardError):
            safe_densify(X, reason="chunk_test", dtype=np.float64, allow_chunk=False)


def test_chunked_row_densify_sparse():
    rng = np.random.default_rng(1)
    X = _make_sparse(rng, 30, 8)
    chunks = list(chunked_row_densify(X, chunk_rows=10, dtype=np.float32))
    assert len(chunks) == 3
    assert chunks[0].dtype == np.float32
    reconstructed = np.vstack(chunks)
    np.testing.assert_allclose(reconstructed, X.toarray(), atol=1e-6)


def test_chunked_row_densify_dense():
    arr = np.arange(40, dtype=np.float64).reshape(8, 5)
    chunks = list(chunked_row_densify(arr, chunk_rows=3, dtype=np.float32))
    assert sum(c.shape[0] for c in chunks) == 8


def test_sparse_groupby_mean_parity():
    rng = np.random.default_rng(2)
    n, g = 60, 15
    X = _make_sparse(rng, n, g, density=0.4)
    labels = np.array(["A"] * 20 + ["B"] * 20 + ["C"] * 20)
    result = sparse_groupby_mean(X, labels)
    dense = X.toarray()
    for grp in ["A", "B", "C"]:
        mask = labels == grp
        expected = dense[mask].mean(axis=0)
        np.testing.assert_allclose(result[grp], expected, atol=1e-12, err_msg=f"group {grp}")


def test_sparse_groupby_mean_single_group():
    X = sp.csr_matrix(np.array([[1.0, 2.0], [3.0, 4.0]]))
    labels = ["A", "A"]
    result = sparse_groupby_mean(X, labels)
    np.testing.assert_allclose(result["A"], [2.0, 3.0], atol=1e-12)


# --- sparse_welch_t parity gate (science gate: 1e-9 tolerance) ---

def _welch_t_reference(X_dense, in_mask, out_mask):
    g1 = X_dense[in_mask]
    g2 = X_dense[out_mask]
    n_genes = X_dense.shape[1]
    t_vals = np.zeros(n_genes)
    p_vals = np.zeros(n_genes)
    for j in range(n_genes):
        t, p = stats.ttest_ind(g1[:, j], g2[:, j], equal_var=False)
        t_vals[j] = t
        p_vals[j] = p
    return t_vals, p_vals


@pytest.mark.parametrize("seed", [10, 42, 99])
def test_sparse_welch_t_parity_vs_scipy(seed):
    rng = np.random.default_rng(seed)
    n_cells, n_genes = 80, 25
    X_dense = rng.random((n_cells, n_genes)) * (rng.random((n_cells, n_genes)) < 0.5)
    # Add some signal to first 5 genes
    X_dense[:40, :5] += rng.random((40, 5)) * 2.0

    in_mask = np.array([True] * 40 + [False] * 40)
    out_mask = ~in_mask

    X_sparse = sp.csr_matrix(X_dense)

    t_sp, p_sp, _, _ = sparse_welch_t(X_sparse, in_mask, out_mask)
    t_ref, p_ref = _welch_t_reference(X_dense, in_mask, out_mask)

    np.testing.assert_allclose(t_sp, t_ref, atol=1e-9,
                               err_msg=f"t-stat parity failed (seed={seed})")
    np.testing.assert_allclose(p_sp, p_ref, atol=1e-9,
                               err_msg=f"p-value parity failed (seed={seed})")


def test_sparse_welch_t_float64_internal():
    """Verify float64 accumulators don't lose precision vs scipy on float64 input."""
    rng = np.random.default_rng(7)
    n_cells, n_genes = 40, 10
    # Small variance to stress numerical precision; keep float64 throughout
    X_dense = rng.random((n_cells, n_genes)) * 1e-3
    X_sparse = sp.csr_matrix(X_dense)  # float64 sparse

    in_mask = np.array([True] * 20 + [False] * 20)
    out_mask = ~in_mask

    t_sp, _, _, _ = sparse_welch_t(X_sparse, in_mask, out_mask)
    t_ref, _ = _welch_t_reference(X_dense, in_mask, out_mask)

    np.testing.assert_allclose(t_sp, t_ref, atol=1e-9,
                               err_msg="float64 internal computation required")


def test_sparse_welch_t_dense_input():
    rng = np.random.default_rng(3)
    n_cells, n_genes = 30, 8
    X_dense = rng.random((n_cells, n_genes))
    in_mask = np.array([True] * 15 + [False] * 15)
    out_mask = ~in_mask
    t_sp, p_sp, _, _ = sparse_welch_t(X_dense, in_mask, out_mask)
    t_ref, p_ref = _welch_t_reference(X_dense, in_mask, out_mask)
    np.testing.assert_allclose(t_sp, t_ref, atol=1e-9)
    np.testing.assert_allclose(p_sp, p_ref, atol=1e-9)
