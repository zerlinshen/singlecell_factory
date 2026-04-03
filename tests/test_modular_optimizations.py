import numpy as np
import pandas as pd
from anndata import AnnData

from workflow.modular.modules.pseudo_velocity import compute_pseudo_velocity_vectorized
from workflow.modular.perf_baseline import baseline_pseudo_velocity


def test_vectorized_pseudo_velocity_matches_baseline():
    rng = np.random.default_rng(42)
    n = 120
    umap = rng.normal(size=(n, 2))
    pseudotime = np.sort(rng.uniform(0, 1, size=n))
    adata = AnnData(np.ones((n, 5)))
    adata.obsm["X_umap"] = umap
    adata.obs["dpt_pseudotime"] = pseudotime

    baseline = baseline_pseudo_velocity(adata, n_neighbors=15)
    optimized = compute_pseudo_velocity_vectorized(umap, pseudotime, n_neighbors=15)
    assert np.allclose(baseline, optimized, atol=1e-10, rtol=1e-7)


def test_annotation_second_score_vectorized():
    """Verify np.partition gives same results as pandas nlargest for second-best score."""
    rng = np.random.default_rng(99)
    vals = rng.standard_normal((100, 5))
    score_mat = pd.DataFrame(vals, columns=["A", "B", "C", "D", "E"])

    # Vectorized (new approach)
    second_vec = np.partition(vals, -2, axis=1)[:, -2]
    # Reference (old approach)
    second_ref = score_mat.apply(lambda row: row.nlargest(2).iloc[-1], axis=1).values

    assert np.allclose(second_vec, second_ref)


def test_bh_fdr_correction():
    """Verify BH FDR correction produces valid adjusted p-values."""
    pvals = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]
    ranks = np.arange(1, n + 1)
    raw_adj = sorted_pvals * n / ranks
    padj_sorted = np.minimum.accumulate(raw_adj[::-1])[::-1]
    padj_sorted = np.clip(padj_sorted, 0, 1.0)

    # Verify monotonicity (adjusted p-values should be non-decreasing)
    assert all(padj_sorted[i] <= padj_sorted[i + 1] + 1e-15 for i in range(len(padj_sorted) - 1))
    # Adjusted p-values should be >= raw p-values
    assert all(padj_sorted[i] >= sorted_pvals[i] - 1e-15 for i in range(n))
    # All should be in [0, 1]
    assert all(0 <= p <= 1 for p in padj_sorted)


def test_cnv_vectorized_smoothing():
    """Verify scipy uniform_filter1d produces correct sliding window means."""
    from scipy.ndimage import uniform_filter1d

    rng = np.random.default_rng(42)
    data = rng.standard_normal((10, 50))
    window = 5

    # Vectorized
    result = uniform_filter1d(data, size=window, axis=1, mode="nearest")

    # Manual reference for center positions
    for col in range(window // 2, data.shape[1] - window // 2):
        expected = data[:, col - window // 2 : col + window // 2 + 1].mean(axis=1)
        assert np.allclose(result[:, col], expected, atol=1e-10)
