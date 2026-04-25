"""DE parity gate: dense fallback path vs sparse Welch-t path (Phase 2).

Phase 1: dense-path test runs unconditionally.
         sparse-path test skips until Phase 2 skeleton is importable.

Thresholds (from master plan):
  - top-50 gene overlap >= 48/50
  - |logFC| Pearson >= 0.999
  - p-value Spearman >= 0.99
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from scipy.stats import pearsonr, spearmanr

from workflow.modular.modules.differential_expression import DifferentialExpressionModule


# ---------------------------------------------------------------------------
# Synthetic data fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def synthetic_adata():
    """5 000 cells × 2 000 genes with planted DE signal."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 5_000, 2_000

    # Two balanced groups
    group = np.array(["A"] * (n_cells // 2) + ["B"] * (n_cells // 2))

    # Plant signal: first 50 genes are truly DE (up in group A)
    base = rng.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(np.float32)
    signal_genes = 50
    base[: n_cells // 2, :signal_genes] *= 4.0  # up-regulated in group A

    X_dense = base.copy()
    X_sparse = sp.csr_matrix(base)

    adata = AnnData(X_sparse)
    adata.obs["leiden"] = pd.Categorical(group)
    adata.var_names = [f"gene_{i:04d}" for i in range(n_genes)]

    adata.uns["_dense_X"] = X_dense  # stashed for round-trip checks
    return adata


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _run_dense_fallback(adata: AnnData, n_genes: int = 50) -> pd.DataFrame:
    return DifferentialExpressionModule._fallback_rank_genes_groups_df(
        adata=adata,
        groupby="leiden",
        n_genes=n_genes,
    )


def _top50_names(df: pd.DataFrame, group: str) -> set[str]:
    sub = df[df["group"] == group].head(50)
    return set(sub["names"].tolist())


# ---------------------------------------------------------------------------
# Phase 1: dense-path correctness
# ---------------------------------------------------------------------------

def test_dense_fallback_runs(synthetic_adata):
    df = _run_dense_fallback(synthetic_adata)
    assert not df.empty
    assert set(df.columns) >= {"group", "names", "scores", "pvals_adj", "logfoldchanges"}


def test_dense_fallback_top50_contains_signal_genes(synthetic_adata):
    """Top-50 genes for group A must include >= 48 of the 50 planted signal genes."""
    df = _run_dense_fallback(synthetic_adata)
    top = _top50_names(df, "A")
    signal = {f"gene_{i:04d}" for i in range(50)}
    overlap = len(top & signal)
    assert overlap >= 48, f"Dense top-50 overlap with signal genes = {overlap}/50 (need >= 48)"


def test_dense_fallback_logfc_reasonable(synthetic_adata):
    df = _run_dense_fallback(synthetic_adata)
    a_rows = df[df["group"] == "A"]
    # Signal genes should have positive logFC
    signal = {f"gene_{i:04d}" for i in range(50)}
    sig_rows = a_rows[a_rows["names"].isin(signal)]
    assert (sig_rows["logfoldchanges"] > 0).mean() >= 0.9, "Most signal genes should be up in A"


# ---------------------------------------------------------------------------
# Phase 2: sparse-path parity (skipped until skeleton lands)
# ---------------------------------------------------------------------------

def _try_import_sparse_welch():
    try:
        from workflow.modular._sparse_utils import sparse_welch_t  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(not _try_import_sparse_welch(), reason="phase 2 not enabled")
def test_sparse_vs_dense_top50_overlap(synthetic_adata):
    """Sparse path top-50 must overlap dense path top-50 by >= 48/50 genes."""
    from workflow.modular._sparse_utils import sparse_welch_t

    adata = synthetic_adata
    dense_df = _run_dense_fallback(adata, n_genes=50)

    # Run sparse path for group A
    groups = pd.Categorical(adata.obs["leiden"].astype(str))
    in_mask = np.asarray(groups == "A")
    out_mask = ~in_mask
    X = adata.X  # CSR sparse

    t_stat, pvals = sparse_welch_t(X, in_mask, out_mask)
    mean_in = np.asarray(X[in_mask].mean(axis=0)).ravel() + 1e-9
    mean_out = np.asarray(X[out_mask].mean(axis=0)).ravel() + 1e-9
    logfc = np.log2(mean_in / mean_out)

    gene_names = np.asarray(adata.var_names)
    order = np.lexsort((-np.abs(logfc), pvals))[:50]
    sparse_top50 = set(gene_names[order].tolist())

    dense_top50 = _top50_names(dense_df, "A")
    overlap = len(sparse_top50 & dense_top50)
    assert overlap >= 48, f"Sparse/dense top-50 overlap = {overlap}/50 (need >= 48)"


@pytest.mark.skipif(not _try_import_sparse_welch(), reason="phase 2 not enabled")
def test_sparse_vs_dense_logfc_pearson(synthetic_adata):
    """|logFC| Pearson correlation between sparse and dense paths >= 0.999."""
    from workflow.modular._sparse_utils import sparse_welch_t

    adata = synthetic_adata
    dense_df = _run_dense_fallback(adata, n_genes=adata.n_vars)

    groups = pd.Categorical(adata.obs["leiden"].astype(str))
    in_mask = np.asarray(groups == "A")
    out_mask = ~in_mask
    X = adata.X

    _, pvals = sparse_welch_t(X, in_mask, out_mask)
    mean_in = np.asarray(X[in_mask].mean(axis=0)).ravel() + 1e-9
    mean_out = np.asarray(X[out_mask].mean(axis=0)).ravel() + 1e-9
    sparse_logfc = np.log2(mean_in / mean_out)

    gene_names = list(adata.var_names)
    dense_a = dense_df[dense_df["group"] == "A"].set_index("names")
    common = [g for g in gene_names if g in dense_a.index]
    dense_logfc = dense_a.loc[common, "logfoldchanges"].values
    sparse_logfc_aligned = np.array([sparse_logfc[gene_names.index(g)] for g in common])

    r, _ = pearsonr(np.abs(dense_logfc), np.abs(sparse_logfc_aligned))
    assert r >= 0.999, f"|logFC| Pearson = {r:.6f} (need >= 0.999)"


@pytest.mark.skipif(not _try_import_sparse_welch(), reason="phase 2 not enabled")
def test_sparse_vs_dense_pval_spearman(synthetic_adata):
    """p-value rank Spearman correlation between sparse and dense paths >= 0.99.

    Both sides use BH-adjusted p-values so the comparison is apples-to-apples.
    """
    from workflow.modular._sparse_utils import sparse_welch_t
    from workflow.modular.modules.differential_expression import DifferentialExpressionModule

    adata = synthetic_adata
    dense_df = _run_dense_fallback(adata, n_genes=adata.n_vars)

    groups = pd.Categorical(adata.obs["leiden"].astype(str))
    in_mask = np.asarray(groups == "A")
    out_mask = ~in_mask
    X = adata.X

    _, sparse_raw_pvals = sparse_welch_t(X, in_mask, out_mask)
    # Apply BH correction to match the dense fallback's pvals_adj
    sparse_pvals_adj = DifferentialExpressionModule._benjamini_hochberg(sparse_raw_pvals)

    gene_names = list(adata.var_names)
    dense_a = dense_df[dense_df["group"] == "A"].set_index("names")
    common = [g for g in gene_names if g in dense_a.index]
    dense_pv = dense_a.loc[common, "pvals_adj"].values
    sparse_pv = np.array([sparse_pvals_adj[gene_names.index(g)] for g in common])

    rho, _ = spearmanr(dense_pv, sparse_pv)
    assert rho >= 0.99, f"p-value Spearman = {rho:.6f} (need >= 0.99)"
