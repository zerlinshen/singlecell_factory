# Sparse math utilities — dormant until Phase 2 imports this module.
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.stats as stats


class MemoryGuardError(MemoryError):
    pass


def estimate_dense_bytes(shape: tuple[int, int], dtype) -> int:
    return int(np.prod(shape)) * np.dtype(dtype).itemsize


def safe_densify(x, *, reason: str, dtype, allow_chunk: bool = False):
    if sp.issparse(x):
        nbytes = estimate_dense_bytes(x.shape, dtype)
        from workflow.modular._mem_guard import MemoryGuard
        decision = MemoryGuard._quick_check(nbytes)
        if decision == "abort":
            raise MemoryGuardError(f"safe_densify blocked ({reason}): would need {nbytes:,} bytes")
        if decision == "chunk" and not allow_chunk:
            raise MemoryGuardError(f"safe_densify requires chunking ({reason}) but allow_chunk=False")
        return np.asarray(x.toarray(), dtype=dtype)
    return np.asarray(x, dtype=dtype)


def chunked_row_densify(x, chunk_rows: int, dtype):
    if not sp.issparse(x):
        for i in range(0, x.shape[0], chunk_rows):
            yield np.asarray(x[i:i + chunk_rows], dtype=dtype)
        return
    for i in range(0, x.shape[0], chunk_rows):
        yield np.asarray(x[i:i + chunk_rows].toarray(), dtype=dtype)


def sparse_groupby_mean(X, group_labels):
    """Compute per-group column means via sparse indicator matrix multiplication."""
    labels = np.asarray(group_labels)
    groups = np.unique(labels)
    if sp.issparse(X):
        X_csr = X.tocsr()
    else:
        X_csr = sp.csr_matrix(X)
    n = X_csr.shape[0]
    result = {}
    for g in groups:
        mask = (labels == g).astype(np.float64)
        count = mask.sum()
        if count == 0:
            result[g] = np.zeros(X_csr.shape[1], dtype=np.float64)
        else:
            result[g] = np.asarray(mask @ X_csr, dtype=np.float64).ravel() / count
    return result


def sparse_welch_t(X, in_mask, out_mask):
    """Closed-form Welch t-statistic computed in float64.

    Returns (t_stats, p_values, mean_in, mean_out) arrays of shape (n_genes,).
    All intermediate accumulators are float64 to ensure numerical stability.
    """
    in_mask = np.asarray(in_mask, dtype=bool)
    out_mask = np.asarray(out_mask, dtype=bool)

    if sp.issparse(X):
        X_csr = X.tocsr()
        X_in = X_csr[in_mask]
        X_out = X_csr[out_mask]
    else:
        X_in = X[in_mask]
        X_out = X[out_mask]

    n1 = int(in_mask.sum())
    n2 = int(out_mask.sum())

    # Compute means in float64
    if sp.issparse(X_in):
        mean1 = np.asarray(X_in.mean(axis=0), dtype=np.float64).ravel()
    else:
        mean1 = np.mean(X_in, axis=0, dtype=np.float64).ravel()

    if sp.issparse(X_out):
        mean2 = np.asarray(X_out.mean(axis=0), dtype=np.float64).ravel()
    else:
        mean2 = np.mean(X_out, axis=0, dtype=np.float64).ravel()

    # Compute variances in float64 via E[X^2] - (E[X])^2
    if sp.issparse(X_in):
        X_in_f64 = X_in.astype(np.float64)
        mean_sq1 = np.asarray(X_in_f64.multiply(X_in_f64).mean(axis=0), dtype=np.float64).ravel()
    else:
        X_in_f64 = np.asarray(X_in, dtype=np.float64)
        mean_sq1 = np.mean(X_in_f64 ** 2, axis=0).ravel()
    var1 = mean_sq1 - mean1 ** 2
    # Bessel correction: n/(n-1)
    var1 = var1 * n1 / max(n1 - 1, 1)

    if sp.issparse(X_out):
        X_out_f64 = X_out.astype(np.float64)
        mean_sq2 = np.asarray(X_out_f64.multiply(X_out_f64).mean(axis=0), dtype=np.float64).ravel()
    else:
        X_out_f64 = np.asarray(X_out, dtype=np.float64)
        mean_sq2 = np.mean(X_out_f64 ** 2, axis=0).ravel()
    var2 = mean_sq2 - mean2 ** 2
    var2 = var2 * n2 / max(n2 - 1, 1)

    # Welch t-statistic
    denom = np.sqrt(np.maximum(var1 / n1 + var2 / n2, 0.0))
    t_stats = np.where(denom > 0, (mean1 - mean2) / denom, 0.0)

    # Welch-Satterthwaite degrees of freedom
    v1 = var1 / n1
    v2 = var2 / n2
    num_df = (v1 + v2) ** 2
    den_df = v1 ** 2 / max(n1 - 1, 1) + v2 ** 2 / max(n2 - 1, 1)
    df = np.where(den_df > 0, num_df / den_df, 1.0)

    p_values = 2.0 * stats.t.sf(np.abs(t_stats), df=df)
    return t_stats, p_values, mean1, mean2
