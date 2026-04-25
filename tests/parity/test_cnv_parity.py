"""CNV parity gate: dense default vs SC_CNV_ENGINE=chunked.

Both paths must produce identical smoothed CNV matrix within max abs diff <= 1e-5.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData


def _try_import_chunked():
    try:
        from workflow.modular._sparse_utils import chunked_row_densify  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.fixture(scope="module")
def synthetic_cnv_adata():
    """Small AnnData mimicking gene_pos-aligned expression."""
    rng = np.random.default_rng(11)
    n_cells, n_genes = 400, 250
    base = rng.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(np.float32)
    # Make first 100 cells "malignant" by boosting genes 0-30
    base[:100, :30] *= 3.0
    X_sparse = sp.csr_matrix(base)
    adata = AnnData(X_sparse)
    adata.var_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    adata.obs["cell_type"] = pd.Categorical(
        ["Fibroblast"] * 200 + ["Tumor"] * 200
    )
    return adata


def _build_gene_pos(adata):
    """Synthetic chromosome positions: 5 chroms, ~50 genes each."""
    n_genes = adata.n_vars
    chrom_size = n_genes // 5
    chroms, starts = [], []
    for c in range(5):
        for i in range(chrom_size):
            chroms.append(f"chr{c+1}")
            starts.append(i * 1000)
    while len(chroms) < n_genes:
        chroms.append(f"chr5")
        starts.append((len(chroms) - 4 * chrom_size) * 1000)
    df = pd.DataFrame(
        {"chromosome": chroms, "start": starts},
        index=adata.var_names[: len(chroms)],
    )
    return df


@pytest.mark.skipif(not _try_import_chunked(), reason="phase 2 not enabled")
def test_cnv_chunked_matches_dense(synthetic_cnv_adata, monkeypatch):
    """SC_CNV_ENGINE=chunked smoothed matrix must match default dense within 1e-5."""
    from scipy.ndimage import uniform_filter1d
    from workflow.modular._sparse_utils import chunked_row_densify
    from workflow.modular.modules.cnv_inference import CNVInferenceModule

    adata = synthetic_cnv_adata
    gene_pos = _build_gene_pos(adata).sort_values(["chromosome", "start"])
    chromosomes = gene_pos["chromosome"].values
    window = 100

    expr_sparse = adata[:, gene_pos.index].X

    # Dense reference path (mirrors the module's default branch)
    expr = np.array(expr_sparse.toarray() if hasattr(expr_sparse, "toarray") else expr_sparse,
                    dtype=np.float32)
    ref_mask = (adata.obs["cell_type"] == "Fibroblast").values
    ref_mean = expr[ref_mask].mean(axis=0)
    centered = expr - ref_mean
    smoothed_dense = np.zeros_like(centered)
    for chrom in np.unique(chromosomes):
        chrom_idx = np.where(chromosomes == chrom)[0]
        if len(chrom_idx) < 3:
            smoothed_dense[:, chrom_idx] = centered[:, chrom_idx]
            continue
        smoothed_dense[:, chrom_idx] = uniform_filter1d(
            centered[:, chrom_idx],
            size=min(window, len(chrom_idx)),
            axis=1,
            mode="nearest",
        )

    # Chunked path via the module's static helper
    smoothed_chunked = CNVInferenceModule._compute_smoothed_chunked(
        expr_sparse, adata, "Fibroblast", chromosomes, window, chunk_rows=64,
    )

    max_diff = np.max(np.abs(smoothed_dense - smoothed_chunked))
    assert max_diff <= 1e-5, f"smoothed matrix max abs diff = {max_diff:.2e} (need <= 1e-5)"


@pytest.mark.skipif(not _try_import_chunked(), reason="phase 2 not enabled")
def test_cnv_chunked_chunk_size_invariance(synthetic_cnv_adata):
    """smoothed matrix is chunk-size invariant: chunk=64 vs chunk=257 must produce identical output."""
    from workflow.modular.modules.cnv_inference import CNVInferenceModule

    adata = synthetic_cnv_adata
    gene_pos = _build_gene_pos(adata).sort_values(["chromosome", "start"])
    chromosomes = gene_pos["chromosome"].values
    expr_sparse = adata[:, gene_pos.index].X

    s1 = CNVInferenceModule._compute_smoothed_chunked(
        expr_sparse, adata, "Fibroblast", chromosomes, 100, chunk_rows=64,
    )
    s2 = CNVInferenceModule._compute_smoothed_chunked(
        expr_sparse, adata, "Fibroblast", chromosomes, 100, chunk_rows=257,
    )
    np.testing.assert_allclose(s1, s2, atol=1e-6,
                               err_msg="chunk-size invariance violated")
