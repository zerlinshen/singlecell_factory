"""Tests for scripts/dev/load_trevino_2021.py — synthetic path only."""
from __future__ import annotations

import numpy as np
import scipy.sparse

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.dev.load_trevino_2021 import load_synthetic


def test_load_synthetic_shapes():
    adata = load_synthetic(n_obs=100, n_rna=200, n_atac=1000, seed=0)
    assert adata.X.shape == (100, 200)
    assert adata.obsm["atac_peaks"].shape == (100, 1000)
    assert isinstance(adata.X, scipy.sparse.csr_matrix)
    assert adata.X.dtype == np.int32
    assert adata.obsm["atac_peaks"].dtype == np.int32


def test_load_synthetic_deterministic():
    a1 = load_synthetic(seed=7)
    a2 = load_synthetic(seed=7)
    diff = (a1.X - a2.X)
    assert diff.nnz == 0
    diff_atac = (a1.obsm["atac_peaks"] - a2.obsm["atac_peaks"])
    assert diff_atac.nnz == 0


def test_load_synthetic_cell_types():
    adata = load_synthetic(n_cell_types=3, seed=42)
    unique = adata.obs["cell_type_synth"].nunique()
    assert unique >= 2


def test_load_synthetic_csr_sparse():
    adata = load_synthetic(seed=42)
    assert isinstance(adata.obsm["atac_peaks"], scipy.sparse.csr_matrix)
