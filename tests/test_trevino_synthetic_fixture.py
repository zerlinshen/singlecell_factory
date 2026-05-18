"""Tests for tests/fixtures/trevino_synthetic.py."""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tests.fixtures.trevino_synthetic import (
    make_atac_dense_variant,
    make_missing_layer_variant,
    trevino_synthetic_adata,
)
from scripts.dev.load_trevino_2021 import load_synthetic


@pytest.fixture(scope="session")
def synth_adata():
    return load_synthetic(seed=42)


def test_round_trip_shapes(synth_adata):
    assert synth_adata.X.shape == (1000, 3000)
    assert synth_adata.obsm["atac_peaks"].shape == (1000, 5000)
    assert isinstance(synth_adata.obsm["atac_peaks"], scipy.sparse.csr_matrix)


def test_cell_type_count(synth_adata):
    assert synth_adata.obs["cell_type_synth"].nunique() >= 2


def test_dense_variant(synth_adata):
    dense = make_atac_dense_variant(synth_adata)
    assert isinstance(dense.obsm["atac_peaks"], np.ndarray)


def test_missing_layer_variant(synth_adata):
    missing = make_missing_layer_variant(synth_adata)
    with pytest.raises(KeyError):
        _ = missing.obsm["atac_peaks"]
