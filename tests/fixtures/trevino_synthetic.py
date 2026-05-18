"""Pytest fixtures and helpers for Trevino 2021 synthetic data.

Schema source: Trevino et al. 2021, Cell (doi:10.1016/j.cell.2021.07.039).
No __references__ here — this is not a pipeline module.
"""
from __future__ import annotations

import numpy as np
import pytest


@pytest.fixture(scope="session")
def trevino_synthetic_adata():
    from scripts.dev.load_trevino_2021 import load_synthetic
    return load_synthetic(seed=42)


def make_atac_dense_variant(adata):
    """Return a copy of adata with obsm['atac_peaks'] densified to ndarray."""
    import anndata
    copy = adata.copy()
    copy.obsm["atac_peaks"] = copy.obsm["atac_peaks"].toarray()
    return copy


def make_missing_layer_variant(adata):
    """Return a copy of adata with obsm['atac_peaks'] removed."""
    copy = adata.copy()
    del copy.obsm["atac_peaks"]
    return copy
