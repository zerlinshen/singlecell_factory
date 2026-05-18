"""Tests for Pearson + permutation FDR peak-to-gene linkage (US-W5-6).

Tests use the synthetic Trevino fixture (1000 cells, 3000 genes, 5000 peaks).
Planted linkage structure (from load_synthetic defaults):
  markers_per_type = max(10, 3000 // (3 * 10)) = 100
  peaks_per_type   = max(10, 5000 // (3 * 10)) = 166
  type_0: genes[0:100]   <-> peaks[0:166]
  type_1: genes[100:200] <-> peaks[166:332]
  type_2: genes[200:300] <-> peaks[332:498]
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse

from tests.fixtures.trevino_synthetic import trevino_synthetic_adata, make_atac_dense_variant
from workflow.modular.modules.peak_to_gene import compute_peak_gene_linkages, __references__
from workflow.modular._contract_violation import ModuleContractError


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _planted_pairs_set(n_cell_types: int = 3, n_rna: int = 3000, n_atac: int = 5000):
    """Reconstruct the set of planted (peak_idx, gene_name) pairs from load_synthetic defaults."""
    markers_per_type = max(10, n_rna // (n_cell_types * 10))
    peaks_per_type = max(10, n_atac // (n_cell_types * 10))
    planted = set()
    for ct in range(n_cell_types):
        gene_start = ct * markers_per_type
        gene_end = gene_start + markers_per_type
        peak_start = ct * peaks_per_type
        peak_end = peak_start + peaks_per_type
        for p in range(peak_start, peak_end):
            for g in range(gene_start, gene_end):
                planted.add((p, f"ENSG{g:011d}"))
    return planted


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_linkage_recovery_on_planted(trevino_synthetic_adata):
    """Top-1000 linkages should be dominated by planted (peak, gene) pairs (>=90%)."""
    adata = trevino_synthetic_adata
    compute_peak_gene_linkages(
        adata,
        peak_block_size=500,  # small blocks for speed on 5000 peaks
        gene_window_bp=500_000,
        n_perms=50,
        fdr_alpha=0.05,
        random_state=42,
    )

    top1000 = adata.uns["peak_to_gene_top1000"]
    assert len(top1000) > 0, "top-1000 table is empty"

    planted = _planted_pairs_set()
    n_top = len(top1000)
    n_planted_in_top = sum(
        1 for _, row in top1000.iterrows()
        if (int(row["peak_idx"]), str(row["gene"])) in planted
    )
    recovery = n_planted_in_top / n_top
    assert recovery >= 0.90, (
        f"Planted linkage recovery in top-{n_top} = {recovery:.3f} < 0.90. "
        f"Planted: {n_planted_in_top}/{n_top}. "
        "Check correlation kernel or planting signal strength."
    )


def test_sparse_axis_assertion_dense_X_raises(trevino_synthetic_adata):
    """Dense adata.X must raise ModuleContractError before any computation."""
    import anndata
    adata = trevino_synthetic_adata.copy()
    adata.X = adata.X.toarray()  # densify X
    with pytest.raises(ModuleContractError, match="adata.X must be a scipy.sparse"):
        compute_peak_gene_linkages(adata, n_perms=2)


def test_sparse_axis_assertion_dense_atac_raises(trevino_synthetic_adata):
    """Dense obsm['atac_peaks'] must raise ModuleContractError before any computation."""
    from tests.fixtures.trevino_synthetic import make_atac_dense_variant
    adata = make_atac_dense_variant(trevino_synthetic_adata)
    with pytest.raises(ModuleContractError, match="atac_peaks.*must be a scipy.sparse|atac_peaks.*dense"):
        compute_peak_gene_linkages(adata, n_perms=2)


def test_sparse_preserved_through_chunks(trevino_synthetic_adata):
    """adata.X and obsm['atac_peaks'] must remain sparse before AND after the call."""
    adata = trevino_synthetic_adata.copy()

    assert scipy.sparse.issparse(adata.X), "X not sparse before call"
    assert scipy.sparse.issparse(adata.obsm["atac_peaks"]), "atac_peaks not sparse before call"

    compute_peak_gene_linkages(
        adata,
        peak_block_size=500,
        n_perms=5,
        fdr_alpha=0.05,
        random_state=0,
    )

    assert scipy.sparse.issparse(adata.X), "X was densified during compute_peak_gene_linkages"
    assert scipy.sparse.issparse(adata.obsm["atac_peaks"]), (
        "obsm['atac_peaks'] was densified during compute_peak_gene_linkages"
    )


def test_references_includes_high_if_new_entries():
    """__references__ must include at least one Cell-journal entry (Trevino 2021 or Ma 2020)."""
    cell_entries = [
        key for key, val in __references__.items()
        if val.get("journal") == "Cell"
    ]
    assert len(cell_entries) >= 1, (
        f"Expected at least one Cell-journal entry in __references__; found: {list(__references__.keys())}"
    )
    # Verify at least one of the two required papers is present
    dois = {val.get("doi", "") for val in __references__.values()}
    has_trevino = any("cell.2021.07.039" in d for d in dois)
    has_ma = any("cell.2020.09.056" in d for d in dois)
    assert has_trevino or has_ma, (
        "Neither Trevino 2021 (doi:10.1016/j.cell.2021.07.039) nor Ma 2020 "
        "(doi:10.1016/j.cell.2020.09.056) found in __references__."
    )
