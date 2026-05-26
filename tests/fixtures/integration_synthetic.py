"""Synthetic 2-batch integration CI fixture.

Generates a small embedding-space dataset for the integration calibration
harness (scripts/bench/integration/score_integration.py). The geometry is
constructed so the *correct* answer is known a-priori, which lets CI assert
that scib-metrics and the native sklearn cross-check agree:

  - K shared biological cell types (same types present in BOTH batches), each
    a well-separated Gaussian blob along the "bio" axes -> bio conservation
    metrics (cLISI / celltype-ASW / ARI / NMI) should be HIGH.
  - A batch-only nuisance axis: a constant additive offset applied to every
    cell of batch 1 along dedicated nuisance dimensions, with NO correlation
    to cell type. A perfect integration removes this offset; an unintegrated
    embedding retains it.

The fixture therefore yields two canonical embeddings:
  - ``X_unintegrated`` : bio blobs + the batch nuisance offset (batches
    separable -> LOW batch mixing).
  - ``X_integrated``   : the same bio blobs with the nuisance offset removed
    (batches mixed -> HIGH batch mixing) while bio structure is preserved.

This is NOT pipeline data and carries no ``__references__`` — it is a CI
fixture, not a pipeline module.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class IntegrationFixture:
    """Container for the synthetic 2-batch fixture.

    Attributes
    ----------
    X_unintegrated, X_integrated:
        (n_cells, n_dims) float embeddings. ``X_integrated`` has the
        batch-only nuisance offset removed.
    batch:
        (n_cells,) string batch labels ("batch0" / "batch1").
    cell_type:
        (n_cells,) string biological cell-type labels (ground truth, shared
        across both batches).
    n_cell_types:
        number of distinct shared cell types.
    """

    X_unintegrated: np.ndarray
    X_integrated: np.ndarray
    batch: np.ndarray
    cell_type: np.ndarray
    n_cell_types: int


def make_integration_fixture(
    *,
    seed: int = 0,
    n_cell_types: int = 4,
    cells_per_type_per_batch: int = 60,
    n_bio_dims: int = 6,
    n_nuisance_dims: int = 2,
    bio_separation: float = 8.0,
    batch_offset: float = 12.0,
    within_blob_sd: float = 1.0,
) -> IntegrationFixture:
    """Build the synthetic 2-batch fixture.

    The bio signal lives in the first ``n_bio_dims`` columns (one blob center
    per cell type, separated by ``bio_separation``). The batch nuisance signal
    is a constant offset (``batch_offset``) applied to batch-1 cells along the
    final ``n_nuisance_dims`` columns; it is uncorrelated with cell type, so
    removing it (``X_integrated``) does not touch the bio geometry.
    """
    rng = np.random.default_rng(seed)
    n_dims = n_bio_dims + n_nuisance_dims

    # One random, well-separated center per cell type in the bio subspace.
    bio_centers = rng.normal(size=(n_cell_types, n_bio_dims)) * bio_separation

    rows_X: list[np.ndarray] = []
    rows_batch: list[str] = []
    rows_ct: list[str] = []

    for batch_idx in (0, 1):
        for ct in range(n_cell_types):
            blob = rng.normal(
                loc=0.0, scale=within_blob_sd,
                size=(cells_per_type_per_batch, n_dims),
            )
            # Place bio blob.
            blob[:, :n_bio_dims] += bio_centers[ct]
            # Apply batch-only nuisance offset to batch 1 along nuisance dims.
            if batch_idx == 1:
                blob[:, n_bio_dims:] += batch_offset
            rows_X.append(blob)
            rows_batch.extend([f"batch{batch_idx}"] * cells_per_type_per_batch)
            rows_ct.extend([f"celltype{ct}"] * cells_per_type_per_batch)

    X_unintegrated = np.vstack(rows_X).astype(np.float64)
    batch = np.asarray(rows_batch, dtype=object)
    cell_type = np.asarray(rows_ct, dtype=object)

    # Integration target: remove the batch-only offset. Bio geometry untouched.
    X_integrated = X_unintegrated.copy()
    is_batch1 = batch == "batch1"
    X_integrated[is_batch1, n_bio_dims:] -= batch_offset

    return IntegrationFixture(
        X_unintegrated=X_unintegrated,
        X_integrated=X_integrated,
        batch=batch,
        cell_type=cell_type,
        n_cell_types=n_cell_types,
    )
