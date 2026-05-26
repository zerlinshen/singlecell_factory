"""Synthetic fixtures for the discovery integration-selection gate (plan §6/COND #6).

THREE fixtures with the construction parameters mandated by COND #6. Each
returns embeddings whose correct discovery verdict is known a-priori, so the
acceptance suite (tests/test_discovery_integration_gate.py) can assert the
metrics fire in the right direction.

Mirrors the parametric style of ``tests/fixtures/integration_synthetic.py`` (a
frozen dataclass + a ``make_*`` builder seeded for determinism). This is a CI
fixture, NOT pipeline data — it carries no ``__references__``.

Geometry convention
--------------------
Every fixture lives in a PCA-like embedding space. Dimensions are split into:
  - "bio" dims    : well-separated Gaussian centers per biological population.
  - a "batch axis" : a dedicated nuisance dim carrying a per-batch offset.
The BASELINE embedding (``X_baseline``) keeps the batch offset (under-corrected:
one population may split into per-batch clusters in this shared X_pca). A GOOD
integration removes the batch offset on the batch axis while preserving the bio
geometry. Transforms ("erase", "fuse", "over_merge") produce the BAD-merge and
rare-erasure embeddings the gates must catch.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class DiscoveryFixture:
    """Container for a discovery-gate synthetic fixture.

    Attributes
    ----------
    embeddings:
        ``{method_name -> (n_cells, n_dims) float embedding}``. ALWAYS contains
        ``"baseline"`` (the shared pre-integration X_pca) plus the named
        candidate / control embeddings the fixture defines.
    batch:
        (n_cells,) string batch labels.
    true_bio:
        (n_cells,) string GROUND-TRUTH biological population labels. For
        ASSERTIONS ONLY — the gate itself is label-agnostic and never consumes
        these.
    meta:
        free-form provenance (construction params, known answers).
    """

    embeddings: dict[str, np.ndarray]
    batch: np.ndarray
    true_bio: np.ndarray
    meta: dict = field(default_factory=dict)


def _blob(rng, n, center, sd):
    out = rng.normal(loc=0.0, scale=sd, size=(n, len(center)))
    return out + np.asarray(center, dtype=np.float64)


# ==========================================================================
# F1 — batch-split-of-one-type
# ==========================================================================

def make_f1_batch_split(
    *,
    seed: int = 0,
    split_type_cells: int = 200,         # per batch -> ~400 total (COND #6)
    other_type_cells: int = 120,         # per batch, per other type
    n_other_types: int = 2,              # >=2 clearly-separated types (COND #6)
    n_bio_dims: int = 6,
    bio_separation: float = 15.0,
    batch_axis_offset: float = 10.0,     # offset along the dedicated batch axis
    within_sd: float = 1.0,
) -> DiscoveryFixture:
    """ONE biological type whose per-batch copies are offset along a dedicated
    BATCH AXIS (splits into 2 per-batch clusters in X_pca) but is ONE
    population; 2 batches; plus ``n_other_types`` clearly-separated types so
    good-merge != collapse-all.

    Embeddings:
      - ``baseline``       : split type offset on the batch axis (under-corrected).
      - ``harmony``        : GOOD integration — batch-axis offset REMOVED for the
                             split type so its two per-batch copies co-localize;
                             other types untouched and distinct.
      - ``neg_control_shuffle_label`` : row-permuted baseline (structure destroyed).
    """
    rng = np.random.default_rng(seed)
    n_dims = n_bio_dims + 1                       # last dim = batch axis
    batch_axis = n_dims - 1

    # Bio centers: split type at index 0, others well separated.
    bio_centers = rng.normal(size=(1 + n_other_types, n_bio_dims)) * bio_separation

    rows_X: list[np.ndarray] = []
    rows_batch: list[str] = []
    rows_bio: list[str] = []

    for b in (0, 1):
        # Split type (ONE population) — present in both batches, offset on the
        # batch axis so it forms TWO per-batch clusters in the baseline X_pca.
        blob = np.zeros((split_type_cells, n_dims), dtype=np.float64)
        blob[:, :n_bio_dims] = _blob(rng, split_type_cells, bio_centers[0], within_sd)
        blob[:, batch_axis] = rng.normal(0, within_sd, split_type_cells) + (
            batch_axis_offset if b == 1 else -batch_axis_offset)
        rows_X.append(blob)
        rows_batch.extend([f"batch{b}"] * split_type_cells)
        rows_bio.extend(["split_type"] * split_type_cells)

        # Other distinct types — also carry the batch offset on the batch axis
        # (so the batch effect is global, not type-specific) but they are NOT
        # the same population (distinct bio centers).
        for t in range(n_other_types):
            ob = np.zeros((other_type_cells, n_dims), dtype=np.float64)
            ob[:, :n_bio_dims] = _blob(rng, other_type_cells, bio_centers[1 + t], within_sd)
            ob[:, batch_axis] = rng.normal(0, within_sd, other_type_cells) + (
                batch_axis_offset if b == 1 else -batch_axis_offset)
            rows_X.append(ob)
            rows_batch.extend([f"batch{b}"] * other_type_cells)
            rows_bio.extend([f"other{t}"] * other_type_cells)

    X_baseline = np.vstack(rows_X)
    batch = np.asarray(rows_batch, dtype=object)
    true_bio = np.asarray(rows_bio, dtype=object)

    # GOOD integration: remove the batch-axis offset everywhere -> the split
    # type's two per-batch copies now co-localize (good merge); other types
    # stay distinct in the bio dims.
    X_harmony = X_baseline.copy()
    X_harmony[:, batch_axis] = rng.normal(0, within_sd, X_harmony.shape[0])

    # Shuffle control: row-permute the baseline (destroys structure).
    perm = np.random.default_rng(seed + 99).permutation(X_baseline.shape[0])
    X_shuffle = X_baseline[perm]

    return DiscoveryFixture(
        embeddings={
            "baseline": X_baseline,
            "harmony": X_harmony,
            "neg_control_shuffle_label": X_shuffle,
        },
        batch=batch,
        true_bio=true_bio,
        meta={
            "fixture": "F1_batch_split_of_one_type",
            "split_type_total_cells": 2 * split_type_cells,
            "n_other_types": n_other_types,
            "batch_axis": batch_axis,
            "known": "split_type forms 2 per-batch clusters in baseline; GOOD "
                     "harmony merges them; other types remain distinct.",
        },
    )


# ==========================================================================
# F2 — rare/novel population
# ==========================================================================

def make_f2_rare_novel(
    *,
    seed: int = 1,
    rare_cells_per_batch: int = 18,      # ~36 total; ~3-5% of batch (COND #6)
    abundant_cells_per_batch: int = 200,
    n_abundant_types: int = 2,
    n_bio_dims: int = 6,
    bio_separation: float = 15.0,
    rare_separation: float = 22.0,       # rare type well-separated on a bio dim
    batch_axis_offset: float = 8.0,
    within_sd: float = 1.0,
) -> DiscoveryFixture:
    """A small rare type (~40 cells, ~3-5% of batch) well-separated on a
    dedicated bio dim, present in >=2 batches, alongside abundant types.

    Embeddings:
      - ``baseline``  : rare type well-isolated (with batch-axis offset).
      - ``harmony``   : GOOD — batch removed, rare type STILL isolated (D1c high).
      - ``erase``     : the rare type's cells are MOVED onto a neighbor abundant
                        type (fused) -> D1c collapses (rare type erased).
      - ``neg_control_shuffle_label`` : row-permuted baseline.
    """
    rng = np.random.default_rng(seed)
    n_dims = n_bio_dims + 1
    batch_axis = n_dims - 1
    rare_dim = 0  # dedicated bio dim for the rare type's separation

    abundant_centers = rng.normal(size=(n_abundant_types, n_bio_dims)) * bio_separation
    rare_center = abundant_centers[0].copy()
    rare_center[rare_dim] += rare_separation  # push rare type far on rare_dim

    rows_X: list[np.ndarray] = []
    rows_batch: list[str] = []
    rows_bio: list[str] = []

    for b in (0, 1):
        for t in range(n_abundant_types):
            ob = np.zeros((abundant_cells_per_batch, n_dims), dtype=np.float64)
            ob[:, :n_bio_dims] = _blob(rng, abundant_cells_per_batch, abundant_centers[t], within_sd)
            ob[:, batch_axis] = rng.normal(0, within_sd, abundant_cells_per_batch) + (
                batch_axis_offset if b == 1 else -batch_axis_offset)
            rows_X.append(ob)
            rows_batch.extend([f"batch{b}"] * abundant_cells_per_batch)
            rows_bio.extend([f"abundant{t}"] * abundant_cells_per_batch)
        # rare type
        rb = np.zeros((rare_cells_per_batch, n_dims), dtype=np.float64)
        rb[:, :n_bio_dims] = _blob(rng, rare_cells_per_batch, rare_center, within_sd)
        rb[:, batch_axis] = rng.normal(0, within_sd, rare_cells_per_batch) + (
            batch_axis_offset if b == 1 else -batch_axis_offset)
        rows_X.append(rb)
        rows_batch.extend([f"batch{b}"] * rare_cells_per_batch)
        rows_bio.extend(["rare_novel"] * rare_cells_per_batch)

    X_baseline = np.vstack(rows_X)
    batch = np.asarray(rows_batch, dtype=object)
    true_bio = np.asarray(rows_bio, dtype=object)

    X_harmony = X_baseline.copy()
    X_harmony[:, batch_axis] = rng.normal(0, within_sd, X_harmony.shape[0])

    # Erase transform: fuse the rare type into abundant0 by PLACING each rare
    # cell at an existing abundant0 cell's exact location (interior points). The
    # rare cells become indistinguishable INTERIOR members of abundant0, so the
    # rare frozen per-batch label gets a strongly NEGATIVE silhouette in E (its
    # cells are closer to abundant0 cells of a different label than to each
    # other) -> D1c collapses below the neg-control floor (TM-2 / TM-3).
    X_erase = X_harmony.copy()
    rare_mask = true_bio == "rare_novel"
    n_rare = int(rare_mask.sum())
    erase_rng = np.random.default_rng(seed + 7)
    abundant0_mask = true_bio == "abundant0"
    donor = X_harmony[abundant0_mask][:, :n_bio_dims]
    donor_idx = erase_rng.integers(0, donor.shape[0], size=n_rare)
    X_erase[rare_mask, :n_bio_dims] = donor[donor_idx]

    perm = np.random.default_rng(seed + 99).permutation(X_baseline.shape[0])
    X_shuffle = X_baseline[perm]

    return DiscoveryFixture(
        embeddings={
            "baseline": X_baseline,
            "harmony": X_harmony,
            "erase": X_erase,
            "neg_control_shuffle_label": X_shuffle,
        },
        batch=batch,
        true_bio=true_bio,
        meta={
            "fixture": "F2_rare_novel",
            "rare_total_cells": 2 * rare_cells_per_batch,
            "rare_fraction_approx": rare_cells_per_batch /
            (rare_cells_per_batch + n_abundant_types * abundant_cells_per_batch),
            "known": "harmony keeps rare type isolated (D1c high); erase fuses "
                     "it into abundant0 (D1c collapses).",
        },
    )


# ==========================================================================
# F3 — single-batch / sample-specific
# ==========================================================================

def make_f3_single_batch(
    *,
    seed: int = 2,
    single_batch_cells: int = 60,        # present in EXACTLY one batch (COND #6)
    shared_cells_per_batch: int = 160,
    n_shared_types: int = 2,
    n_bio_dims: int = 6,
    bio_separation: float = 15.0,
    single_separation: float = 22.0,
    batch_axis_offset: float = 8.0,
    within_sd: float = 1.0,
) -> DiscoveryFixture:
    """A distinct type (~60 cells) present in EXACTLY ONE batch (COND #5),
    well-separated; plus shared types in both batches.

    Used to assert: the single-batch type is EXCLUDED from D1a's denominator
    (no mutual cross-batch match) and is held isolated by D1c.

    Embeddings:
      - ``baseline`` / ``harmony`` as in the other fixtures.
      - ``neg_control_shuffle_label`` : row-permuted baseline.
    """
    rng = np.random.default_rng(seed)
    n_dims = n_bio_dims + 1
    batch_axis = n_dims - 1

    shared_centers = rng.normal(size=(n_shared_types, n_bio_dims)) * bio_separation
    single_center = shared_centers[0].copy()
    single_center[1] += single_separation  # separate on a dedicated bio dim

    rows_X: list[np.ndarray] = []
    rows_batch: list[str] = []
    rows_bio: list[str] = []

    for b in (0, 1):
        for t in range(n_shared_types):
            ob = np.zeros((shared_cells_per_batch, n_dims), dtype=np.float64)
            ob[:, :n_bio_dims] = _blob(rng, shared_cells_per_batch, shared_centers[t], within_sd)
            ob[:, batch_axis] = rng.normal(0, within_sd, shared_cells_per_batch) + (
                batch_axis_offset if b == 1 else -batch_axis_offset)
            rows_X.append(ob)
            rows_batch.extend([f"batch{b}"] * shared_cells_per_batch)
            rows_bio.extend([f"shared{t}"] * shared_cells_per_batch)
        # single-batch-only type — ONLY in batch0
        if b == 0:
            sb = np.zeros((single_batch_cells, n_dims), dtype=np.float64)
            sb[:, :n_bio_dims] = _blob(rng, single_batch_cells, single_center, within_sd)
            sb[:, batch_axis] = rng.normal(0, within_sd, single_batch_cells) - batch_axis_offset
            rows_X.append(sb)
            rows_batch.extend([f"batch{b}"] * single_batch_cells)
            rows_bio.extend(["single_batch_only"] * single_batch_cells)

    X_baseline = np.vstack(rows_X)
    batch = np.asarray(rows_batch, dtype=object)
    true_bio = np.asarray(rows_bio, dtype=object)

    X_harmony = X_baseline.copy()
    X_harmony[:, batch_axis] = rng.normal(0, within_sd, X_harmony.shape[0])

    perm = np.random.default_rng(seed + 99).permutation(X_baseline.shape[0])
    X_shuffle = X_baseline[perm]

    return DiscoveryFixture(
        embeddings={
            "baseline": X_baseline,
            "harmony": X_harmony,
            "neg_control_shuffle_label": X_shuffle,
        },
        batch=batch,
        true_bio=true_bio,
        meta={
            "fixture": "F3_single_batch",
            "single_batch_cells": single_batch_cells,
            "single_batch_label": "single_batch_only",
            "present_in_batch": "batch0",
            "known": "single_batch_only has NO cross-batch match -> excluded "
                     "from D1a denominator; held isolated by D1c.",
        },
    )
