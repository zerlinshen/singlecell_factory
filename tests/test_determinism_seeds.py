"""Determinism / seed-wiring tests for composition and evolution modules.

Covers two reproducibility/perf fixes:

* #20 — ``composition._fallback_test`` and ``evolution._cluster_clones`` used to
  hardcode the seed (``np.random.default_rng(42)`` / ``np.random.RandomState(42)``),
  making their stochastic output INSENSITIVE to ``--random-state``. They now thread
  the canonical ``ctx.random_state``. We assert: same seed -> identical output, two
  different seeds -> different output (proving the literal is gone and the seed is
  wired), at both the unit level and through ``Module.run`` (ctx integration).

* #24 — ``evolution._cluster_clones`` reassigned the >10k-cell tail to clones via a
  per-centroid broadcast loop building k full dense temporaries. It now uses
  ``scipy.spatial.distance.cdist``. We assert NUMERIC EQUIVALENCE of the produced
  labels (old broadcast-argmin vs new cdist-argmin) on a realistic clustered fixture.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.spatial.distance import cdist

from workflow.modular.config import BatchConfig, CellRangerConfig, PipelineConfig
from workflow.modular.context import PipelineContext
from workflow.modular.modules.composition import CompositionModule
from workflow.modular.modules.evolution import EvolutionModule


def _ctx(tmp_path, random_state: int) -> PipelineContext:
    """Smallest PipelineContext exposing ``ctx.random_state`` for module.run()."""
    run_dir = tmp_path / f"run_{random_state}"
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg = PipelineConfig(
        project="det",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(
            sample_root=tmp_path / "ds", outs_dir=tmp_path / "outs"
        ),
        optional_modules=[],
        random_state=random_state,
    )
    return PipelineContext(
        cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir
    )


# --------------------------------------------------------------------------- #
# Composition (#20): single-group permutation enrichment z-scores
# --------------------------------------------------------------------------- #

def _single_group_counts() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.DataFrame(
        [[120, 80, 60, 40, 30, 20, 18, 12]],
        index=["S1"],
        columns=[f"ct{i}" for i in range(8)],
    )
    props = counts.div(counts.sum(axis=1), axis=0)
    return counts, props


def test_composition_fallback_seed_wired_unit():
    """_fallback_test: same seed identical, different seeds differ (literal gone)."""
    counts, props = _single_group_counts()

    a = CompositionModule._fallback_test(counts, props, 11)
    b = CompositionModule._fallback_test(counts, props, 11)
    c = CompositionModule._fallback_test(counts, props, 99)

    # z_score is the stochastic (permutation-derived) column.
    assert a["z_score"].equals(b["z_score"]), "same random_state must be deterministic"
    assert not a["z_score"].equals(c["z_score"]), (
        "different random_state must change permutation z-scores "
        "(proves the hardcoded seed=42 is no longer used)"
    )


def test_composition_random_state_threaded_through_run(tmp_path, monkeypatch):
    """ctx.random_state reaches _fallback_test via CompositionModule.run."""
    captured: list[int] = []
    orig = CompositionModule.__dict__["_fallback_test"].__func__

    def spy(count_df, prop_df, random_state=42):
        captured.append(random_state)
        return orig(count_df, prop_df, random_state)

    monkeypatch.setattr(
        CompositionModule, "_fallback_test", staticmethod(spy)
    )
    # Force the scipy fallback path (skip optional pertpy/scCODA).
    monkeypatch.setattr(
        CompositionModule, "_try_pertpy", staticmethod(lambda adata, key: None)
    )
    # Silence figure side-effects; we only care about the seed threading.
    monkeypatch.setattr(
        CompositionModule, "_plot_barplot", staticmethod(lambda prop_df, ctx: None)
    )
    monkeypatch.setattr(
        CompositionModule, "_plot_boxplot", staticmethod(lambda prop_df, ctx: None)
    )

    # Minimal single-group AnnData: no multi-sample column and a single leiden
    # cluster -> _find_group_key falls back to "leiden", count_df has one row, so
    # _fallback_test takes the n_groups<2 permutation path (the seeded branch).
    n = 40
    rng = np.random.default_rng(0)
    adata = AnnData(rng.random((n, 3)).astype(np.float32))
    adata.obs["cell_type"] = pd.Categorical(
        ["A"] * 10 + ["B"] * 10 + ["C"] * 12 + ["D"] * 8
    )
    adata.obs["leiden"] = pd.Categorical(["0"] * n)  # single group -> permutation path

    ctx = _ctx(tmp_path, random_state=123)
    ctx.cfg.batch = BatchConfig(batch_key="sample")  # absent -> leiden fallback
    ctx.adata = adata

    CompositionModule().run(ctx)

    assert captured == [123], (
        f"_fallback_test must receive ctx.random_state (got {captured})"
    )


# --------------------------------------------------------------------------- #
# Evolution (#20): >10k-cell CNV subsample drives clone clustering
# --------------------------------------------------------------------------- #

def _noise_cnv_adata(n: int = 12000, n_feat: int = 30, seed: int = 0) -> AnnData:
    """Unstructured CNV matrix: the partition is genuinely subsample-dependent,
    so a change in the seed changes clone labels (proving the seed is wired)."""
    rng = np.random.RandomState(seed)
    adata = AnnData(np.zeros((n, 2), dtype=np.float32))
    adata.obsm["X_cnv"] = rng.standard_normal((n, n_feat)).astype(np.float32)
    return adata


def test_evolution_cluster_clones_seed_wired_unit():
    """_cluster_clones: same seed identical, different seeds differ (literal gone)."""
    adata = _noise_cnv_adata()

    a, _ = EvolutionModule._cluster_clones(adata, 11)
    b, _ = EvolutionModule._cluster_clones(adata, 11)
    c, _ = EvolutionModule._cluster_clones(adata, 99)

    assert np.array_equal(a, b), "same random_state must be deterministic"
    assert not np.array_equal(a, c), (
        "different random_state must change clone assignment "
        "(proves the hardcoded RandomState(42) is no longer used)"
    )


def test_evolution_random_state_threaded_through_run(tmp_path, monkeypatch):
    """ctx.random_state reaches _cluster_clones via EvolutionModule.run."""
    captured: list[int] = []
    orig = EvolutionModule.__dict__["_cluster_clones"].__func__

    def spy(adata, random_state=42):
        captured.append(random_state)
        return orig(adata, random_state)

    monkeypatch.setattr(
        EvolutionModule, "_cluster_clones", staticmethod(spy)
    )
    # Stub all downstream tables/plots; only the seed threading is under test.
    for meth in (
        "_compute_clone_stats",
        "_clone_markers",
        "_plot_clone_umap",
        "_plot_clone_composition",
        "_plot_evolution_timeline",
        "_plot_phylo_dendrogram",
        "_plot_cnv_by_clone",
    ):
        monkeypatch.setattr(
            EvolutionModule, meth, staticmethod(lambda *a, **k: None)
        )
    # _compute_clone_stats must still return a DataFrame for .to_csv.
    monkeypatch.setattr(
        EvolutionModule,
        "_compute_clone_stats",
        staticmethod(lambda adata, has_pt: pd.DataFrame()),
    )

    adata = _noise_cnv_adata(n=12000, n_feat=10)
    ctx = _ctx(tmp_path, random_state=123)
    ctx.adata = adata

    EvolutionModule().run(ctx)

    assert captured == [123], (
        f"_cluster_clones must receive ctx.random_state (got {captured})"
    )


# --------------------------------------------------------------------------- #
# Evolution (#24): cdist reassignment is numerically equivalent to the old loop
# --------------------------------------------------------------------------- #

def test_evolution_cdist_reassignment_equivalent_to_broadcast_loop():
    """New cdist nearest-centroid labels == old per-centroid broadcast-argmin.

    Mirrors the >10k tail-reassignment block of _cluster_clones on a realistic
    clustered CNV fixture (distinct cluster means -> no measure-zero exact ties),
    computing labels both the OLD and NEW way and asserting array equality.
    """
    rng = np.random.RandomState(7)
    n, n_feat, k = 12000, 50, 6
    true = rng.randint(0, k, size=n)
    centers = rng.standard_normal((k, n_feat)) * 3.0
    cnv_mat = (centers[true] + rng.standard_normal((n, n_feat))).astype(np.float32)

    sub_idx = rng.choice(n, 5000, replace=False)
    sub_labels = true[sub_idx] + 1  # labels in 1..k
    centroids = np.array(
        [cnv_mat[sub_idx[sub_labels == c]].mean(axis=0) for c in range(1, k + 1)]
    )
    other_idx = np.setdiff1d(np.arange(n), sub_idx)

    # OLD: stack per-centroid distance rows -> (k, n_other), argmin axis=0
    old_dists = np.array(
        [np.linalg.norm(cnv_mat[other_idx] - c, axis=1) for c in centroids]
    )
    old_labels = old_dists.argmin(axis=0) + 1

    # NEW: cdist -> (n_other, k), argmin axis=1
    new_dists = cdist(cnv_mat[other_idx], centroids, metric="euclidean")
    new_labels = new_dists.argmin(axis=1) + 1

    assert np.array_equal(old_labels, new_labels), (
        "cdist reassignment must produce identical clone labels to the "
        "original broadcast euclidean-argmin"
    )
    assert int(new_labels.min()) >= 1 and int(new_labels.max()) <= k
