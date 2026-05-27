"""WS1 science-correctness regressions: #4 ComBat + #5 doublet counts.

Covers two landed fixes that depend on batch_correction.py / doublet_detection.py:
  #4 ComBat correction no-op -> X_pca_combat recompute + use_rep wiring.
  #5 Scrublet fed fractional DecontX counts -> raw-counts priority resolver +
     try/finally .X swap with restore.

The #1 annotation/batch_correction leiden-race regression lives in
test_ws1_leiden_race.py (it depends only on the self-contained #1 surface).

Plan: .omc/plans/2026-05-24-pipeline-governance-remediation.md (WS1 #4/#5).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import pytest

from workflow.modular.modules.doublet_detection import DoubletDetectionModule


# ---------------------------------------------------------------------------
# #5 doublet raw-counts resolver + .X swap
# ---------------------------------------------------------------------------


def _counts_adata(n_obs=120, n_vars=80, integer_X=True, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(1.0, size=(n_obs, n_vars)).astype(np.float32)
    adata = ad.AnnData(counts.copy())
    adata.obs_names = [f"c{i}" for i in range(n_obs)]
    adata.var_names = [f"g{i}" for i in range(n_vars)]
    if not integer_X:
        # Simulate fractional ambient-corrected .X.
        adata.X = counts * 0.137
    return adata, counts


def test_resolver_prefers_layers_counts():
    adata, counts = _counts_adata()
    adata.layers["counts"] = counts.copy()
    adata.layers["counts_raw_pre_decontx"] = counts.copy()
    mat, src = DoubletDetectionModule._resolve_raw_counts_matrix(adata)
    assert src == "layers.counts"


def test_resolver_uses_pre_decontx_layer_when_no_counts():
    adata, counts = _counts_adata(integer_X=False)
    adata.layers["counts_raw_pre_decontx"] = counts.copy()
    mat, src = DoubletDetectionModule._resolve_raw_counts_matrix(adata)
    assert src == "layers.counts_raw_pre_decontx"
    assert DoubletDetectionModule._matrix_looks_integer_like(mat)


def test_resolver_uses_integer_X_when_no_layers():
    adata, counts = _counts_adata(integer_X=True)
    mat, src = DoubletDetectionModule._resolve_raw_counts_matrix(adata)
    assert src == "adata.X_integer_like"


def test_resolver_raises_on_fractional_X_without_raw_layer():
    """Fractional .X + no raw layer -> loud RAISE (no silent fallback)."""
    adata, _ = _counts_adata(integer_X=False)
    with pytest.raises(RuntimeError, match="no raw integer-UMI count source"):
        DoubletDetectionModule._resolve_raw_counts_matrix(adata)


class _DCfg:
    expected_doublet_rate = 0.06
    remove_doublets = False
    backend = "scrublet"


class _DCtx:
    def __init__(self, adata, tmp_path):
        self.adata = adata
        self.metadata: dict = {}
        self.module_status: list = []
        self.figure_dir = tmp_path
        self.table_dir = tmp_path

        class _C:
            doublet = _DCfg()
            doublet_strategy = "whole"
            random_state = 0
        self.cfg = _C()

    @property
    def random_state(self):
        return 0

    def status(self, *a, **kw):
        pass


def _patch_scrublet_capture(monkeypatch):
    """Patch CPU scrublet to capture the matrix it is fed and return singlets."""
    captured = {}
    import scrublet as scr

    class _FakeScrublet:
        def __init__(self, matrix, **kw):
            captured["matrix"] = matrix

        def scrub_doublets(self, **kw):
            n = captured["matrix"].shape[0]
            return np.zeros(n, dtype=float), np.zeros(n, dtype=bool)

        threshold_ = 0.5

    monkeypatch.setattr(scr, "Scrublet", _FakeScrublet)
    # Force CPU path (rsc unavailable in CI).
    monkeypatch.setattr(
        "workflow.modular.modules.doublet_detection._RSC_AVAILABLE", False
    )
    return captured


@pytest.mark.parametrize("state", ["ambient_fired", "ambient_skipped", "pseudobulk"])
def test_doublet_feeds_integer_matrix_to_scrublet(monkeypatch, tmp_path, state):
    """Scrublet receives integer counts in all ambient/pseudobulk states."""
    captured = _patch_scrublet_capture(monkeypatch)

    if state == "ambient_fired":
        adata, counts = _counts_adata(integer_X=False)
        adata.layers["counts_raw_pre_decontx"] = counts.copy()
    elif state == "pseudobulk":
        adata, counts = _counts_adata(integer_X=False)
        adata.layers["counts"] = counts.copy()
    else:  # ambient_skipped: integer .X, no raw layer
        adata, counts = _counts_adata(integer_X=True)

    ctx = _DCtx(adata, tmp_path)
    import unittest.mock as mock
    import workflow.modular.modules.doublet_detection as dd_mod
    with mock.patch.object(dd_mod.plt, "savefig", return_value=None), \
         mock.patch.object(dd_mod.plt, "subplots",
                           return_value=(mock.MagicMock(), mock.MagicMock())), \
         mock.patch.object(dd_mod.plt, "close", return_value=None):
        DoubletDetectionModule().run(ctx)

    fed = captured["matrix"]
    assert DoubletDetectionModule._matrix_looks_integer_like(fed), (
        f"{state}: scrublet was fed non-integer matrix"
    )


def test_doublet_restores_X_after_scrublet(monkeypatch, tmp_path):
    """The swapped .X is restored on the original object after running."""
    captured = _patch_scrublet_capture(monkeypatch)
    adata, counts = _counts_adata(integer_X=False)
    fractional = adata.X.copy()
    adata.layers["counts_raw_pre_decontx"] = counts.copy()

    ctx = _DCtx(adata, tmp_path)
    import unittest.mock as mock
    import workflow.modular.modules.doublet_detection as dd_mod
    with mock.patch.object(dd_mod.plt, "savefig", return_value=None), \
         mock.patch.object(dd_mod.plt, "subplots",
                           return_value=(mock.MagicMock(), mock.MagicMock())), \
         mock.patch.object(dd_mod.plt, "close", return_value=None):
        DoubletDetectionModule().run(ctx)

    # .X restored to the fractional matrix (scrublet saw integer counts).
    assert np.allclose(np.asarray(ctx.adata.X), fractional)
    assert ctx.metadata["doublet_input_source"] == "layers.counts_raw_pre_decontx"


def test_doublet_restores_X_after_forced_scrublet_exception(monkeypatch, tmp_path):
    """try/finally restores .X even when scrublet raises (opt-in null fallback)."""
    import scrublet as scr

    class _BoomScrublet:
        def __init__(self, matrix, **kw):
            pass

        def scrub_doublets(self, **kw):
            raise RuntimeError("forced scrublet failure")

    monkeypatch.setattr(scr, "Scrublet", _BoomScrublet)
    monkeypatch.setattr(
        "workflow.modular.modules.doublet_detection._RSC_AVAILABLE", False
    )
    # Acknowledge the null fallback so the failure path completes (and we can
    # assert .X restoration after the exception is handled).
    monkeypatch.setenv("SC_ALLOW_DOUBLET_NULL_FALLBACK", "1")

    adata, counts = _counts_adata(integer_X=False)
    fractional = adata.X.copy()
    adata.layers["counts_raw_pre_decontx"] = counts.copy()

    ctx = _DCtx(adata, tmp_path)
    import unittest.mock as mock
    import workflow.modular.modules.doublet_detection as dd_mod
    with mock.patch.object(dd_mod.plt, "savefig", return_value=None), \
         mock.patch.object(dd_mod.plt, "subplots",
                           return_value=(mock.MagicMock(), mock.MagicMock())), \
         mock.patch.object(dd_mod.plt, "close", return_value=None):
        DoubletDetectionModule().run(ctx)

    assert np.allclose(np.asarray(ctx.adata.X), fractional)


# ---------------------------------------------------------------------------
# #4 ComBat correction is no longer a no-op
# ---------------------------------------------------------------------------


def _batch_adata(seed=3):
    """2-batch log-normalized fixture with PCA + UMAP + leiden already set."""
    import scanpy as sc
    rng = np.random.default_rng(seed)
    n_per = 60
    n_genes = 50
    # Two batches with a deterministic additive batch offset on raw expression.
    a = rng.normal(2.0, 1.0, size=(n_per, n_genes)).astype(np.float32)
    b = rng.normal(2.0, 1.0, size=(n_per, n_genes)).astype(np.float32) + 3.0
    X = np.vstack([a, b]).astype(np.float32)
    obs = pd.DataFrame(
        {"sample": ["A"] * n_per + ["B"] * n_per},
        index=[f"cell_{i}" for i in range(2 * n_per)],
    )
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))
    sc.pp.pca(adata, n_comps=10)
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor="igraph", directed=False, random_state=0)
    return adata


class _ClusterCfg:
    n_pcs = 10
    n_neighbors = 15
    leiden_resolution = 1.0
    random_state = 0
    leiden_resolution_sweep = ()


class _BatchCfg:
    batch_key = "sample"
    method = "combat"


class _BCtx:
    def __init__(self, adata, tmp_path):
        self.adata = adata
        self.metadata: dict = {}
        self.module_status: list = []
        self.figure_dir = tmp_path
        self.table_dir = tmp_path

        class _C:
            batch = _BatchCfg()
            clustering = _ClusterCfg()
            gpu_mode = "off"
        self.cfg = _C()

    @property
    def random_state(self):
        return 0

    def status(self, *a, **kw):
        pass


def test_combat_produces_x_pca_combat_and_uses_it(tmp_path):
    """ComBat recomputes PCA into X_pca_combat, sets use_rep, mutates leiden."""
    import unittest.mock as mock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod

    adata = _batch_adata()
    pre_leiden = adata.obs["leiden"].astype(str).to_numpy().copy()
    ctx = _BCtx(adata, tmp_path)

    with mock.patch.object(bc_mod.plt, "savefig", return_value=None), \
         mock.patch.object(bc_mod.plt, "subplots",
                           return_value=(mock.MagicMock(), [mock.MagicMock(), mock.MagicMock()])), \
         mock.patch.object(bc_mod.plt, "tight_layout", return_value=None), \
         mock.patch.object(bc_mod.plt, "close", return_value=None), \
         mock.patch.object(bc_mod.sc.pl, "umap", return_value=None):
        BatchCorrectionModule().run(ctx)

    assert "X_pca_combat" in adata.obsm, "ComBat did not write X_pca_combat"
    assert ctx.metadata["batch_correction_use_rep"] == "X_pca_combat"
    assert ctx.metadata["batch_correction_status"] == "completed"
    # leiden_corrected snapshot written (non-destructive audit).
    assert "leiden_corrected" in adata.obs
    # Mixing 'after' metric computed on the corrected representation.
    import json
    audit = json.loads((tmp_path / "batch_mixing_metrics.json").read_text())
    assert audit["after"]["use_rep"] == "X_pca_combat"
    # Embedding actually changed clustering relative to the pre-combat leiden
    # (the no-op bug would have left leiden identical because X_pca was stale).
    post_leiden = adata.obs["leiden"].astype(str).to_numpy()
    assert not np.array_equal(pre_leiden, post_leiden), (
        "ComBat left leiden unchanged — embedding recompute did not take effect"
    )


def test_combat_run_combat_helper_sets_combat_rep(tmp_path):
    """_run_combat alone populates X_pca_combat and preserves pre-combat X_pca."""
    from workflow.modular.modules.batch_correction import BatchCorrectionModule

    adata = _batch_adata()
    prev_pca = adata.obsm["X_pca"].copy()
    ctx = _BCtx(adata, tmp_path)
    BatchCorrectionModule._run_combat(adata, "sample", ctx)
    assert "X_pca_combat" in adata.obsm
    # Pre-combat X_pca preserved for the 'before' baseline.
    assert np.allclose(adata.obsm["X_pca"], prev_pca)
    # Corrected embedding differs from the pre-combat one.
    assert not np.allclose(adata.obsm["X_pca_combat"], prev_pca)


# ---------------------------------------------------------------------------
# W14 batch_correction post-processing failure must not be swallowed
# ---------------------------------------------------------------------------


def _run_batch_correction_with_failing_postproc(adata, ctx):
    """Run batch_correction with CPU post-processing forced to fail.

    The corrected representation (X_pca_combat via _run_combat) is still
    produced, but neighbors/UMAP/Leiden on it raises — exactly the W14
    scenario where leiden/X_umap would be left STALE (pre-correction).
    """
    import unittest.mock as mock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod

    def _boom(*a, **kw):
        raise RuntimeError("forced post-processing failure")

    with mock.patch.object(bc_mod.plt, "savefig", return_value=None), \
         mock.patch.object(bc_mod.plt, "subplots",
                           return_value=(mock.MagicMock(),
                                         [mock.MagicMock(), mock.MagicMock()])), \
         mock.patch.object(bc_mod.plt, "tight_layout", return_value=None), \
         mock.patch.object(bc_mod.plt, "close", return_value=None), \
         mock.patch.object(bc_mod.sc.pl, "umap", return_value=None), \
         mock.patch.object(bc_mod, "get_or_compute_neighbors", _boom):
        BatchCorrectionModule().run(ctx)


def test_w14_postproc_failure_raises_loud_by_default(tmp_path, monkeypatch):
    """W14: CPU post-processing failure FAILS LOUD by default (no swallow).

    Previously the module logged a warning, kept the stale pre-correction
    clustering, and still reported status='completed'. The fix raises unless the
    explicit SC_ALLOW_BATCH_BACKEND_SKIP opt-in is set.
    """
    monkeypatch.delenv("SC_ALLOW_BATCH_BACKEND_SKIP", raising=False)
    adata = _batch_adata()
    ctx = _BCtx(adata, tmp_path)
    with pytest.raises(RuntimeError, match="SC_ALLOW_BATCH_BACKEND_SKIP"):
        _run_batch_correction_with_failing_postproc(adata, ctx)
    # Must NOT have been marked completed on the loud-failure path.
    assert ctx.metadata.get("batch_correction_status") != "completed"


def test_w14_postproc_failure_optin_uses_distinct_status_and_records_staleness(
    tmp_path, monkeypatch
):
    """W14: with the audited opt-in, the module does NOT report 'completed'.

    It uses a distinct status and records that the clustering is stale relative
    to the corrected representation.
    """
    monkeypatch.setenv("SC_ALLOW_BATCH_BACKEND_SKIP", "1")
    adata = _batch_adata()
    ctx = _BCtx(adata, tmp_path)
    _run_batch_correction_with_failing_postproc(adata, ctx)

    assert ctx.metadata["batch_correction_status"] == (
        "completed_with_stale_clustering_opt_in"
    )
    assert ctx.metadata["batch_correction_status"] != "completed"
    assert ctx.metadata["batch_post_processing_stale_clustering"] is True
    assert ctx.metadata["batch_post_processing_skip_opt_in_acknowledged"] is True
    assert "batch_post_error" in ctx.metadata
