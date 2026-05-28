"""Wiring + e2e acceptance for the integration_select gate + Harmony-direct backend.

Covers TASK W1 (working harmonypy-direct backend in batch_correction) and TASK
W2/W3 (integration_select module + registration), per plan TM-10. The discovery
ENGINE itself (discovery_metrics / select_integration / methods) is exercised by
tests/test_discovery_integration_gate.py; this suite asserts the PRODUCTION
wiring: the gate sets cfg.batch.method, routes Harmony to the working direct
backend, never edits clustering's use_rep, and batch_correction re-clusters on
the corrected representation.

scib-metrics + leidenalg + harmonypy are required; tests skip (not silently
pass) if the engine is unavailable.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

pytest.importorskip("scib_metrics")
pytest.importorskip("leidenalg")
pytest.importorskip("scanpy")
pytest.importorskip("harmonypy")

import anndata as ad  # noqa: E402

from workflow.modular.config import CellRangerConfig, PipelineConfig  # noqa: E402
from workflow.modular.context import PipelineContext  # noqa: E402
from workflow.modular.modules.batch_correction import BatchCorrectionModule  # noqa: E402
from workflow.modular.modules.integration_select import IntegrationSelectModule  # noqa: E402


def _load_fixtures():
    path = ROOT / "tests" / "fixtures" / "discovery_integration_synthetic.py"
    name = "discovery_integration_synthetic_fixture"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


fixtures = _load_fixtures()


def _make_ctx(adata, tmp_path) -> PipelineContext:
    cfg = PipelineConfig(
        project="ints",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(
            sample_root=tmp_path / "ds", outs_dir=tmp_path / "ds" / "outs"
        ),
        optional_modules=[],
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    return PipelineContext(
        cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir, adata=adata
    )


def _adata_from_fixture(fx, *, n_pcs: int = 7) -> ad.AnnData:
    """Build an AnnData with X_pca (baseline), a leiden col, and a 2-level batch.

    The fixture baseline IS the shared pre-integration PCA; we install it as
    X_pca. A leiden column is attached so the gate's "clustering ran first"
    contract is visibly satisfied (the gate keys on X_pca, not leiden).
    """
    base = np.asarray(fx.embeddings["baseline"], dtype=np.float32)
    n = base.shape[0]
    counts = np.random.default_rng(0).poisson(1.0, size=(n, 12)).astype(np.float32)
    a = ad.AnnData(counts)
    a.layers["counts"] = counts.copy()
    a.obsm["X_pca"] = base[:, :n_pcs].copy()
    a.obs["sample"] = [str(b) for b in fx.batch]
    # A leiden label so downstream re-clustering has something to overwrite.
    a.obs["leiden"] = (np.arange(n) % 3).astype(str)
    a.obs["leiden"] = a.obs["leiden"].astype("category")
    return a


# ==========================================================================
# Unit: _run_harmony_direct produces (n_cells, n_pcs) — guards the prod shape bug
# ==========================================================================

def test_run_harmony_direct_shape_is_n_cells_by_n_pcs(tmp_path):
    """Direct backend MUST store X_pca_harmony as (n_cells, n_pcs), never (n_pcs,).

    This directly guards against the production scanpy/harmonypy 0.2.0 bug that
    stored the embedding with the wrong shape.
    """
    rng = np.random.default_rng(0)
    n, n_pcs = 90, 15
    counts = rng.poisson(1.0, size=(n, 20)).astype(np.float32)
    a = ad.AnnData(counts)
    a.obs["sample"] = (["b0"] * 45 + ["b1"] * 45)
    a.obsm["X_pca"] = rng.normal(size=(n, n_pcs)).astype(np.float32)
    ctx = _make_ctx(a, tmp_path)

    BatchCorrectionModule._run_harmony_direct(a, "sample", ctx)

    assert "X_pca_harmony" in a.obsm
    Z = a.obsm["X_pca_harmony"]
    assert Z.ndim == 2
    assert Z.shape == (n, n_pcs), f"expected ({n},{n_pcs}), got {Z.shape}"
    assert Z.dtype == np.float32
    assert ctx.metadata["harmony_backend"] == "direct"
    assert ctx.metadata["harmony_converged"] is True


def test_resolve_harmony_backend_accepts_direct(tmp_path, monkeypatch):
    """`direct` is an explicit, valid backend choice; `auto` now routes to `direct`.

    DOCUMENTED ROUTING CHANGE (2026-05-27): both gpu (rapids CUBLAS) and cpu
    (scanpy 1.12 wrapper mis-stores harmonypy Z_corr) are broken on this env, so
    `auto` resolves to the proven-working `direct` path. The escape hatches stay
    intact: explicit gpu still raises without rsc, explicit cpu is unchanged, and
    SC_HARMONY_BACKEND still overrides config.
    """
    monkeypatch.delenv("SC_HARMONY_BACKEND", raising=False)
    a = ad.AnnData(np.zeros((4, 2), dtype=np.float32))
    ctx = _make_ctx(a, tmp_path)
    ctx.cfg.batch.harmony_backend = "direct"
    assert BatchCorrectionModule._resolve_harmony_backend(ctx) == "direct"
    # auto now resolves TO direct (documented routing change on this env).
    ctx.cfg.batch.harmony_backend = "auto"
    assert BatchCorrectionModule._resolve_harmony_backend(ctx) == "direct"
    # Escape hatch: explicit cpu is honored unchanged.
    ctx.cfg.batch.harmony_backend = "cpu"
    assert BatchCorrectionModule._resolve_harmony_backend(ctx) == "cpu"
    # Escape hatch: SC_HARMONY_BACKEND env override wins over config.
    ctx.cfg.batch.harmony_backend = "auto"
    monkeypatch.setenv("SC_HARMONY_BACKEND", "cpu")
    assert BatchCorrectionModule._resolve_harmony_backend(ctx) == "cpu"


# ==========================================================================
# TM-10 (e2e): gate -> batch_correction; method set, direct backend, shapes,
# clustering untouched, re-cluster on corrected rep.
# ==========================================================================

def test_tm10_gate_then_batch_correction_e2e(tmp_path, monkeypatch):
    fx = fixtures.make_f1_batch_split()
    adata = _adata_from_fixture(fx)
    n = adata.n_obs
    n_pcs = adata.obsm["X_pca"].shape[1]
    ctx = _make_ctx(adata, tmp_path)

    # Snapshot clustering's hardcoded use_rep BEFORE the gate runs.
    use_rep_before = ctx.cfg.clustering  # config object; we assert the module never edits it
    pca_before = adata.obsm["X_pca"].copy()

    # Inject the fixture's known-answer embeddings so the gate's scoring +
    # recommendation chain runs deterministically and fast (no real scVI train).
    # The fixture's GOOD harmony embedding makes the gate pick "harmony".
    import scripts.bench.integration.methods as M
    import scripts.bench.integration.select_integration as sel

    def _fake_embeddings(adata_arg, batch_key, scvi_seeds, **kwargs):
        emb = {
            M.BASELINE_METHOD: np.asarray(fx.embeddings["baseline"], dtype=np.float64),
            M.HARMONY_METHOD: np.asarray(fx.embeddings["harmony"], dtype=np.float64),
            sel.SHUFFLE_CONTROL_METHOD: np.asarray(
                fx.embeddings["neg_control_shuffle_label"], dtype=np.float64
            ),
        }
        for seed in scvi_seeds:
            emb[f"{M.SCVI_METHOD}_seed{int(seed)}"] = np.asarray(
                fx.embeddings["baseline"], dtype=np.float64
            )
        return emb, []

    monkeypatch.setattr(
        IntegrationSelectModule, "_compute_embeddings", staticmethod(_fake_embeddings)
    )

    gate = IntegrationSelectModule()
    gate.run(ctx)

    # (a) the gate SET cfg.batch.method.
    assert ctx.cfg.batch.method == "harmony"
    md = ctx.metadata["integration_select"]
    assert md["chosen_method"] == "harmony"

    # (b) harmony chosen -> backend routed to the working direct path.
    assert ctx.cfg.batch.harmony_backend == "direct"

    # (c) clustering's hardcoded use_rep="X_pca" is UNTOUCHED by the gate.
    assert ctx.cfg.clustering is use_rep_before
    assert np.array_equal(adata.obsm["X_pca"], pca_before)
    # The gate must not invent its own clustering rep.
    assert "X_pca_harmony" not in adata.obsm  # only batch_correction writes it

    # Now run batch_correction with the gate-applied method/backend.
    bc = BatchCorrectionModule()
    bc.run(ctx)

    # (b cont.) batch_correction produced X_pca_harmony with (n_cells, n_pcs).
    assert "X_pca_harmony" in adata.obsm
    assert adata.obsm["X_pca_harmony"].shape == (n, n_pcs)
    assert ctx.metadata["harmony_backend"] == "direct"

    # (d) batch_correction re-clustered on the corrected representation.
    assert ctx.metadata["batch_correction_use_rep"] == "X_pca_harmony"
    assert ctx.metadata.get("batch_correction_status") == "completed"


# ==========================================================================
# Integration: 2-batch picks non-none + sets cfg; single-batch skips gracefully.
# ==========================================================================

def test_two_batch_selects_non_none_and_sets_cfg(tmp_path, monkeypatch):
    fx = fixtures.make_f1_batch_split()
    adata = _adata_from_fixture(fx)
    ctx = _make_ctx(adata, tmp_path)
    # Start from a non-harmony method to prove the gate actively sets it.
    ctx.cfg.batch.method = "combat"

    import scripts.bench.integration.methods as M
    import scripts.bench.integration.select_integration as sel

    def _fake_embeddings(adata_arg, batch_key, scvi_seeds, **kwargs):
        emb = {
            M.BASELINE_METHOD: np.asarray(fx.embeddings["baseline"], dtype=np.float64),
            M.HARMONY_METHOD: np.asarray(fx.embeddings["harmony"], dtype=np.float64),
            sel.SHUFFLE_CONTROL_METHOD: np.asarray(
                fx.embeddings["neg_control_shuffle_label"], dtype=np.float64
            ),
        }
        for seed in scvi_seeds:
            emb[f"{M.SCVI_METHOD}_seed{int(seed)}"] = np.asarray(
                fx.embeddings["baseline"], dtype=np.float64
            )
        return emb, []

    monkeypatch.setattr(
        IntegrationSelectModule, "_compute_embeddings", staticmethod(_fake_embeddings)
    )

    IntegrationSelectModule().run(ctx)
    # Strong batch effect -> gate must pick a real integration, never "none".
    assert ctx.cfg.batch.method != "none"
    assert ctx.cfg.batch.method == "harmony"
    assert ctx.metadata["integration_select_status"] == "completed"
    # Recommendation artifacts written under the (non-factory) module dir.
    out = Path(ctx.table_dir)
    assert (out / "integration_recommendation.json").is_file()
    assert (out / "integration_scoreboard.csv").is_file()
    assert (out / "integration_audit.md").is_file()


def test_single_batch_skips_gracefully(tmp_path):
    fx = fixtures.make_f1_batch_split()
    adata = _adata_from_fixture(fx)
    # Collapse to a single batch level.
    adata.obs["sample"] = "only_one"
    ctx = _make_ctx(adata, tmp_path)
    ctx.cfg.batch.method = "harmony"

    IntegrationSelectModule().run(ctx)

    # Skips without failing; leaves the method untouched.
    assert ctx.metadata["integration_select_status"] == "skipped_single_batch"
    assert ctx.cfg.batch.method == "harmony"
    statuses = [e for e in ctx.module_status if e["module"] == "integration_select"]
    assert statuses and statuses[-1]["status"] == "skipped"


def test_missing_batch_key_skips_gracefully(tmp_path):
    fx = fixtures.make_f1_batch_split()
    adata = _adata_from_fixture(fx)
    del adata.obs["sample"]
    ctx = _make_ctx(adata, tmp_path)
    ctx.cfg.batch.batch_key = "sample"

    IntegrationSelectModule().run(ctx)
    assert ctx.metadata["integration_select_status"] == "skipped_missing_batch_key"


def test_missing_pca_fails_loud(tmp_path):
    fx = fixtures.make_f1_batch_split()
    adata = _adata_from_fixture(fx)
    del adata.obsm["X_pca"]
    ctx = _make_ctx(adata, tmp_path)

    with pytest.raises(ValueError, match="X_pca"):
        IntegrationSelectModule().run(ctx)


def test_cache_hit_replays_recommendation(tmp_path, monkeypatch):
    """A second run on identical X_pca/batch HITs the cache (no recompute)."""
    fx = fixtures.make_f1_batch_split()
    import scripts.bench.integration.methods as M
    import scripts.bench.integration.select_integration as sel

    def _fake_embeddings(adata_arg, batch_key, scvi_seeds, **kwargs):
        emb = {
            M.BASELINE_METHOD: np.asarray(fx.embeddings["baseline"], dtype=np.float64),
            M.HARMONY_METHOD: np.asarray(fx.embeddings["harmony"], dtype=np.float64),
            sel.SHUFFLE_CONTROL_METHOD: np.asarray(
                fx.embeddings["neg_control_shuffle_label"], dtype=np.float64
            ),
        }
        for seed in scvi_seeds:
            emb[f"{M.SCVI_METHOD}_seed{int(seed)}"] = np.asarray(
                fx.embeddings["baseline"], dtype=np.float64
            )
        return emb, []

    monkeypatch.setattr(
        IntegrationSelectModule, "_compute_embeddings", staticmethod(_fake_embeddings)
    )

    # First run: compute + store. Use the SAME run_dir so the cache persists.
    adata1 = _adata_from_fixture(fx)
    ctx1 = _make_ctx(adata1, tmp_path)
    IntegrationSelectModule().run(ctx1)
    assert ctx1.metadata["integration_select_status"] == "completed"

    adata2 = _adata_from_fixture(fx)
    ctx2 = _make_ctx(adata2, tmp_path)  # same tmp_path/run -> same cache dir
    IntegrationSelectModule().run(ctx2)
    assert ctx2.metadata["integration_select_status"] == "completed_cache_hit"
    assert ctx2.cfg.batch.method == ctx1.cfg.batch.method
    # W13: the HIT path re-asserts the engine pins and records the resolved map.
    assert "pins_reasserted_on_hit" in ctx2.metadata["integration_select"]


def test_w13_cache_hit_reasserts_engine_pins_and_fails_loud_on_drift(
    tmp_path, monkeypatch
):
    """W13: a cache HIT must re-run assert_pinned_discovery_engine() and FAIL
    LOUD on engine-version drift, instead of silently serving a stored
    recommendation computed under different (now-drifted) pins.

    First populate the cache (clean pins). Then, on the second run, force a
    Leiden-engine version mismatch and assert the HIT path raises (proving the
    pin assertion now runs on HIT, not only on MISS).
    """
    fx = fixtures.make_f1_batch_split()
    import scripts.bench.integration.methods as M
    import scripts.bench.integration.select_integration as sel
    from scripts.bench.integration import discovery_metrics as dm

    def _fake_embeddings(adata_arg, batch_key, scvi_seeds, **kwargs):
        emb = {
            M.BASELINE_METHOD: np.asarray(fx.embeddings["baseline"], dtype=np.float64),
            M.HARMONY_METHOD: np.asarray(fx.embeddings["harmony"], dtype=np.float64),
            sel.SHUFFLE_CONTROL_METHOD: np.asarray(
                fx.embeddings["neg_control_shuffle_label"], dtype=np.float64
            ),
        }
        for seed in scvi_seeds:
            emb[f"{M.SCVI_METHOD}_seed{int(seed)}"] = np.asarray(
                fx.embeddings["baseline"], dtype=np.float64
            )
        return emb, []

    monkeypatch.setattr(
        IntegrationSelectModule, "_compute_embeddings", staticmethod(_fake_embeddings)
    )

    # First run with clean pins: compute + store in the cache.
    monkeypatch.delenv(dm.SC_ALLOW_SCIB_VERSION_DRIFT, raising=False)
    adata1 = _adata_from_fixture(fx)
    ctx1 = _make_ctx(adata1, tmp_path)
    IntegrationSelectModule().run(ctx1)
    assert ctx1.metadata["integration_select_status"] == "completed"

    # Second run on identical X_pca/batch -> cache HIT. Now drift the Leiden pin.
    real = dm._installed_version_any

    def fake(aliases):
        if "leidenalg" in aliases:
            return "0.0.0-bogus"
        return real(aliases)

    monkeypatch.setattr(dm, "_installed_version_any", fake)

    adata2 = _adata_from_fixture(fx)
    ctx2 = _make_ctx(adata2, tmp_path)  # same cache dir -> HIT
    with pytest.raises(RuntimeError, match="Leiden-engine version gate FAILED"):
        IntegrationSelectModule().run(ctx2)


def test_cache_hit_rejects_malformed_payload(tmp_path):
    """A malformed cache entry must not silently become method='none'."""
    fx = fixtures.make_f1_batch_split()
    adata = _adata_from_fixture(fx)
    ctx = _make_ctx(adata, tmp_path)

    import scripts.bench.integration.integration_cache as ic

    batch = adata.obs["sample"].astype(str).to_numpy()
    key = ic.make_cache_key(
        np.asarray(adata.obsm["X_pca"], dtype=np.float64),
        "sample",
        scvi_seeds=ctx.cfg.batch.integration_scvi_seeds,
        margin_mix=ctx.cfg.batch.integration_margin_mix,
        harmony_theta=ctx.cfg.batch.harmony_theta,
        harmony_sigma=ctx.cfg.batch.harmony_sigma,
        harmony_max_iter=ctx.cfg.batch.harmony_max_iter,
        scvi_max_epochs=ctx.cfg.batch.scvi_max_epochs,
        scvi_n_latent=ctx.cfg.batch.scvi_n_latent,
        scvi_early_stopping=ctx.cfg.batch.scvi_early_stopping,
        n_neighbors=ctx.cfg.clustering.n_neighbors,
        batch_labels=batch,
        obs_names=adata.obs_names,
    )
    cache_dir = tmp_path / "run" / "integration_select" / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / f"{key.digest()}.json").write_text(
        json.dumps({"schema_version": "discovery-integration-gate-v1"}) + "\n"
    )

    with pytest.raises(RuntimeError, match="missing chosen_method"):
        IntegrationSelectModule().run(ctx)


def test_degraded_required_candidates_fail_loud(tmp_path, monkeypatch):
    """Harmony/scVI/shuffle are the required candidate set for a completed gate."""
    fx = fixtures.make_f1_batch_split()
    adata = _adata_from_fixture(fx)
    ctx = _make_ctx(adata, tmp_path)

    import scripts.bench.integration.methods as M

    def _fake_embeddings(adata_arg, batch_key, scvi_seeds, **kwargs):
        return ({
            M.BASELINE_METHOD: np.asarray(fx.embeddings["baseline"], dtype=np.float64),
            M.HARMONY_METHOD: np.asarray(fx.embeddings["harmony"], dtype=np.float64),
        }, [{"method": "scvi_seed0", "error": "boom"}])

    monkeypatch.setattr(
        IntegrationSelectModule, "_compute_embeddings", staticmethod(_fake_embeddings)
    )

    with pytest.raises(RuntimeError, match="degraded candidate set"):
        IntegrationSelectModule().run(ctx)
    assert ctx.metadata["integration_select_status"] == "failed_degraded_candidates"


def test_shuffle_control_must_fire_before_recommendation(tmp_path, monkeypatch):
    """A present but non-firing shuffle control invalidates the gate."""
    fx = fixtures.make_f1_batch_split()
    adata = _adata_from_fixture(fx)
    ctx = _make_ctx(adata, tmp_path)

    import scripts.bench.integration.methods as M
    import scripts.bench.integration.select_integration as sel

    def _fake_embeddings(adata_arg, batch_key, scvi_seeds, **kwargs):
        emb = {
            M.BASELINE_METHOD: np.asarray(fx.embeddings["baseline"], dtype=np.float64),
            M.HARMONY_METHOD: np.asarray(fx.embeddings["harmony"], dtype=np.float64),
            # Non-firing control: identical to baseline instead of structure-destroying.
            sel.SHUFFLE_CONTROL_METHOD: np.asarray(
                fx.embeddings["baseline"], dtype=np.float64
            ),
        }
        for seed in scvi_seeds:
            emb[f"{M.SCVI_METHOD}_seed{int(seed)}"] = np.asarray(
                fx.embeddings["baseline"], dtype=np.float64
            )
        return emb, []

    monkeypatch.setattr(
        IntegrationSelectModule, "_compute_embeddings", staticmethod(_fake_embeddings)
    )

    with pytest.raises(RuntimeError, match="shuffle falsifiability control"):
        IntegrationSelectModule().run(ctx)
    assert ctx.metadata["integration_select_status"] == "failed_shuffle_control"


def test_compute_embeddings_propagates_gate_params(monkeypatch):
    """Params included in the cache key must also reach the actual computation."""
    import anndata as ad
    import scripts.bench.integration.methods as M
    import scripts.bench.integration.select_integration as sel

    rng = np.random.default_rng(0)
    adata = ad.AnnData(rng.poisson(1.0, size=(12, 5)).astype(np.float32))
    adata.layers["counts"] = adata.X.copy()
    adata.obs["sample"] = ["a"] * 6 + ["b"] * 6
    adata.obsm["X_pca"] = rng.normal(size=(12, 4)).astype(np.float32)
    captured = {}

    def fake_harmony(adata_arg, batch_key, **kwargs):
        captured["harmony"] = kwargs
        return SimpleNamespace(method=M.HARMONY_METHOD, X=adata_arg.obsm["X_pca"], error=None)

    def fake_scvi(adata_arg, batch_key, *, seeds, **kwargs):
        captured["scvi"] = kwargs
        return [
            SimpleNamespace(method=f"{M.SCVI_METHOD}_seed{int(seed)}",
                            X=adata_arg.obsm["X_pca"], error=None)
            for seed in seeds
        ]

    def fake_neg(adata_arg, batch_key):
        return SimpleNamespace(method=M.NEG_CONTROL_METHOD, X=None, error="singular")

    def fake_shuffle(adata_arg, batch_key, seed=0):
        return SimpleNamespace(method=sel.SHUFFLE_CONTROL_METHOD,
                               X=adata_arg.obsm["X_pca"], error=None)

    monkeypatch.setattr(M, "compute_harmonypy_direct", fake_harmony)
    monkeypatch.setattr(M, "compute_scvi_seed_sweep", fake_scvi)
    monkeypatch.setattr(M, "compute_neg_control", fake_neg)
    monkeypatch.setattr(M, "compute_shuffle_label_control", fake_shuffle)

    IntegrationSelectModule._compute_embeddings(
        adata, "sample", (0, 1, 2),
        harmony_theta=4.0,
        harmony_sigma=0.25,
        harmony_max_iter=99,
        scvi_max_epochs=42,
        scvi_n_latent=17,
        scvi_early_stopping=False,
    )

    assert captured["harmony"]["theta"] == 4.0
    assert captured["harmony"]["sigma"] == 0.25
    assert captured["harmony"]["max_iter_harmony"] == 99
    assert captured["scvi"]["scvi_max_epochs"] == 42
    assert captured["scvi"]["scvi_n_latent"] == 17
    assert captured["scvi"]["scvi_early_stopping"] is False
