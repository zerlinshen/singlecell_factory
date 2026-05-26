"""Acceptance suite for the discovery integration-selection gate (plan §7, TM-1..TM-16).

Implements the full DELIBERATE test matrix as the binding acceptance criteria.
scib-metrics + leidenalg are required; tests are skipped (not silently passed)
if the engine is unavailable. The four reviewed files
(score_integration.py / methods.py / test_integration_scoreboard.py /
integration_synthetic.py) are NOT touched — this suite layers on the NEW
discovery modules only.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

pytest.importorskip("scib_metrics")
pytest.importorskip("leidenalg")
pytest.importorskip("scanpy")

from scripts.bench.integration import discovery_metrics as dm  # noqa: E402
from scripts.bench.integration import integration_cache as ic  # noqa: E402
from scripts.bench.integration import select_integration as sel  # noqa: E402


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


# Module-scoped scored fixtures (per-batch reference + scored embeddings) so the
# (relatively expensive) Leiden passes run once per fixture.

@pytest.fixture(scope="module")
def f1():
    f = fixtures.make_f1_batch_split()
    results, ref = dm.score_all_embeddings_discovery(
        f.embeddings, f.embeddings["baseline"], f.batch)
    return {"fx": f, "results": results, "ref": ref,
            "by": {r.method: r for r in results}}


@pytest.fixture(scope="module")
def f2():
    f = fixtures.make_f2_rare_novel()
    results, ref = dm.score_all_embeddings_discovery(
        f.embeddings, f.embeddings["baseline"], f.batch)
    return {"fx": f, "results": results, "ref": ref,
            "by": {r.method: r for r in results}}


@pytest.fixture(scope="module")
def f3():
    f = fixtures.make_f3_single_batch()
    results, ref = dm.score_all_embeddings_discovery(
        f.embeddings, f.embeddings["baseline"], f.batch)
    return {"fx": f, "results": results, "ref": ref,
            "by": {r.method: r for r in results}}


# ==========================================================================
# TM-1 — correspondence deterministic + matches the known split
# ==========================================================================

def test_tm1_good_merge_high_and_undercorrection_not_winning(f1):
    """GOOD merge of the split type scores HIGH on D1a; under-correction
    (baseline, no merge) does NOT win the gate."""
    by = f1["by"]
    # GOOD integration (harmony removes batch axis) co-localizes the split type.
    assert by["harmony"].good_merge_fraction == pytest.approx(1.0)
    # Under-correction baseline: split type stays split -> low D1a.
    assert by["baseline"].good_merge_fraction == pytest.approx(0.0)
    # Gate: baseline is the reference (never a gate-survivor candidate), harmony
    # is the chosen integration -> under-correction does NOT win the gate.
    rec = sel.recommend(f1["results"])
    assert "baseline" not in rec.gate_results  # baseline is the reference, not a candidate
    assert rec.chosen_method == "harmony"


def test_tm1_correspondence_matches_known_split(f1):
    """The correspondence map links the split type's per-batch copies."""
    ref = f1["ref"]
    fx = f1["fx"]
    bstr = np.asarray([str(b) for b in fx.batch])
    # Identify, for each corresponded pair, the dominant true_bio on each side.
    def dominant_bio(key):
        b, u = key.split("::")
        idx = np.where(bstr == b)[0]
        rows = idx[ref.labels[b] == int(u)]
        vals, counts = np.unique(fx.true_bio[rows], return_counts=True)
        return str(vals[int(np.argmax(counts))])
    matched_bios = {tuple(sorted((dominant_bio(a), dominant_bio(b))))
                    for a, b in ref.correspondence}
    # The split type must appear matched to itself across batches.
    assert ("split_type", "split_type") in matched_bios


def test_tm16_determinism_bit_identical(f1):
    """TM-16: per-batch Leiden labels + correspondence map bit-identical across runs."""
    fx = f1["fx"]
    r1 = dm.compute_per_batch_reference(fx.embeddings["baseline"], fx.batch)
    r2 = dm.compute_per_batch_reference(fx.embeddings["baseline"], fx.batch)
    for b in r1.batches:
        assert np.array_equal(r1.labels[b], r2.labels[b])
    assert r1.correspondence == r2.correspondence


# ==========================================================================
# TM-2 — rare/novel preservation; erase collapses D1c
# ==========================================================================

def test_tm2_rare_preserved_high_erase_collapses(f2):
    """Good integration keeps D1c HIGH; the erase transform collapses D1c."""
    by = f2["by"]
    harmony_d1c = by["harmony"].within_batch_isolation
    erase_d1c = by["erase"].within_batch_isolation
    assert harmony_d1c is not None and erase_d1c is not None
    # Erasing the rare type drops D1c clearly below the good integration.
    assert erase_d1c < harmony_d1c - 0.04
    # And the erase D1c lands at/below the disqualifying isolation floor.
    rec = sel.recommend(f2["results"])
    iso_floor = rec.rule_constants["isolation_floor"]
    assert erase_d1c < iso_floor
    assert harmony_d1c >= iso_floor


# ==========================================================================
# TM-3 — bad-merge transform fuses distinct populations -> D1b/gate disqualified
# ==========================================================================

def test_tm3_bad_merge_disqualified(f2):
    """A bad-merge transform that fuses distinct populations is disqualified.

    The erase transform (rare fused into abundant0) is a bad merge of distinct
    populations; it must NOT survive the gates (its D1c drops below floor)."""
    rec = sel.recommend(f2["results"])
    assert rec.gate_results["erase"]["survived_gates"] is False
    # And it is specifically the isolation/distinctness floor that fails.
    g = rec.gate_results["erase"]
    assert (g["isolation_floor"] is False) or (g["distinctness_floor"] is False)


def test_tm3_severe_overmerge_drops_distinctness_below_floor(f1):
    """A severe over-merge that destroys bio separation collapses D1b below floor."""
    fx = f1["fx"]
    ref = f1["ref"]
    rng = np.random.default_rng(5)
    E = fx.embeddings["harmony"].copy()
    E[:, :6] = rng.normal(0, 2.0, (E.shape[0], 6))  # destroy bio separation
    overmerge = dm.score_embedding_discovery("severe_overmerge", E, ref, fx.batch)
    emb = {
        "baseline": fx.embeddings["baseline"],
        "severe_overmerge": E,
        "neg_control_shuffle_label": fx.embeddings["neg_control_shuffle_label"],
    }
    results, _ = dm.score_all_embeddings_discovery(emb, fx.embeddings["baseline"], fx.batch)
    rec = sel.recommend(results)
    g = rec.gate_results["severe_overmerge"]
    # Severe over-merge must be disqualified by the distinctness or isolation floor.
    assert g["survived_gates"] is False
    assert (g["distinctness_floor"] is False) or (g["isolation_floor"] is False)


# ==========================================================================
# TM-4 (3a) — strong batch effect: baseline can NEVER be recommended
# ==========================================================================

@pytest.mark.parametrize("maker", [
    fixtures.make_f1_batch_split,
    fixtures.make_f2_rare_novel,
    fixtures.make_f3_single_batch,
])
def test_tm4_strong_batch_baseline_never_recommended(maker):
    f = maker()
    results, _ = dm.score_all_embeddings_discovery(
        f.embeddings, f.embeddings["baseline"], f.batch)
    rec = sel.recommend(results)
    # baseline mixing is ~0 under a strong batch effect -> fails the mixing floor.
    by = {r.method: r for r in results}
    assert by["baseline"].batch_mixing_score < 0.1
    # The recommendation is a real integration, never baseline ("none").
    assert rec.chosen_method != "none"
    assert rec.reason != "baseline_fallback_no_survivor"


# ==========================================================================
# TM-5 (3b) — no real batch effect: baseline SHOULD win
# ==========================================================================

def test_tm5_no_batch_effect_baseline_wins():
    f = fixtures.make_f1_batch_split(batch_axis_offset=0.0)
    results, _ = dm.score_all_embeddings_discovery(
        f.embeddings, f.embeddings["baseline"], f.batch)
    rec = sel.recommend(results)
    by = {r.method: r for r in results}
    # baseline is already well-mixed (no confound).
    assert by["baseline"].batch_mixing_score > 0.5
    # No candidate beats baseline by MARGIN_MIX -> baseline (no integration).
    assert rec.chosen_method == "none"
    assert rec.reason == "baseline_fallback_no_survivor"


# ==========================================================================
# TM-6 — extreme-theta control REPORTED, NOT required to fire
# TM-7 — shuffle control REQUIRED anchor, MUST FIRE
# ==========================================================================

def test_tm7_shuffle_control_fires(f1):
    """The shuffle control MUST collapse the discovery metrics (required anchor)."""
    ctrl = sel.check_controls_fire(f1["results"])
    assert ctrl["shuffle_control_present"] is True
    assert ctrl["shuffle_fired"] is True
    assert ctrl["at_least_one_control_fires"] is True


def test_tm6_extreme_theta_reported_not_gating(f1):
    """If no extreme-theta control is present, gating still keys on shuffle only.

    The extreme-theta control is REPORTED-NOT-GATING: its (non-)firing never
    drives at_least_one_control_fires. Here the fixture has no extreme-theta
    embedding, so it is simply absent -> firing stays driven by the shuffle."""
    ctrl = sel.check_controls_fire(f1["results"])
    # extreme-theta absent in the synthetic fixture; logged-only, never gating.
    assert ctrl["extreme_theta_control_present"] is False
    assert ctrl["extreme_theta_fired"] is None
    # at_least_one_control_fires is driven by the shuffle ONLY.
    assert ctrl["at_least_one_control_fires"] == ctrl["shuffle_fired"]


def test_tm6_extreme_theta_nonfiring_does_not_block(f1):
    """A present-but-non-firing extreme-theta control must NOT block firing.

    Synthesize an extreme-theta result that does NOT collapse (mimics harmonypy
    robustness on real data): at_least_one_control_fires stays True via shuffle."""
    results = list(f1["results"])
    # A non-collapsing extreme-theta clone of harmony (robust, doesn't fire).
    theta = dm.DiscoveryMetrics(method=sel.NEG_CONTROL_METHOD)
    h = f1["by"]["harmony"]
    theta.within_batch_distinctness = h.within_batch_distinctness
    theta.within_batch_isolation = h.within_batch_isolation
    theta.batch_mixing_score = h.batch_mixing_score
    results.append(theta)
    ctrl = sel.check_controls_fire(results)
    assert ctrl["extreme_theta_control_present"] is True
    assert ctrl["extreme_theta_fired"] is False  # robust -> did not collapse
    # Still fires overall because the shuffle anchor fired.
    assert ctrl["at_least_one_control_fires"] is True


# ==========================================================================
# TM-8 — CLI fail-loud guards + output emission
# ==========================================================================

def test_tm8_missing_batch_key_fails_loud(tmp_path):
    import anndata as ad
    f = fixtures.make_f1_batch_split()
    a = ad.AnnData(np.zeros((f.batch.shape[0], 1), dtype=np.float32))
    a.obsm["X_baseline"] = f.embeddings["baseline"].astype(np.float32)
    a.obs["sample"] = [str(b) for b in f.batch]
    h5 = tmp_path / "emb.h5ad"
    a.write_h5ad(h5)
    with pytest.raises(ValueError, match="absent from obs"):
        sel._load_embeddings_from_h5ad(h5, "NOT_A_COLUMN")


def test_tm8_missing_baseline_fails_loud(tmp_path):
    import anndata as ad
    f = fixtures.make_f1_batch_split()
    a = ad.AnnData(np.zeros((f.batch.shape[0], 1), dtype=np.float32))
    a.obsm["X_harmony"] = f.embeddings["harmony"].astype(np.float32)  # no baseline
    a.obs["sample"] = [str(b) for b in f.batch]
    h5 = tmp_path / "emb.h5ad"
    a.write_h5ad(h5)
    with pytest.raises(ValueError, match="no baseline embedding"):
        sel._load_embeddings_from_h5ad(h5, "sample")


def test_tm8_run_gate_missing_baseline_fails_loud(f1):
    emb = {"harmony": f1["fx"].embeddings["harmony"]}  # no baseline
    with pytest.raises(ValueError, match="baseline embedding"):
        sel.run_gate(emb, f1["fx"].batch)


def test_tm8_cli_emits_outputs(tmp_path):
    import anndata as ad
    f = fixtures.make_f1_batch_split()
    a = ad.AnnData(np.zeros((f.batch.shape[0], 1), dtype=np.float32))
    a.obsm["X_baseline"] = f.embeddings["baseline"].astype(np.float32)
    a.obsm["X_harmony"] = f.embeddings["harmony"].astype(np.float32)
    a.obsm["X_neg_control_shuffle_label"] = f.embeddings["neg_control_shuffle_label"].astype(np.float32)
    a.obs["sample"] = [str(b) for b in f.batch]
    h5 = tmp_path / "emb.h5ad"
    a.write_h5ad(h5)
    project_root = tmp_path / "proj"
    project_root.mkdir()
    rc = sel.main([
        "--project-root", str(project_root),
        "--batch-key", "sample",
        "--embeddings-h5ad", str(h5),
    ])
    assert rc == 0
    out = project_root / "discovery_gate"
    assert (out / "integration_recommendation.json").is_file()
    assert (out / "integration_scoreboard.csv").is_file()
    assert (out / "integration_audit.md").is_file()
    assert (out / "integration_audit.json").is_file()


def test_tm8_project_root_inside_factory_rejected():
    factory_root = Path(__file__).resolve().parents[1]
    with pytest.raises(ValueError, match="inside the factory"):
        sel.resolve_out_dir(factory_root)


# ==========================================================================
# TM-9 — cache layer: hit returns prior recommendation; miss triggers fresh eval
# ==========================================================================

def test_tm9_cache_hit_miss(tmp_path):
    baseline = np.random.default_rng(0).normal(size=(50, 4)).astype(np.float32)
    cache = ic.IntegrationCache(tmp_path / "cache")
    key = ic.make_cache_key(baseline, "sample")
    # Miss first.
    res = cache.lookup(key)
    assert res.hit is False
    # Store + hit.
    payload = {"cfg.batch.method": "harmony", "scvi_retrained": True}
    cache.store(key, payload)
    res2 = cache.lookup(key)
    assert res2.hit is True
    assert res2.payload["cfg.batch.method"] == "harmony"
    # force-refresh -> miss even with an entry present.
    res3 = cache.lookup(key, force_refresh=True)
    assert res3.hit is False


def test_tm9_cache_key_distinguishes_fast_and_skip(tmp_path):
    baseline = np.random.default_rng(0).normal(size=(50, 4)).astype(np.float32)
    full = ic.make_cache_key(baseline, "sample")
    fast = ic.make_cache_key(baseline, "sample", fast=True)
    skip = ic.make_cache_key(baseline, "sample", skip_scvi=True)
    digests = {full.digest(), fast.digest(), skip.digest()}
    assert len(digests) == 3  # a FAST/skip result never aliases a full one


def test_tm9_cache_key_changes_with_data(tmp_path):
    rng = np.random.default_rng(0)
    a = rng.normal(size=(50, 4)).astype(np.float32)
    b = rng.normal(size=(50, 4)).astype(np.float32)
    assert ic.make_cache_key(a, "sample").digest() != ic.make_cache_key(b, "sample").digest()
    # Same data, different batch_key -> different key.
    assert ic.make_cache_key(a, "sample").digest() != ic.make_cache_key(a, "donor").digest()


def test_tm9_loud_mode_banner_warns():
    notes = ic.loud_mode_banner(fast=True, skip_scvi=True, force_refresh=True)
    assert any("FAST" in n and "NOT FOR FINAL CLAIM" in n for n in notes)
    assert any("skip-scvi" in n for n in notes)
    assert any("force-refresh" in n for n in notes)


# ==========================================================================
# TM-11 — audit JSON records provenance
# ==========================================================================

def test_tm11_audit_records_provenance(f1):
    rec = sel.recommend(f1["results"])
    resolved = dm.assert_pinned_discovery_engine()
    payload = sel.build_recommendation_json(
        rec, f1["ref"], resolved, batch_key="sample", cache_key="abc:sample")
    assert payload["verdict_tier"] == "relative_ranking_only"
    assert "MARGIN_MIX" in payload["rule_constants"]
    assert "distinctness_floor" in payload["rule_constants"]
    assert "isolation_floor" in payload["rule_constants"]
    assert payload["leiden_params"]["resolution"] == 1.0
    assert payload["leiden_params"]["flavor"] == "leidenalg"
    # leidenalg + igraph versions present (COND #4 G1 extension).
    assert "leidenalg" in payload["pinned_versions"]
    assert "igraph" in payload["pinned_versions"]
    assert "harmonypy_direct_caveat" in payload
    assert payload["cache_key"] == "abc:sample"
    assert payload["correspondence_map"] is not None


def test_tm11_g1_extension_asserts_leiden_versions(monkeypatch):
    """COND #4: GATE G1 extension asserts leidenalg/igraph; fail-loud on drift."""
    real = dm._installed_version_any

    def fake(aliases):
        if "leidenalg" in aliases:
            return "0.0.0-bogus"
        return real(aliases)

    monkeypatch.setattr(dm, "_installed_version_any", fake)
    monkeypatch.delenv(dm.SC_ALLOW_SCIB_VERSION_DRIFT, raising=False)
    with pytest.raises(RuntimeError, match="Leiden-engine version gate FAILED"):
        dm.assert_pinned_discovery_engine()


def test_tm11_g1_extension_drift_opt_in(monkeypatch):
    real = dm._installed_version_any

    def fake(aliases):
        if "leidenalg" in aliases:
            return "0.0.0-bogus"
        return real(aliases)

    monkeypatch.setattr(dm, "_installed_version_any", fake)
    monkeypatch.setenv(dm.SC_ALLOW_SCIB_VERSION_DRIFT, "1")
    resolved = dm.assert_pinned_discovery_engine()
    assert resolved["leidenalg"] == "0.0.0-bogus"


# ==========================================================================
# TM-12 — existing 22 tests still pass (regression) — verified in a separate
# pytest invocation; here we assert the reviewed engine imports are intact.
# ==========================================================================

def test_tm12_reviewed_engine_imports_intact():
    from scripts.bench.integration import methods, score_integration
    # The discovery layer reused these symbols; their presence guards against an
    # accidental edit to the reviewed files.
    assert hasattr(score_integration, "band_aware_rank") is False  # lives in methods
    assert hasattr(methods, "band_aware_rank")
    assert hasattr(methods, "compute_harmonypy_direct")
    assert hasattr(methods, "compute_scvi_seed_sweep")
    assert hasattr(methods, "compute_shuffle_label_control")
    assert hasattr(methods, "compute_neg_control")
    assert hasattr(score_integration, "pre_register_over_correcting_threshold")
    assert hasattr(score_integration, "assert_pinned_metric_engine")


# ==========================================================================
# TM-13 (COND #5) — single-batch type excluded from D1a denominator
# ==========================================================================

def test_tm13_single_batch_excluded_from_d1a_denominator(f3):
    """The one-batch-only type has NO mutual cross-batch match -> excluded from
    D1a's denominator, and is held isolated by D1c."""
    ref = f3["ref"]
    fx = f3["fx"]
    bstr = np.asarray([str(b) for b in fx.batch])

    # Find the global cluster keys whose dominant true_bio is single_batch_only.
    single_keys = []
    for b in ref.batches:
        idx = np.where(bstr == b)[0]
        lab = ref.labels[b]
        for u in sorted(set(int(x) for x in lab)):
            rows = idx[lab == u]
            vals, counts = np.unique(fx.true_bio[rows], return_counts=True)
            if str(vals[int(np.argmax(counts))]) == "single_batch_only":
                single_keys.append(f"{b}::{u}")

    assert single_keys, "fixture must contain a single-batch-only cluster"
    # COND #5: none of these participate in any mutual cross-batch match.
    for k in single_keys:
        assert k not in ref.matched_keys
    # And they appear in NO correspondence pair (excluded from D1a denominator).
    flat = {k for pair in ref.correspondence for k in pair}
    for k in single_keys:
        assert k not in flat


def test_tm13_single_batch_held_isolated_by_d1c(f3):
    """The single-batch type's cluster stays well-isolated in a good integration."""
    by = f3["by"]
    # D1c (worst per-batch isolation) on the good integration stays above floor.
    rec = sel.recommend(f3["results"])
    iso_floor = rec.rule_constants["isolation_floor"]
    assert by["harmony"].within_batch_isolation >= iso_floor


# ==========================================================================
# TM-14 (COND #2) — swapping E does NOT change the frozen reference
# ==========================================================================

def test_tm14_reference_invariant_to_embedding(f1):
    """The frozen per-batch reference depends ONLY on the baseline X_pca; scoring
    a different embedding E never mutates it (anti-circularity)."""
    fx = f1["fx"]
    ref = dm.compute_per_batch_reference(fx.embeddings["baseline"], fx.batch)
    labels_before = {b: ref.labels[b].copy() for b in ref.batches}
    corr_before = list(ref.correspondence)
    # Score TWO different embeddings against the SAME ref object.
    _ = dm.score_embedding_discovery("harmony", fx.embeddings["harmony"], ref, fx.batch)
    _ = dm.score_embedding_discovery(
        "shuffle", fx.embeddings["neg_control_shuffle_label"], ref, fx.batch)
    # Reference partition + correspondence are unchanged.
    for b in ref.batches:
        assert np.array_equal(ref.labels[b], labels_before[b])
    assert ref.correspondence == corr_before


# ==========================================================================
# TM-15 (COND #3) — mild-over-merger (highest D1a, D1c<floor) REJECTED
# ==========================================================================

def test_tm15_mild_over_merger_rejected_for_aggressive_correct(f2):
    """A mild-over-merger with the HIGHEST D1a among candidates but D1c below
    floor is REJECTED in favor of an aggressive-correct method with lower D1a but
    D1c above floor (COND #3 disqualifying gate dominates the D1a bonus)."""
    fx = f2["fx"]
    ref = f2["ref"]
    bstr = np.asarray([str(b) for b in fx.batch])

    # mild_over_merger: good cross-batch merge (high D1a) BUT erases the rare
    # type (D1c collapses below floor) -> the F2 'erase' embedding is exactly this.
    mild = fx.embeddings["erase"].copy()

    # aggressive_correct: preserve the rare type (D1c high) but degrade ONE
    # corresponded cross-batch merge so its D1a is lower than the mild merger.
    aggr = fx.embeddings["harmony"].copy()
    ka, _kb = ref.correspondence[0]
    b, u = ka.split("::")
    idx = np.where(bstr == b)[0]
    rows = idx[ref.labels[b] == int(u)]
    aggr[rows, 2] += 8.0  # shift one cluster so it no longer co-localizes

    emb = {
        "baseline": fx.embeddings["baseline"],
        "mild_over_merger": mild,
        "aggressive_correct": aggr,
        "neg_control_shuffle_label": fx.embeddings["neg_control_shuffle_label"],
    }
    results, _ = dm.score_all_embeddings_discovery(emb, fx.embeddings["baseline"], fx.batch)
    by = {r.method: r for r in results}
    rec = sel.recommend(results)

    # mild_over_merger has the higher D1a of the two real candidates...
    assert by["mild_over_merger"].good_merge_fraction > by["aggressive_correct"].good_merge_fraction
    # ...yet it is DISQUALIFIED by the D1c isolation floor.
    assert rec.gate_results["mild_over_merger"]["survived_gates"] is False
    assert rec.gate_results["mild_over_merger"]["isolation_floor"] is False
    # ...and the aggressive-correct method (lower D1a, D1c ok) is chosen.
    assert rec.gate_results["aggressive_correct"]["survived_gates"] is True
    assert rec.chosen_method == "aggressive_correct"


# ==========================================================================
# Extra coverage: scoreboard CSV has the new per-batch columns
# ==========================================================================

def test_scoreboard_has_new_columns(f1, tmp_path):
    out = tmp_path / "sb.csv"
    sel.write_discovery_scoreboard(f1["results"], out)
    header = out.read_text().splitlines()[0].split(",")
    for col in ("good_merge_fraction", "within_batch_distinctness",
                "within_batch_isolation", "n_corresponded_pairs"):
        assert col in header
