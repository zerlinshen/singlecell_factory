"""CI tests for the integration calibration harness (T1 + T3).

Covers GATE G1 (scib-metrics default engine + version-assert fail-loud), the
synthetic-fixture native-vs-scib agreement, GATE G3 (verdict_tier), the
OVER-CORRECTING flag, and GATE G2 fail-loud label/join guards. scib-metrics is
required; tests are skipped (not silently passed) if the engine is unavailable.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

scib_metrics = pytest.importorskip("scib_metrics")

from scripts.bench.integration import score_integration as si  # noqa: E402

# ``tests`` is a pytest namespace dir (no top-level __init__.py), so importing
# ``tests.fixtures.integration_synthetic`` is unreliable across invocation
# modes. Load the generator by file path (same approach the harness uses).
def _load_make_fixture():
    import importlib.util

    fixture_path = ROOT / "tests" / "fixtures" / "integration_synthetic.py"
    mod_name = "integration_synthetic_fixture"
    if mod_name in sys.modules:
        return sys.modules[mod_name].make_integration_fixture
    spec = importlib.util.spec_from_file_location(mod_name, fixture_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod.make_integration_fixture


make_integration_fixture = _load_make_fixture()


@pytest.fixture(scope="module")
def fixture():
    return make_integration_fixture(seed=0)


@pytest.fixture(scope="module")
def scored(fixture):
    embeddings = {
        "unintegrated": fixture.X_unintegrated,
        "integrated": fixture.X_integrated,
    }
    results = si.score_embeddings(
        embeddings, fixture.batch, fixture.cell_type)
    return {r.method: r for r in results}


# --------------------------------------------------------------------------
# GATE G1: pinned engine + version-assert fail-loud
# --------------------------------------------------------------------------

def test_pinned_versions_match_installed():
    resolved = si.assert_pinned_metric_engine()
    for dist, pinned in si.PINNED_VERSIONS.items():
        assert resolved[dist] == pinned


def test_version_mismatch_fails_loud(monkeypatch):
    """A version mismatch must raise RuntimeError by default (GATE G1)."""
    real = si._installed_version

    def fake(dist):
        if dist == "scib-metrics":
            return "0.0.0-bogus"
        return real(dist)

    monkeypatch.setattr(si, "_installed_version", fake)
    monkeypatch.delenv(si.SC_ALLOW_SCIB_VERSION_DRIFT, raising=False)
    with pytest.raises(RuntimeError, match="version gate FAILED"):
        si.assert_pinned_metric_engine()


def test_version_drift_opt_in_does_not_raise(monkeypatch):
    """SC_ALLOW_SCIB_VERSION_DRIFT=1 downgrades the mismatch to a loud log."""
    real = si._installed_version

    def fake(dist):
        if dist == "scib-metrics":
            return "0.0.0-bogus"
        return real(dist)

    monkeypatch.setattr(si, "_installed_version", fake)
    monkeypatch.setenv(si.SC_ALLOW_SCIB_VERSION_DRIFT, "1")
    resolved = si.assert_pinned_metric_engine()
    assert resolved["scib-metrics"] == "0.0.0-bogus"


def test_runtime_versions_recorded_not_asserted():
    rec = si.record_resolved_runtime_versions()
    for key in ("jax", "numpy", "scikit-learn", "python"):
        assert key in rec
        assert rec[key] != ""


# --------------------------------------------------------------------------
# Synthetic fixture: directionality + bio preservation
# --------------------------------------------------------------------------

def test_integrated_outmixes_unintegrated(scored):
    """Integration should raise batch-mixing without sacrificing bio."""
    uni = scored["unintegrated"]
    intg = scored["integrated"]
    assert intg.batch_mixing_score > uni.batch_mixing_score
    # Bio conservation must be preserved (not collapsed by integration).
    assert intg.bio_conservation_score >= uni.bio_conservation_score - 0.05
    # All bio aggregates should remain high (well-separated blobs).
    assert uni.bio_conservation_score > 0.6
    assert intg.bio_conservation_score > 0.6


def test_all_metrics_populated(scored):
    """No metric should be missing on the clean synthetic fixture."""
    for r in scored.values():
        for attr in ("ilisi", "kbet_acceptance", "batch_asw",
                     "clisi", "celltype_asw", "ari", "nmi"):
            assert getattr(r, attr) is not None, f"{r.method}.{attr} missing"
        assert not r.notes, f"unexpected metric notes for {r.method}: {r.notes}"


# --------------------------------------------------------------------------
# GATE G1: native (sklearn) cross-check AGREES with scib-metrics default
# --------------------------------------------------------------------------

def test_native_scib_agreement_within_tolerance(scored):
    """ASW agrees to ~1e-6; ARI/NMI agree exactly on the fixture.

    The native sklearn cross-check is NOT the default engine — this test
    confirms the default (scib-metrics) and the independent native estimate
    converge on the same numbers, so the default is trustworthy.
    """
    for r in scored.values():
        agree = si.native_scib_agreement(r)
        assert agree.get("celltype_asw_abs_diff", 0.0) < 1e-4
        assert agree.get("batch_asw_abs_diff", 0.0) < 1e-4
        assert agree.get("ari_abs_diff", 1.0) < 1e-6
        assert agree.get("nmi_abs_diff", 1.0) < 1e-6


def test_self_test_writes_outputs(tmp_path):
    """Full self-test produces scoreboard CSV + audit json/md."""
    audit = si.run_self_test(tmp_path, seed=0)
    assert (tmp_path / "integration_scoreboard.csv").is_file()
    assert (tmp_path / "integration_audit.json").is_file()
    assert (tmp_path / "integration_audit.md").is_file()
    assert audit["lane"] == "calibration_relative_only"
    assert audit["engine"]["default"] == "scib-metrics"
    # GATE G3: every audit carries a verdict_tier (default contrast => weakest).
    assert audit["verdict_tier"] == "relative_ranking_only"
    # CSV header shape is stable.
    header = (tmp_path / "integration_scoreboard.csv").read_text().splitlines()[0]
    assert header.split(",") == si._SCOREBOARD_COLUMNS


# ==========================================================================
# GATE G3: verdict_tier mapping (enforced)
# ==========================================================================

@pytest.mark.parametrize("contrast_type", [
    "global", "sample_level", "pipeline_derived_labels",
])
def test_sample_level_is_relative_ranking_only(contrast_type):
    spec = si.ContrastSpec(contrast_type=contrast_type, labels_source="author")
    assert si.classify_verdict_tier(spec) == "relative_ranking_only"


def test_within_age_pair_is_conditioned_concordance():
    spec = si.ContrastSpec(
        contrast_type="within_age_replicate_pair", labels_source="author")
    assert si.classify_verdict_tier(spec) == "conditioned_concordance"


def test_within_age_pair_with_pipeline_labels_downgrades():
    """Concordance with your OWN clustering output is not concordance."""
    spec = si.ContrastSpec(
        contrast_type="within_age_replicate_pair", labels_source="pipeline")
    assert si.classify_verdict_tier(spec) == "relative_ranking_only"


def test_unknown_contrast_type_fails_loud():
    with pytest.raises(ValueError, match="unknown contrast_type"):
        si.classify_verdict_tier(si.ContrastSpec(contrast_type="bogus"))


def test_no_tier_exceeds_conditioned_concordance():
    """No mapping value may claim public concordance (GATE G6 deferral)."""
    assert set(si.VERDICT_TIER_MAP.values()) <= {
        "relative_ranking_only", "conditioned_concordance"}


# ==========================================================================
# OVER-CORRECTING flag (conservation floor, pre-registered T)
# ==========================================================================

def _mk(method, bio, mixing=0.9):
    r = si.MetricResult(method=method)
    r.bio_conservation_score = bio
    r.batch_mixing_score = mixing
    return r


def test_over_correcting_threshold_from_baseline_neg_control_spread():
    """T is pre-registered from baseline - neg_control bio spread."""
    T, prov = si.pre_register_over_correcting_threshold(0.90, 0.40)
    assert T == pytest.approx(0.50)
    assert "baseline_bio" in prov


def test_over_correcting_threshold_floor_without_neg_control():
    T, prov = si.pre_register_over_correcting_threshold(0.90, None)
    assert T == si.OVER_CORRECTING_THRESHOLD_FLOOR
    assert "floor" in prov


def test_over_correcting_flag_fires_on_bio_collapse():
    """A high-mixing method that shreds biology must be flagged."""
    results = [
        _mk("baseline", bio=0.90, mixing=0.10),
        _mk("neg_control", bio=0.40, mixing=0.99),
        # corrected method mixes great (0.99) but bio crashes to 0.30:
        # drop = 0.90 - 0.30 = 0.60 > T(=0.50). Must flag regardless of mixing.
        _mk("aggressive", bio=0.30, mixing=0.99),
        # well-behaved method: small bio drop, NOT flagged.
        _mk("good", bio=0.88, mixing=0.85),
    ]
    flags, audit = si.flag_over_correcting(
        results, baseline_method="baseline",
        neg_control_method="neg_control")
    assert audit["over_correcting_threshold"] == pytest.approx(0.50)
    assert any("aggressive" in f for f in flags)
    assert not any("good" in f for f in flags)


def test_over_correcting_on_synthetic_fixture_does_not_overflag(scored):
    """Real fixture: integrated PRESERVES bio, so no over-correcting flag."""
    results = [
        si.MetricResult(method="baseline",
                        bio_conservation_score=scored["unintegrated"].bio_conservation_score,
                        batch_mixing_score=scored["unintegrated"].batch_mixing_score),
        si.MetricResult(method="integrated",
                        bio_conservation_score=scored["integrated"].bio_conservation_score,
                        batch_mixing_score=scored["integrated"].batch_mixing_score),
    ]
    flags, _ = si.flag_over_correcting(results, baseline_method="baseline")
    assert flags == []


# ==========================================================================
# GATE G2: fail-loud on missing label column / dead-symlink join
# ==========================================================================

def test_g2_missing_label_column_raises():
    with pytest.raises(si.BioConservationLabelError, match="absent"):
        si.require_bio_conservation_labels(
            obs_columns=["leiden", "batch"], label_column="cell_type")


def test_g2_present_label_column_ok():
    si.require_bio_conservation_labels(
        obs_columns=["leiden", "batch", "cell_type"], label_column="cell_type")


def test_g2_dead_symlink_join_raises(tmp_path):
    dead = tmp_path / "cluster_names.tsv"
    dead.symlink_to(tmp_path / "does_not_exist.tsv")
    assert dead.is_symlink() and not dead.exists()
    with pytest.raises(si.BioConservationLabelError, match="DEAD SYMLINK"):
        si.require_bio_conservation_labels(
            obs_columns=["cell_type"], label_column="cell_type",
            cluster_names_path=dead)


def test_g2_empty_join_map_raises():
    with pytest.raises(si.BioConservationLabelError, match="empty"):
        si.require_bio_conservation_labels(
            obs_columns=["cell_type"], label_column="cell_type",
            cluster_names_join={})


def test_g2_missing_annotation_file_raises(tmp_path):
    with pytest.raises(si.BioConservationLabelError, match="does not exist"):
        si.require_bio_conservation_labels(
            obs_columns=["cell_type"], label_column="cell_type",
            cluster_names_path=tmp_path / "nope.tsv")
