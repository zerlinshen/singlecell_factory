"""Contract tests for plan-time multi-batch detection and claim downgrade.

Guards the defect fixed on 2026-08-03: the default optional-module set carries
no integration step, so a default run on multi-batch input produced batch-driven
clusters and said nothing about it. Every test here is written so that deleting
the feature makes it fail — the negative cases assert the WARNING and the
DOWNGRADE, not merely that the code runs.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from workflow.modular.batch_risk import (
    BATCH_CANDIDATE_COLUMNS,
    BATCH_STRATEGY_ACCEPT_UNCORRECTED,
    BATCH_STRATEGY_AUTO,
    BATCH_STRATEGY_INTEGRATE,
    BATCH_STRATEGY_SINGLE_BATCH,
    CLAIM_STATUS_CLAIMABLE,
    CLAIM_STATUS_EXPLORATORY,
    CLAIM_STATUS_UNDETERMINED,
    DETECTION_MULTI_BATCH,
    DETECTION_SINGLE_BATCH,
    DETECTION_UNAVAILABLE,
    OBSERVATION_MULTI_BATCH,
    OBSERVATION_SINGLE_BATCH,
    OBSERVATION_UNOBSERVED,
    BatchStrategyConflict,
    announce_batch_risk,
    batch_risk_warning,
    detect_batch_structure,
    observe_batch_structure,
    reconcile_batch_risk,
    reconciliation_warning,
    resolve_batch_risk,
    resolve_obs_source,
)
from workflow.modular.module_catalog import (
    ANALYSIS_PROFILES,
    DEFAULT_OPTIONAL_MODULES,
    MANDATORY_MODULES,
    analysis_profile,
    analysis_profile_names,
    optional_modules_for_profile,
)
from workflow.modular.pipeline import _resolve_execution_order

DEFAULT_PLAN = list(DEFAULT_OPTIONAL_MODULES)


def _write_h5ad(path: Path, obs: pd.DataFrame) -> Path:
    adata = ad.AnnData(
        X=np.zeros((len(obs), 4), dtype="float32"),
        obs=obs,
        var=pd.DataFrame(index=[f"g{i}" for i in range(4)]),
    )
    adata.write_h5ad(path)
    return path


def _multi_batch_obs(n: int = 60, n_samples: int = 6) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample": pd.Categorical([f"S{i % n_samples}" for i in range(n)]),
            "dataset": [f"D{i % 3}" for i in range(n)],
        },
        index=[f"c{i}" for i in range(n)],
    )


def _single_batch_obs(n: int = 60) -> pd.DataFrame:
    return pd.DataFrame(
        {"sample": pd.Categorical(["S0"] * n)},
        index=[f"c{i}" for i in range(n)],
    )


# ---------------------------------------------------------------------------
# detection
# ---------------------------------------------------------------------------


def test_detects_multi_batch_from_obs_only(tmp_path):
    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs(n_samples=6))
    detection = detect_batch_structure(path, preferred_key="sample")
    assert detection["detection_status"] == DETECTION_MULTI_BATCH
    assert detection["detected_batches"] == 6
    assert detection["key"] == "sample"
    assert detection["candidate_batch_columns"] == {"dataset": 3, "sample": 6}


def test_single_level_column_is_not_a_batch(tmp_path):
    path = _write_h5ad(tmp_path / "single.h5ad", _single_batch_obs())
    detection = detect_batch_structure(path, preferred_key="sample")
    assert detection["detection_status"] == DETECTION_SINGLE_BATCH
    assert detection["detected_batches"] == 1
    assert detection["candidate_batch_columns"] == {}


def test_unused_categories_are_not_counted_as_batches(tmp_path):
    """A subsetted cohort keeps unused categories; those are not present batches."""
    obs = pd.DataFrame(
        {
            "sample": pd.Categorical(
                ["S0"] * 30, categories=["S0", "S1", "S2", "S3"]
            )
        },
        index=[f"c{i}" for i in range(30)],
    )
    path = _write_h5ad(tmp_path / "subset.h5ad", obs)
    detection = detect_batch_structure(path, preferred_key="sample")
    assert detection["detected_batches"] == 1, (
        "declared-but-unused categories must not be reported as batches"
    )
    assert detection["detection_status"] == DETECTION_SINGLE_BATCH


def test_preferred_key_wins_over_wider_column(tmp_path):
    obs = _multi_batch_obs(n=60, n_samples=3)
    obs["donor_id"] = [f"P{i % 10}" for i in range(60)]
    path = _write_h5ad(tmp_path / "pref.h5ad", obs)
    assert detect_batch_structure(path, preferred_key="sample")["key"] == "sample"
    # Without the preference, the widest split is reported.
    assert detect_batch_structure(path, preferred_key="")["key"] == "donor_id"


def test_detection_unavailable_without_prepared_input(tmp_path):
    detection = detect_batch_structure(resolve_obs_source(sample_root=tmp_path))
    assert detection["detection_status"] == DETECTION_UNAVAILABLE
    assert detection["detected_batches"] == 0


def test_resolve_obs_source_matches_cellranger_resolution(tmp_path):
    sample_root = tmp_path / "cohort"
    sample_root.mkdir()
    assert resolve_obs_source(sample_root=sample_root) is None
    prepared = _write_h5ad(sample_root / "prepared_input.h5ad", _single_batch_obs())
    assert resolve_obs_source(sample_root=sample_root) == prepared
    direct = _write_h5ad(tmp_path / "explicit.h5ad", _single_batch_obs())
    assert resolve_obs_source(input_h5ad=direct, sample_root=sample_root) == direct


def test_detection_never_raises_on_a_corrupt_input(tmp_path):
    broken = tmp_path / "broken.h5ad"
    broken.write_bytes(b"not an hdf5 file")
    detection = detect_batch_structure(broken)
    assert detection["detection_status"] == DETECTION_UNAVAILABLE
    assert detection["detection_reason"]


# ---------------------------------------------------------------------------
# detection FAILURE must never become an affirmative single-batch claim
# ---------------------------------------------------------------------------


def _break_categorical_column(path: Path, column: str) -> None:
    """Make one obs column present-but-unreadable, without mocking anything.

    Deletes the ``codes`` child of an AnnData categorical, which is what a
    truncated file, an unsupported encoding, or a partially-written object looks
    like to the reader: the column is there, and it cannot be read.
    """
    import h5py

    with h5py.File(str(path), "r+") as handle:
        del handle["obs"][column]["codes"]


def test_unreadable_batch_column_is_not_reported_as_single_batch(tmp_path):
    """HIGH-1: the reader KNOWS it failed; it must not claim success.

    A genuine 6-batch cohort whose `sample` column cannot be read used to
    collapse into `single_batch`/`detected_batches: 1`, granting a silent
    `claimable` with no warning at all.
    """
    path = _write_h5ad(tmp_path / "broken.h5ad", _multi_batch_obs(n_samples=6)[["sample"]])
    _break_categorical_column(path, "sample")
    detection = detect_batch_structure(path, preferred_key="sample")

    assert detection["detection_status"] != DETECTION_SINGLE_BATCH, (
        "an unreadable batch column was reported as an affirmative single batch"
    )
    assert detection["detection_status"] == DETECTION_UNAVAILABLE
    assert detection["unreadable_columns"] == ["sample"]

    risk = resolve_batch_risk(
        detection=detection,
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    assert risk["clustering_claim_status"] != CLAIM_STATUS_CLAIMABLE
    assert risk["clustering_claim_reason"] != "single_batch_input"
    assert batch_risk_warning(risk) != "", "a detection failure must not be silent"


def test_read_error_on_a_batch_column_is_not_single_batch(tmp_path, monkeypatch):
    """The raising branch of the same defect (NFS/filter/IO errors)."""
    import workflow.modular.batch_risk as batch_risk_mod

    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs(n_samples=6))

    def exploding_count(dataset, *, drop_negative):
        raise OSError("Can't read data (filter returned failure during read)")

    monkeypatch.setattr(batch_risk_mod, "_count_levels", exploding_count)
    detection = detect_batch_structure(path, preferred_key="sample")
    assert detection["detection_status"] == DETECTION_UNAVAILABLE
    assert sorted(detection["unreadable_columns"]) == ["dataset", "sample"]


def test_partially_unreadable_candidates_do_not_yield_single_batch(tmp_path):
    """One column read fine and is single-level, another failed to read.

    The failed one could have been the multi-level axis, so an affirmative
    single-batch answer is not supported by what was actually read.
    """
    obs = pd.DataFrame(
        {
            "sample": pd.Categorical(["S0"] * 30),
            "donor_id": pd.Categorical([f"D{i % 4}" for i in range(30)]),
        },
        index=[f"c{i}" for i in range(30)],
    )
    path = _write_h5ad(tmp_path / "partial.h5ad", obs)
    _break_categorical_column(path, "donor_id")
    detection = detect_batch_structure(path, preferred_key="sample")
    assert detection["detection_status"] == DETECTION_UNAVAILABLE
    assert detection["unreadable_columns"] == ["donor_id"]


def test_no_candidate_column_at_all_is_not_single_batch(tmp_path):
    """HIGH-2: a batch axis the heuristic does not know about.

    A 5-level `library_id` with --batch-key left at the default `sample` used to
    read as `single_batch` -> `claimable`. A heuristic may have false negatives;
    it may not report one as an affirmative positive claim.
    """
    obs = pd.DataFrame(
        {"library_id": pd.Categorical([f"L{i % 5}" for i in range(40)])},
        index=[f"c{i}" for i in range(40)],
    )
    path = _write_h5ad(tmp_path / "library.h5ad", obs)
    detection = detect_batch_structure(path, preferred_key="sample")
    assert detection["detection_status"] != DETECTION_SINGLE_BATCH
    assert detection["detection_status"] == DETECTION_UNAVAILABLE
    assert detection["candidate_columns_read"] == []

    risk = resolve_batch_risk(
        detection=detection,
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    assert risk["clustering_claim_status"] != CLAIM_STATUS_CLAIMABLE


def test_orig_ident_is_a_recognised_batch_column(tmp_path):
    """Seurat's standard sample column, and this suite's own R default.

    r_multiomics_factory/R/pipeline_steps.R:18-19 resolves `orig.ident` as the
    batch column when present; plotting_factory/r/general/qc_plots.R:62 uses it
    as the default group.by. The Python detector was blind to it.
    """
    assert "orig.ident" in BATCH_CANDIDATE_COLUMNS
    obs = pd.DataFrame(
        {"orig.ident": pd.Categorical([f"P{i % 5}" for i in range(40)])},
        index=[f"c{i}" for i in range(40)],
    )
    path = _write_h5ad(tmp_path / "seurat.h5ad", obs)
    detection = detect_batch_structure(path, preferred_key="sample")
    assert detection["detection_status"] == DETECTION_MULTI_BATCH
    assert detection["detected_batches"] == 5
    assert detection["key"] == "orig.ident"


def test_batch_column_matching_is_case_insensitive(tmp_path):
    obs = pd.DataFrame(
        {"Sample": pd.Categorical([f"S{i % 4}" for i in range(40)])},
        index=[f"c{i}" for i in range(40)],
    )
    path = _write_h5ad(tmp_path / "capital.h5ad", obs)
    detection = detect_batch_structure(path, preferred_key="sample")
    assert detection["detection_status"] == DETECTION_MULTI_BATCH
    assert detection["detected_batches"] == 4
    assert detection["key"] == "Sample"


def test_runtime_observation_shares_the_same_failure_policy():
    """observe_batch_structure must not grant single_batch on zero candidates.

    Plan-time and runtime deliberately share BATCH_CANDIDATE_COLUMNS, so a blind
    spot in one is a blind spot in both; reconciliation would otherwise confirm
    a wrong answer with a second wrong answer.
    """
    obs = pd.DataFrame(
        {"library_id": pd.Categorical([f"L{i % 5}" for i in range(20)])},
        index=[f"c{i}" for i in range(20)],
    )
    observation = observe_batch_structure(obs)
    assert observation["observation_status"] != OBSERVATION_SINGLE_BATCH
    assert observation["observation_status"] == OBSERVATION_UNOBSERVED


def test_runtime_observation_sees_orig_ident_and_mixed_case():
    for column in ("orig.ident", "Batch"):
        obs = pd.DataFrame(
            {column: pd.Categorical([f"X{i % 3}" for i in range(20)])},
            index=[f"c{i}" for i in range(20)],
        )
        observation = observe_batch_structure(obs)
        assert observation["observation_status"] == OBSERVATION_MULTI_BATCH, column
        assert observation["observed_batches"] == 3, column


def test_genuine_single_batch_still_requires_a_readable_candidate(tmp_path):
    """The one case that IS a real single batch stays claimable."""
    path = _write_h5ad(tmp_path / "single.h5ad", _single_batch_obs())
    detection = detect_batch_structure(path, preferred_key="sample")
    assert detection["detection_status"] == DETECTION_SINGLE_BATCH
    assert detection["candidate_columns_read"] == ["sample"]
    risk = resolve_batch_risk(
        detection=detection,
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    assert risk["clustering_claim_status"] == CLAIM_STATUS_CLAIMABLE


def test_known_read_failure_blocks_a_declared_single_batch_claim(tmp_path):
    """"We could not read it" must not be upgraded by an assertion.

    Distinct from "there is no batch annotation at all", where an affirmative
    declaration is still the operator's to make.
    """
    path = _write_h5ad(tmp_path / "broken.h5ad", _multi_batch_obs(n_samples=6)[["sample"]])
    _break_categorical_column(path, "sample")
    risk = resolve_batch_risk(
        detection=detect_batch_structure(path, preferred_key="sample"),
        declared_strategy=BATCH_STRATEGY_SINGLE_BATCH,
        planned_modules=DEFAULT_PLAN,
    )
    assert risk["clustering_claim_status"] == CLAIM_STATUS_UNDETERMINED
    assert risk["clustering_claim_reason"] == "batch_structure_unreadable_before_run"


# ---------------------------------------------------------------------------
# claim resolution — the negative cases
# ---------------------------------------------------------------------------


def test_multi_batch_default_run_is_downgraded_to_exploratory(tmp_path):
    """THE defect: defaults on multi-batch input must not read as reliable."""
    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs())
    risk = resolve_batch_risk(
        detection=detect_batch_structure(path, preferred_key="sample"),
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    assert risk["clustering_claim_status"] == CLAIM_STATUS_EXPLORATORY
    assert risk["clustering_claim_reason"] == (
        "exploratory_nonclaimable_unintegrated_multi_batch"
    )
    assert risk["strategy_declared"] is False
    assert risk["integration_modules_planned"] == []
    assert risk["downgraded_modules"] == [
        "clustering",
        "differential_expression",
        "annotation",
    ]


def test_multi_batch_default_run_emits_an_actionable_warning(tmp_path):
    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs())
    risk = resolve_batch_risk(
        detection=detect_batch_structure(path, preferred_key="sample"),
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    warning = batch_risk_warning(risk)
    assert warning.startswith("BATCH_RISK:")
    assert "6 batches" in warning
    assert "batch_correction" in warning, "warning must name the remedy"
    assert "EXPLORATORY" in warning


def test_announce_writes_to_stderr_exactly_once(tmp_path, capsys):
    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs())
    risk = resolve_batch_risk(
        detection=detect_batch_structure(path, preferred_key="sample"),
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    announce_batch_risk(risk)
    captured = capsys.readouterr()
    assert captured.err.count("BATCH_RISK:") == 1
    assert "BATCH_RISK:" not in captured.out


def test_single_batch_input_is_claimable_and_silent(tmp_path):
    path = _write_h5ad(tmp_path / "single.h5ad", _single_batch_obs())
    risk = resolve_batch_risk(
        detection=detect_batch_structure(path, preferred_key="sample"),
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    assert risk["clustering_claim_status"] == CLAIM_STATUS_CLAIMABLE
    assert risk["clustering_claim_reason"] == "single_batch_input"
    assert risk["downgraded_modules"] == []
    assert batch_risk_warning(risk) == ""


def test_planned_batch_correction_clears_the_downgrade(tmp_path):
    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs())
    risk = resolve_batch_risk(
        detection=detect_batch_structure(path, preferred_key="sample"),
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=[*DEFAULT_PLAN, "batch_correction"],
    )
    assert risk["clustering_claim_status"] == CLAIM_STATUS_CLAIMABLE
    assert risk["strategy"] == BATCH_STRATEGY_INTEGRATE
    assert risk["strategy_declared"] is True
    assert risk["strategy_source"] == "planned_modules"
    assert batch_risk_warning(risk) == ""


def test_integration_select_alone_does_not_count_as_integration(tmp_path):
    """integration_select only SCORES candidates; it corrects nothing."""
    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs())
    risk = resolve_batch_risk(
        detection=detect_batch_structure(path, preferred_key="sample"),
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=[*DEFAULT_PLAN, "integration_select"],
    )
    assert risk["clustering_claim_status"] == CLAIM_STATUS_EXPLORATORY


def test_accept_uncorrected_declares_but_stays_exploratory(tmp_path):
    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs())
    risk = resolve_batch_risk(
        detection=detect_batch_structure(path, preferred_key="sample"),
        declared_strategy=BATCH_STRATEGY_ACCEPT_UNCORRECTED,
        planned_modules=DEFAULT_PLAN,
    )
    assert risk["strategy_declared"] is True
    assert risk["clustering_claim_status"] == CLAIM_STATUS_EXPLORATORY, (
        "declaring an uncorrected multi-batch run must not make it claimable"
    )
    assert batch_risk_warning(risk).startswith("BATCH_RISK:")


def test_false_single_batch_declaration_fails_loudly(tmp_path):
    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs())
    with pytest.raises(BatchStrategyConflict) as excinfo:
        resolve_batch_risk(
            detection=detect_batch_structure(path, preferred_key="sample"),
            declared_strategy=BATCH_STRATEGY_SINGLE_BATCH,
            planned_modules=DEFAULT_PLAN,
        )
    assert "single-batch" in str(excinfo.value)
    assert "6 levels" in str(excinfo.value)


def test_truthful_single_batch_declaration_is_accepted(tmp_path):
    path = _write_h5ad(tmp_path / "single.h5ad", _single_batch_obs())
    risk = resolve_batch_risk(
        detection=detect_batch_structure(path, preferred_key="sample"),
        declared_strategy=BATCH_STRATEGY_SINGLE_BATCH,
        planned_modules=DEFAULT_PLAN,
    )
    assert risk["clustering_claim_status"] == CLAIM_STATUS_CLAIMABLE
    assert risk["strategy_declared"] is True


def test_integrate_declaration_without_an_integration_module_fails(tmp_path):
    path = _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs())
    with pytest.raises(BatchStrategyConflict) as excinfo:
        resolve_batch_risk(
            detection=detect_batch_structure(path, preferred_key="sample"),
            declared_strategy=BATCH_STRATEGY_INTEGRATE,
            planned_modules=DEFAULT_PLAN,
        )
    assert "batch_correction" in str(excinfo.value)


def test_undetectable_input_is_undetermined_not_claimable():
    risk = resolve_batch_risk(
        detection=detect_batch_structure(None),
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    assert risk["clustering_claim_status"] == CLAIM_STATUS_UNDETERMINED, (
        "an unreadable input must not silently produce a claimable run"
    )
    assert "undetermined" in batch_risk_warning(risk)


def test_unknown_strategy_is_rejected():
    with pytest.raises(BatchStrategyConflict):
        resolve_batch_risk(
            detection=detect_batch_structure(None),
            declared_strategy="whatever",
            planned_modules=DEFAULT_PLAN,
        )


# ---------------------------------------------------------------------------
# catalog profiles + ordering
# ---------------------------------------------------------------------------


def test_named_profiles_exist_and_declare_a_strategy():
    assert set(analysis_profile_names()) == {"single_batch", "multi_batch_harmony"}
    assert analysis_profile("single_batch").batch_strategy == (
        BATCH_STRATEGY_SINGLE_BATCH
    )
    assert analysis_profile("multi_batch_harmony").batch_strategy == (
        BATCH_STRATEGY_INTEGRATE
    )
    with pytest.raises(ValueError):
        analysis_profile("no_such_profile")


def test_single_batch_profile_is_the_canonical_default_set():
    assert optional_modules_for_profile("single_batch") == DEFAULT_OPTIONAL_MODULES


def test_multi_batch_profile_is_defaults_plus_batch_correction():
    modules = optional_modules_for_profile("multi_batch_harmony")
    assert set(modules) == set(DEFAULT_OPTIONAL_MODULES) | {"batch_correction"}


def test_every_profile_names_only_catalogued_modules():
    from workflow.modular.module_catalog import MODULE_SPECS

    for profile in ANALYSIS_PROFILES.values():
        unknown = set(profile.optional_modules) - set(MODULE_SPECS)
        assert not unknown, f"{profile.name} names uncatalogued modules: {unknown}"


def test_leiden_consumers_declare_runs_after_batch_correction():
    """batch_correction declares `provides obs.leiden` and overwrites the labels.

    Every leiden consumer must therefore carry the ordering hint, or its markers
    and labels describe clusters that no longer exist in final_adata. This is
    the declared contract; the resolved-order tests below exercise it.
    """
    from workflow.modular.module_catalog import MODULE_SPECS, module_runs_after

    for name in ("differential_expression", "annotation"):
        assert "batch_correction" in MODULE_SPECS[name].runs_after, (
            f"{name} consumes obs['leiden'] but does not declare runs_after "
            "batch_correction"
        )
        assert "batch_correction" in module_runs_after()[name]


def test_de_is_sequenced_after_batch_correction_under_select_integration():
    """The plan where the hint is load-bearing rather than incidental.

    With --select-integration in the plan, batch_correction gains an inbound
    ordering edge from integration_select, which delays it past the point where
    a hint-less differential_expression is released. Measured on 2026-08-03
    without the hint: clustering, differential_expression, integration_select,
    batch_correction, annotation — markers computed on pre-correction leiden,
    then silently overwritten. Deleting the hint makes this test fail.
    """
    order = _resolve_execution_order(
        list(MANDATORY_MODULES),
        [
            "clustering",
            "integration_select",
            "batch_correction",
            "differential_expression",
            "annotation",
        ],
    )
    correction = order.index("batch_correction")
    assert order.index("differential_expression") > correction
    assert order.index("annotation") > correction


def test_multi_batch_profile_order_keeps_leiden_consumers_last():
    order = _resolve_execution_order(
        list(MANDATORY_MODULES),
        list(optional_modules_for_profile("multi_batch_harmony")),
    )
    correction = order.index("batch_correction")
    assert order.index("differential_expression") > correction
    assert order.index("annotation") > correction


def test_runs_after_never_pulls_batch_correction_into_a_plan():
    order = _resolve_execution_order(list(MANDATORY_MODULES), DEFAULT_PLAN)
    assert "batch_correction" not in order


def test_clustering_shares_the_candidate_columns_with_plan_time_detection():
    from workflow.modular.modules.clustering import ClusteringModule

    assert ClusteringModule._BATCH_CANDIDATE_COLUMNS is BATCH_CANDIDATE_COLUMNS


# ---------------------------------------------------------------------------
# manifest wiring
# ---------------------------------------------------------------------------


def test_batch_risk_reaches_the_governed_manifest(tmp_path):
    from workflow.modular.manifest_writer import write_manifest

    risk = resolve_batch_risk(
        detection=detect_batch_structure(
            _write_h5ad(tmp_path / "multi.h5ad", _multi_batch_obs()),
            preferred_key="sample",
        ),
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    out = write_manifest(
        tmp_path / "run",
        project_id="p",
        run_id="2026-08-03T0000Z-abcdef0",
        modules_run=["clustering"],
        batch_risk=risk,
    )
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["batch_risk"]["clustering_claim_status"] == (
        CLAIM_STATUS_EXPLORATORY
    )
    assert payload["batch_risk"]["detected_batches"] == 6
    assert payload["batch_risk"]["strategy_declared"] is False


def test_manifest_batch_risk_is_always_an_object(tmp_path):
    from workflow.modular.manifest_writer import write_manifest

    out = write_manifest(
        tmp_path / "run",
        project_id="p",
        run_id="2026-08-03T0000Z-abcdef0",
        modules_run=["clustering"],
    )
    assert json.loads(out.read_text(encoding="utf-8"))["batch_risk"] == {}


def _run_pipeline_with_stub_modules(tmp_path, *, obs, batch_strategy, optional_modules,
                                    module_statuses=(), metadata=None):
    """Drive run_pipeline with compute stubbed out, installing ``obs`` as the run's
    AnnData. The sample root deliberately has NO prepared_input.h5ad, so plan-time
    detection is UNAVAILABLE — the normal state of any run that starts from raw.
    """
    from workflow.modular import pipeline as pipeline_mod
    from workflow.modular.config import CellRangerConfig, PipelineConfig

    sample_root = tmp_path / "cohort"
    sample_root.mkdir(exist_ok=True)
    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(
            sample_root=sample_root, outs_dir=sample_root / "outs"
        ),
        optional_modules=list(optional_modules),
        batch_strategy=batch_strategy,
    )

    def fake_run_sequential(modules, registry, ctx, mandatory):
        ctx.adata = ad.AnnData(
            X=np.zeros((len(obs), 3), dtype="float32"),
            obs=obs,
            var=pd.DataFrame(index=["g0", "g1", "g2"]),
        )
        ctx.metadata.update(metadata or {})
        for name, status in module_statuses:
            ctx.status(name, status, "stub")

    original = pipeline_mod._run_sequential
    pipeline_mod._run_sequential = fake_run_sequential
    try:
        manifest_path = pipeline_mod.run_pipeline(cfg)
    finally:
        pipeline_mod._run_sequential = original
    return json.loads(Path(manifest_path).read_text(encoding="utf-8"))


def test_falsified_single_batch_declaration_is_downgraded_at_runtime(tmp_path):
    """A declaration that plan-time could not check must not survive the data.

    Detection is UNAVAILABLE by construction for any run starting from a raw
    sample root, so `--batch-strategy single-batch` was granted `claimable` on
    the operator's word alone. If the loaded object then carries several batch
    levels, the authoritative claim must be downgraded, not left claimable.
    """
    payload = _run_pipeline_with_stub_modules(
        tmp_path,
        obs=_multi_batch_obs(n=40, n_samples=6),
        batch_strategy=BATCH_STRATEGY_SINGLE_BATCH,
        optional_modules=DEFAULT_PLAN,
    )
    risk = payload["metadata"]["batch_risk"]
    assert risk["clustering_claim_status"] != CLAIM_STATUS_CLAIMABLE, (
        "an unverified single-batch declaration was falsified by the loaded "
        "data, but the authoritative claim still reads claimable"
    )
    assert risk["clustering_claim_status"] == CLAIM_STATUS_EXPLORATORY
    assert risk["clustering_claim_reason"] == (
        "exploratory_nonclaimable_single_batch_declaration_falsified"
    )
    assert risk["observed_batches"] == 6
    assert risk["downgraded_modules"] == [
        "clustering", "differential_expression", "annotation",
    ]
    # The plan-time reading is preserved rather than overwritten: the
    # disagreement between the two timepoints is itself evidence.
    assert risk["plan_time_claim"]["clustering_claim_status"] == (
        CLAIM_STATUS_CLAIMABLE
    )
    assert risk["claim_basis"] == "runtime_observation"


def test_verified_single_batch_declaration_survives_reconciliation(tmp_path):
    payload = _run_pipeline_with_stub_modules(
        tmp_path,
        obs=_single_batch_obs(n=40),
        batch_strategy=BATCH_STRATEGY_SINGLE_BATCH,
        optional_modules=DEFAULT_PLAN,
    )
    risk = payload["metadata"]["batch_risk"]
    assert risk["clustering_claim_status"] == CLAIM_STATUS_CLAIMABLE
    assert risk["clustering_claim_reason"] == "single_batch_input_observed"
    assert risk["claim_basis"] == "runtime_observation"
    assert risk["downgraded_modules"] == []


def test_planned_but_incomplete_integration_is_not_claimable(tmp_path):
    """`integrate` was also granted on an unverified basis: planned != ran.

    Plan-time cannot know whether batch_correction will complete. When it does
    not, the claim must fall back to exploratory rather than resting on the
    module having merely been listed.
    """
    payload = _run_pipeline_with_stub_modules(
        tmp_path,
        obs=_multi_batch_obs(n=40, n_samples=6),
        batch_strategy=BATCH_STRATEGY_INTEGRATE,
        optional_modules=[*DEFAULT_PLAN, "batch_correction"],
        module_statuses=[("batch_correction", "skipped")],
        metadata={"batch_correction_status": "skipped_harmony_unavailable_or_failed_opt_in"},
    )
    risk = payload["metadata"]["batch_risk"]
    assert risk["plan_time_claim"]["clustering_claim_status"] == CLAIM_STATUS_CLAIMABLE
    assert risk["clustering_claim_status"] == CLAIM_STATUS_EXPLORATORY
    assert risk["clustering_claim_reason"] == (
        "exploratory_nonclaimable_integration_planned_but_not_completed"
    )
    assert risk["integration_completed"] is False


def test_stale_clustering_opt_in_does_not_count_as_integration(tmp_path):
    """batch_correction's own `completed_with_stale_clustering_opt_in` state.

    The corrected representation exists but leiden is the PRE-correction one, so
    the clusters are exactly what this claim field must refuse to bless.
    """
    payload = _run_pipeline_with_stub_modules(
        tmp_path,
        obs=_multi_batch_obs(n=40, n_samples=6),
        batch_strategy=BATCH_STRATEGY_INTEGRATE,
        optional_modules=[*DEFAULT_PLAN, "batch_correction"],
        module_statuses=[("batch_correction", "ok")],
        metadata={
            "batch_correction_status": "completed_with_stale_clustering_opt_in"
        },
    )
    risk = payload["metadata"]["batch_risk"]
    assert risk["integration_completed"] is False
    assert risk["clustering_claim_status"] == CLAIM_STATUS_EXPLORATORY


def test_genuinely_completed_integration_is_claimable(tmp_path):
    payload = _run_pipeline_with_stub_modules(
        tmp_path,
        obs=_multi_batch_obs(n=40, n_samples=6),
        batch_strategy=BATCH_STRATEGY_INTEGRATE,
        optional_modules=[*DEFAULT_PLAN, "batch_correction"],
        module_statuses=[("batch_correction", "ok")],
        metadata={"batch_correction_status": "completed"},
    )
    risk = payload["metadata"]["batch_risk"]
    assert risk["integration_completed"] is True
    assert risk["clustering_claim_status"] == CLAIM_STATUS_CLAIMABLE
    assert risk["clustering_claim_reason"] == "multi_batch_integrated"


def test_declared_uncorrected_run_keeps_its_declared_reason(tmp_path):
    payload = _run_pipeline_with_stub_modules(
        tmp_path,
        obs=_multi_batch_obs(n=40, n_samples=6),
        batch_strategy=BATCH_STRATEGY_ACCEPT_UNCORRECTED,
        optional_modules=DEFAULT_PLAN,
    )
    risk = payload["metadata"]["batch_risk"]
    assert risk["clustering_claim_status"] == CLAIM_STATUS_EXPLORATORY
    assert risk["clustering_claim_reason"] == (
        "exploratory_nonclaimable_unintegrated_multi_batch_declared"
    )


def test_reconciliation_downgrade_warns_loudly(capsys):
    plan_time = resolve_batch_risk(
        detection=detect_batch_structure(None),
        declared_strategy=BATCH_STRATEGY_SINGLE_BATCH,
        planned_modules=DEFAULT_PLAN,
    )
    assert plan_time["clustering_claim_status"] == CLAIM_STATUS_CLAIMABLE
    reconciled = reconcile_batch_risk(
        plan_time,
        observation=observe_batch_structure(_multi_batch_obs(n=20, n_samples=4)),
        integration_completed=False,
    )
    warning = reconciliation_warning(reconciled)
    assert warning.startswith("BATCH_RISK_FALSIFIED:")
    assert "single-batch" in warning
    assert "DOWNGRADED" in warning
    assert "artifacts are kept" in warning


def test_reconciliation_downgrade_reaches_stderr_exactly_once(tmp_path, capsys):
    _run_pipeline_with_stub_modules(
        tmp_path,
        obs=_multi_batch_obs(n=40, n_samples=6),
        batch_strategy=BATCH_STRATEGY_SINGLE_BATCH,
        optional_modules=DEFAULT_PLAN,
    )
    captured = capsys.readouterr()
    assert captured.err.count("BATCH_RISK_FALSIFIED:") == 1
    assert "BATCH_RISK_FALSIFIED:" not in captured.out


def test_reconciliation_is_silent_when_the_claim_did_not_worsen():
    plan_time = resolve_batch_risk(
        detection=detect_batch_structure(None),
        declared_strategy=BATCH_STRATEGY_SINGLE_BATCH,
        planned_modules=DEFAULT_PLAN,
    )
    reconciled = reconcile_batch_risk(
        plan_time,
        observation=observe_batch_structure(_single_batch_obs(n=20)),
        integration_completed=False,
    )
    assert reconciliation_warning(reconciled) == ""


def test_unobserved_run_keeps_the_plan_time_claim_and_says_so():
    plan_time = resolve_batch_risk(
        detection=detect_batch_structure(None),
        declared_strategy=BATCH_STRATEGY_AUTO,
        planned_modules=DEFAULT_PLAN,
    )
    reconciled = reconcile_batch_risk(
        plan_time, observation=observe_batch_structure(None), integration_completed=False
    )
    assert reconciled["observation_status"] == "unobserved"
    assert reconciled["claim_basis"] == "plan_time_detection"
    assert reconciled["clustering_claim_status"] == CLAIM_STATUS_UNDETERMINED


def test_reconciliation_never_upgrades_a_downgraded_claim_on_declaration_alone():
    """No declaration can make observed, uncorrected multi-batch claimable."""
    for strategy in (
        BATCH_STRATEGY_AUTO,
        BATCH_STRATEGY_ACCEPT_UNCORRECTED,
        BATCH_STRATEGY_INTEGRATE,
    ):
        plan_time = resolve_batch_risk(
            detection=detect_batch_structure(None),
            declared_strategy=strategy,
            planned_modules=(
                [*DEFAULT_PLAN, "batch_correction"]
                if strategy == BATCH_STRATEGY_INTEGRATE
                else DEFAULT_PLAN
            ),
        )
        reconciled = reconcile_batch_risk(
            plan_time,
            observation=observe_batch_structure(_multi_batch_obs(n=20, n_samples=4)),
            integration_completed=False,
        )
        assert reconciled["clustering_claim_status"] == CLAIM_STATUS_EXPLORATORY, (
            f"strategy {strategy!r} produced a claimable run on observed, "
            "uncorrected multi-batch data"
        )


def test_observation_and_clustering_diagnostic_use_one_reading(tmp_path):
    """The two records cannot disagree: they share observe_batch_structure."""
    from workflow.modular.modules.clustering import ClusteringModule

    obs = _multi_batch_obs(n=40, n_samples=6)
    adata = ad.AnnData(X=np.zeros((40, 3), dtype="float32"), obs=obs)
    ctx = type(
        "_Ctx", (), {"cfg": type("_Cfg", (), {"optional_modules": []})(), "metadata": {}}
    )()
    ClusteringModule._record_batch_confounding_risk(ClusteringModule(), adata, ctx)
    diagnostic = ctx.metadata["batch_confounding_risk"]["candidate_batch_columns"]
    assert diagnostic == observe_batch_structure(obs)[
        "observed_candidate_batch_columns"
    ]


def test_governed_manifest_carries_the_reconciled_claim_not_the_plan_time_one(
    tmp_path, monkeypatch
):
    """The field a downstream consumer reads must be the reconciled one.

    The CLI used to copy its own plan-time envelope into manifest.json, so even
    a correct producer-side downgrade would not have reached the governed
    manifest.
    """
    import os

    from workflow.modular import cli as cli_mod
    from workflow.modular import manifest_writer

    run_id = "2026-08-03T0000Z-bbbbbbb"
    project_root = tmp_path / "project"
    sample_root = tmp_path / "sample"
    sample_root.mkdir()
    monkeypatch.setattr(
        manifest_writer,
        "factory_git_state",
        lambda root: {"sha": "bbbbbbb", "dirty": False, "diff_sha256": ""},
    )

    def producer_that_reconciles(cfg, ledger=None):
        # Plan-time granted claimable on the unverifiable declaration.
        assert cfg.batch_risk["clustering_claim_status"] == CLAIM_STATUS_CLAIMABLE
        reconciled = reconcile_batch_risk(
            cfg.batch_risk,
            observation=observe_batch_structure(_multi_batch_obs(n=20, n_samples=4)),
            integration_completed=False,
        )
        producer_dir = Path(cfg.output_dir) / "p_20260803_000000"
        producer_dir.mkdir(parents=True, exist_ok=True)
        producer = producer_dir / "run_manifest.json"
        producer.write_text(
            json.dumps({
                "modules_run": ["clustering"],
                "overall_status": "complete",
                "metadata": {"batch_risk": reconciled},
            }),
            encoding="utf-8",
        )
        return producer

    monkeypatch.setattr(cli_mod, "run_pipeline", producer_that_reconciles)
    monkeypatch.setattr(
        "sys.argv",
        [
            "cli", "--project", "p", "--sample-root", str(sample_root),
            "--project-root", str(project_root), "--run-id", run_id,
            "--allow-dirty", "--batch-strategy", "single-batch",
        ],
    )
    monkeypatch.delenv("SC_REQUIRE_PROJECT_ROOT", raising=False)

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        cli_mod.main()
    finally:
        os.chdir(cwd)

    governed = json.loads(
        (project_root / "runs" / run_id / "manifest.json").read_text(encoding="utf-8")
    )
    risk = governed["batch_risk"]
    assert risk["clustering_claim_status"] == CLAIM_STATUS_EXPLORATORY, (
        "the governed manifest carried the plan-time claim instead of the "
        "producer's reconciled one"
    )
    assert risk["clustering_claim_reason"] == (
        "exploratory_nonclaimable_single_batch_declaration_falsified"
    )
    assert risk["claim_basis"] == "runtime_observation"


def test_run_pipeline_records_batch_risk_for_programmatic_callers(tmp_path):
    """A caller that never touches the CLI still gets the accounting."""
    from workflow.modular import pipeline as pipeline_mod
    from workflow.modular.config import CellRangerConfig, PipelineConfig

    sample_root = tmp_path / "cohort"
    sample_root.mkdir()
    _write_h5ad(sample_root / "prepared_input.h5ad", _multi_batch_obs())
    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(
            sample_root=sample_root, outs_dir=sample_root / "outs"
        ),
        optional_modules=list(DEFAULT_OPTIONAL_MODULES),
    )
    captured: dict = {}

    def fake_run_sequential(modules, registry, ctx, mandatory):
        captured["batch_risk"] = ctx.metadata.get("batch_risk")

    monkey_target = pipeline_mod._run_sequential
    pipeline_mod._run_sequential = fake_run_sequential
    try:
        pipeline_mod.run_pipeline(cfg)
    finally:
        pipeline_mod._run_sequential = monkey_target

    assert captured["batch_risk"] is not None
    assert captured["batch_risk"]["clustering_claim_status"] == (
        CLAIM_STATUS_EXPLORATORY
    )
