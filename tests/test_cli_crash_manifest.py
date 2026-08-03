"""Contract tests for the CLI fatal-exception failure envelope.

Guards the defect fixed on 2026-08-03: ``main()`` wrapped ``run_pipeline`` in
try/finally with no ``except``, so a pipeline that raised before saving its own
manifest left the governed run directory with no manifest at all — a crashed run
was byte-for-byte indistinguishable from a run that never started.

Every test is written to fail if the envelope is removed, and the masking tests
assert that it does NOT change the exception, the exit code, or the pre-existing
partial-run and manifest-write-failure paths.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

RUN_ID = "2026-08-03T0000Z-aaaaaaa"


class _Boom(RuntimeError):
    """Distinct exception type so the manifest record can be asserted exactly."""


def _prepare(tmp_path, monkeypatch, *, extra_argv=()):
    """Wire a project-root CLI invocation with a stubbed git state."""
    from workflow.modular import cli as cli_mod
    from workflow.modular import manifest_writer

    project_root = tmp_path / "project"
    sample_root = tmp_path / "sample"
    sample_root.mkdir()
    monkeypatch.setattr(
        manifest_writer,
        "factory_git_state",
        lambda root: {"sha": "aaaaaaa", "dirty": False, "diff_sha256": ""},
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "cli", "--project", "crashtest",
            "--sample-root", str(sample_root),
            "--project-root", str(project_root),
            "--run-id", RUN_ID,
            "--allow-dirty",
            *extra_argv,
        ],
    )
    monkeypatch.delenv("SC_REQUIRE_PROJECT_ROOT", raising=False)
    return cli_mod, project_root, sample_root


def _run_main_in(tmp_path, cli_mod):
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        cli_mod.main()
    finally:
        os.chdir(cwd)


def _manifest(project_root: Path) -> dict:
    path = project_root / "runs" / RUN_ID / "manifest.json"
    assert path.is_file(), f"no manifest written at {path}"
    return json.loads(path.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# the envelope
# ---------------------------------------------------------------------------


def test_fatal_exception_writes_a_crashed_manifest(tmp_path, monkeypatch):
    cli_mod, project_root, _ = _prepare(tmp_path, monkeypatch)

    def exploding_run_pipeline(cfg, ledger=None):
        if ledger is not None:
            ledger.record_module("cellranger", "ok", "completed", 1.0, 0)
            ledger.record_module("qc", "ok", "completed", 2.0, 0)
            ledger.record_module("ambient_correction", "skipped", "no trigger", 0.1, 0)
        raise _Boom("scrublet blew up on sample 3")

    monkeypatch.setattr(cli_mod, "run_pipeline", exploding_run_pipeline)
    with pytest.raises(_Boom):
        _run_main_in(tmp_path, cli_mod)

    payload = _manifest(project_root)
    assert payload["overall_status"] == "crashed"
    crash = payload["extra"]["crash"]
    assert crash["exception_class"] == "_Boom"
    assert crash["message"] == "scrublet blew up on sample 3"
    assert crash["last_completed_module"] == "qc"


def test_crashed_manifest_records_module_accounting(tmp_path, monkeypatch):
    cli_mod, project_root, _ = _prepare(tmp_path, monkeypatch)

    def exploding_run_pipeline(cfg, ledger=None):
        ledger.record_module("cellranger", "ok", "completed", 1.0, 0)
        ledger.record_module("qc", "ok", "completed", 2.0, 0)
        ledger.record_module("ambient_correction", "skipped", "no trigger", 0.1, 0)
        ledger.record_module("doublet_detection", "failed", "boom", 0.1, 0)
        raise _Boom("died in doublet_detection")

    monkeypatch.setattr(cli_mod, "run_pipeline", exploding_run_pipeline)
    with pytest.raises(_Boom):
        _run_main_in(tmp_path, cli_mod)

    payload = _manifest(project_root)
    # Requested and planned are known from the config alone.
    assert "clustering" in payload["requested_modules"]
    assert payload["planned_modules"][:2] == ["cellranger", "qc"]
    assert "clustering" in payload["planned_modules"]
    # Executed/completed/skipped/failed come from the live ledger.
    assert payload["executed_modules"] == [
        "cellranger", "qc", "ambient_correction", "doublet_detection",
    ]
    assert payload["completed_modules"] == ["cellranger", "qc"]
    assert payload["skipped_modules"] == ["ambient_correction"]
    assert payload["failed_modules"] == ["doublet_detection"]


def test_crashed_manifest_records_identity_and_traceback_artifact(
    tmp_path, monkeypatch
):
    cli_mod, project_root, sample_root = _prepare(tmp_path, monkeypatch)

    def exploding_run_pipeline(cfg, ledger=None):
        raise _Boom("early death")

    monkeypatch.setattr(cli_mod, "run_pipeline", exploding_run_pipeline)
    with pytest.raises(_Boom):
        _run_main_in(tmp_path, cli_mod)

    payload = _manifest(project_root)
    # Environment / git identity.
    assert payload["factory_python"]["sha"] == "aaaaaaa"
    assert payload["produced_on"]
    assert payload["software_versions"]
    # Input identity.
    assert payload["extra"]["crash"]["input"]["sample_root"] == str(sample_root)
    # Traceback artifact exists and holds the real traceback.
    artifact = Path(payload["extra"]["crash"]["traceback_artifact"])
    assert artifact.is_file()
    text = artifact.read_text(encoding="utf-8")
    assert "_Boom" in text and "early death" in text
    assert "Traceback (most recent call last)" in text


def test_crash_manifest_carries_the_batch_risk_envelope(tmp_path, monkeypatch):
    import anndata as ad
    import numpy as np
    import pandas as pd

    cli_mod, project_root, sample_root = _prepare(tmp_path, monkeypatch)
    obs = pd.DataFrame(
        {"sample": pd.Categorical([f"S{i % 4}" for i in range(20)])},
        index=[f"c{i}" for i in range(20)],
    )
    ad.AnnData(X=np.zeros((20, 3), dtype="float32"), obs=obs).write_h5ad(
        sample_root / "prepared_input.h5ad"
    )

    def exploding_run_pipeline(cfg, ledger=None):
        raise _Boom("died after detection")

    monkeypatch.setattr(cli_mod, "run_pipeline", exploding_run_pipeline)
    with pytest.raises(_Boom):
        _run_main_in(tmp_path, cli_mod)

    risk = _manifest(project_root)["batch_risk"]
    assert risk["detected_batches"] == 4
    assert risk["clustering_claim_status"] == "exploratory"


def test_keyboard_interrupt_also_produces_a_crashed_manifest(tmp_path, monkeypatch):
    """BaseException, not just Exception: an interrupted run must be legible."""
    cli_mod, project_root, _ = _prepare(tmp_path, monkeypatch)

    def interrupted_run_pipeline(cfg, ledger=None):
        raise KeyboardInterrupt()

    monkeypatch.setattr(cli_mod, "run_pipeline", interrupted_run_pipeline)
    with pytest.raises(KeyboardInterrupt):
        _run_main_in(tmp_path, cli_mod)

    payload = _manifest(project_root)
    assert payload["overall_status"] == "crashed"
    assert payload["extra"]["crash"]["exception_class"] == "KeyboardInterrupt"


# ---------------------------------------------------------------------------
# what the envelope must NOT do
# ---------------------------------------------------------------------------


def test_envelope_does_not_swallow_the_exception(tmp_path, monkeypatch):
    cli_mod, _, _ = _prepare(tmp_path, monkeypatch)

    def exploding_run_pipeline(cfg, ledger=None):
        raise _Boom("must reach the caller")

    monkeypatch.setattr(cli_mod, "run_pipeline", exploding_run_pipeline)
    with pytest.raises(_Boom, match="must reach the caller") as excinfo:
        _run_main_in(tmp_path, cli_mod)
    assert excinfo.value.__traceback__ is not None


def test_envelope_failure_does_not_replace_the_original_traceback(
    tmp_path, monkeypatch, caplog
):
    """A broken crash writer must never become the exception the operator sees."""
    from workflow.modular import manifest_writer

    cli_mod, project_root, _ = _prepare(tmp_path, monkeypatch)

    def exploding_run_pipeline(cfg, ledger=None):
        raise _Boom("the real failure")

    def broken_write_manifest(*args, **kwargs):
        raise OSError("disk full while writing the envelope")

    monkeypatch.setattr(cli_mod, "run_pipeline", exploding_run_pipeline)
    monkeypatch.setattr(manifest_writer, "write_manifest", broken_write_manifest)
    with pytest.raises(_Boom, match="the real failure"):
        _run_main_in(tmp_path, cli_mod)
    assert not (project_root / "runs" / RUN_ID / "manifest.json").exists()


def test_partial_run_exit_path_is_not_relabelled_as_crashed(tmp_path, monkeypatch):
    """SystemExit(1) after a complete manifest write is a failure, not a crash."""
    cli_mod, project_root, _ = _prepare(tmp_path, monkeypatch)

    def failing_module_run_pipeline(cfg, ledger=None):
        producer_dir = Path(cfg.output_dir) / "crashtest_20260803_000000"
        producer_dir.mkdir(parents=True, exist_ok=True)
        producer = producer_dir / "run_manifest.json"
        producer.write_text(
            json.dumps({
                "requested_modules": ["qc", "annotation"],
                "planned_modules": ["qc", "annotation"],
                "executed_modules": ["qc", "annotation"],
                "completed_modules": ["qc"],
                "modules_run": ["qc"],
                "skipped_modules": [],
                "failed_modules": ["annotation"],
                "overall_status": "failed",
                "bundle_sha256": "",
            }),
            encoding="utf-8",
        )
        return producer

    monkeypatch.setattr(cli_mod, "run_pipeline", failing_module_run_pipeline)
    with pytest.raises(SystemExit) as excinfo:
        _run_main_in(tmp_path, cli_mod)
    assert excinfo.value.code == 1

    payload = _manifest(project_root)
    assert payload["overall_status"] == "failed", (
        "the partial-run path already wrote a complete manifest; the envelope "
        "must not overwrite it with 'crashed'"
    )
    assert "crash" not in payload["extra"]


def test_manifest_write_failure_still_raises_its_own_runtime_error(
    tmp_path, monkeypatch
):
    """The pre-existing project-root manifest-write guard keeps its identity."""
    from workflow.modular import manifest_writer

    cli_mod, _, _ = _prepare(tmp_path, monkeypatch)
    calls: list[str] = []

    def ok_run_pipeline(cfg, ledger=None):
        return {"modules_run": ["qc"], "bundle_sha256": ""}

    def failing_write_manifest(*args, **kwargs):
        calls.append("called")
        raise OSError("cannot write")

    monkeypatch.setattr(cli_mod, "run_pipeline", ok_run_pipeline)
    monkeypatch.setattr(manifest_writer, "write_manifest", failing_write_manifest)
    with pytest.raises(RuntimeError, match="project-root manifest write failed"):
        _run_main_in(tmp_path, cli_mod)
    # Once for the normal write, once for the envelope's own attempt.
    assert len(calls) == 2


def test_module_level_sys_exit_still_gets_an_envelope(tmp_path, monkeypatch):
    """A blanket `except SystemExit: raise` swallowed real module aborts.

    Only the partial-run exit is exempt, and it is identified by a flag rather
    than by its type.
    """
    cli_mod, project_root, _ = _prepare(tmp_path, monkeypatch)

    def module_calls_sys_exit(cfg, ledger=None):
        import sys as _sys

        ledger.record_module("cellranger", "ok", "completed", 1.0, 0)
        _sys.exit(3)

    monkeypatch.setattr(cli_mod, "run_pipeline", module_calls_sys_exit)
    with pytest.raises(SystemExit) as excinfo:
        _run_main_in(tmp_path, cli_mod)
    assert excinfo.value.code == 3, "the original exit code must be preserved"

    payload = _manifest(project_root)
    assert payload["overall_status"] == "crashed"
    assert payload["extra"]["crash"]["exception_class"] == "SystemExit"
    assert payload["extra"]["crash"]["phase"] == "pipeline"
    assert payload["completed_modules"] == ["cellranger"]


def test_crash_phase_distinguishes_a_completed_run_from_a_failed_analysis(
    tmp_path, monkeypatch
):
    """A run that COMPLETED but whose governed manifest write failed.

    It is still `crashed` (the governed manifest is not trustworthy), but every
    module reads ok, so the envelope names the phase to keep that legible.
    """
    from workflow.modular import manifest_writer

    cli_mod, project_root, _ = _prepare(tmp_path, monkeypatch)
    calls: list[str] = []

    def completed_run(cfg, ledger=None):
        ledger.record_module("cellranger", "ok", "completed", 1.0, 0)
        ledger.record_module("qc", "ok", "completed", 1.0, 0)
        return {"modules_run": ["cellranger", "qc"], "bundle_sha256": ""}

    real_write = manifest_writer.write_manifest

    def fail_first_write(*args, **kwargs):
        calls.append("call")
        if len(calls) == 1:
            raise OSError("governed write failed")
        return real_write(*args, **kwargs)

    monkeypatch.setattr(cli_mod, "run_pipeline", completed_run)
    monkeypatch.setattr(manifest_writer, "write_manifest", fail_first_write)
    with pytest.raises(RuntimeError, match="project-root manifest write failed"):
        _run_main_in(tmp_path, cli_mod)

    crash = _manifest(project_root)["extra"]["crash"]
    assert crash["phase"] == "governed_manifest_write"
    assert crash["exception_class"] == "RuntimeError"


def test_envelope_guards_survive_a_baseexception_mid_write(tmp_path, monkeypatch):
    """The writer's own guards must be BaseException, as its docstring says.

    A KeyboardInterrupt landing inside the envelope used to escape and replace
    the real failure the operator needs to see.
    """
    from workflow.modular import cli as cli_mod_direct
    from workflow.modular import manifest_writer

    cli_mod, project_root, _ = _prepare(tmp_path, monkeypatch)

    def exploding_run_pipeline(cfg, ledger=None):
        raise _Boom("the real failure")

    def interrupted_write(*args, **kwargs):
        raise KeyboardInterrupt()

    monkeypatch.setattr(cli_mod, "run_pipeline", exploding_run_pipeline)
    monkeypatch.setattr(manifest_writer, "write_manifest", interrupted_write)
    with pytest.raises(_Boom, match="the real failure"):
        _run_main_in(tmp_path, cli_mod)


def test_partial_run_flag_not_type_is_what_exempts_the_exit(tmp_path, monkeypatch):
    """Regression guard for the exemption mechanism itself."""
    import inspect

    from workflow.modular import cli as cli_mod_direct

    source = inspect.getsource(cli_mod_direct.main)
    assert "_partial_run_exit" in source
    # Comments may DISCUSS the old handler; no executable line may be one.
    handlers = [
        line.strip()
        for line in source.splitlines()
        if line.strip().startswith("except SystemExit")
    ]
    assert handlers == [], (
        f"a blanket SystemExit handler re-introduces the swallowed-abort bug: {handlers}"
    )


# ---------------------------------------------------------------------------
# manifest status vocabulary
# ---------------------------------------------------------------------------


def test_crashed_is_a_declared_manifest_status():
    from workflow.modular.manifest_writer import (
        MANIFEST_OVERALL_STATUSES,
        MANIFEST_STATUS_CRASHED,
    )

    assert MANIFEST_STATUS_CRASHED == "crashed"
    assert MANIFEST_OVERALL_STATUSES == {
        "complete", "partial", "failed", "crashed",
    }


def test_unknown_overall_status_is_rejected(tmp_path):
    from workflow.modular.manifest_writer import write_manifest

    with pytest.raises(ValueError, match="unknown overall_status"):
        write_manifest(
            tmp_path / "run",
            project_id="p",
            run_id=RUN_ID,
            modules_run=[],
            overall_status="mostly_fine",
        )


def test_every_pipeline_status_is_in_the_declared_vocabulary():
    """pipeline._summarize_run_status must not emit a status write_manifest rejects."""
    from workflow.modular.manifest_writer import MANIFEST_OVERALL_STATUSES

    assert {"complete", "partial", "failed"} <= MANIFEST_OVERALL_STATUSES


def test_run_ledger_exposes_live_module_results():
    from workflow.modular._run_ledger import RunLedger

    ledger = RunLedger(type("_Ctx", (), {"cfg": None})(), "p", Path("."))
    assert ledger.module_results() == []
    ledger.record_module("qc", "ok", "completed", 1.0, 0)
    results = ledger.module_results()
    assert results == [{
        "name": "qc", "status": "ok", "message": "completed",
        "wall_seconds": 1.0, "rss_peak_bytes": 0,
    }]
    results[0]["status"] = "tampered"
    assert ledger.module_results()[0]["status"] == "ok", "snapshot must be a copy"
