from __future__ import annotations

import copy
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path


FACTORY_ROOT = Path(__file__).resolve().parents[1]


def test_legacy_output_contract_blocks_early_retirement() -> None:
    from workflow.modular.legacy_output import load_legacy_output_contract

    contract = load_legacy_output_contract()

    assert contract["schema_version"] == "legacy_factory_output_deprecation/1.0"
    assert contract["status"] == "deprecated_compatibility"
    assert contract["replacement"] == "--project-root"
    assert contract["retirement"]["access_window_days"] == 14
    assert contract["retirement"]["monitoring_started_at"] == "2026-08-02T00:00:00Z"
    assert contract["retirement"]["project_root_migration_complete"] is False
    assert contract["retirement"]["owner_approval_recorded"] is False
    assert contract["protection"]["remove_before_eligible"] == "forbidden"
    assert set(contract["legacy_surface"]["entrypoints"]) == {
        "python -m workflow.modular.cli",
        "scripts/scfactory.py run",
        "scripts/export_singlecell_r_bundle.py",
        "scripts/pack_run_for_mac.sh",
    }


def test_legacy_access_is_persisted_as_jsonl_outside_factory(tmp_path: Path) -> None:
    from workflow.modular.legacy_output import record_legacy_output_access

    log_path = tmp_path / "state" / "legacy-output.jsonl"
    output_dir = FACTORY_ROOT / "results"
    recorded = record_legacy_output_access(
        output_dir=output_dir,
        project="legacy-smoke",
        environ={"SC_LEGACY_OUTPUT_ACCESS_LOG": str(log_path)},
        now=datetime(2026, 8, 2, 4, 5, tzinfo=timezone.utc),
    )

    assert recorded == log_path.resolve()
    payload = json.loads(log_path.read_text(encoding="utf-8"))
    assert payload == {
        "contract": "legacy_factory_output_deprecation/1.0",
        "event": "legacy_factory_output_access",
        "factory_local": True,
        "output_dir": str(output_dir.resolve()),
        "project": "legacy-smoke",
        "recorded_at": "2026-08-02T04:05:00Z",
        "replacement": "--project-root",
    }
    assert FACTORY_ROOT not in recorded.parents


def test_recent_access_keeps_retirement_ineligible(tmp_path: Path) -> None:
    from workflow.modular.legacy_output import evaluate_legacy_output_retirement

    contract = {
        "retirement": {
            "access_window_days": 14,
            "project_root_migration_complete": True,
            "owner_approval_recorded": True,
        }
    }
    telemetry = tmp_path / "legacy-output.jsonl"
    telemetry.write_text(
        json.dumps({"recorded_at": "2026-07-25T00:00:00Z"}) + "\n",
        encoding="utf-8",
    )

    status = evaluate_legacy_output_retirement(
        contract,
        telemetry,
        now=datetime(2026, 8, 2, tzinfo=timezone.utc),
    )

    assert status["eligible"] is False
    assert status["recent_access_count"] == 1
    assert "14-day access window is not silent" in status["blockers"]


def test_silent_window_alone_cannot_override_contract_flags(tmp_path: Path) -> None:
    from workflow.modular.legacy_output import (
        evaluate_legacy_output_retirement,
        load_legacy_output_contract,
    )

    contract = copy.deepcopy(load_legacy_output_contract())
    telemetry = tmp_path / "missing-is-zero-events.jsonl"
    status = evaluate_legacy_output_retirement(
        contract,
        telemetry,
        now=datetime(2026, 8, 2, tzinfo=timezone.utc),
    )

    assert status["recent_access_count"] == 0
    assert status["eligible"] is False
    assert "project-root migration is not complete" in status["blockers"]
    assert "owner approval is not recorded" in status["blockers"]


def test_monitoring_must_cover_the_whole_window(tmp_path: Path) -> None:
    from workflow.modular.legacy_output import evaluate_legacy_output_retirement

    contract = {
        "retirement": {
            "access_window_days": 14,
            "monitoring_started_at": "2026-07-25T00:00:00Z",
            "project_root_migration_complete": True,
            "owner_approval_recorded": True,
        }
    }

    status = evaluate_legacy_output_retirement(
        contract,
        tmp_path / "empty.jsonl",
        now=datetime(2026, 8, 2, tzinfo=timezone.utc),
    )

    assert status["eligible"] is False
    assert "14-day monitoring window has not elapsed" in status["blockers"]


def test_retirement_becomes_eligible_only_after_every_condition(tmp_path: Path) -> None:
    from workflow.modular.legacy_output import evaluate_legacy_output_retirement

    contract = {
        "retirement": {
            "access_window_days": 14,
            "monitoring_started_at": "2026-07-01T00:00:00Z",
            "project_root_migration_complete": True,
            "owner_approval_recorded": True,
        }
    }

    status = evaluate_legacy_output_retirement(
        contract,
        tmp_path / "empty.jsonl",
        now=datetime(2026, 8, 2, tzinfo=timezone.utc),
    )

    assert status["eligible"] is True
    assert status["blockers"] == []


def test_cli_and_scfactory_share_one_legacy_default() -> None:
    from unittest import mock

    from workflow.modular.cli import parse_args
    from workflow.modular.legacy_output import LEGACY_DEFAULT_OUTPUT_DIR

    with mock.patch(
        "sys.argv",
        ["cli", "--project", "test", "--sample-root", "/tmp/input"],
    ):
        args = parse_args()

    assert Path(args.output_dir) == LEGACY_DEFAULT_OUTPUT_DIR
    source = (FACTORY_ROOT / "scripts" / "scfactory.py").read_text(encoding="utf-8")
    assert "LEGACY_DEFAULT_OUTPUT_DIR" in source
    assert 'REPO_ROOT / "results"' not in source


def test_warning_names_contract_window_and_telemetry(tmp_path: Path) -> None:
    from workflow.modular.legacy_output import legacy_output_warning_message

    message = legacy_output_warning_message(tmp_path / "legacy.jsonl")

    assert "legacy_factory_output_deprecation/1.0" in message
    assert "14-day silent access window" in message
    assert str((tmp_path / "legacy.jsonl").resolve()) in message


def test_legacy_bundle_export_uses_the_same_contract(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts import export_singlecell_r_bundle as exporter

    telemetry = tmp_path / "bundle-export-access.jsonl"
    monkeypatch.setenv("SC_LEGACY_OUTPUT_ACCESS_LOG", str(telemetry))
    monkeypatch.delenv("SC_REQUIRE_PROJECT_ROOT", raising=False)
    monkeypatch.setattr(
        exporter,
        "export_bundle",
        lambda config: {"schema_version": "singlecell_r_bundle_v2.1"},
    )

    rc = exporter.main(
        [
            "--input",
            str(tmp_path / "source.h5ad"),
            "--output",
            str(tmp_path / "legacy-bundle"),
        ]
    )

    assert rc == 0
    event = json.loads(telemetry.read_text(encoding="utf-8"))
    assert event["event"] == "legacy_factory_output_access"
    assert event["project"] == "bundle_export"
    assert event["output_dir"] == str((tmp_path / "legacy-bundle").resolve())


def test_legacy_pack_wrapper_records_access_before_writing(tmp_path: Path) -> None:
    run_dir = tmp_path / "legacy-run"
    run_dir.mkdir()
    telemetry = tmp_path / "pack-access.jsonl"
    env = dict(os.environ)
    env.pop("SC_REQUIRE_PROJECT_ROOT", None)
    env["SC_LEGACY_OUTPUT_ACCESS_LOG"] = str(telemetry)

    completed = subprocess.run(
        ["bash", str(FACTORY_ROOT / "scripts" / "pack_run_for_mac.sh"), str(run_dir)],
        cwd=FACTORY_ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr or completed.stdout
    assert "legacy_factory_output_deprecation/1.0" in completed.stderr
    event = json.loads(telemetry.read_text(encoding="utf-8"))
    assert event["project"] == "pack_run_for_mac"
    assert event["output_dir"] == str(run_dir.resolve())
