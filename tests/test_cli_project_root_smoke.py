"""Smoke tests for --project-root / --run-id / --allow-dirty CLI flags."""
import json
import subprocess
import sys
from pathlib import Path

import pytest


def _get_parser():
    # Import parse_args directly so we can inspect the argparse namespace.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from workflow.modular.cli import parse_args
    return parse_args


def test_allow_dirty_flag_exists():
    from workflow.modular.cli import parse_args
    import unittest.mock as mock
    with mock.patch(
        "sys.argv",
        [
            "cli",
            "--project", "test",
            "--sample-root", "/tmp/x",
            "--project-root", "/tmp/proj",
            "--allow-dirty",
        ],
    ):
        args = parse_args()
    assert args.allow_dirty is True


def test_allow_partial_run_flag_exists():
    from workflow.modular.cli import parse_args
    import unittest.mock as mock
    with mock.patch(
        "sys.argv",
        [
            "cli", "--project", "test", "--sample-root", "/tmp/x",
            "--allow-partial-run",
        ],
    ):
        args = parse_args()
    assert args.allow_partial_run is True


def test_project_root_flag_exists():
    from workflow.modular.cli import parse_args
    import unittest.mock as mock
    with mock.patch(
        "sys.argv",
        [
            "cli",
            "--project", "test",
            "--sample-root", "/tmp/x",
            "--project-root", "/tmp/myproject",
        ],
    ):
        args = parse_args()
    assert args.project_root == "/tmp/myproject"
    assert args.allow_dirty is False


def test_run_id_flag_exists():
    from workflow.modular.cli import parse_args
    import unittest.mock as mock
    with mock.patch(
        "sys.argv",
        [
            "cli",
            "--project", "test",
            "--sample-root", "/tmp/x",
            "--project-root", "/tmp/myproject",
            "--run-id", "2026-05-14T0001Z-aaaaaaa",
        ],
    ):
        args = parse_args()
    assert args.run_id == "2026-05-14T0001Z-aaaaaaa"


def test_legacy_invocation_no_project_root():
    """Legacy invocations without --project-root must parse successfully."""
    from workflow.modular.cli import parse_args
    import unittest.mock as mock
    with mock.patch(
        "sys.argv",
        [
            "cli",
            "--project", "test",
            "--sample-root", "/tmp/x",
        ],
    ):
        args = parse_args()
    assert args.project_root is None
    assert args.run_id is None
    assert args.allow_dirty is False


def test_sc_require_project_root_unset_legacy_warns(tmp_path):
    """When SC_REQUIRE_PROJECT_ROOT is unset, missing --project-root emits DeprecationWarning."""
    import unittest.mock as mock
    import warnings

    # Ensure env var is not set
    with mock.patch.dict(
        "os.environ",
        {"SC_LEGACY_OUTPUT_ACCESS_LOG": str(tmp_path / "legacy-access.jsonl")},
        clear=False,
    ):
        # Remove the key if present
        import os
        os.environ.pop("SC_REQUIRE_PROJECT_ROOT", None)

        from workflow.modular import cli as cli_mod
        import importlib
        importlib.reload(cli_mod)

        with mock.patch(
            "sys.argv",
            [
                "cli",
                "--project", "test",
                "--sample-root", "/tmp/x",
                "--project-root", "/tmp/proj",
            ],
        ):
            # parse_args with --project-root set: no warning expected
            args = cli_mod.parse_args()
        assert args.project_root == "/tmp/proj"

    # Now test without --project-root: the main() function emits the warning
    # We test by calling main() and checking it issues DeprecationWarning before
    # hitting pipeline work (which will fail, so we catch SystemExit too)
    with mock.patch.dict(
        "os.environ",
        {"SC_LEGACY_OUTPUT_ACCESS_LOG": str(tmp_path / "legacy-access.jsonl")},
        clear=False,
    ):
        import os
        os.environ.pop("SC_REQUIRE_PROJECT_ROOT", None)

        with mock.patch(
            "sys.argv",
            [
                "cli",
                "--project", "test",
                "--sample-root", "/tmp/x",
                "--output-dir", str(tmp_path / "legacy-output"),
            ],
        ):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                try:
                    cli_mod.main()
                except (SystemExit, Exception):
                    pass
            deprecation_warnings = [
                w for w in caught
                if issubclass(w.category, DeprecationWarning)
                and "project-root" in str(w.message).lower()
            ]
            assert len(deprecation_warnings) >= 1, (
                f"Expected DeprecationWarning about --project-root, got: {[str(w.message) for w in caught]}"
            )
    event = json.loads(
        (tmp_path / "legacy-access.jsonl").read_text(encoding="utf-8")
    )
    assert event["event"] == "legacy_factory_output_access"
    assert event["project"] == "test"
    assert event["output_dir"] == str((tmp_path / "legacy-output").resolve())


def test_sc_require_project_root_set_hard_errors():
    """When SC_REQUIRE_PROJECT_ROOT=1, missing --project-root must exit with code 2."""
    import unittest.mock as mock
    with mock.patch.dict("os.environ", {"SC_REQUIRE_PROJECT_ROOT": "1"}):
        from workflow.modular import cli as cli_mod
        import importlib
        importlib.reload(cli_mod)

        with mock.patch(
            "sys.argv",
            [
                "cli",
                "--project", "test",
                "--sample-root", "/tmp/x",
            ],
        ):
            with pytest.raises(SystemExit) as exc_info:
                cli_mod.main()
            assert exc_info.value.code == 2, (
                f"Expected exit code 2, got {exc_info.value.code}"
            )


def test_project_root_run_ledger_writes_under_run_dir(tmp_path, monkeypatch):
    """With --project-root, RunLedger must live under <run>/ops, not output/ops."""
    import importlib
    import os
    from workflow.modular import cli as cli_mod
    from workflow.modular import manifest_writer

    importlib.reload(cli_mod)
    run_id = "2026-05-14T0001Z-aaaaaaa"
    project_root = tmp_path / "project"
    sample_root = tmp_path / "sample"
    sample_root.mkdir()

    monkeypatch.setattr(
        manifest_writer,
        "factory_git_state",
        lambda root: {"sha": "aaaaaaa", "dirty": False, "branch": "test"},
    )

    def fake_run_pipeline(cfg, ledger=None):
        if ledger is not None:
            ledger.record_module("smoke", "ok", "stub pipeline", 0.0, 0)
            ledger.record_end(None)
            ledger.write()
        return {"modules_run": ["smoke"], "bundle_sha256": "abc123"}

    monkeypatch.setattr(cli_mod, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        "sys.argv",
        [
            "cli",
            "--project", "test",
            "--sample-root", str(sample_root),
            "--project-root", str(project_root),
            "--run-id", run_id,
            "--allow-dirty",
        ],
    )
    monkeypatch.delenv("SC_REQUIRE_PROJECT_ROOT", raising=False)

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        cli_mod.main()
    finally:
        os.chdir(cwd)

    run_dir = project_root / "runs" / run_id
    ledgers = list((run_dir / "ops" / "run_ledger").glob("*.json"))
    assert ledgers, "expected RunLedger JSON under project-root run dir"
    assert not (tmp_path / "output" / "ops" / "run_ledger").exists()


@pytest.mark.parametrize("allow_partial, expected_exit", [(False, 1), (True, 0)])
def test_real_path_manifest_propagates_failure_and_controls_exit(
    tmp_path, monkeypatch, allow_partial, expected_exit
):
    import os
    from workflow.modular import cli as cli_mod
    from workflow.modular import manifest_writer

    project_root = tmp_path / "project"
    sample_root = tmp_path / "sample"
    sample_root.mkdir()
    run_id = "2026-05-14T0001Z-aaaaaaa"

    monkeypatch.setattr(
        manifest_writer,
        "factory_git_state",
        lambda root: {"sha": "aaaaaaa", "dirty": False, "diff_sha256": ""},
    )

    def fake_run_pipeline(cfg, ledger=None):
        producer_dir = Path(cfg.output_dir) / "test_20260514_000100"
        producer_dir.mkdir(parents=True, exist_ok=True)
        producer = producer_dir / "run_manifest.json"
        producer.write_text(
            json.dumps(
                {
                    "requested_modules": ["qc", "annotation"],
                    "planned_modules": ["qc", "annotation"],
                    "executed_modules": ["qc", "annotation"],
                    "completed_modules": ["qc"],
                    "modules_run": ["qc"],
                    "skipped_modules": [],
                    "failed_modules": ["annotation"],
                    "overall_status": "failed",
                    "bundle_sha256": "abc123",
                    "factory_python": {
                        "sha": "aaaaaaa", "dirty": False, "diff_sha256": ""
                    },
                }
            ),
            encoding="utf-8",
        )
        return producer

    monkeypatch.setattr(cli_mod, "run_pipeline", fake_run_pipeline)
    argv = [
        "cli", "--project", "test", "--sample-root", str(sample_root),
        "--project-root", str(project_root), "--run-id", run_id, "--allow-dirty",
    ]
    if allow_partial:
        argv.append("--allow-partial-run")
    monkeypatch.setattr("sys.argv", argv)

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        if expected_exit:
            with pytest.raises(SystemExit) as exc_info:
                cli_mod.main()
            assert exc_info.value.code == expected_exit
        else:
            cli_mod.main()
    finally:
        os.chdir(cwd)

    outer = json.loads((project_root / "runs" / run_id / "manifest.json").read_text())
    assert outer["modules_run"] == ["qc"]
    assert outer["requested_modules"] == ["qc", "annotation"]
    assert outer["planned_modules"] == ["qc", "annotation"]
    assert outer["executed_modules"] == ["qc", "annotation"]
    assert outer["completed_modules"] == ["qc"]
    assert outer["failed_modules"] == ["annotation"]
    assert outer["overall_status"] == "failed"
    assert outer["bundle_sha256"] == "abc123"
    assert outer["factory_python"] == {
        "sha": "aaaaaaa", "dirty": False, "diff_sha256": ""
    }
    assert outer["producer_manifest"].endswith("run_manifest.json")
    assert outer["extra"]["allow_partial_run"] is allow_partial


def test_pack_run_for_mac_project_root_tars_python_bundle(tmp_path):
    """Project-root mode packages <run>/python/bundle, not legacy r_bundle."""
    root = Path(__file__).resolve().parent.parent
    run_id = "2026-05-14T0001Z-aaaaaaa"
    project_root = tmp_path / "project"
    run_dir = project_root / "runs" / run_id / "python"
    bundle_dir = run_dir / "bundle"
    bundle_dir.mkdir(parents=True)
    (bundle_dir / "bundle_manifest.json").write_text("{}\n", encoding="utf-8")

    proc = subprocess.run(
        [
            "bash",
            str(root / "scripts" / "pack_run_for_mac.sh"),
            "--project-root",
            str(project_root),
            "--run-id",
            run_id,
        ],
        cwd=root,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert proc.returncode == 0, (
        f"pack_run_for_mac failed\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    )

    tarball = run_dir / "mac_assets" / "bundle.tgz"
    assert tarball.exists()
    listing = subprocess.check_output(["tar", "-tzf", str(tarball)], text=True)
    assert "bundle/bundle_manifest.json" in listing.splitlines()
    assert "r_bundle/bundle_manifest.json" not in listing.splitlines()
