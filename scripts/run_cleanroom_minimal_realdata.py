#!/usr/bin/env python3
"""Run a bounded clean-room real-data gate for the three-factory pipeline.

The gate does not re-download raw data and does not mutate project roots.  It
creates a temporary isolated workspace, links the current real-data project
roots as read-only inputs, verifies required artifacts, and optionally executes
existing no-rerun validators from that workspace with explicit absolute paths.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import validate_nc2024_architecture_contract as nc_validator
import validate_two_realdata_final as two_validator


DEFAULT_REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    DEFAULT_REPO_ROOT
    / "ops"
    / "governance_records"
    / "2026-05-18-two-real-dataset-final-validation"
    / "cleanroom_minimal_realdata.json"
)


def _path_status(path: Path, label: str, *, expect_dir: bool | None = None) -> dict[str, Any]:
    exists = path.exists()
    ok_type = exists and (expect_dir is None or (path.is_dir() if expect_dir else path.is_file()))
    readable = bool(ok_type and os.access(path, os.R_OK))
    payload: dict[str, Any] = {
        "label": label,
        "path": str(path),
        "exists": exists,
        "readable": readable,
        "ok": readable,
    }
    if exists and path.is_file():
        payload["bytes"] = path.stat().st_size
    if exists and path.is_dir():
        payload["entries"] = sum(1 for _ in path.iterdir())
    if not exists:
        payload["reason"] = "missing"
    elif not ok_type:
        payload["reason"] = "wrong_type"
    elif not readable:
        payload["reason"] = "not_readable"
    return payload


def _link_input(target: Path, workspace_inputs: Path, name: str) -> dict[str, Any]:
    link = workspace_inputs / name
    target = target.resolve()
    try:
        link.symlink_to(target, target_is_directory=target.is_dir())
        mode = "symlink"
    except OSError:
        # Some CI filesystems disable symlinks.  The fallback remains read-only:
        # we record the absolute path and never copy/mutate large artifacts.
        link.write_text(str(target) + "\n", encoding="utf-8")
        mode = "path_pointer"
    return {"name": name, "mode": mode, "workspace_path": str(link), "target": str(target)}


def _run_command(cmd: list[str], *, cwd: Path, timeout: int) -> dict[str, Any]:
    started = datetime.now(timezone.utc).isoformat()
    proc = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, timeout=timeout)
    return {
        "cmd": cmd,
        "cwd": str(cwd),
        "started_at_utc": started,
        "returncode": proc.returncode,
        "ok": proc.returncode == 0,
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
    }


def _required_artifacts(
    *,
    nc_project_root: Path,
    nc_run_id: str,
    cell_project_root: Path,
    cell_pipeline_run_id: str,
    cell_evidence_run_id: str,
) -> list[dict[str, Any]]:
    nc_run = nc_project_root / "runs" / nc_run_id
    cell_run = cell_project_root / "runs" / cell_pipeline_run_id
    cell_evidence = cell_project_root / "runs" / cell_evidence_run_id
    checks = [
        _path_status(nc_project_root, "nc_project_root", expect_dir=True),
        _path_status(nc_run, "nc_run_dir", expect_dir=True),
        _path_status(nc_run / "manifest.json", "nc_root_manifest", expect_dir=False),
        _path_status(
            nc_run / "python" / "bundle" / "bundle_manifest.json",
            "nc_bundle_manifest",
            expect_dir=False,
        ),
        _path_status(
            nc_run / "python" / "bundle" / "provenance.json",
            "nc_bundle_provenance",
            expect_dir=False,
        ),
        _path_status(nc_run / "r", "nc_r_output_dir", expect_dir=True),
        _path_status(cell_project_root, "cell_project_root", expect_dir=True),
        _path_status(cell_run, "cell_pipeline_run_dir", expect_dir=True),
        _path_status(
            cell_run / "python" / "bundle" / "bundle_manifest.json",
            "cell_bundle_manifest",
            expect_dir=False,
        ),
        _path_status(
            cell_run / "python" / "bundle" / "provenance.json",
            "cell_bundle_provenance",
            expect_dir=False,
        ),
        _path_status(cell_run / "r", "cell_r_output_dir", expect_dir=True),
        _path_status(cell_evidence, "cell_evidence_run_dir", expect_dir=True),
        _path_status(
            cell_evidence
            / "python"
            / "figure_reproduction_evidence"
            / "all_figures_final_quality_gate_review.json",
            "cell_figure_quality_gate",
            expect_dir=False,
        ),
    ]
    return checks


def run_cleanroom_gate(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = args.repo_root.resolve()
    output_dir = args.output.parent if args.output else None
    output_dir and output_dir.mkdir(parents=True, exist_ok=True)

    workspace_parent = args.workspace_parent.resolve() if args.workspace_parent else None
    temp_kwargs: dict[str, Any] = {"prefix": "scf-cleanroom-"}
    if workspace_parent:
        workspace_parent.mkdir(parents=True, exist_ok=True)
        temp_kwargs["dir"] = str(workspace_parent)

    with tempfile.TemporaryDirectory(**temp_kwargs) as tmp:
        workspace = Path(tmp).resolve()
        inputs = workspace / "inputs"
        inputs.mkdir()
        linked_inputs = [
            _link_input(args.nc_project_root, inputs, "nc_project_root"),
            _link_input(args.cell_project_root, inputs, "cell_project_root"),
            _link_input(repo_root, inputs, "singlecell_factory"),
        ]
        artifacts = _required_artifacts(
            nc_project_root=args.nc_project_root,
            nc_run_id=args.nc_run_id,
            cell_project_root=args.cell_project_root,
            cell_pipeline_run_id=args.cell_pipeline_run_id,
            cell_evidence_run_id=args.cell_evidence_run_id,
        )
        commands: list[dict[str, Any]] = []
        if args.run_validators:
            common_env_note = {
                "cwd_is_cleanroom": str(workspace),
                "repo_root": str(repo_root),
                "boundary": "no-rerun validators only; no raw-data download or project-root writes",
            }
            commands.append({"note": common_env_note})
            commands.append(
                _run_command(
                    [
                        sys.executable,
                        str(repo_root / "scripts" / "validate_nc2024_architecture_contract.py"),
                        "--project-root",
                        str(args.nc_project_root),
                        "--run-id",
                        args.nc_run_id,
                        "--verify-sha",
                    ],
                    cwd=workspace,
                    timeout=args.timeout,
                )
            )
            commands.append(
                _run_command(
                    [
                        sys.executable,
                        str(repo_root / "scripts" / "validate_two_realdata_final.py"),
                        "--nc-project-root",
                        str(args.nc_project_root),
                        "--nc-run-id",
                        args.nc_run_id,
                        "--cell-project-root",
                        str(args.cell_project_root),
                        "--cell-pipeline-run-id",
                        args.cell_pipeline_run_id,
                        "--cell-evidence-run-id",
                        args.cell_evidence_run_id,
                        "--verify-sha",
                    ],
                    cwd=workspace,
                    timeout=args.timeout,
                )
            )

        missing = [c for c in artifacts if not c["ok"]]
        failed_commands = [c for c in commands if c.get("ok") is False]
        if not missing and not failed_commands:
            verdict = "pass"
        elif (
            args.allow_missing_inputs
            and not args.run_validators
            and missing
            and not failed_commands
        ):
            verdict = "conditional"
        else:
            verdict = "fail"
        report = {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "verdict": verdict,
            "execution_mode": "cleanroom_minimal_realdata_no_raw_rerun",
            "scientific_boundary": (
                "This validates the minimal real-data control plane from an isolated workspace; "
                "it is not a full raw FASTQ/fragments/BPNet rerun."
            ),
            "workspace": str(workspace),
            "workspace_preserved": bool(args.preserve_workspace),
            "linked_inputs": linked_inputs,
            "environment": {
                "python": sys.version,
                "executable": sys.executable,
                "platform": platform.platform(),
            },
            "artifact_checks": artifacts,
            "commands": commands,
        }
        if args.preserve_workspace:
            preserved = Path(tempfile.mkdtemp(prefix="scf-cleanroom-preserved-"))
            shutil.copytree(workspace, preserved, dirs_exist_ok=True)
            report["preserved_workspace"] = str(preserved)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(report, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT)
    parser.add_argument("--nc-project-root", type=Path, default=nc_validator.DEFAULT_PROJECT_ROOT)
    parser.add_argument("--nc-run-id", default=nc_validator.DEFAULT_RUN_ID)
    parser.add_argument(
        "--cell-project-root",
        type=Path,
        default=two_validator.DEFAULT_CELL_PROJECT_ROOT,
    )
    parser.add_argument(
        "--cell-pipeline-run-id",
        default=two_validator.DEFAULT_CELL_PIPELINE_RUN_ID,
    )
    parser.add_argument(
        "--cell-evidence-run-id",
        default=two_validator.DEFAULT_CELL_EVIDENCE_RUN_ID,
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workspace-parent", type=Path, default=None)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--no-run-validators", dest="run_validators", action="store_false")
    parser.set_defaults(run_validators=True)
    parser.add_argument(
        "--allow-missing-inputs",
        action="store_true",
        help=(
            "Return a conditional data-free report when project roots are "
            "unavailable and validators are disabled."
        ),
    )
    parser.add_argument("--preserve-workspace", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = run_cleanroom_gate(args)
    print(json.dumps(report, indent=2, sort_keys=True))
    allowed = report["verdict"] == "pass" or (
        report["verdict"] == "conditional" and args.allow_missing_inputs
    )
    return 0 if allowed else 1


if __name__ == "__main__":
    raise SystemExit(main())
