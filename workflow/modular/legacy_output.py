"""Govern the deprecated factory-local output compatibility path.

The legacy path is deliberately still callable, but every real launch must leave
an external JSONL access event.  Its removal is a separate, explicit decision:
telemetry silence alone is insufficient without completed migration and owner
approval recorded by the contract.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Mapping

import yaml

from ..factory_paths import SINGLECELL_FACTORY_ROOT


CONTRACT_PATH = (
    SINGLECELL_FACTORY_ROOT / "contracts" / "legacy_factory_output_deprecation.yaml"
)
LEGACY_DEFAULT_OUTPUT_DIR = SINGLECELL_FACTORY_ROOT / "results"


def load_legacy_output_contract(path: Path | None = None) -> dict:
    """Load and minimally validate the deprecation contract."""
    contract_path = Path(path or CONTRACT_PATH)
    payload = yaml.safe_load(contract_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"legacy output contract must be a mapping: {contract_path}")

    expected = {
        "schema_version": "legacy_factory_output_deprecation/1.0",
        "status": "deprecated_compatibility",
        "replacement": "--project-root",
    }
    for key, value in expected.items():
        if payload.get(key) != value:
            raise ValueError(
                f"legacy output contract {key!r} must be {value!r}: {contract_path}"
            )

    telemetry = payload.get("telemetry")
    retirement = payload.get("retirement")
    protection = payload.get("protection")
    if not isinstance(telemetry, dict) or telemetry.get("format") != "jsonl":
        raise ValueError("legacy output contract requires JSONL telemetry")
    if not isinstance(retirement, dict) or retirement.get("access_window_days") != 14:
        raise ValueError("legacy output contract requires a 14-day access window")
    if (
        not isinstance(protection, dict)
        or protection.get("remove_before_eligible") != "forbidden"
    ):
        raise ValueError("legacy output contract must forbid early removal")
    return payload


def resolve_legacy_access_log(
    environ: Mapping[str, str] | None = None,
) -> Path:
    """Resolve the operational telemetry path outside the factory checkout."""
    env = os.environ if environ is None else environ
    configured = env.get("SC_LEGACY_OUTPUT_ACCESS_LOG", "").strip()
    if configured:
        path = Path(configured).expanduser().resolve()
    else:
        state_root = Path(
            env.get("XDG_STATE_HOME", str(Path.home() / ".local" / "state"))
        ).expanduser()
        path = (
            state_root / "singlecell_factory" / "legacy_output_access.jsonl"
        ).resolve()

    try:
        path.relative_to(SINGLECELL_FACTORY_ROOT.resolve())
    except ValueError:
        return path
    raise ValueError(
        "legacy access telemetry must live outside the singlecell_factory checkout"
    )


def _utc_timestamp(now: datetime | None = None) -> str:
    instant = now or datetime.now(timezone.utc)
    if instant.tzinfo is None:
        instant = instant.replace(tzinfo=timezone.utc)
    return instant.astimezone(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )


def record_legacy_output_access(
    *,
    output_dir: Path,
    project: str,
    environ: Mapping[str, str] | None = None,
    now: datetime | None = None,
) -> Path:
    """Append one privacy-bounded access event and return its JSONL path."""
    contract = load_legacy_output_contract()
    telemetry_path = resolve_legacy_access_log(environ)
    resolved_output = Path(output_dir).expanduser().resolve()
    try:
        resolved_output.relative_to(SINGLECELL_FACTORY_ROOT.resolve())
        factory_local = True
    except ValueError:
        factory_local = False

    event = {
        "contract": contract["schema_version"],
        "event": contract["telemetry"]["event"],
        "factory_local": factory_local,
        "output_dir": str(resolved_output),
        "project": str(project),
        "recorded_at": _utc_timestamp(now),
        "replacement": contract["replacement"],
    }
    telemetry_path.parent.mkdir(parents=True, exist_ok=True)
    with telemetry_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, sort_keys=True) + "\n")
    return telemetry_path


def legacy_output_warning_message(telemetry_path: Path) -> str:
    """Return the one canonical warning shared by launch surfaces."""
    contract = load_legacy_output_contract()
    days = contract["retirement"]["access_window_days"]
    return (
        f"legacy output is governed by {contract['schema_version']}; migrate to "
        f"{contract['replacement']}. Removal remains forbidden until migration and "
        f"owner approval are recorded and a {days}-day silent access window passes. "
        f"This access was recorded at {Path(telemetry_path).resolve()}"
    )


def _parse_recorded_at(value: object) -> datetime:
    if not isinstance(value, str):
        raise ValueError("recorded_at must be an ISO-8601 string")
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def evaluate_legacy_output_retirement(
    contract: Mapping[str, object],
    telemetry_path: Path,
    *,
    now: datetime | None = None,
) -> dict[str, object]:
    """Evaluate, but never mutate, the compatibility-path retirement state."""
    retirement = contract.get("retirement")
    if not isinstance(retirement, Mapping):
        raise ValueError("contract retirement block must be a mapping")
    days = int(retirement.get("access_window_days", 0))
    if days <= 0:
        raise ValueError("access_window_days must be positive")

    instant = now or datetime.now(timezone.utc)
    if instant.tzinfo is None:
        instant = instant.replace(tzinfo=timezone.utc)
    cutoff = instant.astimezone(timezone.utc) - timedelta(days=days)
    recent_access_count = 0
    invalid_records = 0
    path = Path(telemetry_path)
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                event = json.loads(line)
                if _parse_recorded_at(event.get("recorded_at")) >= cutoff:
                    recent_access_count += 1
            except (json.JSONDecodeError, TypeError, ValueError):
                invalid_records += 1

    blockers: list[str] = []
    monitoring_started_at = retirement.get("monitoring_started_at")
    if monitoring_started_at is None:
        blockers.append("telemetry monitoring start is not recorded")
    else:
        try:
            monitoring_start = _parse_recorded_at(monitoring_started_at)
        except ValueError:
            blockers.append("telemetry monitoring start is invalid")
        else:
            if instant.astimezone(timezone.utc) < monitoring_start + timedelta(days=days):
                blockers.append(f"{days}-day monitoring window has not elapsed")
    if not retirement.get("project_root_migration_complete", False):
        blockers.append("project-root migration is not complete")
    if not retirement.get("owner_approval_recorded", False):
        blockers.append("owner approval is not recorded")
    if recent_access_count:
        blockers.append(f"{days}-day access window is not silent")
    if invalid_records:
        blockers.append("telemetry contains invalid records")

    return {
        "eligible": not blockers,
        "recent_access_count": recent_access_count,
        "invalid_record_count": invalid_records,
        "access_window_days": days,
        "blockers": blockers,
    }


def main(argv: list[str] | None = None) -> int:
    """Provide a dependency-light recorder for shell compatibility wrappers."""
    parser = argparse.ArgumentParser(
        description="Record a governed legacy-output compatibility access."
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project", required=True)
    args = parser.parse_args(argv)
    telemetry_path = record_legacy_output_access(
        output_dir=Path(args.output_dir),
        project=args.project,
    )
    print(
        f"DeprecationWarning: {legacy_output_warning_message(telemetry_path)}",
        file=sys.stderr,
    )
    print(telemetry_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
