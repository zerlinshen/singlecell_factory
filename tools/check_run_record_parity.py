"""Wave-5.5 v5.0.4 Stage B5 — factory ↔ ~/projects/ run-record parity walker.

# STATUS: one-off; promote to ops/ when Stream B graduates from warn-now to fail-after.

Per plan v5.0.4 §5 Stage B0+B5 + AC-V5-PARITY-1/2/3/WALKER-1:

- Walks ledger entries under ``ops/run_ledger/`` and corresponding run-record stubs under
  ``ops/run_records/<project-id>/<run-id>/``; reports parity gaps.
- Dir-vs-file class distinction:
    DIR-class:  run output directories (under PROJECTS_ROOT/<id>/runs/<run-id>/)
    FILE-class: ledger metadata files (under singlecell_factory/ops/run_ledger/)
- Asymmetric mode (per Q9 + Architect-iter-2 Fix #5 + AC-V5-PARITY-WALKER-1):
    factory side:  --mode=warn (initial; 3-clean-runs ratchet → fail-after)
    projects side: --mode=warn-only (always)
- Structured JSON output to ``.omc/logs/wave5_parity_walker_<date>.jsonl``.
- Warning-trend monotonically-non-increasing gate enforced by inspecting prior runs in
  the same .jsonl file (Architect-iter-2 Fix #5 / plan §3.5 ratchet).
- Exit codes:
    0 = parity OK or warn-mode green
    1 = parity warnings emitted (warn mode)
    2 = parity violations (fail mode, only when ratchet has graduated)
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
LEDGER_DIR = REPO_ROOT / "ops" / "run_ledger"
RUN_RECORDS_DIR = REPO_ROOT / "ops" / "run_records"
LOGS_DIR = REPO_ROOT / ".omc" / "logs"

# PROJECTS_ROOT — match singlecell_factory factory-project separation v2
PROJECTS_ROOT = Path("/home/zerlinshen/projects")


def _ts() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def enumerate_ledgers() -> list[dict]:
    """FILE-class: enumerate ledger JSON files."""
    entries = []
    for p in sorted(LEDGER_DIR.glob("*.json")):
        if p.parent.name != "run_ledger":
            continue
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        entries.append({
            "class": "FILE",
            "ledger_path": str(p),
            "ledger_name": p.name,
            "wave": data.get("wave"),
            "run_id": data.get("run_id"),
            "project_root": data.get("project_root"),
            "plan_revision": data.get("plan_revision"),
        })
    return entries


def enumerate_run_records() -> list[dict]:
    """Enumerate run-record stub dirs under ops/run_records/<project>/<run-id>/."""
    out = []
    if not RUN_RECORDS_DIR.exists():
        return out
    for project_dir in sorted(RUN_RECORDS_DIR.iterdir()):
        if not project_dir.is_dir():
            continue
        for run_dir in sorted(project_dir.iterdir()):
            if not run_dir.is_dir():
                continue
            out.append({
                "class": "DIR",
                "record_path": str(run_dir),
                "project": project_dir.name,
                "run_id": run_dir.name,
            })
    return out


def enumerate_project_runs() -> list[dict]:
    """Enumerate <PROJECTS_ROOT>/<project>/runs/<run-id>/ dirs (the real outputs)."""
    out = []
    if not PROJECTS_ROOT.exists():
        return out
    for project_dir in sorted(PROJECTS_ROOT.iterdir()):
        runs_dir = project_dir / "runs"
        if not runs_dir.exists():
            continue
        for run_dir in sorted(runs_dir.iterdir()):
            if not run_dir.is_dir():
                continue
            out.append({
                "class": "DIR",
                "project_run_path": str(run_dir),
                "project": project_dir.name,
                "run_id": run_dir.name,
            })
    return out


def compute_parity(ledgers, run_records, project_runs) -> dict:
    """Compute parity status with dir vs file class breakdown."""
    ledger_run_ids = {((l.get("project_root") or "").rstrip("/").split("/")[-1] or "?",
                       l.get("run_id") or "?"): l for l in ledgers}
    record_run_ids = {(r["project"], r["run_id"]): r for r in run_records}
    project_run_ids = {(r["project"], r["run_id"]): r for r in project_runs}

    # ledgers without record stubs (FACTORY → PROJECT parity)
    ledgers_missing_records = []
    for k, l in ledger_run_ids.items():
        if k not in record_run_ids:
            ledgers_missing_records.append({"key": k, "ledger": l["ledger_name"]})

    # project runs without ledger
    project_runs_missing_ledger = []
    for k, r in project_run_ids.items():
        if k not in ledger_run_ids:
            project_runs_missing_ledger.append({"key": k, "project_run_path": r["project_run_path"]})

    # record stubs without backing project run
    records_missing_project_run = []
    for k, r in record_run_ids.items():
        if k not in project_run_ids:
            records_missing_project_run.append({"key": k, "record_path": r["record_path"]})

    return {
        "summary": {
            "n_ledgers": len(ledgers),
            "n_record_stubs": len(run_records),
            "n_project_runs": len(project_runs),
            "n_ledgers_missing_records": len(ledgers_missing_records),
            "n_project_runs_missing_ledger": len(project_runs_missing_ledger),
            "n_records_missing_project_run": len(records_missing_project_run),
        },
        "ledgers_missing_records": ledgers_missing_records,
        "project_runs_missing_ledger": project_runs_missing_ledger,
        "records_missing_project_run": records_missing_project_run,
    }


def emit_jsonl(parity: dict, *, factory_mode: str, projects_mode: str) -> Path:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOGS_DIR / f"wave5_parity_walker_{datetime.now(timezone.utc).strftime('%Y-%m-%d')}.jsonl"
    entry = {
        "ts": _ts(),
        "factory_mode": factory_mode,
        "projects_mode": projects_mode,
        **parity,
    }
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry) + "\n")
    return log_path


def warning_trend_check(log_path: Path) -> tuple[bool, int]:
    """Architect-iter-2 Fix #5: warning-trend monotonically-non-increasing gate.

    Reads the last up-to-3 entries from the .jsonl; returns (trend_ok, current_warning_count).
    """
    if not log_path.exists():
        return True, 0
    lines = [l for l in log_path.read_text().splitlines() if l.strip()]
    last = lines[-3:] if len(lines) >= 1 else []
    counts = []
    for l in last:
        try:
            e = json.loads(l)
            n = e.get("summary", {}).get("n_ledgers_missing_records", 0)
            counts.append(n)
        except json.JSONDecodeError:
            continue
    if not counts:
        return True, 0
    cur = counts[-1]
    # monotonically non-increasing: every step must be <= previous
    trend_ok = all(counts[i] <= counts[i - 1] for i in range(1, len(counts)))
    return trend_ok, cur


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mode", choices=["warn", "fail"], default="warn",
                    help="Factory-side enforcement (default warn; 3-clean-runs ratchet → fail)")
    ap.add_argument("--projects-mode", choices=["warn-only"], default="warn-only",
                    help="Projects-side enforcement (always warn-only)")
    ap.add_argument("--no-log", action="store_true", help="Skip writing the .jsonl log entry")
    args = ap.parse_args(argv)

    ledgers = enumerate_ledgers()
    records = enumerate_run_records()
    project_runs = enumerate_project_runs()
    parity = compute_parity(ledgers, records, project_runs)

    if not args.no_log:
        log_path = emit_jsonl(parity, factory_mode=args.mode, projects_mode=args.projects_mode)
    else:
        log_path = LOGS_DIR / "wave5_parity_walker_<date>.jsonl"

    trend_ok, cur_warn = warning_trend_check(log_path)

    summary = parity["summary"]
    print(json.dumps({
        "factory_mode": args.mode,
        "projects_mode": args.projects_mode,
        "summary": summary,
        "warning_trend_ok": trend_ok,
        "log_path": str(log_path),
    }, indent=2))

    n_factory_warnings = summary["n_ledgers_missing_records"]
    n_projects_warnings = summary["n_project_runs_missing_ledger"]

    if n_factory_warnings == 0 and n_projects_warnings == 0:
        return 0
    if args.mode == "warn":
        if not trend_ok:
            print(f"WARN-TREND VIOLATION: warning count {cur_warn} > prior run; ratchet broken",
                  file=sys.stderr)
            return 1
        return 1 if (n_factory_warnings or n_projects_warnings) else 0
    # fail mode
    if n_factory_warnings:
        print(f"FAIL: {n_factory_warnings} ledger(s) missing factory run-record stub", file=sys.stderr)
        return 2
    return 1 if n_projects_warnings else 0


if __name__ == "__main__":
    raise SystemExit(main())
