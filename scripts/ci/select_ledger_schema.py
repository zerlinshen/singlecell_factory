#!/usr/bin/env python3
"""Canonical ledger schema dispatcher (Arch-2 Fix #1, plan v5.0.4 §3.3; v5.1 Stage 0.1 three-way dispatch).

Selects the appropriate JSON schema for a Wave-5 ledger record based on its
`plan_revision` field, and exposes a single entry-point for downstream gates +
tests so the choice of schema is not duplicated across call sites.

Schema mapping:
  plan_revision in {v3, v4, v4.1, v4.2}  -> ops/run_ledger/schema/wave5_v4_2.schema.json
  plan_revision in {v5.0}                 -> ops/run_ledger/schema/wave5_v5_0.schema.json
  plan_revision in {v5.1}                 -> ops/run_ledger/schema/wave5_v5_1.schema.json

Records missing `plan_revision` are NOT dispatched; all schemas now require
`plan_revision` (B1 mutual-rejection patch). Caller receives ``ValueError``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_DIR = REPO_ROOT / "ops" / "run_ledger" / "schema"

V42_SCHEMA = SCHEMA_DIR / "wave5_v4_2.schema.json"
V50_SCHEMA = SCHEMA_DIR / "wave5_v5_0.schema.json"
V51_SCHEMA = SCHEMA_DIR / "wave5_v5_1.schema.json"

V42_REVISIONS = frozenset({"v3", "v4", "v4.1", "v4.2"})
V50_REVISIONS = frozenset({"v5.0"})
V51_REVISIONS = frozenset({"v5.1"})


def select_schema_path(record: dict) -> Path:
    """Return absolute Path to the schema that governs ``record``.

    Raises ValueError if plan_revision is missing or unknown.
    """
    rev = record.get("plan_revision")
    if rev is None:
        raise ValueError(
            "ledger record missing required 'plan_revision' field; "
            "cannot dispatch schema (all schemas require plan_revision post-B1 patch)"
        )
    if rev in V42_REVISIONS:
        return V42_SCHEMA
    if rev in V50_REVISIONS:
        return V50_SCHEMA
    if rev in V51_REVISIONS:
        return V51_SCHEMA
    raise ValueError(f"unknown plan_revision={rev!r}; no schema mapping")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("ledger_path", type=Path, help="Path to ledger JSON file")
    p.add_argument("--print-schema-path-only", action="store_true",
                   help="Print only the schema path (stdout); exit 0 on success")
    args = p.parse_args(argv)
    try:
        record = json.loads(args.ledger_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: failed to read/parse ledger: {exc}", file=sys.stderr)
        return 2
    try:
        schema_path = select_schema_path(record)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(schema_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
