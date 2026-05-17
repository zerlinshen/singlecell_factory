"""Cross-schema mutual rejection tests (plan v5.0.4 §3.3 + §3.8).

Verifies the B1 fix invariant: a v4.2-revision ledger record MUST fail validation
against `wave5_v5_0.schema.json`, and a v5.0-revision ledger record MUST fail
validation against `wave5_v4_2.schema.json`. Mutual rejection is enforced
impossible-by-construction via:
  - both schemas list `plan_revision` in `required`
  - each schema's `plan_revision.enum` excludes the other family's values
"""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_DIR = REPO_ROOT / "ops" / "run_ledger" / "schema"
V42_SCHEMA_PATH = SCHEMA_DIR / "wave5_v4_2.schema.json"
V50_SCHEMA_PATH = SCHEMA_DIR / "wave5_v5_0.schema.json"

V42_LEDGER_PATH = REPO_ROOT / "ops" / "run_ledger" / "wave5_trevino_20260516T0931Z-d192836f1bb0.v4.2.json"
V50_FIXTURE_PATH = REPO_ROOT / "tests" / "fixtures" / "wave5" / "v5_0_valid_minimal.json"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def v42_schema() -> dict:
    return _load_json(V42_SCHEMA_PATH)


@pytest.fixture(scope="module")
def v50_schema() -> dict:
    return _load_json(V50_SCHEMA_PATH)


@pytest.fixture(scope="module")
def v42_record() -> dict:
    return _load_json(V42_LEDGER_PATH)


@pytest.fixture(scope="module")
def v50_record() -> dict:
    return _load_json(V50_FIXTURE_PATH)


def test_v42_record_validates_against_v42_schema(v42_schema, v42_record):
    """Sanity: B1 patch (adding plan_revision to required) does not break the existing v4.2 canonical ledger."""
    jsonschema.Draft202012Validator(v42_schema).validate(v42_record)


def test_v50_record_validates_against_v50_schema(v50_schema, v50_record):
    jsonschema.Draft202012Validator(v50_schema).validate(v50_record)


def test_v42_record_REJECTED_by_v50_schema(v50_schema, v42_record):
    """A v4.2-revision record MUST fail v5.0 schema validation (mutual rejection)."""
    validator = jsonschema.Draft202012Validator(v50_schema)
    errors = list(validator.iter_errors(v42_record))
    assert errors, "v4.2 ledger should fail v5.0 schema validation (plan_revision enum mismatch + missing required v5 fields)"
    # plan_revision='v4.2' must violate the v5.0 enum constraint
    enum_errors = [e for e in errors if "plan_revision" in [str(p) for p in e.path] and e.validator == "enum"]
    assert enum_errors, f"expected plan_revision enum violation; got errors: {[(e.validator, list(e.path), e.message) for e in errors]}"


def test_v50_record_REJECTED_by_v42_schema(v42_schema, v50_record):
    """A v5.0-revision record MUST fail v4.2 schema validation (mutual rejection)."""
    validator = jsonschema.Draft202012Validator(v42_schema)
    errors = list(validator.iter_errors(v50_record))
    assert errors, "v5.0 ledger should fail v4.2 schema validation (plan_revision enum mismatch)"
    enum_errors = [e for e in errors if "plan_revision" in [str(p) for p in e.path] and e.validator == "enum"]
    assert enum_errors, f"expected plan_revision enum violation; got errors: {[(e.validator, list(e.path), e.message) for e in errors]}"


def test_missing_plan_revision_REJECTED_by_v42_schema(v42_schema):
    """B1 patch: a record missing plan_revision MUST fail v4.2 schema (newly added to required)."""
    rec = {
        "wave": "5",
        "run_id": "20260516T0931Z-d192836f1bb0",
        "project_root": "/tmp",
        "run_manifest_path": "/tmp/run_manifest.json",
        "module_status_path": "/tmp/module_status.csv",
        "conda_env_hash": "0" * 64,
        "geo_sha256s": {},
        "metrics": {},
        "regression_gate_exit_code": 0,
        "dois_by_module": {},
        "created_at": "2026-05-17T08:00:00Z",
    }
    validator = jsonschema.Draft202012Validator(v42_schema)
    errors = list(validator.iter_errors(rec))
    required_errors = [e for e in errors if e.validator == "required" and "plan_revision" in e.message]
    assert required_errors, f"expected required:plan_revision error; got: {[(e.validator, e.message) for e in errors]}"


def test_missing_plan_revision_REJECTED_by_v50_schema(v50_schema):
    """A record missing plan_revision MUST fail v5.0 schema (required)."""
    rec = {
        "wave": "5",
        "run_id": "20260516T0931Z-d192836f1bb0",
        "project_root": "/tmp",
        "run_manifest_path": "/tmp/run_manifest.json",
        "conda_env_hash": "0" * 64,
        "created_at": "2026-05-17T08:00:00Z",
        "parent_canonical_run_id": "20260516T0931Z-d192836f1bb0",
        "parent_v42_ledger_sha": "0" * 64,
        "figure_outputs": {},
    }
    validator = jsonschema.Draft202012Validator(v50_schema)
    errors = list(validator.iter_errors(rec))
    required_errors = [e for e in errors if e.validator == "required" and "plan_revision" in e.message]
    assert required_errors, f"expected required:plan_revision error; got: {[(e.validator, e.message) for e in errors]}"


def test_dispatcher_selects_v42_for_v42_record():
    from scripts.ci.select_ledger_schema import select_schema_path, V42_SCHEMA, V50_SCHEMA
    rec = {"plan_revision": "v4.2"}
    assert select_schema_path(rec) == V42_SCHEMA


def test_dispatcher_selects_v50_for_v50_record():
    from scripts.ci.select_ledger_schema import select_schema_path, V42_SCHEMA, V50_SCHEMA
    rec = {"plan_revision": "v5.0"}
    assert select_schema_path(rec) == V50_SCHEMA


def test_dispatcher_rejects_missing_revision():
    from scripts.ci.select_ledger_schema import select_schema_path
    with pytest.raises(ValueError, match="plan_revision"):
        select_schema_path({})


def test_dispatcher_rejects_unknown_revision():
    from scripts.ci.select_ledger_schema import select_schema_path
    with pytest.raises(ValueError, match="unknown plan_revision"):
        select_schema_path({"plan_revision": "v9.9"})
