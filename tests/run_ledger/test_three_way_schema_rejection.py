"""Three-way mutual schema rejection tests (plan v5.1 §4.2 AC-V51-SCHEMA-1).

Verifies the three-way invariant:
  - A v4.2-revision record MUST fail v5.0 schema validation
  - A v4.2-revision record MUST fail v5.1 schema validation
  - A v5.0-revision record MUST fail v4.2 schema validation
  - A v5.0-revision record MUST fail v5.1 schema validation
  - A v5.1-revision record MUST fail v4.2 schema validation
  - A v5.1-revision record MUST fail v5.0 schema validation
  - Each schema MUST pass its own native record (sanity)
  - Dispatcher routes each plan_revision to the correct schema

Total: 12 test cases (3 schemas x 4 cross-checks per AC-V51-SCHEMA-1).
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
V51_SCHEMA_PATH = SCHEMA_DIR / "wave5_v5_1.schema.json"

V42_LEDGER_PATH = REPO_ROOT / "ops" / "run_ledger" / "wave5_trevino_20260516T0931Z-d192836f1bb0.v4.2.json"
V50_FIXTURE_PATH = REPO_ROOT / "tests" / "fixtures" / "wave5" / "v5_0_valid_minimal.json"
V51_FIXTURE_PATH = REPO_ROOT / "tests" / "fixtures" / "wave5" / "v5_1_valid_minimal.json"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def v42_schema() -> dict:
    return _load_json(V42_SCHEMA_PATH)


@pytest.fixture(scope="module")
def v50_schema() -> dict:
    return _load_json(V50_SCHEMA_PATH)


@pytest.fixture(scope="module")
def v51_schema() -> dict:
    return _load_json(V51_SCHEMA_PATH)


@pytest.fixture(scope="module")
def v42_record() -> dict:
    return _load_json(V42_LEDGER_PATH)


@pytest.fixture(scope="module")
def v50_record() -> dict:
    return _load_json(V50_FIXTURE_PATH)


@pytest.fixture(scope="module")
def v51_record() -> dict:
    return _load_json(V51_FIXTURE_PATH)


# ── Sanity: each schema accepts its own native record ──────────────────────

def test_v42_record_validates_against_v42_schema(v42_schema, v42_record):
    """Sanity: v4.2 canonical ledger passes v4.2 schema."""
    jsonschema.Draft202012Validator(v42_schema).validate(v42_record)


def test_v50_record_validates_against_v50_schema(v50_schema, v50_record):
    """Sanity: v5.0 minimal fixture passes v5.0 schema."""
    jsonschema.Draft202012Validator(v50_schema).validate(v50_record)


def test_v51_record_validates_against_v51_schema(v51_schema, v51_record):
    """Sanity: v5.1 minimal fixture passes v5.1 schema."""
    jsonschema.Draft202012Validator(v51_schema).validate(v51_record)


# ── v4.2 record rejected by v5.0 and v5.1 schemas ─────────────────────────

def test_v42_record_REJECTED_by_v50_schema(v50_schema, v42_record):
    """v4.2 record MUST fail v5.0 schema (plan_revision enum mismatch)."""
    validator = jsonschema.Draft202012Validator(v50_schema)
    errors = list(validator.iter_errors(v42_record))
    assert errors, "v4.2 ledger should fail v5.0 schema validation"
    enum_errors = [e for e in errors if "plan_revision" in [str(p) for p in e.path] and e.validator == "enum"]
    assert enum_errors, f"expected plan_revision enum violation; got: {[(e.validator, list(e.path), e.message) for e in errors]}"


def test_v42_record_REJECTED_by_v51_schema(v51_schema, v42_record):
    """v4.2 record MUST fail v5.1 schema (plan_revision enum mismatch)."""
    validator = jsonschema.Draft202012Validator(v51_schema)
    errors = list(validator.iter_errors(v42_record))
    assert errors, "v4.2 ledger should fail v5.1 schema validation"
    enum_errors = [e for e in errors if "plan_revision" in [str(p) for p in e.path] and e.validator == "enum"]
    assert enum_errors, f"expected plan_revision enum violation; got: {[(e.validator, list(e.path), e.message) for e in errors]}"


# ── v5.0 record rejected by v4.2 and v5.1 schemas ─────────────────────────

def test_v50_record_REJECTED_by_v42_schema(v42_schema, v50_record):
    """v5.0 record MUST fail v4.2 schema (plan_revision enum mismatch)."""
    validator = jsonschema.Draft202012Validator(v42_schema)
    errors = list(validator.iter_errors(v50_record))
    assert errors, "v5.0 ledger should fail v4.2 schema validation"
    enum_errors = [e for e in errors if "plan_revision" in [str(p) for p in e.path] and e.validator == "enum"]
    assert enum_errors, f"expected plan_revision enum violation; got: {[(e.validator, list(e.path), e.message) for e in errors]}"


def test_v50_record_REJECTED_by_v51_schema(v51_schema, v50_record):
    """v5.0 record MUST fail v5.1 schema (plan_revision enum excludes v5.0)."""
    validator = jsonschema.Draft202012Validator(v51_schema)
    errors = list(validator.iter_errors(v50_record))
    assert errors, "v5.0 ledger should fail v5.1 schema validation"
    enum_errors = [e for e in errors if "plan_revision" in [str(p) for p in e.path] and e.validator == "enum"]
    assert enum_errors, f"expected plan_revision enum violation; got: {[(e.validator, list(e.path), e.message) for e in errors]}"


# ── v5.1 record rejected by v4.2 and v5.0 schemas ─────────────────────────

def test_v51_record_REJECTED_by_v42_schema(v42_schema, v51_record):
    """v5.1 record MUST fail v4.2 schema (plan_revision enum mismatch)."""
    validator = jsonschema.Draft202012Validator(v42_schema)
    errors = list(validator.iter_errors(v51_record))
    assert errors, "v5.1 ledger should fail v4.2 schema validation"
    enum_errors = [e for e in errors if "plan_revision" in [str(p) for p in e.path] and e.validator == "enum"]
    assert enum_errors, f"expected plan_revision enum violation; got: {[(e.validator, list(e.path), e.message) for e in errors]}"


def test_v51_record_REJECTED_by_v50_schema(v50_schema, v51_record):
    """v5.1 record MUST fail v5.0 schema (plan_revision mutually rejected via not-clause)."""
    validator = jsonschema.Draft202012Validator(v50_schema)
    errors = list(validator.iter_errors(v51_record))
    assert errors, "v5.1 ledger should fail v5.0 schema validation"
    # v5.0 schema rejects v5.1 via a `not` clause (plan_revision enum still lists v5.1 for
    # backwards-compat read; mutual rejection is enforced by the allOf[not] guard instead).
    rejection_errors = [
        e for e in errors
        if e.validator in ("enum", "not")
    ]
    assert rejection_errors, f"expected plan_revision rejection; got: {[(e.validator, list(e.path), e.message) for e in errors]}"


# ── Dispatcher routing tests ───────────────────────────────────────────────

def test_dispatcher_selects_v42_for_v42_record():
    from scripts.ci.select_ledger_schema import select_schema_path, V42_SCHEMA
    assert select_schema_path({"plan_revision": "v4.2"}) == V42_SCHEMA


def test_dispatcher_selects_v50_for_v50_record():
    from scripts.ci.select_ledger_schema import select_schema_path, V50_SCHEMA
    assert select_schema_path({"plan_revision": "v5.0"}) == V50_SCHEMA


def test_dispatcher_selects_v51_for_v51_record():
    from scripts.ci.select_ledger_schema import select_schema_path, V51_SCHEMA
    assert select_schema_path({"plan_revision": "v5.1"}) == V51_SCHEMA


def test_dispatcher_rejects_missing_revision():
    from scripts.ci.select_ledger_schema import select_schema_path
    with pytest.raises(ValueError, match="plan_revision"):
        select_schema_path({})


def test_dispatcher_rejects_unknown_revision():
    from scripts.ci.select_ledger_schema import select_schema_path
    with pytest.raises(ValueError, match="unknown plan_revision"):
        select_schema_path({"plan_revision": "v9.9"})
