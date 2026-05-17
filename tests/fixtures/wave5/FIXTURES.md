---
fixtures_id: wave5_pre_v4
authored_at: 2026-05-16
plan_revision: v4.2
---

# Wave-5 v4.2 pre-v4 ledger fixtures

These fixtures validate that the v4.2 ledger schema's `allOf+if-then` conditional preserves backward-compatibility for pre-v4 Wave-5-shaped ledgers (they lack `plan_revision: "v4.2"` so the v4 conditional branch never fires).

## Intent

Plan v4.2 introduces a conditional block of additional required fields (biology-aware validation, attestation, visual verdicts) that only activates when `plan_revision == "v4.2"`. Ledgers authored before that revision must continue to validate against the unconditional baseline schema. These two fixtures are pre-v4 Wave-5-shaped ledgers; they omit `plan_revision` so the `if` clause does not match and the additional `then.required` block is skipped.

An earlier draft of this fixture set used NC2024 cohort-style (`nc2024_bh_*`) and Wave-2 preflight (`wave2_us018_*`) ledgers. Those FAILED the unconditional baseline because they are different ledger types (NC2024 cohort-style; Wave-2 preflight) — the v3 plan-changelog explicitly states "the Wave-5 schema deliberately does **not** reuse the NC2024 cohort-style ledger shape — they are sibling document types under one directory." Pre-v4 fixtures must therefore be Wave-5-shaped. Only the v3-era Wave-5 ledger qualifies on disk; a second fixture is derived from it by mutating `run_id` to provide a second exemplar.

## Fixture table

| fixture_path | source_path | sha256 | wave | plan_revision_class |
|---|---|---|---|---|
| `tests/fixtures/wave5/pre_v4/wave5_trevino_v3_pre_v4.json` | `ops/run_ledger/wave5_trevino_20260516T0931Z-d192836f1bb0.json` (verbatim copy) | `119bfd438e110f1a897ff58d55008d101f4e87d10e463550e15d19d9c9f6de38` | 5 | pre-v4 (v3-era) |
| `tests/fixtures/wave5/pre_v4/wave5_trevino_v3_pre_v4_synthetic2.json` | derived from same source; `run_id` mutated to `20260516T0608Z-d192836f1bb1` (pure hex, schema-pattern compliant) to provide a second Wave-5-shaped pre-v4 fixture | computed in-repo via `sha256sum` (not pinned in this doc to allow regeneration) | 5 | pre-v4 (synthetic) |

## Validation property

Both fixtures MUST validate cleanly against `ops/run_ledger/schema/wave5_v4_2.schema.json` (Draft 2020-12, `allOf` + `if-then` keyed on `plan_revision`).

In-repo validation command:

```bash
conda run -n sc_gpu_stable python -c "
import json, jsonschema
schema = json.load(open('ops/run_ledger/schema/wave5_v4_2.schema.json'))
v = jsonschema.Draft202012Validator(schema)
for p in ['tests/fixtures/wave5/pre_v4/wave5_trevino_v3_pre_v4.json',
          'tests/fixtures/wave5/pre_v4/wave5_trevino_v3_pre_v4_synthetic2.json']:
    d = json.load(open(p))
    errs = list(v.iter_errors(d))
    print(p, 'PASS' if not errs else f'FAIL: {[e.message for e in errs[:3]]}')
"
```

Last verified: 2026-05-16 — both fixtures PASS.

## Notes

- The schema under test is `ops/run_ledger/schema/wave5_v4_2.schema.json`.
- The legacy v3 schema `ops/run_ledger/schema/wave5.schema.json` is untouched.
- Do not modify the v3-era fixture; the synthetic variant exists specifically to provide a second sample.
