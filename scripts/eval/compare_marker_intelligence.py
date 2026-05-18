#!/usr/bin/env python3
"""
compare_marker_intelligence.py — T8 / P1A.S7c delta computation + AC-1/AC-2 eval.

Consumes baseline and candidate run directories, validates both input JSON files
against the inline schema (fail-fast on mismatch), computes macro-F1 delta,
and writes the unified marker_intelligence_eval.json.

AC-1: macro_f1_delta_pp >= 5.0 (candidate beats baseline by 5 percentage points).
AC-2: condition_specific_substates non-empty AND each substate has a matching
      DE CSV with at least one row where fdr < 0.05.
AC-10: regression gate — invoked against baseline run outputs via
       scripts/ci/nc_regression_gate.sh (called externally; this script produces
       the eval JSON that gate consumes).

Usage:
    python scripts/eval/compare_marker_intelligence.py \\
        --baseline-run-dir <path>  \\
        --candidate-run-dir <path> \\
        --output-eval-json <path>  \\
        [--ac1-threshold 5.0]      \\
        [--ac2-fdr-threshold 0.05]

Exit codes:
    0 — eval JSON written; AC checks logged (non-blocking).
    1 — schema validation failure or missing required inputs.
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Inline JSON schema for marker_intelligence_eval.json (plan P1A.S3 imp #11).
# Validated with jsonschema if available; logged-only fallback otherwise.
# ---------------------------------------------------------------------------
_EVAL_SCHEMA = {
    "type": "object",
    "required": [
        "baseline_run_id",
        "candidate_run_id",
        "macro_f1_baseline",
        "macro_f1_candidate",
        "macro_f1_delta_pp",
        "per_cluster",
        "condition_specific_substates",
    ],
    "properties": {
        "baseline_run_id":            {"type": "string"},
        "candidate_run_id":           {"type": "string"},
        "macro_f1_baseline":          {"type": "number"},
        "macro_f1_candidate":         {"type": "number"},
        "macro_f1_delta_pp":          {"type": "number"},
        "per_cluster": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["cluster_id", "f1_baseline", "f1_candidate"],
                "properties": {
                    "cluster_id":   {"type": "string"},
                    "f1_baseline":  {"type": "number"},
                    "f1_candidate": {"type": "number"},
                },
            },
        },
        "condition_specific_substates": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["name", "parent_cell_type", "defining_markers", "fdr"],
                "properties": {
                    "name":             {"type": "string"},
                    "parent_cell_type": {"type": "string"},
                    "defining_markers": {"type": "array", "items": {"type": "string"}},
                    "fdr":              {"type": "number"},
                },
            },
        },
    },
    "additionalProperties": True,
}

_CANDIDATE_SCHEMA = {
    "type": "object",
    "required": [
        "candidate_run_id",
        "macro_f1_candidate",
        "per_cluster",
        "condition_specific_substates",
    ],
    "properties": {
        "candidate_run_id":           {"type": "string"},
        "macro_f1_candidate":         {"type": "number"},
        "per_cluster": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["cluster_id", "f1_candidate"],
                "properties": {
                    "cluster_id":   {"type": "string"},
                    "f1_candidate": {"type": "number"},
                },
            },
        },
        "condition_specific_substates": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["name", "parent_cell_type", "defining_markers", "fdr"],
            },
        },
    },
    "additionalProperties": True,
}

_BASELINE_SCHEMA = {
    "type": "object",
    "required": ["baseline_run_id", "macro_f1_baseline", "per_cluster"],
    "properties": {
        "baseline_run_id":   {"type": "string"},
        "macro_f1_baseline": {"type": "number"},
        "per_cluster": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["cluster_id", "f1_baseline"],
                "properties": {
                    "cluster_id":  {"type": "string"},
                    "f1_baseline": {"type": "number"},
                },
            },
        },
    },
    "additionalProperties": True,
}


def _validate_json(data: Any, schema: dict, label: str) -> None:
    """Validate data against schema; fail-fast with non-zero exit on mismatch or missing dep."""
    try:
        import jsonschema
    except ImportError:
        print(
            "[compare_marker_intelligence] ERROR: jsonschema is required but not installed. "
            "Install with: pip install jsonschema",
            file=sys.stderr,
        )
        sys.exit(1)
    try:
        jsonschema.validate(instance=data, schema=schema)
    except jsonschema.ValidationError as exc:
        print(
            f"[compare_marker_intelligence] SCHEMA MISMATCH in {label}: {exc.message}",
            file=sys.stderr,
        )
        sys.exit(1)


def _load_json(path: Path, label: str) -> dict:
    if not path.exists():
        print(f"[compare_marker_intelligence] ERROR: {label} not found: {path}", file=sys.stderr)
        sys.exit(1)
    with open(path) as fh:
        try:
            return json.load(fh)
        except json.JSONDecodeError as exc:
            print(
                f"[compare_marker_intelligence] ERROR: cannot parse {label} ({path}): {exc}",
                file=sys.stderr,
            )
            sys.exit(1)


def _run_id_from_dir(run_dir: Path) -> str:
    """Derive run_id from run directory name."""
    return run_dir.name


def _check_ac1(delta_pp: float, threshold: float) -> bool:
    passed = delta_pp >= threshold
    status = "PASS" if passed else "FAIL"
    print(
        f"[AC-1] macro_f1_delta_pp={delta_pp:.2f}pp vs threshold={threshold:.1f}pp → {status}",
        file=sys.stderr,
    )
    return passed


def _check_ac2(substates: list, candidate_run_dir: Path, fdr_threshold: float) -> bool:
    """AC-2: substates non-empty AND each has a DE CSV with fdr < threshold."""
    if not substates:
        print("[AC-2] FAIL: condition_specific_substates is empty.", file=sys.stderr)
        return False

    de_dir = candidate_run_dir / "differential_expression" / "substates"
    all_pass = True
    for substate in substates:
        name = substate.get("name", "")
        fdr_reported = substate.get("fdr", 1.0)
        csv_path = de_dir / f"{name}.csv"
        if not csv_path.exists():
            print(
                f"[AC-2] FAIL: DE CSV missing for substate '{name}': {csv_path}",
                file=sys.stderr,
            )
            all_pass = False
            continue
        # Check at least one row has fdr < threshold
        found_sig = False
        with open(csv_path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                try:
                    if float(row.get("fdr", 1.0)) < fdr_threshold:
                        found_sig = True
                        break
                except (ValueError, TypeError):
                    continue
        if found_sig:
            print(f"[AC-2] PASS: substate '{name}' has FDR < {fdr_threshold} rows.", file=sys.stderr)
        else:
            print(
                f"[AC-2] FAIL: substate '{name}' DE CSV has no rows with FDR < {fdr_threshold}.",
                file=sys.stderr,
            )
            all_pass = False

    return all_pass


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute marker intelligence delta (AC-1/AC-2) from baseline + candidate runs."
    )
    parser.add_argument(
        "--baseline-run-dir",
        required=True,
        type=Path,
        help="Path to baseline run directory (contains metrics/marker_intelligence_baseline.json).",
    )
    parser.add_argument(
        "--candidate-run-dir",
        required=True,
        type=Path,
        help="Path to candidate run directory (contains metrics/marker_intelligence_candidate.json).",
    )
    parser.add_argument(
        "--output-eval-json",
        required=True,
        type=Path,
        help="Output path for marker_intelligence_eval.json.",
    )
    parser.add_argument(
        "--ac1-threshold",
        type=float,
        default=5.0,
        help="Minimum macro-F1 delta in percentage points for AC-1 pass (default: 5.0).",
    )
    parser.add_argument(
        "--ac2-fdr-threshold",
        type=float,
        default=0.05,
        help="FDR threshold for AC-2 per-substate significance check (default: 0.05).",
    )
    args = parser.parse_args()

    baseline_dir = args.baseline_run_dir.resolve()
    candidate_dir = args.candidate_run_dir.resolve()

    # Load and validate baseline JSON
    baseline_json_path = baseline_dir / "metrics" / "marker_intelligence_baseline.json"
    baseline = _load_json(baseline_json_path, "baseline marker_intelligence_baseline.json")
    _validate_json(baseline, _BASELINE_SCHEMA, "baseline")

    # Load and validate candidate JSON
    candidate_json_path = candidate_dir / "metrics" / "marker_intelligence_candidate.json"
    candidate = _load_json(candidate_json_path, "candidate marker_intelligence_candidate.json")
    _validate_json(candidate, _CANDIDATE_SCHEMA, "candidate")

    # Derive run IDs — prefer embedded field, fall back to directory name
    baseline_run_id = baseline.get("baseline_run_id") or _run_id_from_dir(baseline_dir)
    candidate_run_id = candidate.get("candidate_run_id") or _run_id_from_dir(candidate_dir)

    macro_f1_baseline = float(baseline["macro_f1_baseline"])
    macro_f1_candidate = float(candidate["macro_f1_candidate"])
    macro_f1_delta_pp = macro_f1_candidate - macro_f1_baseline

    # Merge per_cluster entries — outer join on cluster_id
    baseline_clusters = {c["cluster_id"]: c for c in baseline.get("per_cluster", [])}
    candidate_clusters = {c["cluster_id"]: c for c in candidate.get("per_cluster", [])}
    all_cluster_ids = sorted(set(baseline_clusters) | set(candidate_clusters))
    per_cluster = []
    for cid in all_cluster_ids:
        per_cluster.append({
            "cluster_id":   cid,
            "f1_baseline":  float(baseline_clusters[cid]["f1_baseline"]) if cid in baseline_clusters else 0.0,
            "f1_candidate": float(candidate_clusters[cid]["f1_candidate"]) if cid in candidate_clusters else 0.0,
        })

    condition_specific_substates = candidate.get("condition_specific_substates", [])

    eval_doc = {
        "baseline_run_id":              baseline_run_id,
        "candidate_run_id":             candidate_run_id,
        "macro_f1_baseline":            macro_f1_baseline,
        "macro_f1_candidate":           macro_f1_candidate,
        "macro_f1_delta_pp":            macro_f1_delta_pp,
        "per_cluster":                  per_cluster,
        "condition_specific_substates": condition_specific_substates,
    }

    # Validate the assembled eval doc against the full output schema (fail-fast)
    _validate_json(eval_doc, _EVAL_SCHEMA, "output marker_intelligence_eval.json")

    # Write output
    args.output_eval_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_eval_json, "w") as fh:
        json.dump(eval_doc, fh, indent=2)
    print(f"[compare_marker_intelligence] Wrote: {args.output_eval_json}", file=sys.stderr)

    # AC checks (logged, non-blocking exit — caller reads the JSON for hard gates)
    _check_ac1(macro_f1_delta_pp, args.ac1_threshold)
    _check_ac2(condition_specific_substates, candidate_dir, args.ac2_fdr_threshold)


if __name__ == "__main__":
    main()
