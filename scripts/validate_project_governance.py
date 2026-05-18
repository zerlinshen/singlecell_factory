#!/usr/bin/env python3
"""Read-only validator for remote factory/project governance.

This is a control-plane validator. It inspects project roots and run roots but
never writes into project directories. Reports may be written to an explicit
factory-side --output-dir.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

try:
    import yaml as _yaml
    _YAML_AVAILABLE = True
except ImportError:
    _YAML_AVAILABLE = False

_IMAGE_SUFFIXES = {".png", ".pdf", ".svg", ".tif", ".tiff"}

C2_PARITY_TAXONOMY = ("exact", "approximate", "proxy", "unsupported", "resource_gap")
C2_GAP_TARGETS = (
    "singlecell_factory",
    "r_multiomics_factory",
    "plotting_factory",
    "paper_specific_script",
    "parameter_change_only",
)
_C2_TOP_LEVEL_FIELDS = (
    "upstream_repository",
    "raw_data_reproduction",
    "data_object_reproduction",
    "figure_reproduction",
    "module_gap_decisions",
    "context_optimization_decisions",
)


def _rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _finding(severity: str, code: str, path: Path, message: str) -> dict:
    return {"severity": severity, "code": code, "path": str(path), "message": message}


def _glob_existing(root: Path, patterns: Iterable[str]) -> list[str]:
    hits: list[Path] = []
    for pattern in patterns:
        hits.extend(p for p in root.glob(pattern) if p.is_file())
    return sorted({_rel(p, root) for p in hits})


def _count_artifacts(run_dir: Path) -> dict:
    figure_count = 0
    evidence_json_count = 0
    for path in run_dir.rglob("*"):
        if not path.is_file():
            continue
        lower = path.name.lower()
        if path.suffix.lower() in _IMAGE_SUFFIXES:
            figure_count += 1
        if path.suffix.lower() == ".json" and "evidence" in lower:
            evidence_json_count += 1
    return {"figures": figure_count, "evidence_json": evidence_json_count}


def _validate_c2_reproduction_fields(project_yaml_data: dict, project_yaml_path: Path) -> list[dict]:
    """Advisory warnings (never errors) for missing or malformed C2 reproduction fields."""
    warnings: list[dict] = []

    for field in _C2_TOP_LEVEL_FIELDS:
        if field not in project_yaml_data:
            warnings.append(_finding(
                "warning",
                f"c2_missing_{field}",
                project_yaml_path,
                f"C2 reproduction field '{field}' not present; advisory only — populate via PAPER_REPRODUCTION_SOP.md.",
            ))

    upstream = project_yaml_data.get("upstream_repository")
    if isinstance(upstream, dict):
        for subfield in ("url", "commit_or_tag", "doi", "license", "forked_at", "fork_path"):
            if subfield not in upstream:
                warnings.append(_finding(
                    "warning",
                    f"c2_upstream_repository_missing_{subfield}",
                    project_yaml_path,
                    f"upstream_repository.{subfield} not present; record before first reproduction run.",
                ))

    for map_field in ("data_object_reproduction", "figure_reproduction"):
        field_data = project_yaml_data.get(map_field)
        if isinstance(field_data, dict):
            for key, value in field_data.items():
                if value not in C2_PARITY_TAXONOMY:
                    warnings.append(_finding(
                        "warning",
                        f"c2_{map_field}_invalid_parity_class",
                        project_yaml_path,
                        f"{map_field}['{key}'] = '{value}' is not a valid parity class; "
                        f"must be one of: {', '.join(C2_PARITY_TAXONOMY)}.",
                    ))

    gap_decisions = project_yaml_data.get("module_gap_decisions")
    if isinstance(gap_decisions, list):
        for i, item in enumerate(gap_decisions):
            if not isinstance(item, dict):
                continue
            target = item.get("target")
            if target is not None and target not in C2_GAP_TARGETS:
                warnings.append(_finding(
                    "warning",
                    "c2_module_gap_decisions_invalid_target",
                    project_yaml_path,
                    f"module_gap_decisions[{i}].target = '{target}' is not a valid routing target; "
                    f"must be one of: {', '.join(C2_GAP_TARGETS)}.",
                ))
            gap_class = item.get("gap_class")
            if gap_class is not None and gap_class not in C2_PARITY_TAXONOMY:
                warnings.append(_finding(
                    "warning",
                    "c2_module_gap_decisions_invalid_gap_class",
                    project_yaml_path,
                    f"module_gap_decisions[{i}].gap_class = '{gap_class}' is not a valid parity class; "
                    f"must be one of: {', '.join(C2_PARITY_TAXONOMY)}.",
                ))

    context_decisions = project_yaml_data.get("context_optimization_decisions")
    if isinstance(context_decisions, list):
        for i, item in enumerate(context_decisions):
            if not isinstance(item, dict):
                continue
            parity_class = item.get("parity_class")
            if parity_class is not None and parity_class not in C2_PARITY_TAXONOMY:
                warnings.append(_finding(
                    "warning",
                    "c2_context_optimization_decisions_invalid_parity_class",
                    project_yaml_path,
                    f"context_optimization_decisions[{i}].parity_class = '{parity_class}' is not a valid parity class; "
                    f"must be one of: {', '.join(C2_PARITY_TAXONOMY)}.",
                ))

    return warnings


def _validate_run(project_root: Path, run_dir: Path) -> dict:
    findings: list[dict] = []
    run_id = run_dir.name

    root_manifest = run_dir / "manifest.json"
    if not root_manifest.exists():
        findings.append(_finding(
            "warning",
            "missing_root_manifest",
            root_manifest,
            "Missing cross-factory run envelope; acceptable as legacy/partial state when substitute provenance exists.",
        ))

    for dirname in ["python", "r", "logs"]:
        d = run_dir / dirname
        if not d.is_dir():
            findings.append(_finding(
                "warning",
                f"missing_{dirname}_dir",
                d,
                "Missing canonical run subdirectory; warning for legacy/partial runs.",
            ))

    python_dir = run_dir / "python"
    r_dir = run_dir / "r"
    producer_native = {
        "python_run_manifests": _glob_existing(python_dir, ["run_manifest.json", "**/run_manifest.json"]) if python_dir.exists() else [],
        "python_module_status": _glob_existing(python_dir, ["module_status.csv", "**/module_status.csv"]) if python_dir.exists() else [],
        "bundle_provenance": _glob_existing(python_dir, ["bundle/provenance.json", "**/bundle/provenance.json"]) if python_dir.exists() else [],
        "r_manifests": _glob_existing(r_dir, ["*manifest*.json", "**/*manifest*.json", "*provenance*.json", "**/*provenance*.json"]) if r_dir.exists() else [],
    }

    if not any(producer_native.values()):
        findings.append(_finding(
            "info",
            "no_producer_native_manifest_found",
            run_dir,
            "No producer-native manifest/module-status files found at expected paths; inspect legacy evidence manually if this run is cited.",
        ))

    artifact_counts = _count_artifacts(run_dir)
    if artifact_counts["figures"] == 0 and artifact_counts["evidence_json"] == 0:
        findings.append(_finding(
            "info",
            "no_review_artifacts_found",
            run_dir,
            "No figure or evidence JSON artifacts found under this run root.",
        ))

    return {
        "run_id": run_id,
        "path": str(run_dir),
        "root_manifest": str(root_manifest) if root_manifest.exists() else None,
        "producer_native_provenance": producer_native,
        "artifact_counts": artifact_counts,
        "findings": findings,
    }


def _load_yaml_safe(path: Path) -> dict | None:
    """Return parsed YAML dict or None if yaml unavailable or parse fails."""
    if not _YAML_AVAILABLE:
        return None
    try:
        data = _yaml.safe_load(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else None
    except Exception:
        return None


def validate_project(project_root: Path, run_id: str | None = None) -> dict:
    project_root = project_root.resolve()
    findings: list[dict] = []
    runs: list[dict] = []

    if not project_root.exists():
        findings.append(_finding("error", "project_root_missing", project_root, "Project root does not exist."))
        return _build_report(project_root, run_id, findings, runs)

    project_yaml = project_root / "project.yaml"
    if not project_yaml.exists():
        findings.append(_finding(
            "warning",
            "missing_project_yaml",
            project_yaml,
            "Missing project.yaml; warning for legacy projects, required for new-governance projects.",
        ))
    else:
        project_yaml_data = _load_yaml_safe(project_yaml)
        if project_yaml_data is not None:
            findings.extend(_validate_c2_reproduction_fields(project_yaml_data, project_yaml))

    runs_dir = project_root / "runs"
    if not runs_dir.is_dir():
        findings.append(_finding("error", "runs_dir_missing", runs_dir, "Project runs/ directory is missing."))
        return _build_report(project_root, run_id, findings, runs)

    if run_id:
        candidate_runs = [runs_dir / run_id]
        if not candidate_runs[0].is_dir():
            findings.append(_finding("error", "run_id_missing", candidate_runs[0], "Requested run-id does not exist."))
            return _build_report(project_root, run_id, findings, runs)
    else:
        candidate_runs = sorted(p for p in runs_dir.iterdir() if p.is_dir() and p.name != "logs")
        if not candidate_runs:
            findings.append(_finding("warning", "no_run_dirs", runs_dir, "No run directories found."))

    for run_dir in candidate_runs:
        runs.append(_validate_run(project_root, run_dir))

    return _build_report(project_root, run_id, findings, runs)


def _build_report(project_root: Path, run_id: str | None, findings: list[dict], runs: list[dict]) -> dict:
    all_findings = list(findings)
    for run in runs:
        all_findings.extend(run["findings"])
    counts = Counter(f["severity"] for f in all_findings)
    return {
        "schema_version": "project-governance-validation/v1",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "control_plane_only": True,
        "project_root": str(project_root),
        "project_id": project_root.name,
        "requested_run_id": run_id,
        "overall_status": "fail" if counts.get("error", 0) else "pass",
        "severity_counts": dict(counts),
        "project_findings": findings,
        "runs": runs,
    }


def render_markdown(report: dict) -> str:
    lines = [
        f"# Project Governance Validation: {report['project_id']}",
        "",
        "Control-plane validation report only; not a scientific output.",
        "",
        f"- Generated: `{report['generated_at']}`",
        f"- Project root: `{report['project_root']}`",
        f"- Requested run: `{report.get('requested_run_id') or 'all discovered runs'}`",
        f"- Overall status: `{report['overall_status']}`",
        f"- Severity counts: `{report['severity_counts']}`",
        "",
        "## Runs",
    ]
    for run in report["runs"]:
        lines.extend([
            "",
            f"### `{run['run_id']}`",
            f"- Root manifest: `{run['root_manifest'] or 'missing_expected_manifest'}`",
            f"- Figures: `{run['artifact_counts']['figures']}`",
            f"- Evidence JSON: `{run['artifact_counts']['evidence_json']}`",
            f"- Python run manifests: `{len(run['producer_native_provenance']['python_run_manifests'])}`",
            f"- Module status files: `{len(run['producer_native_provenance']['python_module_status'])}`",
            f"- Bundle provenance files: `{len(run['producer_native_provenance']['bundle_provenance'])}`",
            f"- R manifests/provenance: `{len(run['producer_native_provenance']['r_manifests'])}`",
        ])
        if run["findings"]:
            lines.append("- Findings:")
            for finding in run["findings"]:
                lines.append(f"  - `{finding['severity']}` `{finding['code']}`: {finding['message']}")
    if report["project_findings"]:
        lines.extend(["", "## Project Findings"])
        for finding in report["project_findings"]:
            lines.append(f"- `{finding['severity']}` `{finding['code']}`: {finding['message']}")
    lines.append("")
    return "\n".join(lines)


def write_report(report: dict, output_dir: Path) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = report.get("requested_run_id") or "all-runs"
    stem = f"{report['project_id']}_{suffix}".replace("/", "_")
    json_path = output_dir / f"{stem}.governance_validation.json"
    md_path = output_dir / f"{stem}.governance_validation.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(render_markdown(report), encoding="utf-8")
    return json_path, md_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--json", action="store_true", help="Print JSON report to stdout instead of Markdown.")
    args = parser.parse_args()

    report = validate_project(args.project_root, run_id=args.run_id)
    if args.output_dir:
        json_path, md_path = write_report(report, args.output_dir)
        print(f"wrote {json_path}")
        print(f"wrote {md_path}")
    else:
        print(json.dumps(report, indent=2) if args.json else render_markdown(report))
    return 1 if report["overall_status"] == "fail" else 0


if __name__ == "__main__":
    raise SystemExit(main())
