#!/usr/bin/env python3
"""CI governance gate for three-factory validation contracts and docs."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SUITE_NAME = "Bioinformatics Research Pipeline"

EXPECTED_R_OUTPUT_FILES = {
    "cell_fraction.png",
    "marker_dot.png",
    "marker_feature.png",
    "qc_scatter.png",
    "qc_violin.png",
    "r.pdf",
    "umap_by_cluster.png",
    "umap_by_group.png",
}
REQUIRED_SINGLECELL_FILES = [
    "scripts/validate_two_realdata_final.py",
    "scripts/validate_nc2024_architecture_contract.py",
    "scripts/run_cleanroom_minimal_realdata.py",
    "scripts/validate_figure_parity_gate.py",
    "scripts/validate_ci_governance.py",
    "tests/test_two_realdata_final_validator.py",
    "tests/test_nc2024_architecture_contract.py",
    "tests/test_pipeline_hardening_gates.py",
    "contracts/bundle_schema.yaml",
    "docs/OWNER_BY_PRIMARY_OUTPUT_GOVERNANCE.md",
    ".github/workflows/factory-governance.yml",
]
REQUIRED_DOC_TOKENS = {
    "README.md": [
        SUITE_NAME,
        "clean-room minimal real-data gate",
        "figure parity gate",
        "CI governance gate",
        "scripts/run_cleanroom_minimal_realdata.py",
        "scripts/validate_figure_parity_gate.py",
        "scripts/validate_ci_governance.py",
        "owner-by-primary-output",
        "historical path",
        "retired",
        "R-heavy/spatial primary scientific truth",
        "presentation-only plotting surface",
        "biological conclusions",
        "docs/OWNER_BY_PRIMARY_OUTPUT_GOVERNANCE.md",
    ],
    "AI_AGENT_PROTOCOL.md": [
        SUITE_NAME,
        "clean-room minimal real-data gate",
        "figure parity gate",
        "CI governance gate",
        "owner-by-primary-output",
        "retired",
        "R-heavy/spatial primary scientific truth",
        "presentation-only plotting surface",
        "biological conclusions",
    ],
    "AGENTS.md": [
        SUITE_NAME,
        "owner-by-primary-output",
        "retired",
        "R-heavy/spatial primary scientific truth",
        "presentation-only plotting surface",
        "biological conclusions",
    ],
    "docs/OWNER_BY_PRIMARY_OUTPUT_GOVERNANCE.md": [
        SUITE_NAME,
        "owner-by-primary-output",
        "retired",
        "R-heavy/spatial primary scientific truth",
        "presentation-only plotting surface",
        "biological conclusions",
    ],
}
DOWNSTREAM_DOC_TOKENS = {
    "r_multiomics_factory": [
        SUITE_NAME,
        "figure parity gate",
        "CI governance gate",
        "owner-by-primary-output",
        "retired",
        "R-heavy/spatial primary scientific truth",
        "presentation-only plotting surface",
        "biological conclusions",
    ],
    "plotting_factory": [
        SUITE_NAME,
        "figure parity gate",
        "CI governance gate",
        "owner-by-primary-output",
        "retired",
        "presentation-only plotting surface",
        "biological conclusions",
    ],
}


def _finding(severity: str, code: str, path: Path, message: str) -> dict[str, str]:
    return {"severity": severity, "code": code, "path": str(path), "message": message}


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def _check_required_files(repo_root: Path) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []
    for rel in REQUIRED_SINGLECELL_FILES:
        path = repo_root / rel
        if not path.exists():
            findings.append(_finding("error", "missing_required_file", path, f"Missing {rel}"))
    return findings


def _check_expected_outputs(repo_root: Path) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []
    canonical = repo_root / "scripts" / "validate_nc2024_architecture_contract.py"
    canonical_content = _text(canonical)
    for filename in sorted(EXPECTED_R_OUTPUT_FILES):
        if filename not in canonical_content:
            findings.append(_finding(
                "error",
                "expected_r_output_not_declared",
                canonical,
                f"Expected R output {filename!r} is not declared in {canonical.name}.",
            ))

    # Downstream validators should reuse the canonical expected-output set rather
    # than duplicating filenames and drifting independently.
    for rel in [
        "scripts/validate_two_realdata_final.py",
        "scripts/validate_figure_parity_gate.py",
    ]:
        path = repo_root / rel
        content = _text(path)
        if "nc_validator.EXPECTED_R_OUTPUT_FILES" not in content:
            findings.append(_finding(
                "error",
                "expected_r_output_source_not_reused",
                path,
                (
                    "Gate must reuse nc_validator.EXPECTED_R_OUTPUT_FILES "
                    "as the canonical R output set."
                ),
            ))
    return findings


def _check_docs(repo_root: Path, r_repo: Path, plotting_repo: Path) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []
    for rel, tokens in REQUIRED_DOC_TOKENS.items():
        path = repo_root / rel
        content = _text(path)
        for token in tokens:
            if token not in content:
                findings.append(_finding(
                    "error",
                    "doc_token_missing",
                    path,
                    f"Missing doc token: {token}",
                ))
    downstream = {"r_multiomics_factory": r_repo, "plotting_factory": plotting_repo}
    for name, tokens in DOWNSTREAM_DOC_TOKENS.items():
        root = downstream[name]
        combined = (
            _text(root / "README.md")
            + "\n"
            + _text(root / "AI_AGENT_PROTOCOL.md")
            + "\n"
            + _text(root / "AGENTS.md")
        )
        for token in tokens:
            if token not in combined:
                findings.append(_finding(
                    "error",
                    "downstream_doc_token_missing",
                    root,
                    f"{name} missing token: {token}",
                ))
    return findings


def _check_workflow(repo_root: Path) -> list[dict[str, str]]:
    path = repo_root / ".github" / "workflows" / "factory-governance.yml"
    content = _text(path)
    findings: list[dict[str, str]] = []
    for token in [
        "validate_ci_governance.py",
        "validate_figure_parity_gate.py",
        "run_cleanroom_minimal_realdata.py",
        "pytest --no-cov",
    ]:
        if token not in content:
            findings.append(_finding(
                "error",
                "workflow_command_missing",
                path,
                f"Workflow missing command token: {token}",
            ))
    return findings


def validate(repo_root: Path, r_repo: Path, plotting_repo: Path) -> dict[str, Any]:
    repo_root = repo_root.resolve()
    findings: list[dict[str, str]] = []
    findings.extend(_check_required_files(repo_root))
    findings.extend(_check_expected_outputs(repo_root))
    findings.extend(_check_docs(repo_root, r_repo.resolve(), plotting_repo.resolve()))
    findings.extend(_check_workflow(repo_root))
    severity_counts: dict[str, int] = {}
    for finding in findings:
        severity_counts[finding["severity"]] = severity_counts.get(finding["severity"], 0) + 1
    verdict = "fail" if severity_counts.get("error", 0) else "pass"
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": verdict,
        "repo_root": str(repo_root),
        "r_multiomics_factory": str(r_repo.resolve()),
        "plotting_factory": str(plotting_repo.resolve()),
        "suite_name": SUITE_NAME,
        "severity_counts": severity_counts,
        "expected_r_output_files": sorted(EXPECTED_R_OUTPUT_FILES),
        "findings": findings,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    suite_root = root.parent
    parser.add_argument("--repo-root", type=Path, default=root)
    parser.add_argument(
        "--r-repo",
        type=Path,
        default=suite_root / "r_multiomics_factory",
    )
    parser.add_argument(
        "--plotting-repo",
        type=Path,
        default=suite_root / "plotting_factory",
    )
    parser.add_argument("--output", type=Path, default=None)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = validate(args.repo_root, args.r_repo, args.plotting_repo)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["verdict"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
