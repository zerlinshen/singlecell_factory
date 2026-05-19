from __future__ import annotations

import importlib.util
import json
import sys
import zlib
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS_DIR / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


cleanroom = _load("run_cleanroom_minimal_realdata")
figure_gate = _load("validate_figure_parity_gate")
ci_gate = _load("validate_ci_governance")


def _png_chunk(kind: bytes, payload: bytes) -> bytes:
    import struct

    return (
        struct.pack(">I", len(payload))
        + kind
        + payload
        + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF)
    )


def _write_png(path: Path, rgb: tuple[int, int, int]) -> None:
    import struct

    header = b"\x89PNG\r\n\x1a\n"
    ihdr = struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0)
    raw = b"\x00" + bytes(rgb)
    path.write_bytes(
        header
        + _png_chunk(b"IHDR", ihdr)
        + _png_chunk(b"IDAT", zlib.compress(raw))
        + _png_chunk(b"IEND", b"")
    )


def test_figure_parity_passes_exact_png(tmp_path: Path) -> None:
    produced = tmp_path / "produced.png"
    reference = tmp_path / "reference.png"
    _write_png(produced, (0, 0, 0))
    reference.write_bytes(produced.read_bytes())
    report = figure_gate.validate_manifest(
        {"figures": [{"id": "umap", "produced": produced.name, "reference": reference.name}]},
        base=tmp_path,
    )
    assert report["verdict"] == "pass"
    assert report["figures"][0]["metric"] == "exact_sha256"


def test_figure_parity_flags_png_drift(tmp_path: Path) -> None:
    produced = tmp_path / "produced.png"
    reference = tmp_path / "reference.png"
    _write_png(produced, (255, 255, 255))
    _write_png(reference, (0, 0, 0))
    report = figure_gate.validate_manifest(
        {
            "figures": [
                {
                    "id": "marker_feature",
                    "produced": produced.name,
                    "reference": reference.name,
                    "max_rms_normalized": 0.0,
                    "max_byte_mismatch_ratio": 0.0,
                }
            ]
        },
        base=tmp_path,
    )
    assert report["verdict"] == "fail"
    assert report["figures"][0]["status"] == "fail"


def test_figure_parity_missing_optional_reference_is_conditional(tmp_path: Path) -> None:
    produced = tmp_path / "produced.png"
    _write_png(produced, (0, 0, 0))
    report = figure_gate.validate_manifest(
        {"figures": [{"id": "qc", "produced": produced.name, "reference": "missing.png"}]},
        base=tmp_path,
    )
    assert report["verdict"] == "conditional"
    assert report["figures"][0]["reason"] == "reference_missing"


def test_figure_parity_checks_pdf_header(tmp_path: Path) -> None:
    produced = tmp_path / "r.pdf"
    produced.write_bytes(b"%PDF-1.4\n%stub\n")
    report = figure_gate.validate_manifest(
        {"figures": [{"id": "r_pdf", "produced": produced.name}]},
        base=tmp_path,
    )
    assert report["verdict"] == "conditional"
    assert report["figures"][0]["kind"] == "pdf"


def test_cleanroom_required_artifacts_reports_missing(tmp_path: Path) -> None:
    nc_project = tmp_path / "nc"
    cell_project = tmp_path / "cell"
    checks = cleanroom._required_artifacts(
        nc_project_root=nc_project,
        nc_run_id="run-a",
        cell_project_root=cell_project,
        cell_pipeline_run_id="run-b",
        cell_evidence_run_id="run-c",
    )
    assert any(not check["ok"] for check in checks)
    labels = {check["label"] for check in checks if not check["ok"]}
    assert "nc_project_root" in labels
    assert "cell_project_root" in labels


def _write_minimal_governance_repo(root: Path, r_repo: Path, plotting_repo: Path) -> None:
    for rel in ci_gate.REQUIRED_SINGLECELL_FILES:
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".py":
            content = "\n".join(sorted(ci_gate.EXPECTED_R_OUTPUT_FILES)) + "\n"
            if path.name in {"validate_two_realdata_final.py", "validate_figure_parity_gate.py"}:
                content += "nc_validator.EXPECTED_R_OUTPUT_FILES\n"
            path.write_text(content, encoding="utf-8")
        else:
            path.write_text("placeholder\n", encoding="utf-8")
    owner_tokens = (
        "Bioinformatics Research Pipeline\n"
        "owner-by-primary-output\n"
        "historical path\n"
        "retired\n"
        "R-heavy/spatial primary scientific truth\n"
        "presentation-only plotting surface\n"
        "biological conclusions\n"
        "docs/OWNER_BY_PRIMARY_OUTPUT_GOVERNANCE.md\n"
    )
    (root / "README.md").write_text(
        "clean-room minimal real-data gate\nfigure parity gate\nCI governance gate\n"
        "scripts/run_cleanroom_minimal_realdata.py\nscripts/validate_figure_parity_gate.py\n"
        "scripts/validate_ci_governance.py\n"
        + owner_tokens,
        encoding="utf-8",
    )
    (root / "AI_AGENT_PROTOCOL.md").write_text(
        "clean-room minimal real-data gate\nfigure parity gate\nCI governance gate\n"
        + owner_tokens,
        encoding="utf-8",
    )
    (root / "AGENTS.md").write_text(owner_tokens, encoding="utf-8")
    (root / "docs" / "OWNER_BY_PRIMARY_OUTPUT_GOVERNANCE.md").parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    (root / "docs" / "OWNER_BY_PRIMARY_OUTPUT_GOVERNANCE.md").write_text(
        owner_tokens,
        encoding="utf-8",
    )
    workflow = root / ".github" / "workflows" / "factory-governance.yml"
    workflow.write_text(
        "validate_ci_governance.py\nvalidate_figure_parity_gate.py\n"
        "run_cleanroom_minimal_realdata.py\npytest --no-cov\n",
        encoding="utf-8",
    )
    for repo in (r_repo, plotting_repo):
        repo.mkdir(parents=True, exist_ok=True)
        if repo == r_repo:
            doc = (
                "Bioinformatics Research Pipeline\n"
                "figure parity gate\nCI governance gate\nowner-by-primary-output\n"
                "retired\n"
                "R-heavy/spatial primary scientific truth\n"
                "presentation-only plotting surface\n"
                "biological conclusions\n"
            )
        else:
            doc = (
                "Bioinformatics Research Pipeline\n"
                "figure parity gate\nCI governance gate\nowner-by-primary-output\n"
                "retired\n"
                "presentation-only plotting surface\n"
                "biological conclusions\n"
            )
        (repo / "README.md").write_text(doc, encoding="utf-8")
        (repo / "AI_AGENT_PROTOCOL.md").write_text("", encoding="utf-8")
        (repo / "AGENTS.md").write_text("", encoding="utf-8")


def test_ci_governance_gate_passes_minimal_repo(tmp_path: Path) -> None:
    root = tmp_path / "singlecell"
    r_repo = tmp_path / "r"
    plotting_repo = tmp_path / "plotting"
    _write_minimal_governance_repo(root, r_repo, plotting_repo)
    report = ci_gate.validate(root, r_repo, plotting_repo)
    assert report["verdict"] == "pass"
    assert report["severity_counts"] == {}


def test_ci_governance_gate_fails_missing_doc_token(tmp_path: Path) -> None:
    root = tmp_path / "singlecell"
    r_repo = tmp_path / "r"
    plotting_repo = tmp_path / "plotting"
    _write_minimal_governance_repo(root, r_repo, plotting_repo)
    (root / "README.md").write_text("figure parity gate\nCI governance gate\n", encoding="utf-8")
    report = ci_gate.validate(root, r_repo, plotting_repo)
    assert report["verdict"] == "fail"
    assert any(f["code"] == "doc_token_missing" for f in report["findings"])


def test_ci_governance_gate_fails_missing_plotting_owner_boundary(tmp_path: Path) -> None:
    root = tmp_path / "singlecell"
    r_repo = tmp_path / "r"
    plotting_repo = tmp_path / "plotting"
    _write_minimal_governance_repo(root, r_repo, plotting_repo)
    (plotting_repo / "README.md").write_text(
        "figure parity gate\nCI governance gate\nowner-by-primary-output\n",
        encoding="utf-8",
    )
    report = ci_gate.validate(root, r_repo, plotting_repo)
    assert report["verdict"] == "fail"
    assert any(f["code"] == "downstream_doc_token_missing" for f in report["findings"])


def test_ci_governance_gate_fails_missing_suite_name(tmp_path: Path) -> None:
    root = tmp_path / "singlecell"
    r_repo = tmp_path / "r"
    plotting_repo = tmp_path / "plotting"
    _write_minimal_governance_repo(root, r_repo, plotting_repo)
    (root / "docs" / "OWNER_BY_PRIMARY_OUTPUT_GOVERNANCE.md").write_text(
        "owner-by-primary-output\n"
        "R-heavy/spatial primary scientific truth\n"
        "presentation-only plotting surface\n"
        "biological conclusions\n",
        encoding="utf-8",
    )
    report = ci_gate.validate(root, r_repo, plotting_repo)
    assert report["verdict"] == "fail"
    assert any("Bioinformatics Research Pipeline" in f["message"] for f in report["findings"])
