from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_nc2024_architecture_contract.py"
SPEC = importlib.util.spec_from_file_location("validate_nc2024_architecture_contract", SCRIPT)
validator = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = validator
SPEC.loader.exec_module(validator)


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


def _plot_bytes(name: str) -> bytes:
    if name.endswith(".png"):
        return b"\x89PNG\r\n\x1a\nstub"
    if name.endswith(".pdf"):
        return b"%PDF-1.4\nstub"
    raise AssertionError(f"unexpected plot suffix: {name}")


def _make_project_run(tmp_path: Path) -> Path:
    project = tmp_path / "nc-reproduction"
    run_id = "2026-05-18T0900Z-13c2c88"
    run = project / "runs" / run_id
    py_run = run / "python" / "NC2024_WAVE3_SUBSAMPLE_P15_T1_20260518_170048"
    bundle = run / "python" / "bundle"
    r_dir = run / "r"
    (run / "ops" / "run_ledger").mkdir(parents=True)
    r_dir.mkdir(parents=True)
    bundle.mkdir(parents=True)
    py_run.mkdir(parents=True)

    _write_json(run / "manifest.json", {"run_id": run_id, "factory_python": {"sha": "13c2c88"}, "factory_r": {"sha": "4ede0ce"}})
    _write_json(run / "ops" / "run_ledger" / "ledger.json", {"ok": True})
    _write_json(
        py_run / "run_manifest.json",
        {
            "metadata": {
                "cells_after_qc": 5281,
                "genes_after_qc": 19504,
                "context_validation_mismatches": 0,
                "sample_id_filter_applied": ["P15_T1"],
                "clustering_backend": "gpu",
                "de_backend": "cpu",
                "batch_correction_status": "skipped_single_batch",
            }
        },
    )
    (py_run / "module_status.csv").write_text(
        "module,status,message\n"
        "cellranger,ok,completed\n"
        "marker_db_loader,ok,completed\n"
        "qc,ok,completed\n"
        "doublet_detection,ok,completed\n"
        "clustering,ok,completed\n"
        "annotation,ok,completed\n"
        "batch_correction,skipped,one batch\n"
        "context_aware_annotation,ok,completed\n"
        "differential_expression,ok,completed\n",
        encoding="utf-8",
    )
    (py_run / "final_adata.h5ad").write_bytes(b"not-opened-by-validator")

    for name, content in {
        "obs.parquet": b"obs",
        "marker_expr.parquet": b"marker",
        "X_umap.parquet": b"umap",
        "X_pca.parquet": b"pca",
    }.items():
        (bundle / name).write_bytes(content)
    _write_json(bundle / "provenance.json", {"bundle_sha256": "stub"})
    _write_json(
        bundle / "bundle_manifest.json",
        {
            "schema_version": "singlecell_r_bundle_v2.1",
            "source": {"n_obs": 5281, "n_vars": 19504},
            "bundle": {
                "n_cells_exported": 5281,
                "required_files": ["obs", "marker_expr", "X_umap", "X_pca"],
                "markers_present": ["CD3E", "LYZ"],
                "marker_format": "parquet",
            },
            "expression": {"claim_guard": validator.CLAIM_GUARD},
            "files": {
                "obs": {"path": "obs.parquet", "bytes": 3},
                "marker_expr": {"path": "marker_expr.parquet", "bytes": 6},
                "X_umap": {"path": "X_umap.parquet", "bytes": 4},
                "X_pca": {"path": "X_pca.parquet", "bytes": 3},
            },
        },
    )
    for name in validator.EXPECTED_R_OUTPUT_FILES:
        (r_dir / name).write_bytes(_plot_bytes(name))
    return project


def test_validate_project_run_accepts_project_root_layout(tmp_path: Path) -> None:
    project = _make_project_run(tmp_path)
    report = validator.validate_project_run(project, "2026-05-18T0900Z-13c2c88")
    assert report["python"]["shape"] == {"n_obs": 5281, "n_vars": 19504}
    assert report["bundle"]["schema_version"] == "singlecell_r_bundle_v2.1"
    assert report["r_outputs"]["file_count"] == len(validator.EXPECTED_R_OUTPUT_FILES)


def test_validator_defaults_no_longer_point_at_legacy_factory_results() -> None:
    assert "singlecell_factory/results" not in str(validator.DEFAULT_PROJECT_ROOT)
    assert validator.DEFAULT_PROJECT_ROOT.name == "nc-reproduction"


def test_validate_project_run_rejects_missing_expected_r_plot(tmp_path: Path) -> None:
    project = _make_project_run(tmp_path)
    missing = project / "runs" / "2026-05-18T0900Z-13c2c88" / "r" / "umap_by_group.png"
    missing.unlink()

    with pytest.raises(SystemExit, match="missing expected R outputs") as excinfo:
        validator.validate_project_run(project, "2026-05-18T0900Z-13c2c88")
    assert "umap_by_group.png" in str(excinfo.value)


def test_validate_project_run_rejects_corrupt_r_plot_header(tmp_path: Path) -> None:
    project = _make_project_run(tmp_path)
    corrupt = project / "runs" / "2026-05-18T0900Z-13c2c88" / "r" / "marker_dot.png"
    corrupt.write_bytes(b"not-a-png")

    with pytest.raises(SystemExit, match="invalid PNG output header") as excinfo:
        validator.validate_project_run(project, "2026-05-18T0900Z-13c2c88")
    assert "marker_dot.png" in str(excinfo.value)


def test_validate_project_run_rejects_unexpected_extra_r_output(tmp_path: Path) -> None:
    project = _make_project_run(tmp_path)
    extra = project / "runs" / "2026-05-18T0900Z-13c2c88" / "r" / "extra_plot.png"
    extra.write_bytes(_plot_bytes("extra_plot.png"))

    with pytest.raises(SystemExit, match="unexpected R outputs") as excinfo:
        validator.validate_project_run(project, "2026-05-18T0900Z-13c2c88")
    assert "extra_plot.png" in str(excinfo.value)
