from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"

NC_SPEC = importlib.util.spec_from_file_location(
    "validate_nc2024_architecture_contract",
    SCRIPTS_DIR / "validate_nc2024_architecture_contract.py",
)
nc_validator = importlib.util.module_from_spec(NC_SPEC)
assert NC_SPEC.loader is not None
sys.modules[NC_SPEC.name] = nc_validator
NC_SPEC.loader.exec_module(nc_validator)

SCRIPT = SCRIPTS_DIR / "validate_two_realdata_final.py"
SPEC = importlib.util.spec_from_file_location("validate_two_realdata_final", SCRIPT)
validator = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = validator
SPEC.loader.exec_module(validator)


def _sha(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _write_bundle_file(bundle: Path, rel_path: str, content: bytes) -> dict[str, object]:
    path = bundle / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return {"path": rel_path, "bytes": len(content), "sha256": _sha(content)}


def test_r_output_validation_rejects_corrupt_expected_plot(tmp_path: Path) -> None:
    run = tmp_path / "runs" / "2026-05-17T2004Z-13c2c88"
    r_dir = run / "r"
    r_dir.mkdir(parents=True)
    for name in nc_validator.EXPECTED_R_OUTPUT_FILES:
        payload = b"\x89PNG\r\n\x1a\nstub" if name.endswith(".png") else b"%PDF-1.4\nstub"
        (r_dir / name).write_bytes(payload)
    (r_dir / "umap_by_group.png").write_bytes(b"not-a-png")

    with pytest.raises(SystemExit, match="invalid PNG output header") as excinfo:
        validator._validate_r_outputs(run)
    assert "umap_by_group.png" in str(excinfo.value)


def test_r_output_validation_rejects_unexpected_extra_file(tmp_path: Path) -> None:
    run = tmp_path / "runs" / "2026-05-17T2004Z-13c2c88"
    r_dir = run / "r"
    r_dir.mkdir(parents=True)
    for name in nc_validator.EXPECTED_R_OUTPUT_FILES:
        payload = b"\x89PNG\r\n\x1a\nstub" if name.endswith(".png") else b"%PDF-1.4\nstub"
        (r_dir / name).write_bytes(payload)
    (r_dir / "extra_plot.png").write_bytes(b"\x89PNG\r\n\x1a\nstub")

    with pytest.raises(SystemExit, match="unexpected R outputs") as excinfo:
        validator._validate_r_outputs(run)
    assert "extra_plot.png" in str(excinfo.value)


def test_cell_bundle_validation_hashes_required_project_root_files(tmp_path: Path) -> None:
    run = tmp_path / "runs" / "2026-05-17T2004Z-13c2c88"
    bundle = run / "python" / "bundle"
    bundle.mkdir(parents=True)
    files = {
        "obs": _write_bundle_file(bundle, "obs.parquet", b"obs"),
        "marker_expr": _write_bundle_file(bundle, "marker_expr.parquet", b"markers"),
        "X_umap": _write_bundle_file(bundle, "obsm/X_umap.parquet", b"umap"),
        "X_pca": _write_bundle_file(bundle, "obsm/X_pca.parquet", b"pca"),
    }
    (bundle / "provenance.json").write_text(json.dumps({"ok": True}), encoding="utf-8")
    (bundle / "bundle_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "singlecell_r_bundle_v2.1",
                "source": {"n_obs": 55653, "n_vars": 25519},
                "bundle": {"n_cells_exported": 55653, "required_files": list(files)},
                "expression": {"claim_guard": nc_validator.CLAIM_GUARD},
                "files": files,
            }
        ),
        encoding="utf-8",
    )

    report = validator._validate_bundle(run, (55653, 25519), verify_sha=True)

    assert report["schema_version"] == "singlecell_r_bundle_v2.1"
    assert report["checked_file_count"] == 4


def test_cell_bundle_validation_rejects_gene_count_mismatch(tmp_path: Path) -> None:
    run = tmp_path / "runs" / "2026-05-17T2004Z-13c2c88"
    bundle = run / "python" / "bundle"
    bundle.mkdir(parents=True)
    files = {
        "obs": _write_bundle_file(bundle, "obs.parquet", b"obs"),
        "marker_expr": _write_bundle_file(bundle, "marker_expr.parquet", b"markers"),
        "X_umap": _write_bundle_file(bundle, "obsm/X_umap.parquet", b"umap"),
        "X_pca": _write_bundle_file(bundle, "obsm/X_pca.parquet", b"pca"),
    }
    (bundle / "provenance.json").write_text(json.dumps({"ok": True}), encoding="utf-8")
    (bundle / "bundle_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "singlecell_r_bundle_v2.1",
                "source": {"n_obs": 55653, "n_vars": 123},
                "bundle": {"n_cells_exported": 55653, "required_files": list(files)},
                "expression": {"claim_guard": nc_validator.CLAIM_GUARD},
                "files": files,
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(SystemExit, match="expected 25519 genes"):
        validator._validate_bundle(run, (55653, 25519), verify_sha=True)


def test_cell_bundle_validation_rejects_sha_mismatch(tmp_path: Path) -> None:
    run = tmp_path / "runs" / "2026-05-17T2004Z-13c2c88"
    bundle = run / "python" / "bundle"
    bundle.mkdir(parents=True)
    files = {
        "obs": _write_bundle_file(bundle, "obs.parquet", b"obs"),
        "marker_expr": _write_bundle_file(bundle, "marker_expr.parquet", b"markers"),
        "X_umap": _write_bundle_file(bundle, "obsm/X_umap.parquet", b"umap"),
        "X_pca": _write_bundle_file(bundle, "obsm/X_pca.parquet", b"pca"),
    }
    files["obs"]["sha256"] = "0" * 64
    (bundle / "provenance.json").write_text(json.dumps({"ok": True}), encoding="utf-8")
    (bundle / "bundle_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "singlecell_r_bundle_v2.1",
                "source": {"n_obs": 55653, "n_vars": 25519},
                "bundle": {"n_cells_exported": 55653, "required_files": list(files)},
                "expression": {"claim_guard": nc_validator.CLAIM_GUARD},
                "files": files,
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(SystemExit, match="sha256 mismatch for obs"):
        validator._validate_bundle(run, (55653, 25519), verify_sha=True)


def test_cell_pipeline_shape_prefers_final_doublet_removed_cells() -> None:
    shape = validator._cell_pipeline_shape({
        "metadata": {
            "cells_after_qc": 56879,
            "cells_after_doublet_removal": 55653,
            "genes_after_qc": 25519,
        }
    })
    assert shape == (55653, 25519)


def test_cell_quality_panel_count_accepts_current_source_key() -> None:
    quality = {"direct_evidence_panel_count": None, "evidence_panel_count_direct": 31}
    direct_panel_count = quality.get("direct_evidence_panel_count")
    if direct_panel_count is None:
        direct_panel_count = quality.get("evidence_panel_count_direct")
    assert direct_panel_count == 31
