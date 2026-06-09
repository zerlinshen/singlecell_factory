"""Hermetic tests for scripts/validate_realdata_parity_gate.py (#33).

Build a tiny synthetic SOT run directory + manifest in tmp_path that mirrors the
canonical manifest shape (inputs.exported_cells / markers_requested,
outputs.{bundle_dir,r_dir,png_count,pdf_count}, verification with a recorded
PNG min-dimension), then:
  * assert a complete, real evidence chain PASSES;
  * mutate the chain (missing PNG / wrong bundle cell count / manifest path
    escape) and assert each FAILS.

No real project path is touched; everything lives under tmp_path.
"""

from __future__ import annotations

import importlib.util
import json
import struct
import zlib
from pathlib import Path

import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "validate_realdata_parity_gate.py"
)
_spec = importlib.util.spec_from_file_location("validate_realdata_parity_gate", MODULE_PATH)
gate = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(gate)

EXPORTED_CELLS = 5000
MARKERS = ["CD3E", "CD4", "CD8A", "NKG7", "MS4A1", "EPCAM", "PTPRC"]
PNG_COUNT = 2
PDF_COUNT = 1
# The synthetic PNGs are 3300x1300 so they clear the recorded >= 3200x1280
# threshold that the gate parses out of the verification strings.
PNG_W, PNG_H = 3300, 1300


def _write_real_png(path: Path, width: int = PNG_W, height: int = PNG_H) -> None:
    """Write a minimal but valid, non-blank PNG of the given dimensions.

    Hand-rolled (no Pillow dependency in the writer) so the test is hermetic;
    the gate itself opportunistically uses Pillow when present to read dims.
    """

    def chunk(tag: bytes, data: bytes) -> bytes:
        return (
            struct.pack(">I", len(data))
            + tag
            + data
            + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
        )

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)  # 8-bit RGB
    # One filter byte per row + RGB pixels; paint a non-uniform band so the
    # image is non-degenerate (not a single flat colour).
    raw = bytearray()
    for y in range(height):
        raw.append(0)  # filter type 0 (None)
        for x in range(width):
            raw += bytes(((x + y) % 256, (x * 2) % 256, (y * 3) % 256))
    idat = zlib.compress(bytes(raw), 6)
    png = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IDAT", idat) + chunk(b"IEND", b"")
    path.write_bytes(png)


def _write_real_pdf(path: Path) -> None:
    """Write a tiny but magic-byte-valid, above-floor PDF."""
    body = b"%PDF-1.4\n" + b"% minimal hermetic pdf for the parity gate\n" + (b"0" * 2000) + b"\n%%EOF\n"
    path.write_bytes(body)


def _build_run(tmp_path: Path, *, cells: int = EXPORTED_CELLS, n_pngs: int = PNG_COUNT) -> Path:
    """Create a complete synthetic SOT run; return the manifest path."""
    run_dir = tmp_path / "runs" / "2099-01-01T0000Z-deadbee"
    bundle_dir = run_dir / "python" / "bundle"
    r_dir = run_dir / "r"
    bundle_dir.mkdir(parents=True)
    r_dir.mkdir(parents=True)

    for i in range(n_pngs):
        _write_real_png(r_dir / f"panel_{i}.png")
    _write_real_pdf(r_dir / "report.pdf")

    bundle_manifest = {
        "schema_version": "singlecell_r_bundle_v2.1",
        "source": {"n_obs": 65669, "n_vars": 17267},
        "bundle": {
            "n_cells_exported": cells,
            "markers_present": list(MARKERS),
            "markers_requested": list(MARKERS),
        },
        "files": {
            "X_pca": {"path": "obsm/X_pca.parquet"},
            "X_umap": {"path": "obsm/X_umap.parquet"},
            "obs": {"path": "obs.parquet"},
        },
    }
    (bundle_dir / "bundle_manifest.json").write_text(json.dumps(bundle_manifest, indent=2))

    manifest = {
        "project_id": "hermetic-parity-test",
        "run_id": "2099-01-01T0000Z-deadbee",
        "inputs": {
            "exported_cells": cells,
            "export_seed": 17,
            "markers_requested": list(MARKERS),
        },
        "outputs": {
            "bundle_dir": "python/bundle",
            "r_dir": "r",
            "evidence_dir": "evidence",
            "png_count": PNG_COUNT,
            "pdf_count": PDF_COUNT,
        },
        "verification": [
            "PIL QA found 2/2 PNGs nonblank with dimensions >= 3200x1280",
            "bundle_manifest.json reports schema singlecell_r_bundle_v2.1, X_umap, X_pca",
        ],
    }
    manifest_path = run_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return manifest_path


def test_complete_evidence_chain_passes(tmp_path: Path) -> None:
    manifest_path = _build_run(tmp_path)
    # The gate's path-safety allowlist accepts /tmp; pytest's tmp_path is under
    # /tmp on this host. Guard the assumption so a misconfigured tmp dir is loud.
    assert gate._is_relative_to(manifest_path, Path("/tmp").resolve()), (
        f"tmp_path {manifest_path} is not under /tmp; gate allowlist would reject it"
    )
    report = gate.validate_run(manifest_path, allow_missing_run=False)
    assert report["verdict"] == "pass", report["errors"]
    assert report["errors"] == []
    assert report["expected"]["exported_cells"] == EXPORTED_CELLS
    assert report["expected"]["marker_count"] == len(MARKERS)
    assert report["png_min_dims_from_manifest"] == [3200, 1280]


def test_missing_png_fails(tmp_path: Path) -> None:
    manifest_path = _build_run(tmp_path)
    # Delete one PNG so the recorded png_count is no longer satisfied.
    pngs = sorted((manifest_path.parent / "r").glob("*.png"))
    pngs[0].unlink()
    report = gate.validate_run(manifest_path, allow_missing_run=False)
    assert report["verdict"] == "fail"
    assert any("png" in err.lower() for err in report["errors"]), report["errors"]


def test_wrong_cell_count_fails(tmp_path: Path) -> None:
    # Bundle manifest reports a different cell count than inputs.exported_cells.
    manifest_path = _build_run(tmp_path)
    bundle_manifest_path = manifest_path.parent / "python" / "bundle" / "bundle_manifest.json"
    bm = json.loads(bundle_manifest_path.read_text())
    bm["bundle"]["n_cells_exported"] = EXPORTED_CELLS + 1
    bundle_manifest_path.write_text(json.dumps(bm))
    report = gate.validate_run(manifest_path, allow_missing_run=False)
    assert report["verdict"] == "fail"
    assert any("cell" in err.lower() for err in report["errors"]), report["errors"]


def test_blank_png_fails(tmp_path: Path) -> None:
    # A truncated 1-byte "PNG" (bad magic + below floor) must be rejected.
    manifest_path = _build_run(tmp_path)
    (manifest_path.parent / "r" / "panel_0.png").write_bytes(b"\x00")
    report = gate.validate_run(manifest_path, allow_missing_run=False)
    assert report["verdict"] == "fail"
    assert any("png" in err.lower() for err in report["errors"]), report["errors"]


def test_path_escape_in_manifest_outputs_is_rejected(tmp_path: Path) -> None:
    # A manifest whose outputs.r_dir tries to escape the run dir must SystemExit.
    manifest_path = _build_run(tmp_path)
    manifest = json.loads(manifest_path.read_text())
    manifest["outputs"]["r_dir"] = "../../../../etc"
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit):
        gate.validate_run(manifest_path, allow_missing_run=False)


def test_missing_manifest_strict_raises(tmp_path: Path) -> None:
    # An absent run manifest is a hard SystemExit in strict mode (real
    # enforcement — no silent skip).
    absent_manifest = tmp_path / "absent" / "manifest.json"  # never created
    with pytest.raises(SystemExit):
        gate.validate_run(absent_manifest, allow_missing_run=False)


def test_allow_missing_run_downgrades_to_structural(tmp_path: Path) -> None:
    # With --allow-missing-run an absent manifest downgrades to a
    # structural-only pass (data-free CI/manual probe), not an error.
    absent_manifest = tmp_path / "absent" / "manifest.json"  # never created
    report = gate.validate_run(absent_manifest, allow_missing_run=True)
    assert report["verdict"] == "structural_only_pass"
    assert report.get("structural_only") is True
