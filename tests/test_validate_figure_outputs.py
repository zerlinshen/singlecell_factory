from __future__ import annotations

import importlib.util
import struct
import zlib
from pathlib import Path


_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_figure_outputs.py"
_SPEC = importlib.util.spec_from_file_location("validate_figure_outputs", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
validator = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(validator)


def _png_chunk(name: bytes, payload: bytes) -> bytes:
    return (
        struct.pack(">I", len(payload))
        + name
        + payload
        + struct.pack(">I", zlib.crc32(name + payload) & 0xFFFFFFFF)
    )


def _write_png(path: Path, width: int = 2, height: int = 3) -> None:
    rows = b"".join(b"\x00" + (b"\xff\x00\x00" * width) for _ in range(height))
    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    phys = struct.pack(">IIB", 11811, 11811, 1)
    data = (
        validator.PNG_SIGNATURE
        + _png_chunk(b"IHDR", ihdr)
        + _png_chunk(b"pHYs", phys)
        + _png_chunk(b"IDAT", zlib.compress(rows))
        + _png_chunk(b"IEND", b"")
    )
    path.write_bytes(data)


def test_validate_run_dir_accepts_valid_png_svg_and_pdf(tmp_path: Path) -> None:
    _write_png(tmp_path / "umap.png")
    (tmp_path / "panel.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 480">'
        '<rect width="640" height="480"/></svg>\n',
        encoding="utf-8",
    )
    (tmp_path / "supplement.pdf").write_bytes(
        b"%PDF-1.4\n1 0 obj\n<< /Type /Page >>\nendobj\n%%EOF\n"
    )

    report = validator.validate_run_dir(tmp_path, require_vector=True)

    assert report["status"] == "pass"
    assert report["summary"]["figure_count"] == 3
    assert report["summary"]["valid_count"] == 3
    assert report["summary"]["vector_count"] == 2
    png = next(record for record in report["files"] if record["path"] == "umap.png")
    assert png["metadata"]["width"] == 2
    assert png["metadata"]["height"] == 3
    assert png["metadata"]["dpi"]["x"] == 300.0
    assert report["publication_gate"]["manual_visual_review_required"] is True


def test_validate_run_dir_fails_empty_directory_by_default(tmp_path: Path) -> None:
    report = validator.validate_run_dir(tmp_path)

    assert report["status"] == "fail"
    assert report["summary"]["figure_count"] == 0
    assert {finding["code"] for finding in report["findings"]} == {"no_figure_outputs"}


def test_validate_run_dir_fails_missing_directory(tmp_path: Path) -> None:
    report = validator.validate_run_dir(tmp_path / "missing")

    assert report["status"] == "fail"
    assert report["summary"]["figure_count"] == 0
    assert {finding["code"] for finding in report["findings"]} == {"run_dir_missing"}


def test_validate_run_dir_allows_empty_when_requested(tmp_path: Path) -> None:
    report = validator.validate_run_dir(tmp_path, allow_empty=True)

    assert report["status"] == "pass"
    assert report["summary"]["figure_count"] == 0


def test_validate_run_dir_rejects_corrupt_png_crc(tmp_path: Path) -> None:
    _write_png(tmp_path / "bad.png")
    data = bytearray((tmp_path / "bad.png").read_bytes())
    data[-5] ^= 0xFF
    (tmp_path / "bad.png").write_bytes(bytes(data))

    report = validator.validate_run_dir(tmp_path)

    assert report["status"] == "fail"
    assert report["files"][0]["status"] == "fail"
    assert "png_bad_crc" in {finding["code"] for finding in report["files"][0]["findings"]}


def test_validate_run_dir_still_parses_tiny_files(tmp_path: Path) -> None:
    (tmp_path / "tiny.png").write_bytes(b"notpng")

    report = validator.validate_run_dir(tmp_path)

    assert report["status"] == "fail"
    codes = {finding["code"] for finding in report["files"][0]["findings"]}
    assert "small_file" in codes
    assert "png_bad_signature" in codes


def test_validate_run_dir_rejects_svg_without_dimensions(tmp_path: Path) -> None:
    (tmp_path / "panel.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg"><text>missing dimensions</text></svg>\n',
        encoding="utf-8",
    )

    report = validator.validate_run_dir(tmp_path)

    assert report["status"] == "fail"
    assert "svg_missing_dimensions" in {
        finding["code"] for finding in report["files"][0]["findings"]
    }
