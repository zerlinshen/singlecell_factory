#!/usr/bin/env python3
"""Validate mechanical integrity of pipeline figure outputs.

This is a control-plane QA gate for run artifacts. It checks that PNG/PDF/SVG
files are present, non-empty, parseable, and have basic dimensions where the
format exposes them cheaply. It does not replace human or visual-review checks
for biological correctness, label overlap, palette quality, or publication
layout.
"""
from __future__ import annotations

import argparse
import json
import re
import struct
import zlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from xml.etree import ElementTree


SUPPORTED_SUFFIXES = {".png", ".pdf", ".svg"}
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
SVG_DIMENSION_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?)")


MANUAL_REVIEW_CHECKLIST = [
    "scientific_correctness",
    "readable_at_manuscript_scale",
    "no_label_or_legend_overlap",
    "consistent_theme_and_palette",
    "axes_legends_and_titles_not_clipped",
    "cell_type_colors_distinguishable_in_print",
    "debug_plots_separated_from_manuscript_panels",
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _finding(severity: str, code: str, message: str) -> dict[str, str]:
    return {"severity": severity, "code": code, "message": message}


def _parse_svg_dimension(value: str | None) -> float | None:
    if value is None:
        return None
    match = SVG_DIMENSION_RE.match(value)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _parse_viewbox(value: str | None) -> tuple[float, float] | None:
    if not value:
        return None
    parts = re.split(r"[\s,]+", value.strip())
    if len(parts) != 4:
        return None
    try:
        width = float(parts[2])
        height = float(parts[3])
    except ValueError:
        return None
    if width <= 0 or height <= 0:
        return None
    return width, height


def _validate_png(path: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    data = path.read_bytes()
    findings: list[dict[str, str]] = []
    meta: dict[str, Any] = {"format": "png", "bytes": len(data), "vector": False}

    if not data.startswith(PNG_SIGNATURE):
        findings.append(_finding("error", "png_bad_signature", "PNG signature is invalid."))
        return meta, findings

    offset = len(PNG_SIGNATURE)
    width: int | None = None
    height: int | None = None
    bit_depth: int | None = None
    color_type: int | None = None
    iend_seen = False
    idat_seen = False

    while offset < len(data):
        if offset + 12 > len(data):
            findings.append(_finding("error", "png_truncated_chunk", "PNG chunk header is truncated."))
            return meta, findings
        length = struct.unpack(">I", data[offset : offset + 4])[0]
        chunk_type = data[offset + 4 : offset + 8]
        chunk_start = offset + 8
        chunk_end = chunk_start + length
        crc_end = chunk_end + 4
        if crc_end > len(data):
            findings.append(_finding("error", "png_truncated_chunk", "PNG chunk payload is truncated."))
            return meta, findings
        chunk_data = data[chunk_start:chunk_end]
        expected_crc = struct.unpack(">I", data[chunk_end:crc_end])[0]
        actual_crc = zlib.crc32(chunk_type + chunk_data) & 0xFFFFFFFF
        if actual_crc != expected_crc:
            name = chunk_type.decode("latin1", errors="replace")
            findings.append(_finding("error", "png_bad_crc", f"PNG chunk {name} has an invalid CRC."))
            return meta, findings

        if chunk_type == b"IHDR":
            if length != 13:
                findings.append(_finding("error", "png_bad_ihdr", "PNG IHDR chunk has invalid length."))
                return meta, findings
            width, height = struct.unpack(">II", chunk_data[:8])
            bit_depth = chunk_data[8]
            color_type = chunk_data[9]
        elif chunk_type == b"pHYs" and length == 9:
            x_ppm, y_ppm, unit = struct.unpack(">IIB", chunk_data)
            if unit == 1 and x_ppm > 0 and y_ppm > 0:
                meta["dpi"] = {
                    "x": round(x_ppm * 0.0254, 2),
                    "y": round(y_ppm * 0.0254, 2),
                }
        elif chunk_type == b"IDAT":
            idat_seen = True
        elif chunk_type == b"IEND":
            iend_seen = True
            break
        offset = crc_end

    meta.update({
        "width": width,
        "height": height,
        "bit_depth": bit_depth,
        "color_type": color_type,
    })

    if width is None or height is None or width <= 0 or height <= 0:
        findings.append(_finding("error", "png_missing_dimensions", "PNG dimensions are missing or invalid."))
    if not idat_seen:
        findings.append(_finding("error", "png_missing_idat", "PNG contains no IDAT image data chunk."))
    if not iend_seen:
        findings.append(_finding("error", "png_missing_iend", "PNG IEND chunk is missing."))
    return meta, findings


def _validate_pdf(path: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    data = path.read_bytes()
    findings: list[dict[str, str]] = []
    meta: dict[str, Any] = {"format": "pdf", "bytes": len(data), "vector": True}
    if not data.startswith(b"%PDF-"):
        findings.append(_finding("error", "pdf_bad_header", "PDF header is missing."))
    if b"%%EOF" not in data[-2048:]:
        findings.append(_finding("error", "pdf_missing_eof", "PDF EOF marker is missing."))
    if b"/Page" not in data:
        findings.append(_finding("warning", "pdf_no_page_marker", "No PDF /Page marker found; inspect manually."))
    return meta, findings


def _validate_svg(path: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    findings: list[dict[str, str]] = []
    meta: dict[str, Any] = {"format": "svg", "bytes": path.stat().st_size, "vector": True}
    try:
        root = ElementTree.parse(path).getroot()
    except ElementTree.ParseError as exc:
        findings.append(_finding("error", "svg_parse_error", f"SVG XML parse failed: {exc}."))
        return meta, findings

    if not root.tag.lower().endswith("svg"):
        findings.append(_finding("error", "svg_bad_root", "SVG root element is not <svg>."))
        return meta, findings

    width = _parse_svg_dimension(root.attrib.get("width"))
    height = _parse_svg_dimension(root.attrib.get("height"))
    viewbox = _parse_viewbox(root.attrib.get("viewBox"))
    if (width is None or height is None) and viewbox is not None:
        width, height = viewbox

    meta.update({"width": width, "height": height, "viewBox": root.attrib.get("viewBox")})
    if width is None or height is None:
        findings.append(_finding(
            "error",
            "svg_missing_dimensions",
            "SVG must include numeric width/height or a valid viewBox.",
        ))
    return meta, findings


def discover_figures(run_dir: Path) -> list[Path]:
    return sorted(
        path for path in run_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in SUPPORTED_SUFFIXES
    )


def validate_figure(path: Path, root: Path, min_file_bytes: int) -> dict[str, Any]:
    findings: list[dict[str, str]] = []
    suffix = path.suffix.lower()
    size = path.stat().st_size
    if size == 0:
        findings.append(_finding("error", "empty_file", "Figure file is empty."))
        meta: dict[str, Any] = {"format": suffix.lstrip("."), "bytes": 0}
    else:
        if size < min_file_bytes:
            findings.append(_finding(
                "warning",
                "small_file",
                f"Figure file is smaller than {min_file_bytes} bytes; inspect manually.",
            ))
        if suffix == ".png":
            meta, parsed_findings = _validate_png(path)
        elif suffix == ".pdf":
            meta, parsed_findings = _validate_pdf(path)
        elif suffix == ".svg":
            meta, parsed_findings = _validate_svg(path)
        else:
            meta = {"format": suffix.lstrip("."), "bytes": size}
            parsed_findings = [_finding("error", "unsupported_format", "Unsupported figure format.")]
        findings.extend(parsed_findings)

    error_count = sum(1 for finding in findings if finding["severity"] == "error")
    warning_count = sum(1 for finding in findings if finding["severity"] == "warning")
    return {
        "path": _rel(path, root),
        "status": "fail" if error_count else "pass",
        "error_count": error_count,
        "warning_count": warning_count,
        "metadata": meta,
        "findings": findings,
    }


def validate_run_dir(
    run_dir: Path,
    *,
    min_file_bytes: int = 24,
    allow_empty: bool = False,
    fail_on_warning: bool = False,
    require_vector: bool = False,
) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    findings: list[dict[str, str]] = []

    if not run_dir.is_dir():
        findings.append(_finding("error", "run_dir_missing", "Run directory does not exist or is not a directory."))
        return {
            "schema_version": "figure-output-validation/v1",
            "generated_at": _utc_now(),
            "run_dir": str(run_dir),
            "status": "fail",
            "summary": {
                "figure_count": 0,
                "valid_count": 0,
                "failed_count": 0,
                "warning_count": 0,
                "vector_count": 0,
                "raster_count": 0,
            },
            "publication_gate": {
                "mechanical_integrity": "fail",
                "manual_visual_review_required": True,
                "manual_review_checklist": MANUAL_REVIEW_CHECKLIST,
            },
            "findings": findings,
            "files": [],
        }

    figures = discover_figures(run_dir)
    records = [validate_figure(path, run_dir, min_file_bytes) for path in figures]

    if not figures and not allow_empty:
        findings.append(_finding(
            "error",
            "no_figure_outputs",
            "No PNG/PDF/SVG figure outputs found under run directory.",
        ))

    vector_count = sum(
        1 for record in records
        if record["metadata"].get("vector") is True and record["status"] == "pass"
    )
    if require_vector and vector_count == 0:
        findings.append(_finding(
            "error",
            "missing_vector_output",
            "No valid PDF/SVG vector figure found, but --require-vector was set.",
        ))

    file_errors = sum(record["error_count"] for record in records)
    file_warnings = sum(record["warning_count"] for record in records)
    top_errors = sum(1 for finding in findings if finding["severity"] == "error")
    top_warnings = sum(1 for finding in findings if finding["severity"] == "warning")
    warning_failure = fail_on_warning and (file_warnings + top_warnings) > 0

    status = "fail" if file_errors or top_errors or warning_failure else "pass"
    return {
        "schema_version": "figure-output-validation/v1",
        "generated_at": _utc_now(),
        "run_dir": str(run_dir),
        "status": status,
        "summary": {
            "figure_count": len(records),
            "valid_count": sum(1 for record in records if record["status"] == "pass"),
            "failed_count": sum(1 for record in records if record["status"] == "fail"),
            "warning_count": file_warnings + top_warnings,
            "vector_count": vector_count,
            "raster_count": sum(
                1 for record in records
                if record["metadata"].get("vector") is False and record["status"] == "pass"
            ),
        },
        "publication_gate": {
            "mechanical_integrity": status,
            "manual_visual_review_required": True,
            "manual_review_checklist": MANUAL_REVIEW_CHECKLIST,
        },
        "findings": findings,
        "files": records,
    }


def render_markdown(report: dict[str, Any]) -> str:
    lines = [
        f"# Figure Output Validation: {Path(report['run_dir']).name}",
        "",
        "Control-plane figure QA report; mechanical checks only.",
        "",
        f"- Generated: `{report['generated_at']}`",
        f"- Run directory: `{report['run_dir']}`",
        f"- Status: `{report['status']}`",
        f"- Summary: `{report['summary']}`",
        f"- Manual visual review required: `{report['publication_gate']['manual_visual_review_required']}`",
        "",
        "## Manual Review Checklist",
    ]
    lines.extend(f"- `{item}`" for item in report["publication_gate"]["manual_review_checklist"])
    if report["findings"]:
        lines.extend(["", "## Run Findings"])
        lines.extend(
            f"- `{finding['severity']}` `{finding['code']}`: {finding['message']}"
            for finding in report["findings"]
        )
    lines.extend(["", "## Files"])
    for record in report["files"]:
        lines.append(
            f"- `{record['status']}` `{record['path']}` "
            f"errors={record['error_count']} warnings={record['warning_count']}"
        )
        for finding in record["findings"]:
            lines.append(f"  - `{finding['severity']}` `{finding['code']}`: {finding['message']}")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--json-out", type=Path, default=None)
    parser.add_argument("--md-out", type=Path, default=None)
    parser.add_argument("--json", action="store_true", help="Print JSON instead of Markdown.")
    parser.add_argument("--allow-empty", action="store_true", help="Pass when no figure outputs are present.")
    parser.add_argument("--fail-on-warning", action="store_true", help="Return non-zero on warnings.")
    parser.add_argument("--require-vector", action="store_true", help="Require at least one valid PDF or SVG.")
    parser.add_argument("--min-file-bytes", type=int, default=24)
    args = parser.parse_args()

    report = validate_run_dir(
        args.run_dir,
        min_file_bytes=args.min_file_bytes,
        allow_empty=args.allow_empty,
        fail_on_warning=args.fail_on_warning,
        require_vector=args.require_vector,
    )
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.md_out:
        args.md_out.parent.mkdir(parents=True, exist_ok=True)
        args.md_out.write_text(render_markdown(report), encoding="utf-8")

    print(json.dumps(report, indent=2) if args.json else render_markdown(report))
    return 1 if report["status"] == "fail" else 0


if __name__ == "__main__":
    raise SystemExit(main())
