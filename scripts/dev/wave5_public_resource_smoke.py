#!/usr/bin/env python3
"""Smoke-check public Trevino resource files without materialising matrices.

The full GSE162170 matrices are too large to load just to validate resource
presence.  This script performs bounded checks:

* manifest status versus local file size,
* gzip magic/first-line readability for complete files,
* format contract checks for the currently available PCW21 multiome files, and
* in-progress ``*.parts`` download accounting.

Outputs are written to the project run directory.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
from pathlib import Path
from typing import Any


def _first_lines(path: Path, n: int = 2) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    mode = "rt"
    lines: list[str] = []
    with opener(path, mode, encoding="utf-8", errors="replace", newline="") as fh:  # type: ignore[arg-type]
        for _ in range(n):
            line = fh.readline()
            if not line:
                break
            lines.append(line.rstrip("\n\r"))
    return lines


def _count_gzip_rows(path: Path, limit: int | None = None) -> int:
    count = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        for count, _ in enumerate(fh, start=1):
            if limit and count >= limit:
                break
    return count


def _shaish_status(row: dict[str, Any], path: Path) -> dict[str, Any]:
    expected = row.get("remote_bytes")
    actual = path.stat().st_size if path.exists() else None
    return {
        "filename": path.name,
        "manifest_status": row.get("status"),
        "expected_bytes": expected,
        "actual_bytes": actual,
        "size_matches_manifest": (actual == expected) if actual is not None and expected else None,
        "exists": path.exists(),
        "is_symlink": path.is_symlink(),
        "sha256_recorded": bool(row.get("sha256")),
    }


def _file_smoke(path: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"filename": path.name, "exists": path.exists()}
    if not path.exists():
        return out
    out["bytes"] = path.stat().st_size
    try:
        lines = _first_lines(path, 2)
        out["readable"] = True
        out["first_line_fields"] = len(lines[0].split("\t")) if lines else 0
        out["second_line_fields"] = len(lines[1].split("\t")) if len(lines) > 1 else None
        out["first_line_prefix"] = lines[0][:160] if lines else ""
    except Exception as e:  # pragma: no cover - diagnostic path
        out["readable"] = False
        out["read_error"] = repr(e)
    return out


def _multiome_contract(input_dir: Path) -> dict[str, Any]:
    rna = input_dir / "GSE162170_multiome_rna_counts.tsv.gz"
    atac = input_dir / "GSE162170_multiome_atac_counts.tsv.gz"
    meta = input_dir / "GSE162170_multiome_cell_metadata.txt.gz"
    clusters = input_dir / "GSE162170_multiome_cluster_names.txt.gz"
    peaks = input_dir / "GSE162170_multiome_atac_consensus_peaks.txt.gz"
    result: dict[str, Any] = {
        "files_present": all(p.exists() for p in [rna, atac, meta, clusters, peaks])
    }
    if not result["files_present"]:
        return result
    rna_lines = _first_lines(rna, 2)
    atac_lines = _first_lines(atac, 1)
    meta_lines = _first_lines(meta, 2)
    cluster_lines = _first_lines(clusters, 3)
    peak_lines = _first_lines(peaks, 2)
    rna_cells = rna_lines[0].split("\t")
    rna_data_fields = len(rna_lines[1].split("\t")) if len(rna_lines) > 1 else None
    atac_fields = len(atac_lines[0].split("\t")) if atac_lines else None
    result.update(
        {
            "rna_header_cells": len(rna_cells),
            "rna_first_gene": rna_lines[1].split("\t", 1)[0] if len(rna_lines) > 1 else None,
            "rna_first_data_fields": rna_data_fields,
            "atac_first_row_fields": atac_fields,
            "rna_data_fields_equals_cells_plus_gene_id": rna_data_fields == len(rna_cells) + 1,
            "atac_fields_equals_rna_cells": atac_fields == len(rna_cells),
            "metadata_header_has_cell_id": bool(meta_lines and meta_lines[0].split("\t")[0] == "Cell.ID"),
            "metadata_rows_including_header": _count_gzip_rows(meta),
            "cluster_header": cluster_lines[0] if cluster_lines else None,
            "peak_header": peak_lines[0] if peak_lines else None,
            "peak_first_id": peak_lines[1].split("\t")[5] if len(peak_lines) > 1 and len(peak_lines[1].split("\t")) > 5 else None,
        }
    )
    result["metadata_rows_minus_header_equals_rna_cells"] = result["metadata_rows_including_header"] - 1 == len(rna_cells)
    result["ok"] = all(
        [
            result["rna_data_fields_equals_cells_plus_gene_id"],
            result["atac_fields_equals_rna_cells"],
            result["metadata_header_has_cell_id"],
            result["metadata_rows_minus_header_equals_rna_cells"],
        ]
    )
    return result


def _parts_status(input_dir: Path) -> list[dict[str, Any]]:
    rows = []
    for part_dir in sorted(input_dir.glob("*.parts")):
        rows.append(
            {
                "parts_dir": str(part_dir),
                "target": part_dir.name.removesuffix(".parts"),
                "part_count": len(list(part_dir.glob("part_*"))),
                "bytes": sum(p.stat().st_size for p in part_dir.glob("part_*")),
            }
        )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--input-dir", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path)
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    manifest_rows = []
    file_smokes = []
    for row in manifest.get("geo_supplementary_files", []):
        path = Path(row.get("project_input_path", ""))
        manifest_rows.append(_shaish_status(row, path))
        if path.exists() and not path.name.endswith(".parts"):
            file_smokes.append(_file_smoke(path))

    payload = {
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "run_dir": str(args.run_dir),
        "input_dir": str(args.input_dir),
        "manifest": str(args.manifest),
        "manifest_rows": manifest_rows,
        "file_smokes": file_smokes,
        "multiome_contract": _multiome_contract(args.input_dir),
        "in_progress_parts": _parts_status(args.input_dir),
    }
    out_dir = args.run_dir / "python" / "resource_smoke"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "resource_smoke.json").write_text(json.dumps(payload, indent=2) + "\n")
    status_counts: dict[str, int] = {}
    for row in manifest.get("geo_supplementary_files", []):
        status_counts[row.get("status", "unknown")] = status_counts.get(row.get("status", "unknown"), 0) + 1
    md = [
        "# Wave-5 Public Resource Smoke Check",
        "",
        f"Created UTC: {payload['created_at_utc']}",
        f"Input dir: `{args.input_dir}`",
        "",
        "## Manifest status counts",
        "",
    ]
    for k in sorted(status_counts):
        md.append(f"- `{k}`: {status_counts[k]}")
    md.extend(
        [
            "",
            "## Multiome contract",
            "",
            "```json",
            json.dumps(payload["multiome_contract"], indent=2),
            "```",
            "",
            "## In-progress parts",
            "",
            "```json",
            json.dumps(payload["in_progress_parts"], indent=2),
            "```",
            "",
        ]
    )
    (out_dir / "resource_smoke.md").write_text("\n".join(md) + "\n")
    print(json.dumps({"out_dir": str(out_dir), "status_counts": status_counts, "multiome_ok": payload["multiome_contract"].get("ok")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
