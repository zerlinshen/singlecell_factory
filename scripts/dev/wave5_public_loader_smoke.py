#!/usr/bin/env python3
"""Bounded loader smoke for Trevino public matrices.

The goal is to verify that the public GSE162170 matrices are structurally
loadable and mutually consistent without materialising full dense matrices.
The checks only read headers, the first data row, and small metadata/cell-name
files.  Large count matrices are not converted to arrays here.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
from pathlib import Path
from typing import Any


def _first_lines(path: Path, n: int = 2) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    lines: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:  # type: ignore[arg-type]
        for _ in range(n):
            line = fh.readline()
            if not line:
                break
            lines.append(line.rstrip("\n\r"))
    return lines


def _count_rows(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    count = 0
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:  # type: ignore[arg-type]
        for count, _ in enumerate(fh, start=1):
            pass
    return count


def _field_counts(line: str) -> dict[str, int]:
    return {"tab": len(line.split("\t")), "whitespace": len(line.split())}


def _smart_fields(line: str) -> list[str]:
    tab = line.split("\t")
    if len(tab) > 1:
        return tab
    return line.split()


def _matrix_head(path: Path, *, header_has_gene_label: bool = False) -> dict[str, Any]:
    lines = _first_lines(path, 2)
    out: dict[str, Any] = {
        "path": str(path),
        "exists": path.exists(),
        "bytes": path.stat().st_size if path.exists() else None,
        "readable": bool(lines),
    }
    if not lines:
        return out
    first_fields = _smart_fields(lines[0])
    out.update(
        {
            "line1_field_counts": _field_counts(lines[0]),
            "line1_smart_fields": len(first_fields),
            "line1_prefix": lines[0][:120],
        }
    )
    if len(lines) > 1:
        second_fields = _smart_fields(lines[1])
        out.update(
            {
                "line2_field_counts": _field_counts(lines[1]),
                "line2_smart_fields": len(second_fields),
                "line2_first_token": second_fields[0] if second_fields else None,
                "line2_prefix": lines[1][:120],
            }
        )
    if header_has_gene_label:
        out["header_cell_count"] = max(len(first_fields) - 1, 0)
        out["first_data_field_count_expected"] = len(first_fields)
    else:
        out["header_cell_count"] = len(first_fields)
        out["first_data_field_count_expected"] = len(first_fields) + 1
    out["first_data_has_gene_plus_cells"] = (
        out.get("line2_smart_fields") == out.get("first_data_field_count_expected")
        if "line2_smart_fields" in out
        else None
    )
    return out


def _metadata_contract(path: Path, id_col: str) -> dict[str, Any]:
    lines = _first_lines(path, 2)
    rows = _count_rows(path) if path.exists() else 0
    header = lines[0].split("\t") if lines else []
    return {
        "path": str(path),
        "exists": path.exists(),
        "rows_including_header": rows,
        "header_first": header[0] if header else None,
        "id_col_ok": bool(header and header[0] == id_col),
        "line2_prefix": lines[1][:120] if len(lines) > 1 else None,
    }


def _cell_names_contract(path: Path) -> dict[str, Any]:
    lines = _first_lines(path, 3)
    return {
        "path": str(path),
        "exists": path.exists(),
        "rows": _count_rows(path) if path.exists() else 0,
        "head": lines,
    }


def _rna_contract(input_dir: Path) -> dict[str, Any]:
    counts = _matrix_head(input_dir / "GSE162170_rna_counts.tsv.gz")
    spliced = _matrix_head(input_dir / "GSE162170_rna_spliced_counts.tsv.gz")
    unspliced = _matrix_head(input_dir / "GSE162170_rna_unspliced_counts.tsv.gz")
    ambiguous = _matrix_head(input_dir / "GSE162170_rna_ambiguous_counts.tsv.gz")
    meta = _metadata_contract(input_dir / "GSE162170_rna_cell_metadata.txt.gz", "Cell.ID")
    names = _cell_names_contract(input_dir / "GSE162170_rna_cell_names.txt.gz")
    cells = counts.get("header_cell_count")
    same_cells = all(
        m.get("header_cell_count") == cells
        for m in [counts, spliced, unspliced, ambiguous]
    )
    meta_ok = meta.get("rows_including_header", 0) - 1 == cells
    names_ok = names.get("rows") == cells
    return {
        "counts": counts,
        "spliced": spliced,
        "unspliced": unspliced,
        "ambiguous": ambiguous,
        "metadata": meta,
        "cell_names": names,
        "cell_count": cells,
        "all_matrices_share_cell_count": same_cells,
        "metadata_rows_match_cells": meta_ok,
        "cell_names_match_cells": names_ok,
        "ok": bool(cells and same_cells and meta_ok and names_ok),
    }


def _atac_contract(input_dir: Path) -> dict[str, Any]:
    counts = _matrix_head(input_dir / "GSE162170_atac_counts.tsv.gz", header_has_gene_label=True)
    gene_activity = _matrix_head(input_dir / "GSE162170_atac_gene_activities.tsv.gz")
    chromvar = _matrix_head(input_dir / "GSE162170_atac_chromVAR_TF_activities.tsv.gz")
    meta = _metadata_contract(input_dir / "GSE162170_atac_cell_metadata.txt.gz", "Cell.ID")
    cells = counts.get("line1_smart_fields")
    # ATAC counts have no header or feature-id column; one line is one peak.
    counts["cell_count"] = cells
    same_cells = all(
        m.get("header_cell_count") == cells
        for m in [gene_activity, chromvar]
    )
    meta_ok = meta.get("rows_including_header", 0) - 1 == cells
    return {
        "counts": counts,
        "gene_activity": gene_activity,
        "chromvar": chromvar,
        "metadata": meta,
        "cell_count": cells,
        "gene_activity_and_chromvar_match_counts_cells": same_cells,
        "metadata_rows_match_cells": meta_ok,
        "ok": bool(cells and same_cells and meta_ok),
    }


def _multiome_contract(input_dir: Path) -> dict[str, Any]:
    rna = _matrix_head(input_dir / "GSE162170_multiome_rna_counts.tsv.gz")
    logcounts = _matrix_head(input_dir / "GSE162170_multiome_rna_logcounts.tsv.gz")
    gene_activity = _matrix_head(input_dir / "GSE162170_multiome_atac_gene_activities.tsv.gz")
    atac = _matrix_head(input_dir / "GSE162170_multiome_atac_counts.tsv.gz", header_has_gene_label=True)
    spliced = _matrix_head(input_dir / "GSE162170_multiome_spliced_rna_counts.tsv.gz", header_has_gene_label=True)
    unspliced = _matrix_head(input_dir / "GSE162170_multiome_unspliced_rna_counts.tsv.gz", header_has_gene_label=True)
    meta = _metadata_contract(input_dir / "GSE162170_multiome_cell_metadata.txt.gz", "Cell.ID")
    cells = rna.get("header_cell_count")
    atac["cell_count"] = atac.get("line1_smart_fields")
    same_cells = (
        logcounts.get("header_cell_count") == cells
        and gene_activity.get("header_cell_count") == cells
        and atac.get("cell_count") == cells
        and spliced.get("header_cell_count") == cells
        and unspliced.get("header_cell_count") == cells
    )
    meta_ok = meta.get("rows_including_header", 0) - 1 == cells
    return {
        "rna": rna,
        "rna_logcounts": logcounts,
        "atac": atac,
        "atac_gene_activity": gene_activity,
        "spliced": spliced,
        "unspliced": unspliced,
        "metadata": meta,
        "cell_count": cells,
        "all_modalities_share_cell_count": same_cells,
        "metadata_rows_match_cells": meta_ok,
        "ok": bool(cells and same_cells and meta_ok),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--input-dir", required=True, type=Path)
    args = ap.parse_args()

    payload = {
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "run_dir": str(args.run_dir),
        "input_dir": str(args.input_dir),
        "policy": "headers and small metadata only; no dense matrix materialisation",
        "rna": _rna_contract(args.input_dir),
        "atac": _atac_contract(args.input_dir),
        "multiome": _multiome_contract(args.input_dir),
    }
    payload["ok"] = bool(payload["rna"]["ok"] and payload["atac"]["ok"] and payload["multiome"]["ok"])
    out_dir = args.run_dir / "python" / "loader_smoke"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "public_loader_smoke.json").write_text(json.dumps(payload, indent=2) + "\n")
    md = [
        "# Wave-5 Public Matrix Loader Smoke",
        "",
        f"Created UTC: {payload['created_at_utc']}",
        f"Input dir: `{args.input_dir}`",
        "",
        f"- RNA cells: `{payload['rna']['cell_count']}`; ok: `{payload['rna']['ok']}`",
        f"- ATAC cells: `{payload['atac']['cell_count']}`; ok: `{payload['atac']['ok']}`",
        f"- Multiome cells: `{payload['multiome']['cell_count']}`; ok: `{payload['multiome']['ok']}`",
        f"- Overall ok: `{payload['ok']}`",
        "",
        "This smoke intentionally avoids dense loading of full matrices.",
        "",
    ]
    (out_dir / "public_loader_smoke.md").write_text("\n".join(md))
    print(
        json.dumps(
            {
                "out_dir": str(out_dir),
                "ok": payload["ok"],
                "rna_cells": payload["rna"]["cell_count"],
                "atac_cells": payload["atac"]["cell_count"],
                "multiome_cells": payload["multiome"]["cell_count"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
