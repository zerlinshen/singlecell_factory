#!/usr/bin/env python3
"""Extract a sheet from a simple XLSX workbook to TSV without extra deps.

This is intentionally small and stdlib-only so the Trevino reproduction can
read local Cell supplementary tables even when ``openpyxl`` is unavailable.
It supports normal shared-string and inline-string cells, preserving formula
cached values when present. It is not a general Excel calculation engine.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import zipfile
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
NS = {"a": MAIN_NS, "r": REL_NS, "pr": PKG_REL_NS}
RID_ATTR = f"{{{REL_NS}}}id"


def column_number(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    n = 0
    for ch in letters.upper():
        n = n * 26 + (ord(ch) - 64)
    return n


def load_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    strings: list[str] = []
    for si in root.findall("a:si", NS):
        parts = [
            t.text or ""
            for t in si.iter(f"{{{MAIN_NS}}}t")
        ]
        strings.append("".join(parts))
    return strings


def workbook_sheets(zf: zipfile.ZipFile) -> dict[str, str]:
    workbook = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rel_map = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in rels.findall("pr:Relationship", NS)
    }
    sheets: dict[str, str] = {}
    sheets_el = workbook.find("a:sheets", NS)
    if sheets_el is None:
        return sheets
    for sheet in sheets_el:
        rid = sheet.attrib.get(RID_ATTR)
        if not rid or rid not in rel_map:
            continue
        target = rel_map[rid]
        if not target.startswith("xl/"):
            target = f"xl/{target}"
        sheets[sheet.attrib["name"]] = target
    return sheets


def cell_value(cell: ET.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    value = cell.find("a:v", NS)
    if cell_type == "s" and value is not None:
        return shared_strings[int(value.text or "0")]
    if cell_type == "inlineStr":
        inline = cell.find("a:is", NS)
        if inline is None:
            return ""
        return "".join(
            t.text or ""
            for t in inline.iter(f"{{{MAIN_NS}}}t")
        )
    if value is not None:
        return value.text or ""
    return ""


def iter_sheet_rows(
    zf: zipfile.ZipFile,
    sheet_path: str,
    shared_strings: list[str],
) -> Iterable[list[str]]:
    root = ET.fromstring(zf.read(sheet_path))
    for row in root.findall(".//a:sheetData/a:row", NS):
        values: dict[int, str] = {}
        max_col = 0
        for cell in row.findall("a:c", NS):
            ref = cell.attrib.get("r", "")
            col = column_number(ref)
            if col <= 0:
                continue
            values[col] = cell_value(cell, shared_strings)
            max_col = max(max_col, col)
        if max_col:
            yield [values.get(i, "") for i in range(1, max_col + 1)]


def trim_leading_title_rows(rows: Iterable[list[str]]) -> Iterable[list[str]]:
    """Drop workbook title rows before the real table header.

    Trevino Cell supplemental sheets put a prose title in row 1 and the real
    tabular header in row 2. Preserve rows unchanged for sheets that already
    start with a table header.
    """

    buffered: list[list[str]] = []
    for row in rows:
        non_empty = [x for x in row if x != ""]
        if not non_empty:
            continue
        buffered.append(row)
        if len(non_empty) > 1:
            yield row
            break
    for row in rows:
        yield row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xlsx", required=True, type=Path)
    parser.add_argument("--sheet", required=True, help="sheet name, e.g. B")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--trim-title-row",
        action="store_true",
        help="drop leading prose/title rows before the real table header",
    )
    args = parser.parse_args()

    if not args.xlsx.exists():
        parser.error(f"missing workbook: {args.xlsx}")

    with zipfile.ZipFile(args.xlsx) as zf:
        sheets = workbook_sheets(zf)
        if args.sheet not in sheets:
            available = ", ".join(sorted(sheets))
            parser.error(f"sheet {args.sheet!r} not found; available: {available}")
        shared_strings = load_shared_strings(zf)
        rows: Iterable[list[str]] = iter_sheet_rows(zf, sheets[args.sheet], shared_strings)
        if args.trim_title_row:
            rows = trim_leading_title_rows(rows)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        n = 0
        max_cols = 0
        with args.out.open("w", newline="") as fh:
            writer = csv.writer(fh, dialect="excel-tab", lineterminator="\n")
            for row in rows:
                writer.writerow(row)
                n += 1
                max_cols = max(max_cols, len(row))
    print(f"WROTE {args.out} rows={n} max_cols={max_cols}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
