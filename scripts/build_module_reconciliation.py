from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


def read_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def read_status(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    return {row["module"]: row for row in rows}


def read_excluded(path: Path) -> list[tuple[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [(row["module"], row["reason"]) for row in rows]


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit(
            "Usage: build_module_reconciliation.py <planned_modules.txt> <excluded_modules.tsv> <module_status.csv> <out.tsv>"
        )

    planned = read_lines(Path(sys.argv[1]))
    excluded = read_excluded(Path(sys.argv[2]))
    status = read_status(Path(sys.argv[3]))
    out = Path(sys.argv[4])

    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["module", "contract_class", "expected", "actual_status", "message_or_reason"])
        for module in planned:
            row = status.get(module)
            writer.writerow(
                [
                    module,
                    "planned",
                    "present",
                    row.get("status", "absent") if row else "absent",
                    row.get("message", "") if row else "missing from module_status.csv",
                ]
            )
        for module, reason in excluded:
            row = status.get(module)
            writer.writerow(
                [
                    module,
                    "excluded",
                    "absent",
                    row.get("status", "absent") if row else "absent",
                    reason if not row else f"unexpectedly present: {row.get('message', '')}",
                ]
            )


if __name__ == "__main__":
    main()
