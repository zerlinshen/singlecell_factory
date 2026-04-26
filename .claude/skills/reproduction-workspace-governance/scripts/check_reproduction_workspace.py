#!/usr/bin/env python3
"""Check the minimal human/agent governance layout for a reproduction workspace."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


REQUIRED_PATHS = [
    "00_HUMAN_START_HERE.md",
    "human_review/README.md",
    "human_review/figures/README.md",
    "human_review/comparisons/README.md",
    "human_review/roi/README.md",
    "agent_runs/README.md",
    "agent_runs/TEMPLATE_RUN_README.md",
    "README.md",
    "AGENTS.md",
    "CODEX_PROFILE.md",
    "SKILL_PLAYBOOK.md",
]

README_KEYWORDS = [
    "00_HUMAN_START_HERE.md",
    "human_review",
    "agent_runs",
    "remote R",
]

AGENT_KEYWORDS = [
    "agent_runs",
    "before-every-run",
]


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("workspace_root", type=Path)
    args = parser.parse_args()

    root = args.workspace_root.expanduser().resolve()
    failures: list[str] = []

    for rel in REQUIRED_PATHS:
        path = root / rel
        if not path.exists():
            failures.append(f"missing: {rel}")

    if (root / "README.md").exists():
        text = read_text(root / "README.md")
        for keyword in README_KEYWORDS:
            if keyword not in text:
                failures.append(f"README.md does not mention: {keyword}")

    if (root / "AGENTS.md").exists():
        text = read_text(root / "AGENTS.md")
        for keyword in AGENT_KEYWORDS:
            if keyword not in text:
                failures.append(f"AGENTS.md does not mention: {keyword}")

    run_dirs = [
        path
        for path in (root / "agent_runs").glob("20*-*")
        if path.is_dir() and (path / "RUN.md").exists()
    ] if (root / "agent_runs").exists() else []
    if not run_dirs:
        failures.append("agent_runs has no dated RUN.md record")

    if failures:
        for failure in failures:
            print(f"FAIL {failure}")
        return 1

    print(f"OK reproduction workspace governance: {root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
