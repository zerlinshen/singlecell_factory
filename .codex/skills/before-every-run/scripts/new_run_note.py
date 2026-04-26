#!/usr/bin/env python3
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import re


TEMPLATE = """# {title}

## Objective

- 

## Starting Context

- 

## What Was Run

- Host:
- Working directory:
- Command or controlling script:
- Run directory:

## Outcome

- `success`

## Evidence

- 

## Problems Encountered

- 

## What Was Changed

- 

## Resolution

- 

## Cautions for the Next Run

- 

## Improvement Ideas

- 

## Classification

- `evidence-only`
"""


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("topic", help="Short run/topic name")
    parser.add_argument("--out-dir", default=".", help="Journal directory")
    args = parser.parse_args()

    date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    slug = slugify(args.topic)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    title = f"{date} - {args.topic}"
    out_path = out_dir / f"{date}-{slug}.md"
    out_path.write_text(TEMPLATE.format(title=title), encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
