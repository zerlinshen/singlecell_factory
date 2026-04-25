from __future__ import annotations

import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.modular.reporting.nc2024_checkpoint_audit import summarize_hierarchy


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize the paper-facing checkpoint hierarchy from a formalized audit directory.",
    )
    parser.add_argument("--audit-dir", required=True)
    parser.add_argument("--output-dir", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = summarize_hierarchy(
        audit_dir=Path(args.audit_dir),
        output_dir=Path(args.output_dir) if args.output_dir else None,
    )
    print(out_dir)


if __name__ == "__main__":
    main()
