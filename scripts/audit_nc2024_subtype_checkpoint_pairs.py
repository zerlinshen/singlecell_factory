from __future__ import annotations

import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.modular.reporting.nc2024_checkpoint_audit import run_audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit paper-relevant subtype checkpoint pairs from existing subtype LIANA outputs.",
    )
    parser.add_argument("--luad-raw-csv", required=True)
    parser.add_argument("--lusc-raw-csv", required=True)
    parser.add_argument("--luad-top20-csv", required=True)
    parser.add_argument("--lusc-top20-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = run_audit(
        luad_raw_csv=Path(args.luad_raw_csv),
        lusc_raw_csv=Path(args.lusc_raw_csv),
        luad_top20_csv=Path(args.luad_top20_csv),
        lusc_top20_csv=Path(args.lusc_top20_csv),
        output_dir=Path(args.output_dir),
        project=args.project,
    )
    print(out_dir)


if __name__ == "__main__":
    main()
