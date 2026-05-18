"""Regenerate prepared zarr inputs after schema version bump.

Usage: python scripts/dev/regenerate_prepared_zarr.py <output_dir>

This is a trivial wrapper that invokes the canonical
prepare_emtab13526_full_cohort_zarr.py with appropriate args.
Use after a CheckpointSchemaMismatch error.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_FACTORY_ROOT = Path(__file__).resolve().parent.parent.parent
_PREPARE_SCRIPT = _FACTORY_ROOT / "scripts" / "prepare_emtab13526_full_cohort_zarr.py"


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.stderr.write(
            "ERROR: <output_dir> argument is required.\n"
            "Example: python scripts/dev/regenerate_prepared_zarr.py "
            "/path/to/project/data/raw/nc2024_nsclc_emtab13526/full_cohort\n"
        )
        sys.exit(1)

    output_dir = Path(sys.argv[1]).resolve()
    print(f"Regenerating prepared zarr in: {output_dir}", flush=True)
    print(f"Running: {_PREPARE_SCRIPT}", flush=True)

    cmd = [sys.executable, str(_PREPARE_SCRIPT)] + sys.argv[2:]
    result = subprocess.run(cmd, check=False)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
