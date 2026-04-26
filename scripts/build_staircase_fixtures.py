#!/usr/bin/env python3
"""Build staircase test fixtures from the NC2024 NSCLC prepared zarr.

Tiers
-----
nano        : 5k cells, 2 samples  — pure-synthetic, no disk dependency
small_real  : ~100k cells, 10 samples  — from prepared_input.zarr; exercises
              auto-grouped-Scrublet trigger threshold (n_obs >= 100k)
medium_real : ~500k cells, 50 samples — from prepared_input.zarr; exercises
              memory + integration behavior at near-production scale
full_real   : all 884 050 cells, 81 samples — symlink / copy of prepared_input.zarr

Only small_real and medium_real are written to tests/data/staircase/ as zarr
subsets.  nano is generated synthetically at test-collection time via conftest.
full_real uses the original zarr in-place (path recorded in tests/data/staircase/full_real.path).

Usage
-----
    python scripts/build_staircase_fixtures.py [--tiers small_real medium_real full_real]
    python scripts/build_staircase_fixtures.py --help
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import zarr

ROOT = Path(__file__).resolve().parents[1]
SRC_ZARR = ROOT / "data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr"
DEST_DIR = ROOT / "tests/data/staircase"

TIER_SPECS: dict[str, dict] = {
    # small_real must exceed the 100k auto-grouped-Scrublet threshold; this is
    # the smallest tier that exercises real-data routing decisions.
    "small_real": {"n_cells": 100_000, "n_samples": 10},
    # medium_real: ~half of full cohort, exercises memory + integration
    # behavior that 100k cannot reach (HVG sparsity / Harmony scaling).
    "medium_real": {"n_cells": 500_000, "n_samples": 50},
    "full_real": {"n_cells": None, "n_samples": None},  # handled specially
}


def _open_src(path: Path) -> zarr.Group:
    if not path.exists():
        print(f"ERROR: source zarr not found: {path}", file=sys.stderr)
        sys.exit(1)
    return zarr.open_group(str(path), mode="r")


def _read_csr(z: zarr.Group) -> sp.csr_matrix:
    data = z["X"]["data"][:]
    indices = z["X"]["indices"][:]
    indptr = z["X"]["indptr"][:]
    n_cells = len(indptr) - 1
    n_genes = int(indices.max()) + 1 if len(indices) > 0 else 0
    return sp.csr_matrix((data, indices, indptr), shape=(n_cells, n_genes))


def _get_sample_codes(z: zarr.Group) -> np.ndarray:
    return z["obs"]["sample"]["codes"][:]


def _get_sample_categories(z: zarr.Group) -> np.ndarray:
    return z["obs"]["sample"]["categories"][:]


def _select_cells(
    z: zarr.Group,
    n_cells: int,
    n_samples: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return sorted cell indices covering >= n_samples samples, total >= n_cells.

    Uses the largest samples first so that the total cell budget is reliably
    met even when individual samples are small. Always includes >= n_samples
    distinct samples; expands beyond when needed to reach the cell budget.
    """
    codes = _get_sample_codes(z)
    cats = _get_sample_categories(z)
    sample_sizes = np.array([int((codes == s).sum()) for s in range(len(cats))])

    # Rank samples by size descending so we always have enough cells available.
    order_desc = np.argsort(sample_sizes)[::-1]
    chosen_samples: list[int] = []
    cumulative = 0
    for s in order_desc:
        size = int(sample_sizes[s])
        if size == 0:
            continue
        chosen_samples.append(int(s))
        cumulative += size
        if len(chosen_samples) >= n_samples and cumulative >= n_cells:
            break

    if cumulative < n_cells:
        print(
            f"  warning: requested {n_cells:,} cells but only {cumulative:,} "
            f"available across {len(chosen_samples)} samples; using all of them.",
            flush=True,
        )

    # Pool all candidate cells from chosen samples
    pool_masks = [np.where(codes == s)[0] for s in chosen_samples]
    pool = np.concatenate(pool_masks) if pool_masks else np.array([], dtype=np.int64)
    if len(pool) <= n_cells:
        idx = pool
    else:
        idx = rng.choice(pool, size=n_cells, replace=False)
    return np.sort(idx)


def _write_anndata_subset(
    src_zarr: Path,
    cell_idx: np.ndarray,
    dest: Path,
    tier_name: str,
) -> None:
    """Read source via anndata.read_zarr, slice, write canonical h5ad.

    Using anndata's writer ensures encoding metadata that newer versions of
    anndata require for read-back. The subset is written as ``<tier>.h5ad``
    next to a ``<tier>.summary.json`` provenance sidecar.
    """
    import anndata as ad  # local import keeps script importable on lean envs

    print(f"  reading source zarr (full cohort) ...", flush=True)
    full = ad.read_zarr(src_zarr)
    print(f"  source: {full.n_obs:,} cells x {full.n_vars:,} genes", flush=True)
    sub = full[cell_idx, :].copy()
    del full
    dest.parent.mkdir(parents=True, exist_ok=True)
    h5ad_path = dest.with_suffix(".h5ad")
    print(f"  writing {sub.n_obs:,} cells -> {h5ad_path}", flush=True)
    sub.write_h5ad(h5ad_path)
    summary = {
        "tier": tier_name,
        "n_cells": int(sub.n_obs),
        "n_genes": int(sub.n_vars),
        "source": str(src_zarr),
        "obs_columns": list(sub.obs.columns),
    }
    (dest.parent / f"{tier_name}.summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(f"  wrote {sub.n_obs:>7,} cells -> {h5ad_path}")


def build_small_real(
    z: zarr.Group, dest_base: Path, rng: np.random.Generator, src_zarr: Path = SRC_ZARR
) -> None:
    spec = TIER_SPECS["small_real"]
    idx = _select_cells(z, spec["n_cells"], spec["n_samples"], rng)
    _write_anndata_subset(src_zarr, idx, dest_base / "small_real", "small_real")


def build_medium_real(
    z: zarr.Group, dest_base: Path, rng: np.random.Generator, src_zarr: Path = SRC_ZARR
) -> None:
    spec = TIER_SPECS["medium_real"]
    idx = _select_cells(z, spec["n_cells"], spec["n_samples"], rng)
    _write_anndata_subset(src_zarr, idx, dest_base / "medium_real", "medium_real")


def build_full_real(dest_base: Path, src_zarr: Path = SRC_ZARR) -> None:
    path_file = dest_base / "full_real.path"
    dest_base.mkdir(parents=True, exist_ok=True)
    path_file.write_text(str(src_zarr) + "\n")
    print(f"  recorded full_real path -> {path_file}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Build staircase test fixtures from NC2024 NSCLC prepared zarr.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Tiers
  small_real  : 2 000 cells, 4 samples  -> tests/data/staircase/small_real/
  medium_real : 20 000 cells, 16 samples -> tests/data/staircase/medium_real/
  full_real   : records path to original zarr (884 050 cells, 81 samples)

After running this script, re-run:
  pytest -m small_real tests/test_staircase_smoke.py
  pytest -m medium_real tests/test_staircase_smoke.py
  pytest -m full_real tests/test_staircase_smoke.py
        """,
    )
    parser.add_argument(
        "--tiers",
        nargs="+",
        choices=list(TIER_SPECS.keys()),
        default=list(TIER_SPECS.keys()),
        help="Which tiers to build (default: all)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible cell selection (default: 42)",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=DEST_DIR,
        help=f"Output directory (default: {DEST_DIR})",
    )
    parser.add_argument(
        "--src-zarr",
        type=Path,
        default=SRC_ZARR,
        help=f"Source prepared_input.zarr (default: {SRC_ZARR})",
    )
    args = parser.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    args.dest.mkdir(parents=True, exist_ok=True)

    src_zarr = args.src_zarr

    need_src = [t for t in args.tiers if t != "full_real"]
    if need_src:
        print(f"Opening source zarr: {src_zarr}")
        z = _open_src(src_zarr)
    else:
        z = None

    for tier in args.tiers:
        print(f"Building tier: {tier}")
        if tier == "small_real":
            build_small_real(z, args.dest, rng)
        elif tier == "medium_real":
            build_medium_real(z, args.dest, rng)
        elif tier == "full_real":
            build_full_real(args.dest, src_zarr=src_zarr)

    print("Done.")


if __name__ == "__main__":
    main()
