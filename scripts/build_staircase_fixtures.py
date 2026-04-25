#!/usr/bin/env python3
"""Build staircase test fixtures from the NC2024 NSCLC prepared zarr.

Tiers
-----
nano        : 200 cells, 2 samples  — pure-synthetic, no disk dependency
small_real  : 2 000 cells, 4 samples  — from prepared_input.zarr
medium_real : 20 000 cells, 16 samples — from prepared_input.zarr
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
    "small_real": {"n_cells": 2_000, "n_samples": 4},
    "medium_real": {"n_cells": 20_000, "n_samples": 16},
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
    """Return sorted cell indices covering exactly n_samples, total ~n_cells."""
    codes = _get_sample_codes(z)
    cats = _get_sample_categories(z)
    all_samples = np.arange(len(cats))
    chosen_samples = rng.choice(all_samples, size=min(n_samples, len(cats)), replace=False)
    chosen_samples.sort()

    per_sample = max(1, n_cells // len(chosen_samples))
    chosen_cells: list[np.ndarray] = []
    for s in chosen_samples:
        mask = np.where(codes == s)[0]
        take = min(per_sample, len(mask))
        chosen_cells.append(rng.choice(mask, size=take, replace=False))

    idx = np.concatenate(chosen_cells)
    # trim to exactly n_cells if we overshot
    if len(idx) > n_cells:
        idx = rng.choice(idx, size=n_cells, replace=False)
    return np.sort(idx)


def _write_zarr_subset(
    z: zarr.Group,
    cell_idx: np.ndarray,
    dest: Path,
    compressor: zarr.codecs.Blosc | None = None,
) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    out = zarr.open_group(str(dest), mode="w")

    # --- X (CSR subset) ---
    indptr_full = z["X"]["indptr"][:]
    starts = indptr_full[cell_idx]
    ends = indptr_full[cell_idx + 1]
    slices = np.concatenate([np.arange(s, e) for s, e in zip(starts, ends)])
    new_data = z["X"]["data"][slices]
    new_indices = z["X"]["indices"][slices]
    new_indptr = np.zeros(len(cell_idx) + 1, dtype=np.int64)
    counts = ends - starts
    np.cumsum(counts, out=new_indptr[1:])

    out.require_group("X")
    out["X"].create_array("data", data=new_data, overwrite=True)
    out["X"].create_array("indices", data=new_indices, overwrite=True)
    out["X"].create_array("indptr", data=new_indptr, overwrite=True)

    # --- obs ---
    out.require_group("obs")
    for key in z["obs"].keys():
        src_item = z["obs"][key]
        if isinstance(src_item, zarr.Group):
            # categorical
            codes_full = src_item["codes"][:]
            cats_full = src_item["categories"][:]
            new_codes = codes_full[cell_idx]
            used = np.unique(new_codes)
            code_remap = np.full(len(cats_full), -1, dtype=np.int8)
            for new_code, old_code in enumerate(used):
                code_remap[old_code] = new_code
            remapped = code_remap[new_codes]
            new_cats = cats_full[used]
            out["obs"].require_group(key)
            out["obs"][key].create_array("codes", data=remapped, overwrite=True)
            out["obs"][key].create_array("categories", data=new_cats, overwrite=True)
        else:
            # plain array
            out["obs"].create_array(key, data=src_item[cell_idx], overwrite=True)

    # --- var (unchanged) ---
    out.require_group("var")
    for key in z["var"].keys():
        src_item = z["var"][key]
        if isinstance(src_item, zarr.Group):
            out["var"].require_group(key)
            for sub in src_item.keys():
                out["var"][key].create_array(sub, data=src_item[sub][:], overwrite=True)
        else:
            out["var"].create_array(key, data=src_item[:], overwrite=True)

    # metadata
    meta = {
        "n_cells": int(len(cell_idx)),
        "n_genes": int(z["X"]["indices"][:].max()) + 1 if len(z["X"]["indices"]) > 0 else 0,
        "source": str(SRC_ZARR),
    }
    (dest / "staircase_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  wrote {len(cell_idx):>7,} cells -> {dest}")


def build_small_real(z: zarr.Group, dest_base: Path, rng: np.random.Generator) -> None:
    spec = TIER_SPECS["small_real"]
    idx = _select_cells(z, spec["n_cells"], spec["n_samples"], rng)
    _write_zarr_subset(z, idx, dest_base / "small_real")


def build_medium_real(z: zarr.Group, dest_base: Path, rng: np.random.Generator) -> None:
    spec = TIER_SPECS["medium_real"]
    idx = _select_cells(z, spec["n_cells"], spec["n_samples"], rng)
    _write_zarr_subset(z, idx, dest_base / "medium_real")


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
