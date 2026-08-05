#!/usr/bin/env python3
"""LUCA patient-stratified LIANA expansion driver (Wave-9).

Expands the Wave-8 patient-level LIANA lane from 3+3 to N+N donors per
histology on the LUCA core atlas (Salcher 2022), reusing the identical
engine (``li.mt.rank_aggregate`` consensus) and output contract
(``cell_communication_liana.csv`` per donor + ``liana_patient_recurrence.py``
summary).

Science / scope
---------------
- Descriptive donor recurrence of panel LR pairs across more patients;
  still **not** the De Zuani paper's CellphoneDB full-cohort multi-condition
  design and not causal inference. Claim class stays ``partial``.
- Donor eligibility is recorded explicitly (min cells / min cell types);
  low-diversity donors are excluded with a recorded reason, never silently
  included (Wave-8 n_types=2 case).

Usage
-----
  python scripts/liana_luca_patient_expansion.py \\
    --atlas /home/zerlinshen/data/raw/luca_core_salcher2022/luca_core_atlas_892k.h5ad \\
    --output-dir <run_dir>/python/patient_liana \\
    --n-per-histology 8 --cell-cap 1000 --seed 47
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import zlib
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

HISTOLOGIES = {
    "LUAD": "lung adenocarcinoma",
    "LUSC": "squamous cell lung carcinoma",
}
CT_COL = "cell_type_predicted"
MIN_CELLS_PER_TYPE = 5
MIN_TYPES_PER_DONOR = 6
MIN_DONOR_CELLS = 800
LIANA_REQUIRED_COLUMNS = {
    "ligand_complex",
    "receptor_complex",
    "source",
    "target",
    "magnitude_rank",
    "specificity_rank",
}


def _decode_strings(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            value.decode("utf-8")
            if isinstance(value, (bytes, np.bytes_))
            else str(value)
            for value in np.asarray(values)
        ],
        dtype=object,
    )


def _read_h5_column(frame: h5py.Group, column: str) -> np.ndarray:
    """Read one AnnData dataframe column without loading unrelated slots."""
    if column not in frame:
        raise KeyError(f"missing H5AD column: {frame.name}/{column}")
    node = frame[column]
    if isinstance(node, h5py.Group):
        encoding = node.attrs.get("encoding-type")
        if isinstance(encoding, bytes):
            encoding = encoding.decode("utf-8")
        if encoding != "categorical":
            raise TypeError(f"unsupported H5AD column encoding: {node.name}={encoding}")
        categories = _decode_strings(node["categories"][:])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        values = np.empty(len(codes), dtype=object)
        values[:] = None
        valid = codes >= 0
        values[valid] = categories[codes[valid]]
        return values
    values = node[:]
    if values.dtype.kind in {"O", "S", "U"}:
        return _decode_strings(values)
    return np.asarray(values)


def read_obs_columns(atlas_path: Path, columns: list[str]) -> pd.DataFrame:
    """Read only the requested observation columns from a large H5AD."""
    with h5py.File(atlas_path, "r") as handle:
        obs = handle["obs"]
        return pd.DataFrame({column: _read_h5_column(obs, column) for column in columns})


def load_expression_subset(atlas_path: Path, rows: np.ndarray):
    """Load selected ``X`` CSR rows while ignoring large atlas graphs/embeddings."""
    import anndata as ad

    rows = np.asarray(rows, dtype=np.int64)
    if rows.ndim != 1 or len(rows) == 0:
        raise ValueError("rows must be a non-empty one-dimensional array")
    if np.any(np.diff(rows) <= 0):
        raise ValueError("rows must be strictly increasing")

    with h5py.File(atlas_path, "r") as handle:
        matrix = handle["X"]
        encoding = matrix.attrs.get("encoding-type")
        if isinstance(encoding, bytes):
            encoding = encoding.decode("utf-8")
        if encoding != "csr_matrix":
            raise TypeError(f"expected CSR-encoded X, found {encoding!r}")

        shape = tuple(int(value) for value in matrix.attrs["shape"])
        if rows[-1] >= shape[0]:
            raise IndexError(f"row {rows[-1]} is outside X with {shape[0]} rows")
        starts = np.asarray(matrix["indptr"][rows], dtype=np.int64)
        stops = np.asarray(matrix["indptr"][rows + 1], dtype=np.int64)
        lengths = stops - starts
        subset_indptr = np.empty(len(rows) + 1, dtype=np.int64)
        subset_indptr[0] = 0
        np.cumsum(lengths, out=subset_indptr[1:])
        subset_data = np.empty(subset_indptr[-1], dtype=matrix["data"].dtype)
        subset_indices = np.empty(subset_indptr[-1], dtype=matrix["indices"].dtype)
        cursor = 0
        for start, stop in zip(starts, stops, strict=True):
            count = int(stop - start)
            subset_data[cursor : cursor + count] = matrix["data"][start:stop]
            subset_indices[cursor : cursor + count] = matrix["indices"][start:stop]
            cursor += count

        feature_names = _read_h5_column(handle["var"], "feature_name")

    expression = sparse.csr_matrix(
        (subset_data, subset_indices, subset_indptr),
        shape=(len(rows), shape[1]),
    )
    obs = pd.DataFrame(index=pd.Index([f"atlas_row_{row}" for row in rows]))
    var = pd.DataFrame(index=pd.Index(feature_names.astype(str)))
    return ad.AnnData(X=expression, obs=obs, var=var)


def donor_table(obs: pd.DataFrame) -> pd.DataFrame:
    """Per-donor tumor_primary cell counts and post-filter type counts."""
    t = obs[obs["origin"].astype(str) == "tumor_primary"].copy()
    t["donor_id"] = t["donor_id"].astype(str)
    t["disease"] = t["disease"].astype(str)
    rows = []
    for (donor, disease), g in t.groupby(["donor_id", "disease"]):
        ct = g[CT_COL].astype(str).value_counts()
        rows.append(
            {
                "donor": donor,
                "disease": disease,
                "n_cells": int(len(g)),
                "n_types_ge5": int((ct >= MIN_CELLS_PER_TYPE).sum()),
            }
        )
    return pd.DataFrame(rows)


def select_donors(dt: pd.DataFrame, histology: str, n: int) -> tuple[list[str], list[dict], list[str]]:
    """Top-n donors by tumor cell count passing eligibility.

    Returns (selected, ineligible, eligible_not_selected): full accounting —
    every same-histology donor lands in exactly one of the three buckets.
    """
    disease = HISTOLOGIES[histology]
    sub = dt[dt["disease"] == disease].sort_values("n_cells", ascending=False)
    eligible_mask = (sub["n_cells"] >= MIN_DONOR_CELLS) & (
        sub["n_types_ge5"] >= MIN_TYPES_PER_DONOR
    )
    eligible = sub[eligible_mask]
    selected = eligible["donor"].head(n).tolist()
    selected_set = set(selected)
    eligible_not_selected = [
        d for d in eligible["donor"].tolist() if d not in selected_set
    ]
    ineligible = [
        {
            "donor": r["donor"],
            "n_cells": int(r["n_cells"]),
            "n_types_ge5": int(r["n_types_ge5"]),
            "reason": (
                f"below threshold: needs >= {MIN_DONOR_CELLS} cells and "
                f">= {MIN_TYPES_PER_DONOR} types (>= {MIN_CELLS_PER_TYPE} cells/type)"
            ),
        }
        for _, r in sub[~eligible_mask].iterrows()
    ]
    return selected, ineligible, eligible_not_selected


def stratified_cap(idx: np.ndarray, types: pd.Series, cap: int, seed_key: str) -> np.ndarray:
    """Proportional per-type downsample to ``cap`` cells, min 1 per type."""
    if len(idx) <= cap:
        return idx
    rng = np.random.default_rng(zlib.crc32(seed_key.encode()) & 0xFFFFFFFF)
    counts = types.value_counts()
    alloc = (counts / counts.sum() * cap).round().astype(int).clip(lower=1)
    # fix rounding drift deterministically on the largest types
    drift = int(alloc.sum() - cap)
    if drift != 0:
        order = counts.sort_values(ascending=False).index
        i = 0
        while drift != 0:
            t = order[i % len(order)]
            step = -1 if drift > 0 else 1
            if alloc[t] + step >= 1:
                alloc[t] += step
                drift += step
            i += 1
    keep: list[np.ndarray] = []
    type_arr = types.to_numpy()
    for t, k in alloc.items():
        pool = idx[type_arr == t]
        keep.append(rng.choice(pool, size=min(int(k), len(pool)), replace=False))
    return np.sort(np.concatenate(keep))


def run_donor(
    atlas_path: Path,
    donor: str,
    out_dir: Path,
    cell_cap: int,
    seed: int,
    obs: pd.DataFrame | None = None,
) -> dict:
    import scanpy as sc

    if obs is None:
        obs = read_obs_columns(
            atlas_path, ["donor_id", "origin", "disease", CT_COL]
        )
    mask = (
        (obs["donor_id"].astype(str) == donor)
        & (obs["origin"].astype(str) == "tumor_primary")
    ).to_numpy()
    idx_all = np.where(mask)[0]
    types_all = obs[CT_COL].astype(str)[mask]
    # >=5 cells/type filter (applied before cap; drives eligibility)
    ct = types_all.value_counts()
    keep_types = set(ct[ct >= MIN_CELLS_PER_TYPE].index)
    keep_mask = types_all.isin(keep_types).to_numpy()
    idx = idx_all[keep_mask]
    types = types_all[keep_mask]
    idx = stratified_cap(idx, types, cell_cap, f"{seed}:{donor}")

    sub = load_expression_subset(atlas_path, idx)
    sub.obs["cell_type"] = obs.iloc[idx][CT_COL].astype(str).to_numpy()

    # Symbol namespace for LIANA (Wave-6/7 HCI contract).
    sub = sub[:, ~sub.var_names.isin(["nan", "None", ""])]
    sub.var_names = pd.Index(sub.var_names.astype(str))
    sub.var_names_make_unique()

    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    import liana as li

    li.mt.rank_aggregate(
        sub,
        groupby="cell_type",
        resource_name="consensus",
        use_raw=False,
        verbose=False,
    )
    res = sub.uns.get("liana_res")
    if res is None or len(res) == 0:
        raise RuntimeError(f"empty LIANA results for {donor}")

    table_dir = out_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    csv_path = table_dir / "cell_communication_liana.csv"
    res.to_csv(csv_path, index=False)
    return {
        "n_cells": int(sub.n_obs),
        "n_types": int(sub.obs["cell_type"].nunique()),
        "n_liana": int(len(res)),
        "csv": str(csv_path),
        "peak_rss_mb_process": round(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1
        ),
        "atlas_io": "direct_h5py_csr_rows_only",
    }


def load_existing_donor_stats(out_dir: Path) -> dict | None:
    summary_path = out_dir / "donor_summary.json"
    csv_path = out_dir / "tables" / "cell_communication_liana.csv"
    if not summary_path.exists() or not csv_path.exists():
        return None
    stats = json.loads(summary_path.read_text(encoding="utf-8"))
    if stats.get("csv") != str(csv_path) or int(stats.get("n_liana", 0)) <= 0:
        return None
    try:
        table = pd.read_csv(csv_path)
    except Exception:
        return None
    if len(table) != int(stats["n_liana"]):
        return None
    if not LIANA_REQUIRED_COLUMNS.issubset(table.columns):
        return None
    stats["resumed"] = True
    return stats


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--atlas", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--n-per-histology", type=int, default=8)
    p.add_argument("--cell-cap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=47)
    p.add_argument(
        "--resume",
        action="store_true",
        help="reuse donors with a valid CSV and donor_summary.json checkpoint",
    )
    args = p.parse_args(argv)

    obs = read_obs_columns(
        args.atlas, ["donor_id", "origin", "disease", CT_COL]
    )

    dt = donor_table(obs)
    meta: dict = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "atlas": str(args.atlas),
        "design": {
            "origin": "tumor_primary",
            "histology_source": "obs.disease",
            "ct_col": CT_COL,
            "min_cells_per_type": MIN_CELLS_PER_TYPE,
            "min_types_per_donor": MIN_TYPES_PER_DONOR,
            "min_donor_cells": MIN_DONOR_CELLS,
            "n_per_histology": args.n_per_histology,
            "cell_cap": args.cell_cap,
            "sampling": "proportional per-type stratified cap, crc32-seeded per donor",
            "seed": args.seed,
            "engine": "liana_consensus (li.mt.rank_aggregate)",
            "gene_namespace": "feature_name symbols",
            "atlas_io": (
                "direct HDF5 obs-column + CSR X-row reads; atlas obsp/obsm/layers "
                "are not materialized"
            ),
            "eligible_donor_pool": {
                h: int(
                    (
                        (dt["disease"] == d)
                        & (dt["n_cells"] >= MIN_DONOR_CELLS)
                        & (dt["n_types_ge5"] >= MIN_TYPES_PER_DONOR)
                    ).sum()
                )
                for h, d in HISTOLOGIES.items()
            },
        },
        "donors": {},
        "excluded_donors": {},
        "errors": [],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    meta_path = args.output_dir / "patient_liana_meta.json"
    donor_args: list[str] = []
    for hist in HISTOLOGIES:
        selected, ineligible, eligible_not_selected = select_donors(
            dt, hist, args.n_per_histology
        )
        meta["excluded_donors"][hist] = ineligible
        meta["design"].setdefault("eligible_not_selected", {})[hist] = (
            eligible_not_selected
        )
        for donor in selected:
            out_dir = args.output_dir / hist / donor
            print(f"=== {hist}:{donor} ===", flush=True)
            try:
                stats = load_existing_donor_stats(out_dir) if args.resume else None
                if stats is None:
                    stats = run_donor(
                        args.atlas,
                        donor,
                        out_dir,
                        args.cell_cap,
                        args.seed,
                        obs=obs,
                    )
                    (out_dir / "donor_summary.json").write_text(
                        json.dumps(stats, indent=2) + "\n", encoding="utf-8"
                    )
            except Exception as exc:
                meta["errors"].append({"histology": hist, "donor": donor, "error": str(exc)})
                print(f"ERROR {donor}: {exc}", file=sys.stderr, flush=True)
                meta_path.write_text(
                    json.dumps(meta, indent=2) + "\n", encoding="utf-8"
                )
                continue
            meta["donors"][f"{hist}:{donor}"] = stats
            donor_args.append(f"{hist}:{donor}={stats['csv']}")
            print(json.dumps(stats, indent=None), flush=True)
            meta_path.write_text(
                json.dumps(meta, indent=2) + "\n", encoding="utf-8"
            )

    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {meta_path}")

    if donor_args:
        from scripts.liana_patient_recurrence import main as recurrence_main

        rc = recurrence_main(
            [x for pair in donor_args for x in ("--donor", pair)]
            + ["--output-dir", str(args.output_dir / "summary")]
        )
        if rc != 0:
            return rc

    n_ok = len(donor_args)
    print(
        json.dumps(
            {
                "donors_ok": n_ok,
                "errors": len(meta["errors"]),
                "per_histology": {
                    h: sum(1 for d in donor_args if d.startswith(h + ":"))
                    for h in HISTOLOGIES
                },
            },
            indent=2,
        )
    )
    if n_ok < 2:
        print("FAIL: fewer than 2 donors completed", file=sys.stderr)
        return 3
    print("LUCA PATIENT LIANA EXPANSION COMPLETE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
