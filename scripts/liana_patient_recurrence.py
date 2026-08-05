#!/usr/bin/env python3
"""Patient-level LIANA recurrence summary (Wave-8).

Given one LIANA consensus CSV per donor (or donor=path pairs), compute
checkpoint / VEGF panel **patient recurrence** — fraction of donors where a
ligand|receptor (or full edge) appears.

Science / scope
---------------
- Descriptive recurrence across donors, not a hierarchical multi-condition
  paper model and not causal network inference.
- Assumes each CSV was produced by LIANA consensus (``liana_consensus``).

Usage
-----
  python scripts/liana_patient_recurrence.py \\
    --donor LUAD:D1=/path/to/d1/liana.csv \\
    --donor LUSC:D2=/path/to/d2/liana.csv \\
    --output-dir /path/to/out/
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.filter_liana_claim_panels import filter_panels  # noqa: E402


def lr_keys(df: pd.DataFrame) -> set[str]:
    return set(
        df["ligand_complex"].astype(str) + "|" + df["receptor_complex"].astype(str)
    )


def edge_keys(df: pd.DataFrame) -> set[str]:
    return set(
        df["ligand_complex"].astype(str)
        + "|"
        + df["receptor_complex"].astype(str)
        + "|"
        + df["source"].astype(str)
        + "->"
        + df["target"].astype(str)
    )


def recurrence(donor_sets: dict[str, set[str]]) -> list[dict]:
    """Rank LR keys by fraction of donors containing them."""
    if not donor_sets:
        return []
    n = len(donor_sets)
    counts: dict[str, int] = defaultdict(int)
    for s in donor_sets.values():
        for k in s:
            counts[k] += 1
    rows = [
        {
            "key": k,
            "n_donors": c,
            "n_donors_total": n,
            "fraction": c / n,
        }
        for k, c in counts.items()
    ]
    rows.sort(key=lambda r: (-r["fraction"], -r["n_donors"], r["key"]))
    return rows


def summarize_histology(
    donors: dict[str, pd.DataFrame],
) -> dict:
    meta = {}
    ckpt_lr: dict[str, set[str]] = {}
    vegf_lr: dict[str, set[str]] = {}
    classic_lr: dict[str, set[str]] = {}
    ckpt_edge: dict[str, set[str]] = {}

    for donor, df in donors.items():
        panels = filter_panels(df)
        meta[donor] = {
            "n_liana": int(len(df)),
            "n_checkpoint_any": int(len(panels["checkpoint_any"])),
            "n_checkpoint_immune": int(len(panels["checkpoint_immune"])),
            "n_vegf_egfr": int(len(panels["vegf_egfr"])),
            "n_classic_vegfr": int(len(panels["classic_vegf_vegfr"])),
        }
        ckpt_lr[donor] = lr_keys(panels["checkpoint_any"])
        vegf_lr[donor] = lr_keys(panels["vegf_egfr"])
        classic_lr[donor] = lr_keys(panels["classic_vegf_vegfr"])
        ckpt_edge[donor] = edge_keys(panels["checkpoint_any"])

    # Flagship pairs of interest for NSCLC communication narratives
    flagship = [
        "CD274|PDCD1",
        "PDCD1LG2|PDCD1",
        "CD80|CTLA4",
        "CD86|CTLA4",
        "VEGFA|KDR",
        "VEGFA|FLT1",
        "VEGFA|EGFR",
        "LGALS9|HAVCR2",
    ]
    n = len(donors) or 1
    flagship_rec = {}
    for pair in flagship:
        present = [d for d, s in ckpt_lr.items() if pair in s]
        present_v = [d for d, s in vegf_lr.items() if pair in s]
        present_c = [d for d, s in classic_lr.items() if pair in s]
        # pair may be checkpoint or VEGF
        donors_hit = sorted(set(present) | set(present_v) | set(present_c))
        flagship_rec[pair] = {
            "n_donors": len(donors_hit),
            "fraction": len(donors_hit) / n,
            "donors": donors_hit,
        }

    return {
        "n_donors": len(donors),
        "per_donor": meta,
        "checkpoint_lr_recurrence_top": recurrence(ckpt_lr)[:40],
        "vegf_lr_recurrence_top": recurrence(vegf_lr)[:40],
        "classic_vegfr_lr_recurrence_top": recurrence(classic_lr)[:20],
        "checkpoint_edge_recurrence_top": recurrence(ckpt_edge)[:30],
        "flagship_pair_recurrence": flagship_rec,
        "mean_n_checkpoint_any": float(
            sum(m["n_checkpoint_any"] for m in meta.values()) / n
        ),
        "mean_n_classic_vegfr": float(
            sum(m["n_classic_vegfr"] for m in meta.values()) / n
        ),
        "donors_with_any_classic_vegfr": sum(
            1 for m in meta.values() if m["n_classic_vegfr"] > 0
        ),
        "donors_with_cd274_pdcd1": flagship_rec["CD274|PDCD1"]["n_donors"],
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--donor",
        action="append",
        required=True,
        metavar="HIST:DONOR=PATH",
        help="e.g. LUAD:He_Fan_2021_P1=/path/to/liana.csv",
    )
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args(argv)

    by_hist: dict[str, dict[str, pd.DataFrame]] = defaultdict(dict)
    for item in args.donor:
        if "=" not in item or ":" not in item.split("=", 1)[0]:
            raise SystemExit(f"bad --donor {item!r}; want HIST:DONOR=PATH")
        left, path = item.split("=", 1)
        hist, donor = left.split(":", 1)
        by_hist[hist][donor] = pd.read_csv(path)

    out: dict = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "patient_level_LIANA_panel_recurrence",
        "scientific_boundaries": [
            "Descriptive donor recurrence of LIANA panel pairs",
            "Not multi-condition paper figure parity",
            "Not causal inference",
        ],
        "histologies": {},
    }
    for hist, donors in by_hist.items():
        out["histologies"][hist] = summarize_histology(donors)

    # Cross-histology: flagship prevalence contrast
    if "LUAD" in out["histologies"] and "LUSC" in out["histologies"]:
        luad_f = out["histologies"]["LUAD"]["flagship_pair_recurrence"]
        lusc_f = out["histologies"]["LUSC"]["flagship_pair_recurrence"]
        contrast = {}
        for pair in sorted(set(luad_f) | set(lusc_f)):
            contrast[pair] = {
                "luad_fraction": luad_f.get(pair, {}).get("fraction"),
                "lusc_fraction": lusc_f.get(pair, {}).get("fraction"),
                "luad_n": luad_f.get(pair, {}).get("n_donors"),
                "lusc_n": lusc_f.get(pair, {}).get("n_donors"),
            }
        out["luad_vs_lusc_flagship_fraction"] = contrast

    args.output_dir.mkdir(parents=True, exist_ok=True)
    path = args.output_dir / "patient_liana_recurrence_summary.json"
    path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2)[:4000])
    print(f"WROTE {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
