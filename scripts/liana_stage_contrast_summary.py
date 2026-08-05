#!/usr/bin/env python3
"""Summarize an already-produced per-stage LIANA directory tree.

HCI helper: given ``stage_I/tables/cell_communication_liana.csv`` etc., write
Jaccard contrasts for VEGF/EGF and checkpoint panels using
``filter_liana_claim_panels.filter_panels``.

Does not re-run LIANA (compute-heavy). For full generation see audit Wave-5
run ``2026-08-05T0726Z-9b69735``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

# Allow running as script from repo root
_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.filter_liana_claim_panels import filter_panels  # noqa: E402


def pair_key(df: pd.DataFrame) -> set[str]:
    return set(
        df["ligand_complex"].astype(str)
        + "|"
        + df["receptor_complex"].astype(str)
        + "|"
        + df["source"].astype(str)
        + "->"
        + df["target"].astype(str)
    )


def jaccard(a: set[str], b: set[str]) -> float | None:
    if not a and not b:
        return None
    u = a | b
    return len(a & b) / len(u) if u else None


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--stage-dir",
        action="append",
        required=True,
        metavar="STAGE=PATH",
        help="e.g. I=/path/to/stage_I/tables/cell_communication_liana.csv",
    )
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args(argv)

    stages: dict[str, pd.DataFrame] = {}
    for item in args.stage_dir:
        if "=" not in item:
            raise SystemExit(f"bad --stage-dir {item!r}; want STAGE=PATH")
        st, path = item.split("=", 1)
        stages[st] = pd.read_csv(path)

    ckpt_sets = {}
    vegf_sets = {}
    meta = {}
    for st, df in stages.items():
        panels = filter_panels(df)
        ckpt_sets[st] = pair_key(panels["checkpoint_immune"])
        vegf_sets[st] = pair_key(panels["vegf_egfr"])
        meta[st] = {
            "n_liana": int(len(df)),
            "n_checkpoint_immune": int(len(panels["checkpoint_immune"])),
            "n_vegf_egfr": int(len(panels["vegf_egfr"])),
            "n_classic_vegfr": int(len(panels["classic_vegf_vegfr"])),
        }

    keys = list(stages)
    out = {"stages": meta, "checkpoint_pair_jaccard": {}, "vegf_pair_jaccard": {}}
    for i, s1 in enumerate(keys):
        for s2 in keys[i + 1 :]:
            out["checkpoint_pair_jaccard"][f"{s1}_vs_{s2}"] = {
                "jaccard": jaccard(ckpt_sets[s1], ckpt_sets[s2]),
                "n1": len(ckpt_sets[s1]),
                "n2": len(ckpt_sets[s2]),
                "intersection": len(ckpt_sets[s1] & ckpt_sets[s2]),
            }
            out["vegf_pair_jaccard"][f"{s1}_vs_{s2}"] = {
                "jaccard": jaccard(vegf_sets[s1], vegf_sets[s2]),
                "n1": len(vegf_sets[s1]),
                "n2": len(vegf_sets[s2]),
                "intersection": len(vegf_sets[s1] & vegf_sets[s2]),
            }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2)[:2000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
