#!/usr/bin/env python3
"""Filter a LIANA consensus table into literature-informed claim panels.

Science / scope
---------------
This is a **post-hoc panel filter**, not a new inference engine. It assumes the
input CSV was produced by ``cell_communication`` with engine ``liana_consensus``
(Dimitrov et al., Nat Commun 2022 — multi-method consensus). Panels are gene-
symbol / complex substring matches against commonly cited VEGF/EGF-axis and
immune-checkpoint families.

Claim honesty
-------------
- Non-empty VEGF/EGF or checkpoint filters support **partial** method-lane
  evidence that those families appear in the LR table.
- They do **not** by themselves prove multi-condition paper figures, classic
  VEGF–KDR dominance, or LUAD-vs-LUSC checkpoint wiring without a designed
  contrast.

Usage
-----
  python scripts/filter_liana_claim_panels.py \\
    --input path/to/cell_communication_liana.csv \\
    --output-dir path/to/out/
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

# VEGF / EGF axis (ligands + primary / co-receptors often present in LR resources)
VEGF_EGFR_TERMS: tuple[str, ...] = (
    "VEGFA",
    "VEGFB",
    "VEGFC",
    "VEGFD",
    "PGF",
    "EGFR",
    "ERBB2",
    "ERBB3",
    "ERBB4",
    "KDR",
    "FLT1",
    "FLT4",
    "NRP1",
    "NRP2",
    "ITGAV",
    "ITGB3",
)

# Immune checkpoint / co-stimulatory families (PD-1, CTLA-4, TIGIT, TIM-3, LAG-3, …)
CHECKPOINT_TERMS: tuple[str, ...] = (
    "PDCD1",
    "CD274",
    "PDCD1LG2",
    "CD80",
    "CD86",
    "CTLA4",
    "TIGIT",
    "PVR",
    "PVRL2",
    "NECTIN2",
    "NECTIN3",
    "HAVCR2",
    "LGALS9",
    "LAG3",
    "FGL1",
    "CD276",
    "VTCN1",
    "HHLA2",
    "CD47",
    "SIRPA",
    "TNFRSF9",
    "TNFSF9",
    "TNFRSF4",
    "TNFSF4",
    "TNFRSF18",
    "TNFSF18",
    "ICOS",
    "ICOSLG",
    "CD28",
    "CD27",
    "CD70",
)

IMMUNE_LABELS = frozenset(
    {
        "T cell",
        "NK cell",
        "B cell",
        "Plasma cell",
        "Myeloid/Macro",
        "Dendritic cell",
        "Mast cell",
    }
)


def _match_panel(df: pd.DataFrame, terms: tuple[str, ...]) -> pd.DataFrame:
    lig = df["ligand_complex"].astype(str).str.upper()
    rec = df["receptor_complex"].astype(str).str.upper()
    mask = False
    for t in terms:
        mask = mask | lig.str.contains(t, regex=False) | rec.str.contains(t, regex=False)
    return df.loc[mask].copy()


def filter_panels(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    required = {"ligand_complex", "receptor_complex", "source", "target"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"LIANA table missing columns: {sorted(missing)}")

    vegf = _match_panel(df, VEGF_EGFR_TERMS)
    classic = vegf[
        vegf["ligand_complex"].astype(str).str.upper().str.contains("VEGF", regex=False)
        & vegf["receptor_complex"]
        .astype(str)
        .str.upper()
        .str.contains("KDR|FLT1|FLT4", regex=True)
    ]
    vegfa_egfr = vegf[
        vegf["ligand_complex"].astype(str).str.upper().str.contains("VEGFA", regex=False)
        & vegf["receptor_complex"].astype(str).str.upper().str.contains("EGFR", regex=False)
    ]
    ckpt = _match_panel(df, CHECKPOINT_TERMS)
    ckpt_immune = ckpt[ckpt["source"].isin(IMMUNE_LABELS) | ckpt["target"].isin(IMMUNE_LABELS)]
    return {
        "vegf_egfr": vegf,
        "classic_vegf_vegfr": classic,
        "vegfa_egfr": vegfa_egfr,
        "checkpoint_any": ckpt,
        "checkpoint_immune": ckpt_immune,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path, help="cell_communication_liana.csv")
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args(argv)

    df = pd.read_csv(args.input)
    panels = filter_panels(df)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    vegf_path = args.output_dir / "communication_filter_vegf_egfr.csv"
    ckpt_path = args.output_dir / "communication_filter_checkpoint.csv"
    panels["vegf_egfr"].to_csv(vegf_path, index=False)
    panels["checkpoint_immune"].to_csv(ckpt_path, index=False)

    def _top(frame: pd.DataFrame, n: int = 15) -> list[dict]:
        if frame.empty or "magnitude_rank" not in frame.columns:
            return []
        cols = [
            c
            for c in (
                "source",
                "target",
                "ligand_complex",
                "receptor_complex",
                "magnitude_rank",
                "specificity_rank",
            )
            if c in frame.columns
        ]
        return frame.nsmallest(n, "magnitude_rank")[cols].to_dict(orient="records")

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_liana_csv": str(args.input.resolve()),
        "n_liana_total": int(len(df)),
        "vegf_egfr": {
            "n_pairs": int(len(panels["vegf_egfr"])),
            "n_classic_vegf_vegfr_kdr_flt": int(len(panels["classic_vegf_vegfr"])),
            "n_vegfa_egfr": int(len(panels["vegfa_egfr"])),
            "top_by_magnitude": _top(panels["vegf_egfr"]),
            "artifact": str(vegf_path),
            "scientific_note": (
                "Panel filter on LIANA consensus. Classic VEGF–KDR/FLT1 may be "
                "absent; alternative pairs (e.g. VEGFA–EGFR, VEGFB–NRP1) are "
                "reported when present."
            ),
        },
        "checkpoint": {
            "n_pairs_any": int(len(panels["checkpoint_any"])),
            "n_pairs_immune_endpoint": int(len(panels["checkpoint_immune"])),
            "top_by_magnitude": _top(panels["checkpoint_immune"]),
            "artifact": str(ckpt_path),
            "scientific_note": (
                "Immune-checkpoint family filter with immune cell endpoint. "
                "Not a LUAD-vs-LUSC differential without contrast design."
            ),
        },
        "claim_recommendations": {
            "NC2024-F2-03": (
                "partial" if len(panels["vegf_egfr"]) > 0 else "unsupported"
            ),
            "NC2024-F2-04": (
                "partial"
                if len(panels["checkpoint_immune"]) > 0
                else "unsupported"
            ),
        },
    }
    sum_path = args.output_dir / "communication_filter_summary.json"
    sum_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["claim_recommendations"] | {
        "vegf_n": summary["vegf_egfr"]["n_pairs"],
        "ckpt_immune_n": summary["checkpoint"]["n_pairs_immune_endpoint"],
        "summary": str(sum_path),
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
