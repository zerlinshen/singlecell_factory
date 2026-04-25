from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


FOCUS_PAIRS: tuple[dict[str, str], ...] = (
    {
        "pair_id": "TIM3_axis",
        "ligand": "LGALS9",
        "receptor": "HAVCR2",
        "paper_note": "Paper highlights LGALS9-HAVCR2 (TIM3) among subtype-specific checkpoint interactions.",
    },
    {
        "pair_id": "DNAM1_axis",
        "ligand": "NECTIN2",
        "receptor": "CD226",
        "paper_note": "Paper highlights NECTIN2-CD226 (DNAM1) in the LUAD/LUSC checkpoint discussion.",
    },
    {
        "pair_id": "TIGIT_axis_nectin2",
        "ligand": "NECTIN2",
        "receptor": "TIGIT",
        "paper_note": "Paper highlights NECTIN2-TIGIT among subtype-specific checkpoint interactions.",
    },
    {
        "pair_id": "TIGIT_axis_nectin3",
        "ligand": "NECTIN3",
        "receptor": "TIGIT",
        "paper_note": "Paper highlights NECTIN3-TIGIT among subtype-specific checkpoint interactions.",
    },
    {
        "pair_id": "CD96_axis",
        "ligand": "CD96",
        "receptor": "NECTIN1",
        "paper_note": "Paper highlights CD96-NECTIN1 as preferentially found in LUSC.",
    },
    {
        "pair_id": "CTLA4_axis_cd80",
        "ligand": "CD80",
        "receptor": "CTLA4",
        "paper_note": "Paper states CD80-CTLA4 is found in both tumour subtypes.",
    },
    {
        "pair_id": "CTLA4_axis_cd86",
        "ligand": "CD86",
        "receptor": "CTLA4",
        "paper_note": "Paper states CD86-CTLA4 is found in both tumour subtypes.",
    },
    {
        "pair_id": "LILRB1_axis",
        "ligand": "HLA-F",
        "receptor": "LILRB1",
        "paper_note": "Paper states HLA-F-LILRB1 is found in both tumour subtypes.",
    },
    {
        "pair_id": "LILRB2_axis",
        "ligand": "HLA-F",
        "receptor": "LILRB2",
        "paper_note": "Paper states HLA-F-LILRB2 is found in both tumour subtypes.",
    },
)


def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _ensure_cols(df: pd.DataFrame, required: tuple[str, ...], *, name: str) -> None:
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{name} missing required columns: {missing}")


def _best_row(subset: pd.DataFrame) -> pd.Series:
    return subset.sort_values(
        ["magnitude_rank", "specificity_rank", "lrscore"],
        ascending=[True, True, False],
    ).iloc[0]


def _surface_status(present_in_raw: bool, present_in_current_top20: bool) -> tuple[str, str]:
    if present_in_current_top20:
        return "visible_in_current_top20", "supported_in_surface"
    if present_in_raw:
        return "raw_only", "supported_raw_only"
    return "missing_in_raw", "unsupported_in_subtype"


def run_audit(
    *,
    luad_raw_csv: Path,
    lusc_raw_csv: Path,
    luad_top20_csv: Path,
    lusc_top20_csv: Path,
    output_dir: Path,
    project: str | None = None,
) -> Path:
    out_dir = output_dir / project if project else output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    subtype_meta = {
        "luad": {
            "disease": "lung adenocarcinoma",
            "raw_csv": Path(luad_raw_csv),
            "top20_csv": Path(luad_top20_csv),
        },
        "lusc": {
            "disease": "lung squamous cell carcinoma",
            "raw_csv": Path(lusc_raw_csv),
            "top20_csv": Path(lusc_top20_csv),
        },
    }

    all_rows: list[pd.DataFrame] = []
    best_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    visibility_rows: list[dict[str, object]] = []

    for subtype, meta in subtype_meta.items():
        raw_df = pd.read_csv(meta["raw_csv"])
        top20_df = pd.read_csv(meta["top20_csv"])
        _ensure_cols(
            raw_df,
            (
                "source",
                "target",
                "ligand_complex",
                "receptor_complex",
                "lrscore",
                "specificity_rank",
                "magnitude_rank",
                "cellphone_pvals",
            ),
            name=f"{subtype} raw csv",
        )
        _ensure_cols(
            top20_df,
            ("ligand", "receptor"),
            name=f"{subtype} top20 csv",
        )
        top20_pairs = set(zip(top20_df["ligand"], top20_df["receptor"]))

        for pair in FOCUS_PAIRS:
            subset = raw_df[
                (raw_df["ligand_complex"] == pair["ligand"])
                & (raw_df["receptor_complex"] == pair["receptor"])
            ].copy()
            present_in_top20 = (pair["ligand"], pair["receptor"]) in top20_pairs

            visibility_rows.append(
                {
                    "subtype": subtype,
                    "disease": meta["disease"],
                    "pair_id": pair["pair_id"],
                    "ligand": pair["ligand"],
                    "receptor": pair["receptor"],
                    "present_in_raw": not subset.empty,
                    "present_in_current_top20": present_in_top20,
                }
            )

            if subset.empty:
                summary_rows.append(
                    {
                        "subtype": subtype,
                        "disease": meta["disease"],
                        "pair_id": pair["pair_id"],
                        "ligand": pair["ligand"],
                        "receptor": pair["receptor"],
                        "paper_note": pair["paper_note"],
                        "present_in_raw": False,
                        "present_in_current_top20": present_in_top20,
                        "n_matching_rows": 0,
                        "best_source": "",
                        "best_target": "",
                        "best_magnitude_rank": None,
                        "best_lrscore": None,
                        "best_specificity_rank": None,
                        "best_cellphone_pvals": None,
                    }
                )
                continue

            subset.insert(0, "subtype", subtype)
            subset.insert(1, "disease", meta["disease"])
            subset.insert(2, "pair_id", pair["pair_id"])
            subset.insert(3, "paper_note", pair["paper_note"])
            best = _best_row(subset)
            all_rows.append(subset)
            best_rows.append(
                {
                    "subtype": subtype,
                    "disease": meta["disease"],
                    "pair_id": pair["pair_id"],
                    "ligand": pair["ligand"],
                    "receptor": pair["receptor"],
                    "paper_note": pair["paper_note"],
                    "source": best["source"],
                    "target": best["target"],
                    "magnitude_rank": best["magnitude_rank"],
                    "specificity_rank": best["specificity_rank"],
                    "lrscore": best["lrscore"],
                    "cellphone_pvals": best["cellphone_pvals"],
                    "present_in_current_top20": present_in_top20,
                }
            )
            summary_rows.append(
                {
                    "subtype": subtype,
                    "disease": meta["disease"],
                    "pair_id": pair["pair_id"],
                    "ligand": pair["ligand"],
                    "receptor": pair["receptor"],
                    "paper_note": pair["paper_note"],
                    "present_in_raw": True,
                    "present_in_current_top20": present_in_top20,
                    "n_matching_rows": int(len(subset)),
                    "best_source": best["source"],
                    "best_target": best["target"],
                    "best_magnitude_rank": float(best["magnitude_rank"]),
                    "best_lrscore": float(best["lrscore"]),
                    "best_specificity_rank": float(best["specificity_rank"]),
                    "best_cellphone_pvals": float(best["cellphone_pvals"]),
                }
            )

    all_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    best_df = pd.DataFrame(best_rows)
    summary_df = pd.DataFrame(summary_rows)
    visibility_df = pd.DataFrame(visibility_rows)

    if not all_df.empty:
        all_df.to_csv(out_dir / "checkpoint_pair_all_rows.csv", index=False)
    best_df.to_csv(out_dir / "checkpoint_pair_best_rows.csv", index=False)
    summary_df.to_csv(out_dir / "checkpoint_pair_summary.csv", index=False)
    visibility_df.to_csv(out_dir / "checkpoint_pair_visibility.csv", index=False)

    presence_matrix = summary_df.pivot_table(
        index=["pair_id", "ligand", "receptor", "paper_note"],
        columns="subtype",
        values=[
            "present_in_raw",
            "present_in_current_top20",
            "best_magnitude_rank",
            "best_lrscore",
        ],
        aggfunc="first",
    )
    presence_matrix.to_csv(out_dir / "checkpoint_pair_presence_matrix.csv")

    manifest = {
        "project": out_dir.name,
        "generated_at": _now_utc(),
        "input_subtype_cellcomm": {
            subtype: str(meta["raw_csv"]) for subtype, meta in subtype_meta.items()
        },
        "input_current_top20": {
            subtype: str(meta["top20_csv"]) for subtype, meta in subtype_meta.items()
        },
        "focus_pairs": list(FOCUS_PAIRS),
        "ranking_note": "Lower magnitude_rank is better. lrscore is retained as secondary context.",
        "n_focus_pairs": len(FOCUS_PAIRS),
        "n_summary_rows": int(len(summary_df)),
    }
    (out_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (out_dir / "module_status.csv").write_text(
        "module,status,message\ncheckpoint_pair_audit,ok,completed\n",
        encoding="utf-8",
    )
    return out_dir


def summarize_hierarchy(*, audit_dir: Path, output_dir: Path | None = None) -> Path:
    audit_dir = Path(audit_dir)
    out_dir = Path(output_dir) if output_dir is not None else audit_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_df = pd.read_csv(audit_dir / "checkpoint_pair_summary.csv")
    report_rows: list[dict[str, object]] = []
    subtype_summary: dict[str, dict[str, object]] = {}

    for _, row in summary_df.iterrows():
        surface_status, claimable_status = _surface_status(
            bool(row["present_in_raw"]),
            bool(row["present_in_current_top20"]),
        )
        best_edge_label = (
            f"{row['best_source']} -> {row['best_target']}"
            if bool(row["present_in_raw"])
            else "not detected in raw subtype LIANA"
        )
        report_rows.append(
            {
                "subtype": row["subtype"],
                "disease": row["disease"],
                "pair_id": row["pair_id"],
                "ligand": row["ligand"],
                "receptor": row["receptor"],
                "paper_note": row["paper_note"],
                "display_label": f"{row['ligand']} -> {row['receptor']}",
                "best_edge_label": best_edge_label,
                "surface_status": surface_status,
                "claimable_status": claimable_status,
                "present_in_raw": bool(row["present_in_raw"]),
                "present_in_current_top20": bool(row["present_in_current_top20"]),
                "best_magnitude_rank": row["best_magnitude_rank"],
                "best_lrscore": row["best_lrscore"],
            }
        )

    report_df = pd.DataFrame(report_rows)
    report_df.to_csv(out_dir / "checkpoint_hierarchy_report.csv", index=False)

    for subtype, group in report_df.groupby("subtype"):
        visible = group[group["surface_status"] == "visible_in_current_top20"]["display_label"].tolist()
        raw_only = group[group["surface_status"] == "raw_only"]["display_label"].tolist()
        missing = group[group["surface_status"] == "missing_in_raw"]["display_label"].tolist()
        subtype_summary[subtype] = {
            "hierarchy_status": "partial",
            "visible_in_current_top20": visible,
            "raw_only": raw_only,
            "missing_in_raw": missing,
            "note": "Raw pair presence does not by itself establish full paper-faithful hierarchy reproduction.",
        }

    payload = {
        "generated_at": _now_utc(),
        "audit_dir": str(audit_dir),
        "report_csv": str(out_dir / "checkpoint_hierarchy_report.csv"),
        "subtypes": subtype_summary,
        "global_note": "The checkpoint hierarchy remains partial until the paper-facing hierarchy itself, not just raw pair presence, is reproduced.",
    }
    (out_dir / "checkpoint_hierarchy_report.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    lines = [
        "# NC2024 checkpoint hierarchy summary",
        "",
        f"- Generated at: `{payload['generated_at']}`",
        f"- Audit dir: `{audit_dir}`",
        "",
        "## Overall conclusion",
        "",
        "- The checkpoint hierarchy remains `partial`.",
        "- Raw subtype LIANA pair presence does not by itself establish full paper-faithful hierarchy reproduction.",
        "",
    ]
    for subtype in ["luad", "lusc"]:
        if subtype not in subtype_summary:
            continue
        info = subtype_summary[subtype]
        lines.extend(
            [
                f"## {subtype.upper()}",
                "",
                f"- Hierarchy status: `{info['hierarchy_status']}`",
                f"- Visible in current top20: {', '.join(info['visible_in_current_top20']) or 'none'}",
                f"- Raw only: {', '.join(info['raw_only']) or 'none'}",
                f"- Missing in raw: {', '.join(info['missing_in_raw']) or 'none'}",
                "",
            ]
        )
    (out_dir / "checkpoint_hierarchy_memo.md").write_text("\n".join(lines), encoding="utf-8")
    return out_dir
