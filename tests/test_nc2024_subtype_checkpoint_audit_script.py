from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from workflow.modular.reporting.nc2024_checkpoint_audit import run_audit, summarize_hierarchy


def write_raw_csv(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def write_top20_csv(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def base_raw_row(ligand: str, receptor: str, source: str, target: str, magnitude_rank: float) -> dict[str, object]:
    return {
        "source": source,
        "target": target,
        "ligand_complex": ligand,
        "receptor_complex": receptor,
        "lr_means": 1.0,
        "cellphone_pvals": 0.0,
        "expr_prod": 1.0,
        "scaled_weight": 1.0,
        "lr_logfc": 1.0,
        "spec_weight": 1.0,
        "lrscore": 0.8,
        "specificity_rank": 0.2,
        "magnitude_rank": magnitude_rank,
    }


def test_run_audit_writes_expected_outputs(tmp_path):
    luad_raw = tmp_path / "luad_raw.csv"
    lusc_raw = tmp_path / "lusc_raw.csv"
    luad_top20 = tmp_path / "luad_top20.csv"
    lusc_top20 = tmp_path / "lusc_top20.csv"

    write_raw_csv(
        luad_raw,
        [
            base_raw_row("LGALS9", "HAVCR2", "Fibroblast", "Tumor epithelial", 0.11),
            base_raw_row("NECTIN2", "CD226", "Unknown", "Dendritic cell", 0.05),
            base_raw_row("NECTIN2", "TIGIT", "Unknown", "Fibroblast", 0.24),
            base_raw_row("CD80", "CTLA4", "Myeloid/Macro", "Fibroblast", 0.36),
            base_raw_row("CD86", "CTLA4", "Dendritic cell", "Fibroblast", 0.42),
            base_raw_row("HLA-F", "LILRB1", "NK cell", "Myeloid/Macro", 0.37),
            base_raw_row("HLA-F", "LILRB2", "NK cell", "Myeloid/Macro", 0.29),
        ],
    )
    write_raw_csv(
        lusc_raw,
        [
            base_raw_row("LGALS9", "HAVCR2", "Unknown", "Dendritic cell", 0.43),
            base_raw_row("NECTIN2", "CD226", "Myeloid/Macro", "Mast cell", 0.46),
            base_raw_row("NECTIN2", "TIGIT", "Myeloid/Macro", "Unknown", 0.94),
            base_raw_row("NECTIN3", "TIGIT", "Tumor epithelial", "Unknown", 0.16),
            base_raw_row("CD96", "NECTIN1", "Myeloid/Macro", "Tumor epithelial", 0.02),
            base_raw_row("CD80", "CTLA4", "Dendritic cell", "Unknown", 0.53),
            base_raw_row("CD86", "CTLA4", "Unknown", "Dendritic cell", 0.16),
            base_raw_row("HLA-F", "LILRB1", "NK cell", "Myeloid/Macro", 0.35),
            base_raw_row("HLA-F", "LILRB2", "NK cell", "Myeloid/Macro", 0.29),
        ],
    )
    write_top20_csv(
        luad_top20,
        [{"ligand": "CD80", "receptor": "CTLA4"}, {"ligand": "CD86", "receptor": "CTLA4"}],
    )
    write_top20_csv(
        lusc_top20,
        [{"ligand": "CD80", "receptor": "CTLA4"}, {"ligand": "CD86", "receptor": "CTLA4"}],
    )

    out_dir = run_audit(
        luad_raw_csv=luad_raw,
        lusc_raw_csv=lusc_raw,
        luad_top20_csv=luad_top20,
        lusc_top20_csv=lusc_top20,
        output_dir=tmp_path / "out",
        project="demo_audit",
    )

    assert out_dir.name == "demo_audit"
    assert (out_dir / "checkpoint_pair_all_rows.csv").exists()
    assert (out_dir / "checkpoint_pair_summary.csv").exists()
    assert (out_dir / "checkpoint_pair_best_rows.csv").exists()
    assert (out_dir / "checkpoint_pair_visibility.csv").exists()
    assert (out_dir / "checkpoint_pair_presence_matrix.csv").exists()
    assert (out_dir / "run_manifest.json").exists()
    assert (out_dir / "module_status.csv").exists()


def test_run_audit_marks_raw_only_pairs_correctly(tmp_path):
    luad_raw = tmp_path / "luad_raw.csv"
    lusc_raw = tmp_path / "lusc_raw.csv"
    luad_top20 = tmp_path / "luad_top20.csv"
    lusc_top20 = tmp_path / "lusc_top20.csv"

    write_raw_csv(luad_raw, [base_raw_row("LGALS9", "HAVCR2", "Fibroblast", "Tumor epithelial", 0.11)])
    write_raw_csv(lusc_raw, [base_raw_row("CD96", "NECTIN1", "Myeloid/Macro", "Tumor epithelial", 0.02)])
    write_top20_csv(luad_top20, [{"ligand": "CD80", "receptor": "CTLA4"}])
    write_top20_csv(lusc_top20, [{"ligand": "CD80", "receptor": "CTLA4"}])

    out_dir = run_audit(
        luad_raw_csv=luad_raw,
        lusc_raw_csv=lusc_raw,
        luad_top20_csv=luad_top20,
        lusc_top20_csv=lusc_top20,
        output_dir=tmp_path / "out",
    )
    summary = pd.read_csv(out_dir / "checkpoint_pair_summary.csv")

    luad_tim3 = summary[(summary["subtype"] == "luad") & (summary["pair_id"] == "TIM3_axis")].iloc[0]
    lusc_cd96 = summary[(summary["subtype"] == "lusc") & (summary["pair_id"] == "CD96_axis")].iloc[0]

    assert bool(luad_tim3["present_in_raw"]) is True
    assert bool(luad_tim3["present_in_current_top20"]) is False
    assert bool(lusc_cd96["present_in_raw"]) is True
    assert bool(lusc_cd96["present_in_current_top20"]) is False


def test_run_audit_marks_missing_pairs_correctly(tmp_path):
    luad_raw = tmp_path / "luad_raw.csv"
    lusc_raw = tmp_path / "lusc_raw.csv"
    luad_top20 = tmp_path / "luad_top20.csv"
    lusc_top20 = tmp_path / "lusc_top20.csv"

    write_raw_csv(luad_raw, [base_raw_row("CD80", "CTLA4", "Myeloid/Macro", "Fibroblast", 0.36)])
    write_raw_csv(lusc_raw, [base_raw_row("CD80", "CTLA4", "Dendritic cell", "Unknown", 0.53)])
    write_top20_csv(luad_top20, [{"ligand": "CD80", "receptor": "CTLA4"}])
    write_top20_csv(lusc_top20, [{"ligand": "CD80", "receptor": "CTLA4"}])

    out_dir = run_audit(
        luad_raw_csv=luad_raw,
        lusc_raw_csv=lusc_raw,
        luad_top20_csv=luad_top20,
        lusc_top20_csv=lusc_top20,
        output_dir=tmp_path / "out",
    )
    summary = pd.read_csv(out_dir / "checkpoint_pair_summary.csv")

    luad_cd96 = summary[(summary["subtype"] == "luad") & (summary["pair_id"] == "CD96_axis")].iloc[0]
    assert bool(luad_cd96["present_in_raw"]) is False
    assert bool(luad_cd96["present_in_current_top20"]) is False
    assert int(luad_cd96["n_matching_rows"]) == 0


def test_run_audit_writes_manifest_with_expected_inputs(tmp_path):
    luad_raw = tmp_path / "luad_raw.csv"
    lusc_raw = tmp_path / "lusc_raw.csv"
    luad_top20 = tmp_path / "luad_top20.csv"
    lusc_top20 = tmp_path / "lusc_top20.csv"

    write_raw_csv(luad_raw, [base_raw_row("CD80", "CTLA4", "Myeloid/Macro", "Fibroblast", 0.36)])
    write_raw_csv(lusc_raw, [base_raw_row("CD80", "CTLA4", "Dendritic cell", "Unknown", 0.53)])
    write_top20_csv(luad_top20, [{"ligand": "CD80", "receptor": "CTLA4"}])
    write_top20_csv(lusc_top20, [{"ligand": "CD80", "receptor": "CTLA4"}])

    out_dir = run_audit(
        luad_raw_csv=luad_raw,
        lusc_raw_csv=lusc_raw,
        luad_top20_csv=luad_top20,
        lusc_top20_csv=lusc_top20,
        output_dir=tmp_path / "out",
        project="manifest_case",
    )
    manifest = json.loads((out_dir / "run_manifest.json").read_text())
    status = (out_dir / "module_status.csv").read_text()

    assert manifest["project"] == "manifest_case"
    assert manifest["input_subtype_cellcomm"]["luad"] == str(luad_raw)
    assert manifest["input_current_top20"]["lusc"] == str(lusc_top20)
    assert manifest["n_focus_pairs"] >= 9
    assert "checkpoint_pair_audit,ok,completed" in status


def test_summarize_hierarchy_writes_report_files_and_statuses(tmp_path):
    luad_raw = tmp_path / "luad_raw.csv"
    lusc_raw = tmp_path / "lusc_raw.csv"
    luad_top20 = tmp_path / "luad_top20.csv"
    lusc_top20 = tmp_path / "lusc_top20.csv"

    write_raw_csv(
        luad_raw,
        [
            base_raw_row("LGALS9", "HAVCR2", "Fibroblast", "Tumor epithelial", 0.11),
            base_raw_row("CD80", "CTLA4", "Myeloid/Macro", "Fibroblast", 0.36),
        ],
    )
    write_raw_csv(
        lusc_raw,
        [
            base_raw_row("CD96", "NECTIN1", "Myeloid/Macro", "Tumor epithelial", 0.02),
            base_raw_row("CD80", "CTLA4", "Dendritic cell", "Unknown", 0.53),
        ],
    )
    write_top20_csv(luad_top20, [{"ligand": "CD80", "receptor": "CTLA4"}])
    write_top20_csv(lusc_top20, [{"ligand": "CD80", "receptor": "CTLA4"}])

    audit_dir = run_audit(
        luad_raw_csv=luad_raw,
        lusc_raw_csv=lusc_raw,
        luad_top20_csv=luad_top20,
        lusc_top20_csv=lusc_top20,
        output_dir=tmp_path / "out",
        project="audit_case",
    )
    summarize_hierarchy(audit_dir=audit_dir)

    report = pd.read_csv(audit_dir / "checkpoint_hierarchy_report.csv")
    payload = json.loads((audit_dir / "checkpoint_hierarchy_report.json").read_text())
    memo = (audit_dir / "checkpoint_hierarchy_memo.md").read_text()

    luad_tim3 = report[(report["subtype"] == "luad") & (report["pair_id"] == "TIM3_axis")].iloc[0]
    lusc_cd96 = report[(report["subtype"] == "lusc") & (report["pair_id"] == "CD96_axis")].iloc[0]

    assert luad_tim3["surface_status"] == "raw_only"
    assert luad_tim3["claimable_status"] == "supported_raw_only"
    assert lusc_cd96["surface_status"] == "raw_only"
    assert payload["subtypes"]["luad"]["hierarchy_status"] == "partial"
    assert "# NC2024 checkpoint hierarchy summary" in memo
    assert "Hierarchy status: `partial`" in memo


def load_verify_script_text() -> str:
    repo_root = Path(__file__).resolve().parents[1]
    return (repo_root / "scripts" / "verify_nc2024_subtype_checkpoint_audit.sh").read_text()


def test_verify_wrapper_uses_sc_gpu_and_runs_both_scripts():
    script = load_verify_script_text()
    assert "conda run -n sc_gpu" in script
    assert "python -m pytest -q -o addopts=''" in script
    assert "scripts/audit_nc2024_subtype_checkpoint_pairs.py" in script
    assert "scripts/summarize_nc2024_checkpoint_hierarchy.py" in script
