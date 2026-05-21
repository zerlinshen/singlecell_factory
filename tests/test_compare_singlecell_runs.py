from __future__ import annotations

import importlib.util
import json
from pathlib import Path


_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "compare_singlecell_runs.py"
_SPEC = importlib.util.spec_from_file_location("compare_singlecell_runs", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
compare = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(compare)


def _write_run(
    run_dir: Path,
    *,
    project: str,
    backend: str,
    method: str,
    doublet_rate: float,
    doublets: int,
    cells: int,
    clusters: int,
    de_genes: int,
    undercall_ratio: float | None,
) -> None:
    run_dir.mkdir(parents=True)
    metadata = {
        "raw_cells": 1200,
        "raw_genes": 400,
        "cells_after_qc": 1000,
        "genes_after_qc": 380,
        "cells_after_doublet_removal": cells,
        "doublet_backend_requested": backend,
        "doublet_backend": backend,
        "doublet_method": method,
        "doublets_detected": doublets,
        "doublet_rate_pct": doublet_rate,
        "doublet_undercall_ratio": undercall_ratio,
        "doublet_undercall_warning": undercall_ratio is not None and undercall_ratio > 5,
        "module_runtime_sec": {
            "doublet_detection": 3.0,
            "ambient_correction": 0.1,
        },
        "ambient_correction_engine": "decontx",
        "ambient_correction_decision": "skipped_no_trigger",
        "n_clusters": clusters,
        "annotation_unknown_pct": 0.0,
        "de_significant_genes": de_genes,
        "grn_top_tfs": ["STAT1", "MYC", "SNAI1"],
        "pipeline_wall_seconds": 20.0,
        "pseudobulk_de_status": "skipped_missing_contrast_contract",
        "pseudobulk_de_mode": "skipped",
    }
    (run_dir / "run_manifest.json").write_text(
        json.dumps({"project": project, "generated_at": "2026-05-22T00:00:00Z", "metadata": metadata}),
        encoding="utf-8",
    )
    (run_dir / "module_status.csv").write_text(
        "module,status,message\n"
        "cellranger,ok,completed\n"
        "qc,ok,completed\n"
        "doublet_detection,ok,completed\n"
        "pseudobulk_de,skipped,No contrast\n",
        encoding="utf-8",
    )
    (run_dir / "pathway_analysis").mkdir()
    (run_dir / "pathway_analysis" / "pathway_enrichment.csv").write_text(
        "cluster,term,overlap,gene_set_size,pval,genes,padj\n"
        "0,HALLMARK_INTERFERON_GAMMA_RESPONSE,4,10,1e-6,IRF1,1e-5\n"
        "1,HALLMARK_TNFA_SIGNALING_VIA_NFKB,3,10,1e-4,NFKB1,1e-3\n",
        encoding="utf-8",
    )
    (run_dir / "cell_communication").mkdir()
    (run_dir / "cell_communication" / "cell_communication_lr.csv").write_text(
        "ligand,receptor,source,target,ligand_expr,receptor_expr,lr_score\n"
        "SPP1,CD44,Myeloid,T cell,1.0,1.2,1.2\n",
        encoding="utf-8",
    )
    (run_dir / "panel.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100"></svg>\n',
        encoding="utf-8",
    )


def test_build_report_selects_conditional_consensus_lane(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline"
    consensus = tmp_path / "consensus"
    _write_run(
        baseline,
        project="baseline",
        backend="scrublet",
        method="scrublet",
        doublet_rate=0.03,
        doublets=3,
        cells=997,
        clusters=12,
        de_genes=500,
        undercall_ratio=100.0,
    )
    _write_run(
        consensus,
        project="consensus",
        backend="consensus",
        method="consensus_or_scrublet_scdblfinder",
        doublet_rate=8.5,
        doublets=85,
        cells=915,
        clusters=11,
        de_genes=470,
        undercall_ratio=None,
    )

    report = compare.build_report([("baseline", baseline), ("consensus", consensus)])

    assert report["evaluation"]["accepted_lane"] == "consensus"
    assert report["evaluation"]["global_default_change"] is False
    assert report["evaluation"]["lane_decisions"]["baseline"]["status"] == "accepted_as_baseline_only"
    assert report["runs"][0]["figure_inventory"]["vector"] == 1
    assert report["comparisons"][0]["pathway_overlap"]["jaccard"] == 1.0


def test_render_markdown_includes_key_decision(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline"
    _write_run(
        baseline,
        project="baseline",
        backend="scrublet",
        method="scrublet",
        doublet_rate=0.5,
        doublets=5,
        cells=995,
        clusters=8,
        de_genes=300,
        undercall_ratio=12.0,
    )

    report = compare.build_report([("baseline", baseline)])
    markdown = compare.render_markdown(report)

    assert "Global default changed: `False`" in markdown
    assert "Keep Scrublet as the global default" in markdown
