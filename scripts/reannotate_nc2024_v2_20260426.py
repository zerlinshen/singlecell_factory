"""A.3: Re-annotate NC2024 tumor + BH cohorts from after_clustering checkpoints.

Loads checkpoint zarr, runs AnnotationModule with cluster_voting, writes:
- results/<cohort>_v2/tables/cluster_majority_cell_type.csv
- results/<cohort>_v2/tables/cluster_score_matrix.csv
- results/<cohort>_v2/final_adata.h5ad
- ops/run_ledger/ entry
- ops/nc2024_methodology_audit/ALIGNMENT_REPORT_v2_2026-04-26.md
"""
from __future__ import annotations

import json
import sys
import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")

from workflow.modular.modules.annotation import AnnotationModule, DEFAULT_MARKERS
from workflow.modular.config import PipelineConfig, CellRangerConfig

# ── Acceptance gate clusters (tumor cohort) ──────────────────────────────────
# NOTE: these cluster IDs are from the after_clustering checkpoint (37 clusters, res=1.0)
ACCEPTANCE_GATE = {
    "7":  "B cell",
    "17": "Mast cell",
    "13": "Plasma cell",
    "14": "Myeloid/Macro",
    "10": "NK cell",
}


class _Ctx:
    def __init__(self, adata, cfg, out_dir):
        self.adata = adata
        self.cfg = cfg
        self.metadata = {}
        self.table_dir = out_dir / "tables"
        self.figure_dir = out_dir / "figures"
        self.table_dir.mkdir(parents=True, exist_ok=True)
        self.figure_dir.mkdir(parents=True, exist_ok=True)


def run_cohort(name: str, checkpoint_zarr: Path, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load checkpoint, annotate, write outputs. Returns (majority, score_matrix)."""
    print(f"\n{'='*60}", flush=True)
    print(f"Cohort: {name}", flush=True)
    print(f"Checkpoint: {checkpoint_zarr}", flush=True)
    print(f"Output dir: {out_dir}", flush=True)

    print("Loading checkpoint...", flush=True)
    adata = ad.read_zarr(checkpoint_zarr)
    print(f"  n_obs={adata.n_obs}, n_vars={adata.n_vars}, leiden clusters={adata.obs['leiden'].nunique()}", flush=True)

    cfg = PipelineConfig(
        project=name,
        output_dir=out_dir,
        cellranger=CellRangerConfig(
            sample_root=ROOT / "data/raw",
            outs_dir=ROOT / "data/raw",
        ),
        annotation_strategy="cluster_voting",
        annotation_confidence_threshold=0.1,
    )

    ctx = _Ctx(adata, cfg, out_dir)
    print("Running AnnotationModule (cluster_voting)...", flush=True)
    AnnotationModule().run(ctx)
    print(f"  Done. Unknown pct: {ctx.metadata.get('annotation_unknown_pct')}%", flush=True)

    # Write final adata
    final_path = out_dir / "final_adata.h5ad"
    print(f"Writing final_adata.h5ad -> {final_path}", flush=True)
    adata.write_h5ad(final_path)

    majority = pd.read_csv(ctx.table_dir / "cluster_majority_cell_type.csv", index_col=0)
    majority.index = majority.index.astype(str)
    score_matrix = pd.read_csv(ctx.table_dir / "cluster_score_matrix.csv", index_col=0)
    score_matrix.index = score_matrix.index.astype(str)

    print("\ncluster_majority_cell_type:", flush=True)
    print(majority.to_string(), flush=True)

    return majority, score_matrix, ctx.metadata


def check_acceptance(majority: pd.DataFrame, cohort_name: str) -> tuple[int, list[str]]:
    """Return (n_passing, failures)."""
    failures = []
    for cluster, expected in ACCEPTANCE_GATE.items():
        if cluster not in majority.index:
            failures.append(f"cluster {cluster}: NOT IN INDEX (expected {expected})")
            continue
        actual = majority.loc[cluster, "majority_cell_type"]
        if actual != expected:
            failures.append(f"cluster {cluster}: got '{actual}', expected '{expected}'")
    n_pass = len(ACCEPTANCE_GATE) - len(failures)
    return n_pass, failures


def write_ledger_entry(cohort_name: str, out_dir: Path, metadata: dict, n_pass: int, failures: list[str]) -> None:
    ledger_dir = ROOT / "ops" / "run_ledger"
    ledger_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    entry_path = ledger_dir / f"{cohort_name}_{ts}.json"
    entry = {
        "cohort": cohort_name,
        "timestamp": ts,
        "strategy": "cluster_voting",
        "output_dir": str(out_dir),
        "annotation_unknown_pct": metadata.get("annotation_unknown_pct"),
        "acceptance_gate": {"n_pass": n_pass, "n_total": len(ACCEPTANCE_GATE), "failures": failures},
        "metadata": {k: v for k, v in metadata.items() if k != "annotation_unknown_pct"},
    }
    entry_path.write_text(json.dumps(entry, indent=2))
    print(f"Ledger entry: {entry_path}", flush=True)


# ── Run tumor cohort ─────────────────────────────────────────────────────────
tumor_checkpoint = ROOT / "results/nc2024_tumor_20260426/nc2024_tumor_20260426_121139/.checkpoints/after_clustering.zarr"
tumor_out = ROOT / "results/nc2024_tumor_20260426_v2"
tumor_majority, tumor_score_matrix, tumor_meta = run_cohort(
    "nc2024_tumor_20260426_v2", tumor_checkpoint, tumor_out
)

tumor_n_pass, tumor_failures = check_acceptance(tumor_majority, "tumor")
print(f"\nTumor acceptance gate: {tumor_n_pass}/5 passing", flush=True)
if tumor_failures:
    print("FAILURES:", flush=True)
    for f in tumor_failures:
        print(f"  {f}", flush=True)

write_ledger_entry("nc2024_tumor_20260426_v2", tumor_out, tumor_meta, tumor_n_pass, tumor_failures)

# ── Run BH cohort ─────────────────────────────────────────────────────────────
bh_checkpoint = ROOT / "results/nc2024_bh_20260426/nc2024_bh_20260426_132950/.checkpoints/after_clustering.zarr"
bh_out = ROOT / "results/nc2024_bh_20260426_v2"
bh_majority, bh_score_matrix, bh_meta = run_cohort(
    "nc2024_bh_20260426_v2", bh_checkpoint, bh_out
)

write_ledger_entry("nc2024_bh_20260426_v2", bh_out, bh_meta, 0, ["BH cohort: no acceptance gate defined"])

# ── Alignment report ──────────────────────────────────────────────────────────
report_dir = ROOT / "ops" / "nc2024_methodology_audit"
report_dir.mkdir(parents=True, exist_ok=True)
report_path = report_dir / "ALIGNMENT_REPORT_v2_2026-04-26.md"

# Load v1 majority for comparison
v1_majority_path = ROOT / "results/nc2024_tumor_20260426/nc2024_tumor_20260426_121139"
v1_majority = None
for candidate in [v1_majority_path / "tables/cluster_majority_cell_type.csv",
                   v1_majority_path / "cluster_majority_cell_type.csv"]:
    if candidate.exists():
        v1_majority = pd.read_csv(candidate, index_col=0)
        v1_majority.index = v1_majority.index.astype(str)
        break

with open(report_path, "w") as f:
    f.write("# NC2024 Annotation Alignment Report v2\n\n")
    f.write(f"**Date:** 2026-04-26\n")
    f.write(f"**Strategy:** cluster_voting (aggregate score_genes per cluster then argmax)\n")
    f.write(f"**Acceptance gate:** {tumor_n_pass}/5 critical clusters correct\n\n")

    f.write("## Acceptance Gate — Tumor Cohort (NC2024, 877k cells)\n\n")
    f.write("| Cluster | Markers | Paper Expected | Tumor v2 (cluster_voting) | Status |\n")
    f.write("|---------|---------|---------------|--------------------------|--------|\n")
    gate_data = [
        ("7",  "MS4A1+CD79A+BANK1+", "B cell"),
        ("17", "TPSAB1+CPA3+",       "Mast cell"),
        ("13", "IGLV3-1+MZB1+",      "Plasma cell"),
        ("14", "CD163+OLR1+LYZ+",    "Myeloid/Macro"),
        ("10", "KLRD1+GNLY+KLRF1+",  "NK cell"),
    ]
    for cluster, markers, expected in gate_data:
        actual = tumor_majority.loc[cluster, "majority_cell_type"] if cluster in tumor_majority.index else "MISSING"
        status = "PASS" if actual == expected else "FAIL"
        f.write(f"| {cluster} | {markers} | {expected} | {actual} | {status} |\n")

    f.write(f"\n## Full Tumor v2 Cluster Assignments\n\n")
    f.write("| Cluster | v2 Cell Type |\n")
    f.write("|---------|-------------|\n")
    for cluster in sorted(tumor_majority.index, key=lambda x: int(x) if str(x).isdigit() else x):
        ct = tumor_majority.loc[cluster, "majority_cell_type"]
        f.write(f"| {cluster} | {ct} |\n")

    f.write(f"\n## BH Cohort v2 Cluster Assignments\n\n")
    f.write("| Cluster | v2 Cell Type |\n")
    f.write("|---------|-------------|\n")
    for cluster in sorted(bh_majority.index, key=lambda x: int(x) if str(x).isdigit() else x):
        ct = bh_majority.loc[cluster, "majority_cell_type"]
        f.write(f"| {cluster} | {ct} |\n")

    if tumor_n_pass == 5:
        f.write("\n## Conclusion\n\nAll 5 critical clusters pass the acceptance gate. ")
        f.write("cluster_voting strategy successfully corrects the v1 annotation drift.\n")
    else:
        f.write(f"\n## Conclusion\n\n{tumor_n_pass}/5 clusters pass. Failures:\n")
        for failure in tumor_failures:
            f.write(f"- {failure}\n")

print(f"\nAlignment report: {report_path}", flush=True)

# ── Final summary ─────────────────────────────────────────────────────────────
print(f"\n{'='*60}", flush=True)
print(f"FINAL: Tumor acceptance gate = {tumor_n_pass}/5", flush=True)
if tumor_n_pass == 5:
    print("ALL 5 CRITICAL CLUSTERS PASS. A.3 COMPLETE.", flush=True)
else:
    print("PARTIAL PASS. Failures:", flush=True)
    for f in tumor_failures:
        print(f"  {f}", flush=True)
    print("Reporting to lead.", flush=True)
    sys.exit(1)
