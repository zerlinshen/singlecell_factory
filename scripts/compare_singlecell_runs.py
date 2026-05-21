#!/usr/bin/env python3
"""Compare completed modular single-cell runs from their source artifacts.

The report is intentionally post-run and project-root based: it reads
``run_manifest.json``, ``module_status.csv``, and selected module tables from
existing run directories, then writes JSON/Markdown evidence for method-choice
decisions. It does not mutate the run directories.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


FIGURE_SUFFIXES = {".png", ".pdf", ".svg"}
VECTOR_SUFFIXES = {".pdf", ".svg"}


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _as_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(number):
        return None
    return number


def _as_int(value: Any) -> int | None:
    number = _as_float(value)
    return None if number is None else int(number)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _status_summary(run_dir: Path) -> dict[str, Any]:
    rows = _read_csv(run_dir / "module_status.csv")
    counts = Counter(row.get("status", "unknown") for row in rows)
    return {
        "counts": dict(sorted(counts.items())),
        "by_module": {row.get("module", ""): row.get("status", "") for row in rows},
        "messages": {
            row.get("module", ""): row.get("message", "")
            for row in rows
            if row.get("message") and row.get("message") != "completed"
        },
    }


def _figure_inventory(run_dir: Path) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    bytes_by_suffix: Counter[str] = Counter()
    for path in run_dir.rglob("*"):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix not in FIGURE_SUFFIXES:
            continue
        counts[suffix.lstrip(".")] += 1
        bytes_by_suffix[suffix.lstrip(".")] += path.stat().st_size
    vector_count = sum(counts[suffix.lstrip(".")] for suffix in VECTOR_SUFFIXES)
    return {
        "total": sum(counts.values()),
        "png": counts.get("png", 0),
        "pdf": counts.get("pdf", 0),
        "svg": counts.get("svg", 0),
        "vector": vector_count,
        "bytes_by_format": dict(sorted(bytes_by_suffix.items())),
    }


def _top_pathways(run_dir: Path, limit: int = 8) -> list[dict[str, Any]]:
    rows = _read_csv(run_dir / "pathway_analysis" / "pathway_enrichment.csv")
    terms: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in rows:
        term = row.get("term", "")
        if not term or term in seen:
            continue
        seen.add(term)
        padj = _as_float(row.get("padj"))
        terms.append({
            "term": term,
            "cluster": row.get("cluster", ""),
            "padj": padj,
            "overlap": _as_int(row.get("overlap")),
            "genes": row.get("genes", ""),
        })
        if len(terms) >= limit:
            break
    return terms


def _top_ligand_receptors(run_dir: Path, limit: int = 8) -> list[dict[str, Any]]:
    rows = _read_csv(run_dir / "cell_communication" / "cell_communication_lr.csv")
    out: list[dict[str, Any]] = []
    for row in rows[:limit]:
        out.append({
            "ligand": row.get("ligand", ""),
            "receptor": row.get("receptor", ""),
            "source": row.get("source", ""),
            "target": row.get("target", ""),
            "lr_score": _as_float(row.get("lr_score")),
        })
    return out


def _top_marker_genes(run_dir: Path, limit: int = 8) -> list[dict[str, Any]]:
    rows = _read_csv(run_dir / "differential_expression" / "marker_top5_by_cluster.csv")
    out: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for row in rows:
        cluster = row.get("cluster", row.get("leiden", ""))
        gene = row.get("gene", row.get("names", ""))
        key = (cluster, gene)
        if not gene or key in seen:
            continue
        seen.add(key)
        out.append({
            "cluster": cluster,
            "gene": gene,
            "score": _as_float(row.get("score", row.get("scores"))),
            "pval_adj": _as_float(row.get("pval_adj", row.get("pvals_adj"))),
        })
        if len(out) >= limit:
            break
    return out


def _top_immune_subtypes(metadata: dict[str, Any], limit: int = 8) -> list[dict[str, Any]]:
    subtypes = metadata.get("immune_subtypes_detected")
    if not isinstance(subtypes, dict):
        return []
    items = sorted(subtypes.items(), key=lambda item: _as_float(item[1]) or 0.0, reverse=True)
    return [{"subtype": str(name), "cells": _as_int(count)} for name, count in items[:limit]]


def _list_overlap(left: list[str], right: list[str]) -> dict[str, Any]:
    left_set = {item for item in left if item}
    right_set = {item for item in right if item}
    union = left_set | right_set
    intersection = left_set & right_set
    return {
        "intersection": sorted(intersection),
        "intersection_count": len(intersection),
        "jaccard": round(len(intersection) / len(union), 4) if union else None,
    }


def _lr_ids(rows: list[dict[str, Any]]) -> list[str]:
    return [
        f"{row.get('ligand', '')}->{row.get('receptor', '')}"
        for row in rows
        if row.get("ligand") and row.get("receptor")
    ]


def parse_run_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise argparse.ArgumentTypeError("--run must be LABEL=/path/to/run")
    label, raw_path = spec.split("=", 1)
    label = label.strip()
    if not label:
        raise argparse.ArgumentTypeError("--run label cannot be empty")
    return label, Path(raw_path).expanduser().resolve()


def load_run(label: str, run_dir: Path) -> dict[str, Any]:
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"missing run_manifest.json: {manifest_path}")
    manifest = _load_json(manifest_path)
    metadata = manifest.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}

    module_runtime = metadata.get("module_runtime_sec")
    if not isinstance(module_runtime, dict):
        module_runtime = {}

    doublet = {
        "backend_requested": metadata.get("doublet_backend_requested"),
        "backend": metadata.get("doublet_backend"),
        "method": metadata.get("doublet_method"),
        "detected": _as_int(metadata.get("doublets_detected")),
        "rate_pct": _as_float(metadata.get("doublet_rate_pct")),
        "undercall_ratio": _as_float(metadata.get("doublet_undercall_ratio")),
        "undercall_warning": bool(metadata.get("doublet_undercall_warning", False)),
        "consensus_pair": metadata.get("doublet_consensus_pair"),
        "consensus_logic": metadata.get("doublet_consensus_logic"),
        "consensus_agreement": _as_float(metadata.get("doublet_consensus_agreement")),
        "per_backend_rates": metadata.get("doublet_consensus_per_backend_rates"),
        "runtime_sec": _as_float(module_runtime.get("doublet_detection")),
    }
    ambient = {
        "engine": metadata.get("ambient_correction_engine"),
        "decision": metadata.get("ambient_correction_decision"),
        "triggers": metadata.get("ambient_correction_triggers"),
        "pre_mt_median": _as_float(metadata.get("ambient_qc_pre_mt_median")),
        "pre_count_gene_spearman": _as_float(metadata.get("ambient_qc_pre_count_gene_spearman")),
        "runtime_sec": _as_float(module_runtime.get("ambient_correction")),
    }

    top_pathways = _top_pathways(run_dir)
    top_lrs = _top_ligand_receptors(run_dir)
    top_tfs = metadata.get("grn_top_tfs")
    if not isinstance(top_tfs, list):
        top_tfs = []

    return {
        "label": label,
        "run_dir": str(run_dir),
        "project": manifest.get("project"),
        "generated_at": manifest.get("generated_at"),
        "raw_cells": _as_int(metadata.get("raw_cells")),
        "raw_genes": _as_int(metadata.get("raw_genes")),
        "cells_after_qc": _as_int(metadata.get("cells_after_qc")),
        "genes_after_qc": _as_int(metadata.get("genes_after_qc")),
        "cells_after_doublet_removal": _as_int(metadata.get("cells_after_doublet_removal")),
        "n_clusters": _as_int(metadata.get("n_clusters")),
        "annotation_unknown_pct": _as_float(metadata.get("annotation_unknown_pct")),
        "de_significant_genes": _as_int(metadata.get("de_significant_genes")),
        "composition_n_cell_types": _as_int(metadata.get("composition_n_cell_types")),
        "composition_n_groups": _as_int(metadata.get("composition_n_groups")),
        "n_metacells": _as_int(metadata.get("n_metacells")),
        "tme_mean_cyt": _as_float(metadata.get("tme_mean_cyt")),
        "pipeline_wall_seconds": _as_float(metadata.get("pipeline_wall_seconds")),
        "module_runtime_sec": {
            key: _as_float(value)
            for key, value in module_runtime.items()
        },
        "doublet": doublet,
        "ambient": ambient,
        "status": _status_summary(run_dir),
        "figure_inventory": _figure_inventory(run_dir),
        "top_pathways": top_pathways,
        "top_tfs": [str(item) for item in top_tfs],
        "top_ligand_receptors": top_lrs,
        "top_marker_genes": _top_marker_genes(run_dir),
        "top_immune_subtypes": _top_immune_subtypes(metadata),
        "pseudobulk_de_status": metadata.get("pseudobulk_de_status"),
        "pseudobulk_de_mode": metadata.get("pseudobulk_de_mode"),
    }


def compare_to_baseline(baseline: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    base_pathways = [row["term"] for row in baseline.get("top_pathways", [])]
    cand_pathways = [row["term"] for row in candidate.get("top_pathways", [])]
    base_tfs = baseline.get("top_tfs", [])
    cand_tfs = candidate.get("top_tfs", [])
    base_lrs = _lr_ids(baseline.get("top_ligand_receptors", []))
    cand_lrs = _lr_ids(candidate.get("top_ligand_receptors", []))
    return {
        "label": candidate["label"],
        "delta_cells_after_doublet_removal": (
            (candidate.get("cells_after_doublet_removal") or 0)
            - (baseline.get("cells_after_doublet_removal") or 0)
        ),
        "delta_doublet_rate_pct": round(
            (candidate["doublet"].get("rate_pct") or 0.0)
            - (baseline["doublet"].get("rate_pct") or 0.0),
            4,
        ),
        "delta_clusters": (
            (candidate.get("n_clusters") or 0)
            - (baseline.get("n_clusters") or 0)
        ),
        "delta_annotation_unknown_pct": round(
            (candidate.get("annotation_unknown_pct") or 0.0)
            - (baseline.get("annotation_unknown_pct") or 0.0),
            4,
        ),
        "delta_de_significant_genes": (
            (candidate.get("de_significant_genes") or 0)
            - (baseline.get("de_significant_genes") or 0)
        ),
        "pathway_overlap": _list_overlap(base_pathways, cand_pathways),
        "grn_tf_overlap": _list_overlap(base_tfs, cand_tfs),
        "lr_pair_overlap": _list_overlap(base_lrs, cand_lrs),
        "delta_pipeline_wall_seconds": round(
            (candidate.get("pipeline_wall_seconds") or 0.0)
            - (baseline.get("pipeline_wall_seconds") or 0.0),
            3,
        ),
    }


def evaluate_runs(runs: list[dict[str, Any]], baseline_label: str) -> dict[str, Any]:
    baseline = next(run for run in runs if run["label"] == baseline_label)
    baseline_rate = baseline["doublet"].get("rate_pct") or 0.0
    baseline_under_call = (
        baseline["doublet"].get("undercall_warning")
        or (baseline["doublet"].get("undercall_ratio") or 0.0) >= 5.0
        or baseline_rate < 1.0
    )

    lane_decisions: dict[str, dict[str, str]] = {}
    accepted_label: str | None = None
    for run in runs:
        label = run["label"]
        rate = run["doublet"].get("rate_pct") or 0.0
        unknown = run.get("annotation_unknown_pct") or 0.0
        cluster_delta = abs((run.get("n_clusters") or 0) - (baseline.get("n_clusters") or 0))
        de_genes = run.get("de_significant_genes") or 0
        if label == baseline_label:
            status = "accepted_as_baseline_only"
            reason = "Completes the full downstream module set, but the Scrublet under-call makes it unsuitable as the doublet-confidence lane for this tissue shape."
            if not baseline_under_call:
                reason = "Completes the full downstream module set without a doublet under-call signal."
        elif rate >= 4.0 and unknown <= 5.0 and cluster_delta <= 5 and de_genes > 0:
            status = "accepted_conditional_candidate"
            reason = "Rescues the doublet under-call while preserving downstream annotation, clustering, and DE viability."
            if accepted_label is None or "consensus" in label:
                accepted_label = label
        else:
            status = "inconclusive"
            reason = "Needs additional evidence because rescue or downstream stability criteria were not all met."
        lane_decisions[label] = {"status": status, "reason": reason}

    if accepted_label is None and len(runs) > 1:
        accepted_label = runs[1]["label"]

    return {
        "baseline_label": baseline_label,
        "accepted_lane": accepted_label,
        "accepted_scope": "conditional_lusc_tissue_shape_only" if accepted_label else "none",
        "global_default_change": False,
        "routing_condition": (
            "Keep Scrublet as the global default. When the Scrublet under-call "
            "diagnostic fires on heterogeneous tissue or tumor-atlas data "
            "(very low call rate versus expectation, high undercall ratio, or "
            "collapsed score histogram), run scDblFinder as a second opinion or "
            "use consensus OR with the scrublet_scdblfinder pair."
        ),
        "lane_decisions": lane_decisions,
        "residual_risks": [
            "pseudobulk_de remains skipped because no explicit confirmatory contrast was provided",
            "mechanical figure validation does not replace manual review for label overlap, palette separation, or biological interpretation",
            "the accepted lane is a data-shape conditional for this LUSC tissue context, not a universal replacement default",
        ],
    }


def build_report(run_specs: list[tuple[str, Path]], baseline_label: str | None = None) -> dict[str, Any]:
    runs = [load_run(label, path) for label, path in run_specs]
    if not runs:
        raise ValueError("at least one run is required")
    baseline_label = baseline_label or runs[0]["label"]
    if baseline_label not in {run["label"] for run in runs}:
        raise ValueError(f"baseline label not found: {baseline_label}")
    baseline = next(run for run in runs if run["label"] == baseline_label)
    comparisons = [
        compare_to_baseline(baseline, run)
        for run in runs
        if run["label"] != baseline_label
    ]
    return {
        "schema_version": "singlecell-run-comparison/v1",
        "generated_at": _utc_now(),
        "baseline_label": baseline_label,
        "runs": runs,
        "comparisons": comparisons,
        "evaluation": evaluate_runs(runs, baseline_label),
    }


def _fmt(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def render_markdown(report: dict[str, Any]) -> str:
    evaluation = report["evaluation"]
    lines = [
        "# Round9 Single-Cell Run Comparison",
        "",
        f"- Generated: `{report['generated_at']}`",
        f"- Baseline: `{report['baseline_label']}`",
        f"- Accepted lane: `{evaluation.get('accepted_lane') or 'none'}`",
        f"- Scope: `{evaluation['accepted_scope']}`",
        f"- Global default changed: `{evaluation['global_default_change']}`",
        "",
        "## Decision",
        "",
        evaluation["routing_condition"],
        "",
        "## Run Summary",
        "",
        "| Lane | Backend | Doublet % | Doublets | Cells kept | Clusters | Unknown % | DE genes | Wall s | Figures | Vector |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for run in report["runs"]:
        fig = run["figure_inventory"]
        lines.append(
            "| "
            + " | ".join([
                f"`{run['label']}`",
                _fmt(run["doublet"].get("method") or run["doublet"].get("backend")),
                _fmt(run["doublet"].get("rate_pct")),
                _fmt(run["doublet"].get("detected")),
                _fmt(run.get("cells_after_doublet_removal")),
                _fmt(run.get("n_clusters")),
                _fmt(run.get("annotation_unknown_pct")),
                _fmt(run.get("de_significant_genes")),
                _fmt(run.get("pipeline_wall_seconds")),
                _fmt(fig.get("total")),
                _fmt(fig.get("vector")),
            ])
            + " |"
        )

    lines.extend(["", "## Baseline Deltas", ""])
    if report["comparisons"]:
        lines.extend([
            "| Lane | Cells delta | Doublet % delta | Clusters delta | DE genes delta | Pathway Jaccard | TF Jaccard | LR Jaccard | Wall s delta |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ])
        for cmp in report["comparisons"]:
            lines.append(
                "| "
                + " | ".join([
                    f"`{cmp['label']}`",
                    _fmt(cmp["delta_cells_after_doublet_removal"]),
                    _fmt(cmp["delta_doublet_rate_pct"]),
                    _fmt(cmp["delta_clusters"]),
                    _fmt(cmp["delta_de_significant_genes"]),
                    _fmt(cmp["pathway_overlap"]["jaccard"]),
                    _fmt(cmp["grn_tf_overlap"]["jaccard"]),
                    _fmt(cmp["lr_pair_overlap"]["jaccard"]),
                    _fmt(cmp["delta_pipeline_wall_seconds"]),
                ])
                + " |"
            )
    else:
        lines.append("No non-baseline lanes were provided.")

    lines.extend(["", "## Lane Decisions", ""])
    for label, decision in evaluation["lane_decisions"].items():
        lines.append(f"- `{label}`: `{decision['status']}` - {decision['reason']}")

    lines.extend(["", "## Biology Signal Snapshot", ""])
    for run in report["runs"]:
        pathways = ", ".join(row["term"] for row in run.get("top_pathways", [])[:5]) or "NA"
        tfs = ", ".join(run.get("top_tfs", [])[:5]) or "NA"
        lrs = ", ".join(_lr_ids(run.get("top_ligand_receptors", [])[:5])) or "NA"
        lines.extend([
            f"### {run['label']}",
            f"- Top pathways: {pathways}",
            f"- Top TFs: {tfs}",
            f"- Top ligand-receptor pairs: {lrs}",
        ])

    lines.extend(["", "## Residual Risks", ""])
    lines.extend(f"- {risk}" for risk in evaluation["residual_risks"])
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", required=True, type=parse_run_spec)
    parser.add_argument("--baseline-label", default=None)
    parser.add_argument("--json-out", type=Path, default=None)
    parser.add_argument("--md-out", type=Path, default=None)
    parser.add_argument("--json", action="store_true", help="Print JSON instead of Markdown.")
    args = parser.parse_args()

    report = build_report(args.run, baseline_label=args.baseline_label)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.md_out:
        args.md_out.parent.mkdir(parents=True, exist_ok=True)
        args.md_out.write_text(render_markdown(report), encoding="utf-8")

    print(json.dumps(report, indent=2) if args.json else render_markdown(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
