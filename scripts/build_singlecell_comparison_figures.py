#!/usr/bin/env python3
"""Build a compact vector figure bundle from a run-comparison JSON report."""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


PALETTE = {
    "baseline": "#4C78A8",
    "scrublet": "#4C78A8",
    "scdblfinder": "#F58518",
    "consensus": "#54A24B",
    "other": "#8E6C8A",
}


def _color_for(label: str, method: str | None = None) -> str:
    text = f"{label} {method or ''}".lower()
    if "consensus" in text:
        return PALETTE["consensus"]
    if "scdbl" in text:
        return PALETTE["scdblfinder"]
    if "scrublet" in text or "baseline" in text:
        return PALETTE["scrublet"]
    return PALETTE["other"]


def _clean_label(label: str) -> str:
    return label.replace("_", "\n")


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return default if math.isnan(number) else number


def _save_all(fig: plt.Figure, out_dir: Path, stem: str) -> list[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    for suffix in (".svg", ".pdf", ".png"):
        path = out_dir / f"{stem}{suffix}"
        fig.savefig(path, bbox_inches="tight")
        written.append(str(path))
    plt.close(fig)
    return written


def _style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)


def _annotate_bars(ax: plt.Axes, values: list[float], fmt: str = "{:.1f}") -> None:
    ymax = max(values) if values else 0.0
    offset = ymax * 0.03 if ymax else 0.1
    for idx, value in enumerate(values):
        ax.text(idx, value + offset, fmt.format(value), ha="center", va="bottom", fontsize=8)


def figure_doublet_summary(report: dict[str, Any], out_dir: Path) -> list[str]:
    runs = report["runs"]
    labels = [run["label"] for run in runs]
    colors = [_color_for(run["label"], run["doublet"].get("method")) for run in runs]
    x = range(len(runs))

    rates = [_as_float(run["doublet"].get("rate_pct")) for run in runs]
    kept = [_as_float(run.get("cells_after_doublet_removal")) / 1000.0 for run in runs]
    runtimes = [_as_float(run["doublet"].get("runtime_sec")) for run in runs]

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2), constrained_layout=True)
    panels = [
        (axes[0], rates, "Doublet rate (%)", "{:.2f}", 6.0),
        (axes[1], kept, "Cells retained (thousand)", "{:.1f}", None),
        (axes[2], runtimes, "Doublet module runtime (s)", "{:.0f}", None),
    ]
    for ax, values, ylabel, fmt, expected in panels:
        ax.bar(x, values, color=colors, width=0.68)
        ax.set_xticks(list(x), [_clean_label(label) for label in labels], fontsize=8)
        ax.set_ylabel(ylabel, fontsize=9)
        if expected is not None:
            ax.axhline(expected, color="#6B7280", linestyle="--", linewidth=1.0)
            ax.text(len(runs) - 0.5, expected, "expected 6%", fontsize=8, va="bottom", ha="right")
        _style_axes(ax)
        _annotate_bars(ax, values, fmt)
    fig.suptitle("Doublet backend comparison on LUSC tissue atlas", fontsize=12, fontweight="bold")
    return _save_all(fig, out_dir, "figure_A_doublet_backend_summary")


def figure_downstream_stability(report: dict[str, Any], out_dir: Path) -> list[str]:
    runs = report["runs"]
    labels = [run["label"] for run in runs]
    colors = [_color_for(run["label"], run["doublet"].get("method")) for run in runs]
    x = range(len(runs))
    comparison_by_label = {cmp["label"]: cmp for cmp in report.get("comparisons", [])}

    clusters = [_as_float(run.get("n_clusters")) for run in runs]
    de_genes = [_as_float(run.get("de_significant_genes")) / 1000.0 for run in runs]
    unknown = [_as_float(run.get("annotation_unknown_pct")) for run in runs]
    pathway_overlap = [
        1.0 if run["label"] == report["baseline_label"]
        else _as_float(comparison_by_label.get(run["label"], {}).get("pathway_overlap", {}).get("jaccard"))
        for run in runs
    ]

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.0), constrained_layout=True)
    panels = [
        (axes[0, 0], clusters, "Leiden clusters", "{:.0f}"),
        (axes[0, 1], de_genes, "Significant DE genes (thousand)", "{:.1f}"),
        (axes[1, 0], unknown, "Annotation unknown (%)", "{:.2f}"),
        (axes[1, 1], pathway_overlap, "Top pathway overlap vs baseline", "{:.2f}"),
    ]
    for ax, values, ylabel, fmt in panels:
        ax.bar(x, values, color=colors, width=0.68)
        ax.set_xticks(list(x), [_clean_label(label) for label in labels], fontsize=8)
        ax.set_ylabel(ylabel, fontsize=9)
        _style_axes(ax)
        _annotate_bars(ax, values, fmt)
        if "overlap" in ylabel:
            ax.set_ylim(0, 1.08)
    fig.suptitle("Downstream stability after doublet rescue", fontsize=12, fontweight="bold")
    return _save_all(fig, out_dir, "figure_B_downstream_stability_summary")


def _accepted_run(report: dict[str, Any]) -> dict[str, Any]:
    accepted = report.get("evaluation", {}).get("accepted_lane")
    if accepted:
        for run in report["runs"]:
            if run["label"] == accepted:
                return run
    return report["runs"][-1]


def figure_biology_snapshot(report: dict[str, Any], out_dir: Path) -> list[str]:
    run = _accepted_run(report)
    pathways = run.get("top_pathways", [])[:6]
    lrs = run.get("top_ligand_receptors", [])[:6]
    tf_counts = Counter(run.get("top_tfs", [])[:8])

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)

    path_labels = [row.get("term", "").replace("HALLMARK_", "").replace("_", " ") for row in pathways]
    path_values = [
        -math.log10(max(_as_float(row.get("padj"), 1.0), 1e-300))
        for row in pathways
    ]
    axes[0].barh(range(len(path_values)), path_values, color="#4C78A8")
    axes[0].set_yticks(range(len(path_labels)), path_labels, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("-log10 adjusted P", fontsize=9)
    axes[0].set_title("Pathway enrichment", fontsize=10)
    _style_axes(axes[0])

    lr_labels = [
        f"{row.get('ligand', '')}->{row.get('receptor', '')}\n{row.get('source', '')} to {row.get('target', '')}"
        for row in lrs
    ]
    lr_values = [_as_float(row.get("lr_score")) for row in lrs]
    axes[1].barh(range(len(lr_values)), lr_values, color="#F58518")
    axes[1].set_yticks(range(len(lr_labels)), lr_labels, fontsize=7)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("LR score", fontsize=9)
    axes[1].set_title("Cell communication", fontsize=10)
    _style_axes(axes[1])

    tf_labels = list(tf_counts.keys()) or run.get("top_tfs", [])[:6]
    tf_values = [tf_counts.get(tf, 1) for tf in tf_labels]
    axes[2].barh(range(len(tf_values)), tf_values, color="#54A24B")
    axes[2].set_yticks(range(len(tf_labels)), tf_labels, fontsize=8)
    axes[2].invert_yaxis()
    axes[2].set_xlabel("Listed in top TF set", fontsize=9)
    axes[2].set_title("GRN top TFs", fontsize=10)
    _style_axes(axes[2])

    fig.suptitle(f"Biology signal snapshot: {run['label']}", fontsize=12, fontweight="bold")
    return _save_all(fig, out_dir, "figure_C_biology_signal_snapshot")


def write_manual_checklist(report: dict[str, Any], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    accepted = report.get("evaluation", {}).get("accepted_lane") or "none"
    path = out_dir / "manual_visual_review_checklist.md"
    lines = [
        "# Manual Visual Review Checklist",
        "",
        f"- Accepted evidence lane: `{accepted}`",
        "- Mechanical QA target: `scripts/validate_figure_outputs.py --require-vector`",
        "- Scope: method-comparison summary figures for Round9; not a final full biological figure set.",
        "",
        "## Required Human Checks",
        "",
        "- [ ] Scientific conclusion matches the comparison report.",
        "- [ ] Labels and legends are readable at manuscript scale.",
        "- [ ] No labels, bars, legends, or panel titles overlap.",
        "- [ ] Palette separates baseline, scDblFinder, and consensus in print.",
        "- [ ] Vector PDF/SVG files open cleanly in an external viewer.",
        "- [ ] Debug PNGs from module runs are not presented as final manuscript panels.",
        "- [ ] Pseudobulk limitation is stated wherever DE claims are made.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def build_bundle(report: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    written: list[str] = []
    written.extend(figure_doublet_summary(report, out_dir))
    written.extend(figure_downstream_stability(report, out_dir))
    written.extend(figure_biology_snapshot(report, out_dir))
    checklist = write_manual_checklist(report, out_dir)
    return {
        "schema_version": "singlecell-comparison-figure-bundle/v1",
        "source_report": report.get("generated_at"),
        "out_dir": str(out_dir),
        "figures": written,
        "manual_checklist": str(checklist),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-json", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--manifest-out", type=Path, default=None)
    parser.add_argument("--json", action="store_true", help="Print bundle manifest as JSON.")
    args = parser.parse_args()

    report = json.loads(args.comparison_json.read_text(encoding="utf-8"))
    manifest = build_bundle(report, args.out_dir)
    if args.manifest_out:
        args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
        args.manifest_out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2) if args.json else "\n".join(manifest["figures"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
