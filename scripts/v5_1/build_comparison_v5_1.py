"""Build a Wave-5.1 Trevino Cell reproduction comparison PDF.

This v5.1 builder is intentionally evidence-driven: it reads the project-run
artifacts produced by HV1/HV2/MV1a/MV3 and renders a reviewer-facing internal
PDF without embedding the Cell source figures by default.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

FOOTER = (
    "INTERNAL USE ONLY — Trevino et al. 2021 reproduction evidence; "
    "do not distribute outside authorized research review"
)

FIGURE_ORDER = [
    "F-1B", "F-1C", "F-1D", "F-2A", "F-2BD", "F-3", "F-4A", "F-4D",
    "F-4EG", "F-5", "F-6", "F-7", "F-RNA-1B", "F-RNA-1C", "F-RNA-4A",
]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _wrap(text: str, width: int = 96) -> str:
    lines: list[str] = []
    for raw in str(text).splitlines() or [""]:
        if not raw.strip():
            lines.append("")
        else:
            lines.extend(textwrap.wrap(raw, width=width, replace_whitespace=False))
    return "\n".join(lines)


def _footer(ax, page: int, total: int) -> None:
    ax.text(0.5, 0.015, f"{FOOTER} · page {page}/{total}",
            ha="center", va="bottom", transform=ax.transAxes,
            fontsize=6.5, color="#666666", style="italic")


def _text_page(pdf: PdfPages, *, title: str, body: str, page: int, total: int,
               monospace: bool = False) -> None:
    fig = plt.figure(figsize=(8.5, 11))
    ax = fig.add_axes([0.07, 0.06, 0.86, 0.88])
    ax.axis("off")
    ax.text(0.0, 1.0, title, ha="left", va="top", fontsize=15, weight="bold",
            transform=ax.transAxes)
    ax.text(0.0, 0.94, _wrap(body, 98 if monospace else 92), ha="left", va="top",
            fontsize=8.2 if monospace else 9.0,
            family="monospace" if monospace else "DejaVu Sans",
            transform=ax.transAxes, linespacing=1.25)
    _footer(ax, page, total)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)


def collect(run_dir: Path) -> dict:
    py = run_dir / "python"
    fig = py / "figures" / "v5"
    qdir = fig / "quantitative_metrics"
    metrics = {}
    for path in sorted(qdir.glob("*.json")):
        data = _load_json(path)
        metrics[data.get("fig_id", path.stem)] = data
    return {
        "run_dir": str(run_dir),
        "manifest": _load_json(py / "run_manifest.json"),
        "module_status_path": str(py / "module_status.csv"),
        "reproduction_rate": _load_json(fig / "reproduction_rate.json"),
        "hv2_selected_config": _load_json(fig / "rna_only" / "selected_config.json"),
        "mv1a_parity": _load_json(py / "parity" / "parity_vs_v4_2.json"),
        "mv3_disposition": _load_json(fig / "mv3" / "mv3_disposition.json"),
        "metrics": metrics,
    }


def _cover_body(bundle: dict, created_at: str) -> str:
    repro = bundle["reproduction_rate"]
    hist = repro.get("historical_reproduction_rate", {})
    hv2 = bundle["hv2_selected_config"]
    parity = bundle["mv1a_parity"]
    mv3 = bundle["mv3_disposition"]
    return f"""
Created: {created_at}
Run: {bundle['run_dir']}
Project: Trevino Cell 2021 PCW21 multiome reproduction (Wave-5.1)

Headline reproduction-rate:
  v5.0.4 baseline: {hist.get('v5_0_4', repro.get('baseline_v5_0_4'))}
  v5.1 measured:   {hist.get('v5_1', repro.get('overall'))}
  delta:           {hist.get('delta')}

Key evidence gates:
  HV1 real Layer-3 metrics: {len(bundle['metrics'])} metric JSON files; overall={repro.get('overall')}
  HV2 WNN sweep/fix: status={hv2.get('wnn_fix_status')}, ARI={hv2.get('ari_vs_seurat_clusters')}, n_eligible={hv2.get('n_eligible')}
  MV1a processed-count replay parity: {parity.get('parity_verdict')}
  MV3 F-3/F-6/F-7 disposition: {mv3.get('posture')}

Important interpretation:
  This PDF is an internal evidence package. It does not embed Cell source panels and must not be distributed externally. Deferred figures are explicit and are not counted as reproduced.
""".strip()


def _toc_body(bundle: dict) -> str:
    lines = [
        "1. Cover and headline metrics",
        "2. Evidence summary (HV1/HV2/MV1a/MV3)",
        "3. Per-figure metric pages in canonical Trevino order",
        "4. Limitations and follow-up gates",
        "",
        "Per-figure pages:",
    ]
    for idx, fig_id in enumerate(FIGURE_ORDER, start=1):
        data = bundle["metrics"].get(fig_id)
        if data is None:
            lines.append(f"  {idx:02d}. {fig_id}: MISSING metric JSON")
        else:
            score = data.get("per_figure_score")
            score_str = f"{score:.4f}" if isinstance(score, (int, float)) else "n/a"
            lines.append(
                f"  {idx:02d}. {fig_id}: status={data.get('status', 'computed')}, "
                f"score={score_str}, metrics={data.get('n_scored_metrics')}/{data.get('n_total_metrics_declared')}"
            )
    return "\n".join(lines)


def _evidence_body(bundle: dict) -> str:
    repro = bundle["reproduction_rate"]
    limitations = "\n".join(f"  - {x}" for x in repro.get("limitations", []))
    return f"""
Run manifest:
{json.dumps({k: bundle['manifest'].get(k) for k in ['run_id', 'project', 'repo_sha_at_manifest_write', 'factory_tree_dirty', 'allow_dirty', 'modules_run']}, indent=2)}

HV2 selected configuration:
{json.dumps(bundle['hv2_selected_config'], indent=2)}

MV1a parity:
{json.dumps(bundle['mv1a_parity'], indent=2)}

MV3 disposition:
{json.dumps(bundle['mv3_disposition'].get('figure_status', {}), indent=2)}

Limitations:
{limitations}
""".strip()


def _figure_body(fig_id: str, data: dict | None) -> str:
    if data is None:
        return f"{fig_id}: metric JSON missing. This is a blocking gap for final ledger closure."
    lines = [
        f"Figure: {fig_id}",
        f"Status: {data.get('status', 'computed')}",
        f"Per-figure score: {data.get('per_figure_score')}",
        f"Metric coverage: {data.get('n_scored_metrics')}/{data.get('n_total_metrics_declared')}",
        "",
        "Raw metrics / scores:",
    ]
    raw = data.get("metrics_raw", {})
    scored = data.get("metrics_scored", {})
    rationale = data.get("rationale", {})
    if raw:
        for key in sorted(raw):
            lines.append(f"  - {key}: raw={raw.get(key)} scored={scored.get(key)}")
            lines.append(f"    rationale: {rationale.get(key)}")
    else:
        lines.append(f"  {data.get('rationale', 'No numeric metrics for this figure.')}")
    source_paths = data.get("source_data_paths")
    if source_paths:
        lines.extend(["", "Source data paths:", json.dumps(source_paths, indent=2)])
    return "\n".join(lines)


def build(run_dir: Path) -> Path:
    bundle = collect(run_dir)
    fig_root = run_dir / "python" / "figures" / "v5"
    created_at = datetime.now(timezone.utc).isoformat()
    total_pages = 3 + len(FIGURE_ORDER) + 1
    tmp = fig_root / "comparison_tmp_v5.1.pdf"
    with PdfPages(tmp) as pdf:
        _text_page(pdf, title="Wave-5.1 Trevino Cell Reproduction Comparison", body=_cover_body(bundle, created_at), page=1, total=total_pages)
        _text_page(pdf, title="Table of Contents", body=_toc_body(bundle), page=2, total=total_pages, monospace=True)
        _text_page(pdf, title="Evidence Summary", body=_evidence_body(bundle), page=3, total=total_pages, monospace=True)
        page = 4
        for fig_id in FIGURE_ORDER:
            _text_page(pdf, title=f"Per-figure Evidence — {fig_id}", body=_figure_body(fig_id, bundle["metrics"].get(fig_id)), page=page, total=total_pages, monospace=True)
            page += 1
        followup = """
Remaining gates after this PDF:
  - Validate and write the v5.1 run ledger.
  - Preserve v4.2/v5.0.4 artifact immutability.
  - If F-3/F-6/F-7 are required as reproduced figures, run a future approved MV3 lane with genome/motif/GWAS resources and focused module tests.
""".strip()
        _text_page(pdf, title="Limitations and Follow-up", body=followup, page=page, total=total_pages)
    digest = _sha256(tmp)
    final = fig_root / f"comparison_{digest[:16]}_v5.1.pdf"
    tmp.replace(final)
    symlink = fig_root / "comparison_v5.1.pdf"
    old_target = os.readlink(symlink) if symlink.is_symlink() else None
    if symlink.exists() or symlink.is_symlink():
        symlink.unlink()
    symlink.symlink_to(final.name)
    manifest = {
        "schema_version": "v5.1-comparison-pdf",
        "created_at": created_at,
        "pdf_path": str(final),
        "pdf_sha256": digest,
        "symlink": str(symlink),
        "old_symlink_target": old_target,
        "run_dir": str(run_dir),
        "source_artifacts": {
            "reproduction_rate": str(fig_root / "reproduction_rate.json"),
            "hv2_selected_config": str(fig_root / "rna_only" / "selected_config.json"),
            "mv1a_parity": str(run_dir / "python" / "parity" / "parity_vs_v4_2.json"),
            "mv3_disposition": str(fig_root / "mv3" / "mv3_disposition.json"),
        },
    }
    (fig_root / "comparison_v5_1_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return final


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, type=Path)
    args = ap.parse_args(argv)
    build(args.run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
