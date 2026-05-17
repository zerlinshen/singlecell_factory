"""Wave-5.5 v5.0.4 Stage C5 — three-layer comparison PDF assembly.

# STATUS: one-off, promote-on-reuse

Per plan v5.0.4 AC-V5-PDF-1/2/3 + AC-V5-COMP-LAYER1/PARAM-1/QUANT-1 + AC-V5-LICENSE-1:

- Cover page: overall reproduction-rate, NOT-RUN list, ack-file pointer, license disclaimer.
- Per-figure page: Layer 1 (Trevino panel | our render) + Layer 2 (parameter divergence table)
                   + Layer 3 (quantitative metrics + per-fig reproduction-rate).
- Footer on every page: 'INTERNAL USE ONLY — Trevino et al. 2021 panels reproduced under
  fair-use research exception; do not distribute'.
- Supplementary F-5 has its own divider page.
- Output: content-addressed `comparison_<sha>.pdf` + `comparison.pdf` symlink + `NOTICE.txt`.
  No page-count / MB cap (per [U-Q5]).

Flags:
- `--embed-cell-figures` (default OFF; placeholder text used otherwise)
  Per AC-V5-LICENSE-1, embedding requires (a) opt-in flag, AND (b) a non-expired
  acknowledgment file at .omc/research/wave5/CELL_FIGURE_USE_ACKNOWLEDGMENT_zerlinshen_2026-05-17.md.

Usage:
    python3 scripts/dev/wave5_assemble_comparison_pdf.py [--embed-cell-figures]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.image import imread

REPO_ROOT = Path(__file__).resolve().parents[2]
CANONICAL_RUN = Path("/home/zerlinshen/projects/wave5-trevino/runs/20260516T0931Z-d192836f1bb0")
FIG_ROOT = CANONICAL_RUN / "python" / "figures" / "v5"
PARAM_DIV_DIR = FIG_ROOT / "parameter_divergence"
QUANT_DIR = FIG_ROOT / "quantitative_metrics"
REPRO_PATH = FIG_ROOT / "reproduction_rate.json"
TREVINO_DIR = REPO_ROOT / "data" / "external" / "trevino_2021" / "figures_extracted"
TREVINO_MAPPING = TREVINO_DIR / "MAPPING.md"
ACK_PATH = REPO_ROOT / ".omc" / "research" / "wave5" / "CELL_FIGURE_USE_ACKNOWLEDGMENT_zerlinshen_2026-05-17.md"
ACK_EXPIRY_ISO = "2026-08-15T08:38:00Z"

FOOTER = ("INTERNAL USE ONLY — Trevino et al. 2021 panels reproduced under fair-use "
          "research exception; do not distribute")

# Figure-ID → (display name, our render path, Trevino source page, NOT-RUN/supp flag)
FIGURES = [
    ("F-1B",   "Figure 1B — UMAP × cell type",                  FIG_ROOT/"main/F-1B/F-1B.png",   TREVINO_DIR/"page-04.png", "main"),
    ("F-1C",   "Figure 1C — Marker gene overlay",               FIG_ROOT/"main/F-1C/F-1C.png",   TREVINO_DIR/"page-04.png", "main"),
    ("F-1D",   "Figure 1D — Cell-type composition × age",       FIG_ROOT/"main/F-1D/F-1D.png",   TREVINO_DIR/"page-04.png", "main"),
    ("F-2A",   "Figure 2A — ATAC LSI",                          FIG_ROOT/"main/F-2A/F-2A.png",   TREVINO_DIR/"page-05.png", "main"),
    ("F-2BD",  "Figure 2B-D — ATAC peak landscape",             FIG_ROOT/"main/F-2BD/F-2BD.png", TREVINO_DIR/"page-05.png", "main"),
    ("F-3",    "Figure 3 — chromVAR TF motif (NOT-RUN)",        FIG_ROOT/"main/F-3/F-3.png",     TREVINO_DIR/"page-07.png", "not_run"),
    ("F-4A",   "Figure 4A — Trajectory inference",              FIG_ROOT/"main/F-4A/F-4A.png",   TREVINO_DIR/"page-09.png", "main"),
    ("F-4D",   "Figure 4D — UMAP × pseudotime",                 FIG_ROOT/"main/F-4D/F-4D.png",   TREVINO_DIR/"page-07.png", "main"),
    ("F-4EG",  "Figure 4E-G — Branch-specific gene dynamics",   FIG_ROOT/"main/F-4EG/F-4EG.png", TREVINO_DIR/"page-09.png", "main"),
    ("F-6",    "Figure 6 — TF-gene network (NOT-RUN)",          FIG_ROOT/"main/F-6/F-6.png",     TREVINO_DIR/"page-12.png", "not_run"),
    ("F-7",    "Figure 7 — Disease GWAS overlay (NOT-RUN)",     FIG_ROOT/"main/F-7/F-7.png",     TREVINO_DIR/"page-13.png", "not_run"),
    ("F-5",    "Figure 5 (supp) — Peak-gene linkage",           FIG_ROOT/"supplementary/F-5/F-5.png", TREVINO_DIR/"page-05.png", "supplementary"),
    ("F-RNA-1B", "Stream A2 RNA-only — UMAP × Leiden",         FIG_ROOT/"rna_only/F-RNA-1B/F-RNA-1B.png", TREVINO_DIR/"page-04.png", "rna_only"),
    ("F-RNA-1C", "Stream A2 RNA-only — Marker overlay",        FIG_ROOT/"rna_only/F-RNA-1C/F-RNA-1C.png", TREVINO_DIR/"page-04.png", "rna_only"),
    ("F-RNA-4A", "Stream A2 RNA-only — Trajectory",            FIG_ROOT/"rna_only/F-RNA-4A/F-RNA-4A.png", TREVINO_DIR/"page-07.png", "rna_only"),
]


def sha256_file(p: Path) -> str:
    h = hashlib.sha256(); h.update(p.read_bytes()); return h.hexdigest()


def check_ack(embed_flag: bool) -> tuple[bool, str]:
    """Return (allow_embed, message)."""
    if not embed_flag:
        return False, "embed disabled (--embed-cell-figures flag NOT set; placeholders used)"
    if not ACK_PATH.exists():
        return False, f"acknowledgment file missing at {ACK_PATH}; placeholders used"
    now = datetime.now(timezone.utc)
    expiry = datetime.fromisoformat(ACK_EXPIRY_ISO.replace("Z", "+00:00"))
    if now > expiry:
        return False, f"acknowledgment expired at {ACK_EXPIRY_ISO}; placeholders used"
    return True, f"ack OK (expires {ACK_EXPIRY_ISO}); Cell figures embedded"


def draw_footer(ax, page_num: int, total_pages_estimate: int) -> None:
    ax.text(0.5, 0.005, f"{FOOTER}  ·  page {page_num} / {total_pages_estimate}",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=6,
            color="#666666", style="italic")


def add_layer1_image_panel(ax, image_path: Path, label: str, *, placeholder_text: str = "",
                           allow_embed: bool = True):
    """Render Layer-1 figure panel (either embedded image or placeholder)."""
    if not allow_embed and "Trevino" in label:
        ax.text(0.5, 0.5,
                f"[PLACEHOLDER — Cell figure not embedded]\n"
                f"{placeholder_text}\n\nRun with --embed-cell-figures + valid ack to embed.",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=9, color="#888888",
                bbox=dict(boxstyle="round", facecolor="#EEEEEE", edgecolor="#AAAAAA"))
        ax.set_title(label, fontsize=10, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        return
    try:
        img = imread(image_path)
        ax.imshow(img)
    except Exception as exc:
        ax.text(0.5, 0.5, f"[image not available: {image_path.name}]\n{exc}",
                transform=ax.transAxes, ha="center", va="center", fontsize=8, color="red")
    ax.set_title(label, fontsize=10, fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([])


def add_layer2_param_table(ax, fig_id: str):
    p = PARAM_DIV_DIR / f"{fig_id}.json"
    if not p.exists():
        ax.text(0.5, 0.5, f"(no parameter divergence file for {fig_id})",
                transform=ax.transAxes, ha="center", va="center", fontsize=8, color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("Layer 2 — Parameter divergence (Trevino Methods vs our v5.0)", fontsize=9)
        return
    data = json.loads(p.read_text())
    if data.get("status") == "NOT-RUN":
        ax.text(0.5, 0.5, f"NOT-RUN: {data.get('_description','')}",
                transform=ax.transAxes, ha="center", va="center", fontsize=9, color="#CC3333")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("Layer 2 — Parameter divergence", fontsize=9)
        return
    rows = data.get("rows", [])
    headers = ["Parameter", "Trevino", "Our v5.0", "Flag"]
    cell_text = []
    cell_colors = []
    for r in rows[:10]:  # cap at 10 rows per page; further rows truncated for readability
        trev = str(r.get("trevino_value", ""))[:35]
        ours = str(r.get("our_v5_value", ""))[:35]
        flag = r.get("flag", "?")
        cell_text.append([r["parameter"][:25], trev, ours, flag])
        if flag == "OK":
            color = "#D4EDDA"
        elif flag == "WITHIN-TOLERANCE":
            color = "#FFF3CD"
        else:
            color = "#F8D7DA"
        cell_colors.append(["#FFFFFF", "#FFFFFF", "#FFFFFF", color])
    ax.axis("off")
    if cell_text:
        tbl = ax.table(cellText=cell_text, colLabels=headers, cellLoc="left",
                       loc="upper left", cellColours=cell_colors)
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(7)
        tbl.scale(1, 1.2)
    ax.set_title("Layer 2 — Parameter divergence (Trevino Methods vs our v5.0)", fontsize=9, loc="left")


def add_layer3_metrics(ax, fig_id: str, repro_data: dict):
    p = QUANT_DIR / f"{fig_id}.json"
    if not p.exists():
        ax.text(0.5, 0.5, f"(no quantitative metrics for {fig_id})",
                transform=ax.transAxes, ha="center", va="center", fontsize=8, color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("Layer 3 — Quantitative similarity + per-figure reproduction-rate", fontsize=9)
        return
    data = json.loads(p.read_text())
    status = data.get("status", "computed")
    ax.axis("off")
    per_fig = data.get("per_figure_score")
    n_scored = data.get("n_scored_metrics", 0)
    n_decl = data.get("n_total_metrics_declared", 0)
    lines = []
    if status == "NOT-RUN":
        lines.append(f"NOT-RUN: {data.get('rationale', '')}")
    else:
        for m, raw in data.get("metrics_raw", {}).items():
            scored = data.get("metrics_scored", {}).get(m)
            rationale = data.get("rationale", {}).get(m, "")
            if raw is None:
                lines.append(f"  • {m}: NOT-COMPUTED  ({rationale[:75]}{'…' if len(rationale)>75 else ''})")
            else:
                lines.append(f"  • {m}: raw={raw:.4f}, scored={scored:.4f}  ({rationale[:60]}{'…' if len(rationale)>60 else ''})")
        score_str = f"{per_fig:.4f}" if per_fig is not None else "None (insufficient metric coverage)"
        lines.append("")
        lines.append(f"PER-FIGURE REPRODUCTION-RATE: {score_str}    (scored {n_scored}/{n_decl} declared metrics)")
    txt = "\n".join(lines)
    ax.text(0.0, 0.95, txt, transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace")
    ax.set_title("Layer 3 — Quantitative similarity + per-figure reproduction-rate", fontsize=9, loc="left")


def render_cover_page(pdf: PdfPages, repro_data: dict, embed_msg: str, total_pages: int):
    fig = plt.figure(figsize=(8.5, 11))
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.05, 0.05, 0.90, 0.92])
    ax.axis("off")
    ax.text(0.5, 0.96, "Wave-5.5 v5.0.4 — Trevino 2021 Reproduction Comparison",
            transform=ax.transAxes, ha="center", va="top", fontsize=16, weight="bold")
    ax.text(0.5, 0.92, "Three-layer comparison: figure-pair + parameter divergence + quantitative similarity",
            transform=ax.transAxes, ha="center", va="top", fontsize=10, style="italic")
    overall = repro_data.get("overall")
    overall_str = f"{overall:.4f}" if overall is not None else "n/a"
    ax.text(0.5, 0.85,
            f"Overall reproduction-rate (Stream A main + Stream A2;\nexcludes supplementary F-5 + NOT-RUN F-3/F-6/F-7 from numerator AND denominator):\n\n{overall_str}",
            transform=ax.transAxes, ha="center", va="top", fontsize=12,
            bbox=dict(boxstyle="round", facecolor="#F0F0F0"))
    # NOT-RUN list
    ax.text(0.05, 0.70, "NOT-RUN figures (out of v4.2 pipeline scope):",
            transform=ax.transAxes, ha="left", va="top", fontsize=11, weight="bold")
    not_run = [
        "F-3 — chromVAR TF motif scoring (Wave-6 scheduled)",
        "F-6 — TF-gene regulatory network (Wave-6 scheduled)",
        "F-7 — Disease GWAS overlay (Wave-6 scheduled)",
    ]
    for i, nr in enumerate(not_run):
        ax.text(0.07, 0.66 - i*0.025, f"• {nr}", transform=ax.transAxes,
                ha="left", va="top", fontsize=9)
    # supplementary
    supp = repro_data.get("supplementary_score")
    supp_str = f"{supp:.4f}" if supp is not None else "n/a"
    ax.text(0.05, 0.55, f"Supplementary F-5 reproduction-rate (separate; not in overall): {supp_str}",
            transform=ax.transAxes, ha="left", va="top", fontsize=10)

    # per-figure breakdown
    ax.text(0.05, 0.50, "Per-figure reproduction-rate (n_scored/n_declared shown):",
            transform=ax.transAxes, ha="left", va="top", fontsize=11, weight="bold")
    per_fig = repro_data.get("per_figure", {})
    per_fig_cov = repro_data.get("per_figure_n_scored_over_n_declared", {})
    y = 0.46
    for fid, score in per_fig.items():
        s = f"{score:.3f}" if score is not None else " n/a "
        cov = per_fig_cov.get(fid, "?/?")
        ax.text(0.07, y, f"  {fid:<10}: {s}   ({cov} metrics)",
                transform=ax.transAxes, ha="left", va="top", fontsize=8, family="monospace")
        y -= 0.018

    # limitations
    ax.text(0.05, 0.18, "Limitations:",
            transform=ax.transAxes, ha="left", va="top", fontsize=10, weight="bold")
    for i, lim in enumerate(repro_data.get("limitations", [])):
        ax.text(0.07, 0.155 - i*0.025, f"• {lim}",
                transform=ax.transAxes, ha="left", va="top", fontsize=7, wrap=True)

    # license + ack pointer + embed status
    ax.text(0.05, 0.07,
            f"Acknowledgment: {ACK_PATH.relative_to(REPO_ROOT)}  (expires {ACK_EXPIRY_ISO})\n"
            f"Embed status: {embed_msg}\n"
            f"License: §107 fair-use research exception — internal scientific review only; not for distribution.",
            transform=ax.transAxes, ha="left", va="top", fontsize=7, family="monospace",
            bbox=dict(boxstyle="round", facecolor="#FFFDDB", edgecolor="#CC3333"))
    draw_footer(ax, 1, total_pages)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)


def render_supp_divider(pdf: PdfPages, page_num: int, total_pages: int):
    fig = plt.figure(figsize=(8.5, 11))
    ax = fig.add_axes([0.05, 0.05, 0.90, 0.92])
    ax.axis("off")
    ax.text(0.5, 0.55,
            "— Supplementary —\nF-5 peak-gene linkage (EVIDENCE-ONLY)\n\n"
            "Per plan v5.0.4 [B4]: F-5 is demoted to supplementary tier.\n"
            "Not included in overall reproduction-rate numerator/denominator.",
            transform=ax.transAxes, ha="center", va="center", fontsize=14,
            bbox=dict(boxstyle="round", facecolor="#E6F0FF", edgecolor="#3366CC"))
    draw_footer(ax, page_num, total_pages)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)


def render_figure_page(pdf: PdfPages, fig_id: str, display_name: str, our_path: Path,
                       trevino_path: Path, classification: str, repro_data: dict,
                       allow_embed: bool, page_num: int, total_pages: int):
    fig = plt.figure(figsize=(8.5, 11))
    fig.patch.set_facecolor("white")
    # Layer 1 (top): two image panels side-by-side (heights ~0.40 of page)
    gs = fig.add_gridspec(3, 2, height_ratios=[2.5, 1.5, 1.5], hspace=0.30, wspace=0.05,
                           top=0.94, bottom=0.04, left=0.05, right=0.95)
    ax_trev = fig.add_subplot(gs[0, 0])
    ax_ours = fig.add_subplot(gs[0, 1])
    add_layer1_image_panel(ax_trev, trevino_path,
                            label=f"Trevino 2021 — source page (mapping {trevino_path.name})",
                            placeholder_text=f"Source: {trevino_path.relative_to(REPO_ROOT)}",
                            allow_embed=allow_embed)
    add_layer1_image_panel(ax_ours, our_path,
                            label=f"Our v5.0 render — {display_name}",
                            placeholder_text=f"Source: {our_path.relative_to(CANONICAL_RUN)}",
                            allow_embed=True)
    # Layer 2 (middle): parameter divergence table
    ax_l2 = fig.add_subplot(gs[1, :])
    add_layer2_param_table(ax_l2, fig_id)
    # Layer 3 (bottom): quantitative metrics
    ax_l3 = fig.add_subplot(gs[2, :])
    add_layer3_metrics(ax_l3, fig_id, repro_data)
    # title bar
    fig.suptitle(f"{fig_id} — {display_name} [{classification}]", fontsize=12, weight="bold", y=0.98)
    # footer (transparent overlay axis)
    foot_ax = fig.add_axes([0, 0, 1, 0.025])
    foot_ax.axis("off")
    draw_footer(foot_ax, page_num, total_pages)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)


def assemble(allow_embed: bool, embed_msg: str) -> Path:
    repro_data = json.loads(REPRO_PATH.read_text())
    # initial total-page estimate (cover + 1 supp divider + n figures); refined re-stamp not implemented
    total_pages = 1 + 1 + len(FIGURES)

    # build to temp file first, then content-address
    tmp_pdf = FIG_ROOT / "comparison_tmp.pdf"
    with PdfPages(tmp_pdf) as pdf:
        render_cover_page(pdf, repro_data, embed_msg, total_pages)
        page_num = 2
        # main figures first (in canonical order), then NOT-RUN, then supp divider + F-5, then RNA-only
        # Use the order in FIGURES which has main → not_run → supp → rna_only
        in_supp_block = False
        for fig_id, display, our_path, trev_path, klass in FIGURES:
            if klass == "supplementary" and not in_supp_block:
                render_supp_divider(pdf, page_num, total_pages)
                page_num += 1
                in_supp_block = True
            render_figure_page(pdf, fig_id, display, our_path, trev_path, klass,
                                repro_data, allow_embed, page_num, total_pages)
            page_num += 1

    # content-address by SHA-256 of the produced PDF
    pdf_sha = sha256_file(tmp_pdf)
    final_pdf = FIG_ROOT / f"comparison_{pdf_sha[:16]}.pdf"
    tmp_pdf.replace(final_pdf)

    # symlink comparison.pdf -> comparison_<sha>.pdf
    symlink = FIG_ROOT / "comparison.pdf"
    if symlink.is_symlink() or symlink.exists():
        old_target = os.readlink(symlink) if symlink.is_symlink() else None
        symlink.unlink()
    else:
        old_target = None
    symlink.symlink_to(final_pdf.name)

    # append to symlink_history file
    hist_path = FIG_ROOT / "comparison_pdf_symlink_history.json"
    if hist_path.exists():
        hist = json.loads(hist_path.read_text())
    else:
        hist = []
    hist.append({
        "old_target": old_target,
        "new_target": final_pdf.name,
        "swap_iso": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pdf_sha256": pdf_sha,
    })
    hist_path.write_text(json.dumps(hist, indent=2) + "\n")

    # NOTICE.txt sibling
    (FIG_ROOT / "NOTICE.txt").write_text(
        FOOTER + "\n\n"
        f"Acknowledgment: {ACK_PATH}\nExpiry: {ACK_EXPIRY_ISO}\n"
        f"PDF SHA-256: {pdf_sha}\n"
        f"Plan: wave5-completion-consensus-2026-05-17-v5.0.4.md\n"
        + ("Cell figures embedded under fair-use opt-in.\n" if allow_embed
           else "Cell figures NOT embedded (placeholders shown — use --embed-cell-figures + valid ack).\n")
    )
    print(f"PDF SHA-256: {pdf_sha}")
    print(f"PDF: {final_pdf}")
    print(f"Symlink: {symlink} -> {final_pdf.name}")
    return final_pdf


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--embed-cell-figures", action="store_true",
                    help="Per AC-V5-LICENSE-1: opt-in to embed Trevino panels (default placeholders)")
    args = ap.parse_args(argv)
    allow_embed, msg = check_ack(args.embed_cell_figures)
    print(f"embed decision: {msg}")
    assemble(allow_embed, msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
