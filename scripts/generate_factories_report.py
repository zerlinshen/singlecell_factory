#!/usr/bin/env python3
"""Generate FACTORIES_OVERVIEW.{md,pdf} — a self-contained technical reference.

Reads project state from disk (module catalog, recipes, README/docs index, test
counts via ``pytest --collect-only``) and emits a deterministic Markdown +
ReportLab PDF describing both ``singlecell_factory`` and
``r_multiomics_factory``.

Re-run with::

    python scripts/generate_factories_report.py

Stdlib + reportlab (>=4.4) + pyyaml (>=6.0). Deterministic output (the only
non-deterministic line is the ``generated_at`` ISO timestamp on the cover).
"""
from __future__ import annotations

import datetime as _dt
import importlib.util
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import yaml
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

REPO_PY = Path("/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory")
REPO_R = Path("/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory")
DOCS_DIR = REPO_PY / "docs"
OUT_PDF = DOCS_DIR / "FACTORIES_OVERVIEW.pdf"
OUT_MD = DOCS_DIR / "FACTORIES_OVERVIEW.md"

DOC_TITLE = "Factories Overview"
DOC_SUBTITLE = "singlecell_factory + r_multiomics_factory technical reference"
FOOTER_TITLE = "Factories Overview"

PAGE_WIDTH, PAGE_HEIGHT = A4
MARGIN = 2 * cm

TEST_FILES = (
    "tests/test_r_bundle_contract.py",
    "tests/test_python_r_parity.py",
    "tests/test_modular.py",
    "tests/test_scfactory.py",
    "tests/test_singlecell_r_bundle_export.py",
    "tests/test_densify_audit.py",
    "tests/test_module_catalog.py",
)

TEST_PURPOSE = {
    "tests/test_r_bundle_contract.py": (
        "Python <-> R bundle v2/v2.1 round-trip subprocess contract; guards CLAIM_GUARD surfacing in R."
    ),
    "tests/test_python_r_parity.py": (
        "Cross-language sentinels (scanpy vs Seurat) — fires when defaults silently diverge."
    ),
    "tests/test_modular.py": (
        "Modular pipeline contracts: dependency resolution, ctx.set_module_dir, manifest/status outputs."
    ),
    "tests/test_scfactory.py": (
        "scfactory CLI wrapper: auto-detection planner, recipe precedence, doctor JSON shape."
    ),
    "tests/test_singlecell_r_bundle_export.py": (
        "Bundle export internals: schema versions, atomic publish, extension manifest entries."
    ),
    "tests/test_densify_audit.py": (
        "Grep-ban: every .toarray()/.todense() in modules must carry a densify-allowed comment."
    ),
    "tests/test_module_catalog.py": (
        "Module catalog stability: mandatory chain, default optional list, layer coverage."
    ),
}


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def _import_module_catalog():
    """Import workflow.modular.module_catalog without importing the package side-effects."""
    name = "_factories_overview_module_catalog"
    spec = importlib.util.spec_from_file_location(
        name,
        REPO_PY / "workflow" / "modular" / "module_catalog.py",
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load module_catalog.py")
    mod = importlib.util.module_from_spec(spec)
    # Register before exec so @dataclass introspection works on Py3.13.
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@dataclass(frozen=True)
class ModuleRow:
    name: str
    layer: str
    modality: str
    depends_on: str
    bridge_ready: str
    status: str
    description: str


def load_module_rows() -> list[ModuleRow]:
    try:
        cat = _import_module_catalog()
    except Exception as exc:  # pragma: no cover
        warn(f"module catalog import failed: {exc}")
        return []
    rows: list[ModuleRow] = []
    for name, spec in cat.MODULE_SPECS.items():
        is_experimental = "[EXPERIMENTAL]" in (spec.description or "")
        rows.append(
            ModuleRow(
                name=name,
                layer=spec.layer,
                modality=spec.modality,
                depends_on=", ".join(spec.depends_on) if spec.depends_on else "-",
                bridge_ready="yes" if spec.bridge_ready else "no",
                status="experimental" if is_experimental else "stable",
                description=spec.description or "",
            )
        )
    return rows


def load_recipes() -> dict[str, dict]:
    out: dict[str, dict] = {}
    for path in sorted((REPO_PY / "recipes").glob("*.yaml")):
        try:
            with path.open() as fh:
                out[path.stem] = yaml.safe_load(fh) or {}
        except Exception as exc:  # pragma: no cover
            warn(f"recipe parse failed for {path}: {exc}")
            out[path.stem] = {"_error": str(exc)}
    return out


def load_scfactory_help() -> dict[str, str]:
    helps: dict[str, str] = {}
    cmd_root = [sys.executable, str(REPO_PY / "scripts" / "scfactory.py")]
    for sub in ("", "run", "doctor", "report"):
        argv = cmd_root + ([sub] if sub else []) + ["--help"]
        try:
            r = subprocess.run(
                argv, capture_output=True, text=True, timeout=30, cwd=REPO_PY
            )
            helps[sub or "_root"] = (r.stdout or "") + (r.stderr or "")
        except Exception as exc:  # pragma: no cover
            warn(f"scfactory --help failed for '{sub}': {exc}")
            helps[sub or "_root"] = f"(extraction failed: {exc})"
    return helps


def collect_test_counts() -> dict[str, int]:
    """Run pytest --collect-only and tally per-file counts."""
    argv = [
        sys.executable,
        "-m",
        "pytest",
        "--collect-only",
        "-q",
        *TEST_FILES,
    ]
    counts: dict[str, int] = {f: 0 for f in TEST_FILES}
    try:
        r = subprocess.run(
            argv, capture_output=True, text=True, timeout=120, cwd=REPO_PY
        )
        text = (r.stdout or "") + "\n" + (r.stderr or "")
    except Exception as exc:  # pragma: no cover
        warn(f"pytest collect-only failed: {exc}")
        return counts

    for line in text.splitlines():
        line = line.strip()
        if "::" not in line:
            continue
        path = line.split("::", 1)[0]
        if path in counts:
            counts[path] += 1
    return counts


def load_claim_matrix() -> tuple[list[list[str]], dict[str, int], str]:
    """Return (rows_for_table, status_counts, raw_summary_text)."""
    p = REPO_R / "docs" / "NSCLC_CLAIM_VALIDATION_MATRIX.md"
    table_rows: list[list[str]] = []
    counts: dict[str, int] = {"direct": 0, "partial": 0, "unsupported": 0}
    raw = ""
    if not p.exists():
        warn(f"claim matrix missing at {p}")
        return table_rows, counts, raw
    try:
        text = p.read_text()
        raw = text
        # Parse markdown table rows: pipe-delimited lines with claim_id.
        for line in text.splitlines():
            line = line.strip()
            if not line.startswith("|") or "claim_id" in line or "---" in line:
                continue
            cells = [c.strip() for c in line.strip("|").split("|")]
            if len(cells) < 9:
                continue
            claim_id = cells[0]
            family = cells[1]
            module_family = cells[3]
            empirical_lane = cells[4]
            support_status = cells[6]
            evidence = cells[7]
            if not claim_id or not claim_id.startswith("NC2024"):
                continue
            counts[support_status] = counts.get(support_status, 0) + 1
            table_rows.append(
                [
                    claim_id,
                    family,
                    module_family,
                    empirical_lane,
                    support_status,
                    _truncate(evidence, 80),
                ]
            )
    except Exception as exc:  # pragma: no cover
        warn(f"claim matrix parse failed: {exc}")
    return table_rows, counts, raw


def list_doc_index() -> list[tuple[str, str]]:
    """Return sorted (relative_path, one_line_summary) tuples for both repos."""
    entries: list[tuple[str, str]] = []
    targets = [
        ("singlecell_factory", REPO_PY),
        ("r_multiomics_factory", REPO_R),
    ]
    for label, root in targets:
        # top-level *.md
        for p in sorted(root.glob("*.md")):
            entries.append(_doc_entry(label, root, p))
        # docs/*.md
        for p in sorted((root / "docs").glob("*.md")):
            entries.append(_doc_entry(label, root, p))
    entries.sort(key=lambda x: x[0])
    return entries


def _doc_entry(label: str, root: Path, p: Path) -> tuple[str, str]:
    rel = f"{label}/{p.relative_to(root)}"
    summary = "(no summary)"
    try:
        with p.open() as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    summary = line.lstrip("# ").strip()
                    # try to also pull the next non-empty non-heading line if heading is just title-ish
                    continue
                if line.startswith("|") or line.startswith("---"):
                    continue
                summary = line
                break
    except Exception as exc:  # pragma: no cover
        warn(f"doc summary failed for {p}: {exc}")
    return rel, _truncate(summary, 110)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


_WARNINGS: list[str] = []


def warn(msg: str) -> None:
    _WARNINGS.append(msg)
    print(f"[generate_factories_report] WARN: {msg}", file=sys.stderr)


def _truncate(text: str, n: int) -> str:
    text = " ".join((text or "").split())
    if len(text) <= n:
        return text
    return text[: n - 1].rstrip() + "..."


def _esc(text: str) -> str:
    """Escape for ReportLab Paragraph (XML-ish)."""
    return (
        (text or "")
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


# ---------------------------------------------------------------------------
# Markdown emitter
# ---------------------------------------------------------------------------


def render_markdown(
    module_rows: list[ModuleRow],
    recipes: dict[str, dict],
    helps: dict[str, str],
    test_counts: dict[str, int],
    claim_rows: list[list[str]],
    claim_counts: dict[str, int],
    doc_index: list[tuple[str, str]],
    generated_at: str,
) -> str:
    out: list[str] = []
    o = out.append
    o(f"# {DOC_TITLE}")
    o("")
    o(f"_{DOC_SUBTITLE}_")
    o("")
    o(f"generated_at: {generated_at}")
    o("")

    # 1. Executive summary
    o("## 1. Executive summary")
    o("")
    o(
        "`singlecell_factory` is a Python/AnnData heavy-compute pipeline organized as 33 "
        "concrete module files across the layers `ingest`, `quality_control`, `latent_structure`, "
        "`markers`, `annotation`, `state_dynamics`, `biology`, `aggregation`, `reproduction`, "
        "`spatial`, `multimodal`, `genomic_optional`, `external_validation`, `covariates`, and "
        "`reporting`. The catalog declares 29 active module specs (3 mandatory + 26 optional) "
        "with explicit dependencies, modality tags, and `bridge_ready` flags. The CLI entrypoint "
        "is `scfactory run`."
    )
    o("")
    o(
        "`r_multiomics_factory` is the R-side downstream visualization workspace. It consumes "
        "the bundle v2.1 contract emitted by `scripts/export_singlecell_r_bundle.py`, reads it "
        "via `R_bundle/io_bundle.R`, and renders Seurat-based plots via the `R/*_module.R` files "
        "(QC, dim, expression, composition, marker, annotation, batch_integration, integration, "
        "preprocessing, protein, spatial, theme_config). Together: 97 collected tests on the "
        "named contract surface, 4 starter recipes, 3 modality extensions (protein/spatial/"
        "multimodal_obsm), and an NC2024 paper-faithful reproduction lane."
    )
    o("")
    o(
        "Use `scfactory run` for routine pipeline runs, the documented bundle v2.1 contract for "
        "agent-mediated handoff to R, and `--recipe nc2024_paper` for paper-faithful "
        "reproduction."
    )
    o("")

    # 2. Architecture (textual fallback in markdown)
    o("## 2. Architecture")
    o("")
    o(
        "Two side-by-side stacks bridged by the bundle v2.1 contract. Python writes a temp dir, "
        "validates manifest + SHA256 sweep, then `os.replace` swaps it into place; R reads via "
        "`io_bundle.R::read_bundle_v2`, surfaces `CLAIM_GUARD` as both a stderr message and a "
        "matrix attribute, and routes to module plotters."
    )
    o("")
    layers_py = _modules_by_layer(module_rows)
    o("### Python (singlecell_factory) layers")
    o("")
    for layer, names in layers_py.items():
        bridge = sum(
            1 for r in module_rows if r.layer == layer and r.bridge_ready == "yes"
        )
        o(f"- **{layer}** ({len(names)} modules, {bridge} bridge_ready): {', '.join(names)}")
    o("")
    o("### R (r_multiomics_factory) modules")
    o("")
    r_modules = sorted(p.stem for p in (REPO_R / "R").glob("*.R"))
    o(", ".join(r_modules))
    o("")

    # 3. Capability matrix
    o("## 3. Capability matrix")
    o("")
    o(
        "| module | layer | modality | depends_on | bridge_ready | status |"
    )
    o("|---|---|---|---|---|---|")
    for r in module_rows:
        o(
            f"| {r.name} | {r.layer} | {r.modality} | {r.depends_on} | "
            f"{r.bridge_ready} | {r.status} |"
        )
    o("")

    # 4. Bundle contract
    o("## 4. Bundle v2.1 contract")
    o("")
    o(
        "- Schema version string: `singlecell_r_bundle_v2.1`. Back-compat literal: "
        "`singlecell_r_bundle_v2`. The exporter's `ExportConfig.schema_version` accepts "
        "`v1`, `v2`, or `v2.1`."
    )
    o(
        "- Top-level manifest fields: `schema_version`, `source`, `bundle`, `cell_alignment`, "
        "`expression`, `files`, `extensions`."
    )
    o("- Known extension keys (`KNOWN_EXTENSION_KEYS`): `protein`, `spatial`, `multimodal_obsm`.")
    o("    - `protein`: CLR-normalized ADT matrix (`protein.parquet`); fields include "
      "`normalization`, optional `isotype_controls`.")
    o("    - `spatial`: (x,y) coordinates + library metadata (`spatial.parquet`); image data is "
      "STRING-pointer only, never serialized.")
    o("    - `multimodal_obsm` (EXPERIMENTAL): one parquet per joint embedding (`X_wnn`, "
      "`X_mofa`); same plotting-only `claim_guard`.")
    o("")
    o(
        "- `CLAIM_GUARD` literal: "
        "`not_for_de_or_new_quantitative_claims_without_full_object_validation`. "
        "Surfaced in R as a `message()` call AND as `attr(expr_sparse, \"claim_guard\")`. "
        "Meaning: bundle data is plotting-only; do not draw new quantitative claims from it "
        "without full-object validation."
    )
    o(
        "- Atomicity: writers stage every file under a `TemporaryDirectory`; SHA256 is computed "
        "for each, recorded in the manifest, and cross-checked against on-disk byte sizes. The "
        "manifest is written last, then `os.replace(temp_dir, final_output_dir)` performs the "
        "atomic swap; any pre-existing target directory is moved to a sibling backup and "
        "restored on failure."
    )
    o("")

    # 5. Workflow walkthroughs
    o("## 5. Workflow walkthroughs")
    o("")
    for stem in ("nc2024_paper", "quick_explore", "cite_seq_full", "visium_neighborhoods"):
        rec = recipes.get(stem) or {}
        o(f"### {stem}")
        o("")
        o(f"_{rec.get('description', '(no description)')}_  ")
        o("")
        o(f"**Command**: `python scripts/scfactory.py run --recipe {stem} <input>`")
        o("")
        opt = rec.get("optional_modules") or []
        if opt:
            o("Optional modules:")
            for m in opt:
                o(f"- {m}")
            o("")
        env = rec.get("env") or {}
        if env:
            o("Env:")
            for k, v in env.items():
                o(f"- `{k}`={v}")
            o("")
        bundle = rec.get("bundle") or {}
        if bundle:
            o("Bundle config:")
            for k, v in bundle.items():
                o(f"- `{k}`: {v}")
            o("")
        notes = (rec.get("notes") or "").strip()
        if notes:
            o("When to use:")
            o("")
            for line in notes.splitlines():
                if line.strip():
                    o(f"> {line.strip()}")
            o("")

    # 6. CLI reference
    o("## 6. scfactory CLI reference")
    o("")
    for sub in ("_root", "run", "doctor", "report"):
        title = "scfactory" if sub == "_root" else f"scfactory {sub}"
        o(f"### {title}")
        o("")
        o("```")
        o(helps.get(sub, "(no output)").rstrip())
        o("```")
        o("")

    # 7. NC2024 claim coverage
    o("## 7. NC2024 claim coverage")
    o("")
    if claim_rows:
        total = sum(claim_counts.values())
        o(
            f"Source: `r_multiomics_factory/docs/NSCLC_CLAIM_VALIDATION_MATRIX.md`. "
            f"{total} claims parsed: "
            f"{claim_counts.get('direct', 0)} direct, "
            f"{claim_counts.get('partial', 0)} partial, "
            f"{claim_counts.get('unsupported', 0)} unsupported."
        )
        o("")
        o("| claim_id | family | module_family | lane | status | evidence |")
        o("|---|---|---|---|---|---|")
        for row in claim_rows:
            o("| " + " | ".join(row) + " |")
    else:
        o(
            "Claim matrix not parsed; refer directly to "
            "`r_multiomics_factory/docs/NSCLC_CLAIM_VALIDATION_MATRIX.md`."
        )
    o("")

    # 8. Tests
    o("## 8. Test coverage and quality gates")
    o("")
    total_tests = sum(test_counts.values())
    o(f"Total collected on the named contract surface: **{total_tests}** tests.")
    o("")
    o("| file | tests | protects against |")
    o("|---|---:|---|")
    for f in TEST_FILES:
        o(f"| `{f}` | {test_counts.get(f, 0)} | {TEST_PURPOSE.get(f, '(no description)')} |")
    o("")

    # 9. Limitations / roadmap
    o("## 9. Known limitations and roadmap")
    o("")
    o(
        "- `multimodal_integration` is **EXPERIMENTAL**: off by default, gated by "
        "`SC_MULTIMODAL_ENGINE`. Emits the bundle v2.1 `multimodal_obsm` extension only when "
        "explicitly enabled."
    )
    o(
        "- Tissue images in the spatial extension are **STRING-pointer only**. Image bytes "
        "are never serialized into the bundle; readers must resolve paths themselves."
    )
    o(
        "- Two bundle-v1 schema tests are pinned to the v1 literal for back-compat coverage and "
        "are intentionally NOT removed."
    )
    o(
        "- No cell-cell communication / pathway scoring / RNA-velocity / trajectory parity "
        "checks against R yet — listed as roadmap."
    )
    o(
        "- `--include-multimodal-obsm` is opt-in only and is never auto-passed even on the "
        "multimodal modality."
    )
    o(
        "- A pre-existing repo-doc-sync hook on `singlecell_factory` blocks commits until "
        "citation/doc references are staged."
    )
    o(
        "- Follow-up: a dedicated \"How to add a new module\" extension guide is on the roadmap "
        "but not part of this overview."
    )
    o("")

    # Appendix
    o("## Appendix: file index")
    o("")
    o("| path | summary |")
    o("|---|---|")
    for path, summary in doc_index:
        o(f"| `{path}` | {summary} |")
    o("")

    return "\n".join(out) + "\n"


def _modules_by_layer(rows: Iterable[ModuleRow]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for r in rows:
        out.setdefault(r.layer, []).append(r.name)
    return out


# ---------------------------------------------------------------------------
# PDF rendering
# ---------------------------------------------------------------------------


def _make_styles():
    base = getSampleStyleSheet()
    styles = {
        "title": ParagraphStyle(
            "title", parent=base["Title"], fontName="Helvetica-Bold", fontSize=22,
            leading=26, spaceAfter=12,
        ),
        "subtitle": ParagraphStyle(
            "subtitle", parent=base["Normal"], fontName="Helvetica", fontSize=12,
            leading=14, textColor=colors.grey, spaceAfter=18,
        ),
        "h1": ParagraphStyle(
            "h1", parent=base["Heading1"], fontName="Helvetica-Bold", fontSize=14,
            leading=18, spaceBefore=14, spaceAfter=8, textColor=colors.HexColor("#222222"),
        ),
        "h2": ParagraphStyle(
            "h2", parent=base["Heading2"], fontName="Helvetica-Bold", fontSize=12,
            leading=15, spaceBefore=10, spaceAfter=6, textColor=colors.HexColor("#333333"),
        ),
        "body": ParagraphStyle(
            "body", parent=base["Normal"], fontName="Helvetica", fontSize=10,
            leading=13, spaceAfter=6,
        ),
        "small": ParagraphStyle(
            "small", parent=base["Normal"], fontName="Helvetica", fontSize=8.5,
            leading=11,
        ),
        "code": ParagraphStyle(
            "code", parent=base["Normal"], fontName="Courier", fontSize=8,
            leading=10, textColor=colors.HexColor("#222222"),
        ),
        "footer": ParagraphStyle(
            "footer", parent=base["Normal"], fontName="Helvetica", fontSize=8,
            leading=10, textColor=colors.grey,
        ),
    }
    return styles


_TABLE_HEADER_STYLE = TableStyle(
    [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#444444")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#888888")),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]
)


def _zebra_style(n_rows: int) -> TableStyle:
    cmds = list(_TABLE_HEADER_STYLE.getCommands())
    for i in range(1, n_rows):
        bg = colors.HexColor("#f5f5f5") if i % 2 == 1 else colors.white
        cmds.append(("BACKGROUND", (0, i), (-1, i), bg))
    return TableStyle(cmds)


def _wrap_cells(rows: list[list[str]], styles, style_name="small") -> list[list]:
    return [[Paragraph(_esc(str(c)), styles[style_name]) for c in row] for row in rows]


# ---------------------------------------------------------------------------
# Document with footer + TOC tracking
# ---------------------------------------------------------------------------


class _Doc(BaseDocTemplate):
    def __init__(self, filename, **kw):
        super().__init__(filename, **kw)
        frame = Frame(
            MARGIN, MARGIN, PAGE_WIDTH - 2 * MARGIN,
            PAGE_HEIGHT - 2 * MARGIN, id="main", showBoundary=0,
        )
        self.addPageTemplates([
            PageTemplate(id="default", frames=[frame], onPage=_draw_footer)
        ])
        # toc_entries: list of (label, page_number) added at flowable time.
        self.toc_entries: list[tuple[str, int]] = []
        self.total_pages = 0

    def afterFlowable(self, flowable):
        if isinstance(flowable, Paragraph):
            sname = flowable.style.name
            if sname == "h1":
                text = flowable.getPlainText()
                self.toc_entries.append((text, self.page))


def _draw_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.grey)
    text = f"{FOOTER_TITLE} - page {doc.page}"
    canvas.drawCentredString(PAGE_WIDTH / 2.0, 1.0 * cm, text)
    canvas.restoreState()


# ---------------------------------------------------------------------------
# PDF flowables
# ---------------------------------------------------------------------------


def _architecture_flowables(module_rows: list[ModuleRow], styles) -> list:
    layers_py = _modules_by_layer(module_rows)
    rows = [["Python layer", "modules", "bridge_ready"]]
    for layer in layers_py:
        names = layers_py[layer]
        bridge = sum(
            1 for r in module_rows if r.layer == layer and r.bridge_ready == "yes"
        )
        rows.append([layer, str(len(names)), f"{bridge} of {len(names)}"])

    r_modules = sorted(p.stem for p in (REPO_R / "R").glob("*.R"))
    r_rows = [["R module"]] + [[m] for m in r_modules]

    py_table = Table(_wrap_cells(rows, styles), colWidths=[5.5 * cm, 2.0 * cm, 3.0 * cm])
    py_table.setStyle(_zebra_style(len(rows)))
    r_table = Table(_wrap_cells(r_rows, styles), colWidths=[6.0 * cm])
    r_table.setStyle(_zebra_style(len(r_rows)))

    bridge_box_rows = [
        ["Bridge: bundle v2.1"],
        ["schema_version: singlecell_r_bundle_v2.1 (back-compat: ...v2)"],
        ["extensions: protein, spatial, multimodal_obsm"],
        ["claim_guard: not_for_de_or_new_quantitative_claims_..."],
        ["atomicity: temp dir + os.replace + sha256 + size cross-check"],
    ]
    bridge_table = Table(
        _wrap_cells(bridge_box_rows, styles),
        colWidths=[PAGE_WIDTH - 2 * MARGIN],
    )
    bridge_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1f4e79")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
                ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#eaf1f8")),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#888888")),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )

    side_by_side = Table(
        [
            [
                Paragraph("<b>Python (singlecell_factory)</b>", styles["body"]),
                Paragraph("<b>R (r_multiomics_factory)</b>", styles["body"]),
            ],
            [py_table, r_table],
        ],
        colWidths=[(PAGE_WIDTH - 2 * MARGIN) / 2, (PAGE_WIDTH - 2 * MARGIN) / 2],
    )
    side_by_side.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )

    caption = Paragraph(
        "Architecture: Python heavy-compute on the left writes the bundle v2.1 manifest in the "
        "middle; R reads it on the right. Bridge_ready Python modules are the ones whose outputs "
        "are intended for R-side consumption. The bundle contract is the only stable agent "
        "handoff surface between the two repos.",
        styles["small"],
    )

    return [side_by_side, Spacer(1, 0.4 * cm), bridge_table, Spacer(1, 0.3 * cm), caption]


def _capability_table(module_rows: list[ModuleRow], styles) -> Table:
    header = ["module", "layer", "modality", "depends_on", "bridge", "status"]
    rows = [header]
    for r in module_rows:
        rows.append(
            [
                r.name,
                r.layer,
                r.modality,
                r.depends_on,
                r.bridge_ready,
                r.status,
            ]
        )
    col_widths = [3.2 * cm, 2.6 * cm, 3.2 * cm, 4.0 * cm, 1.5 * cm, 2.0 * cm]
    t = Table(_wrap_cells(rows, styles), colWidths=col_widths, repeatRows=1)
    t.setStyle(_zebra_style(len(rows)))
    return t


def _recipe_flowables(stem: str, rec: dict, styles) -> list:
    items: list = []
    items.append(Paragraph(f"5.{stem}", styles["h2"]))
    items.append(
        Paragraph(_esc(rec.get("description", "(no description)")), styles["body"])
    )
    items.append(
        Paragraph(
            f"<b>Command:</b><br/>python scripts/scfactory.py run --recipe {stem} &lt;input&gt;",
            styles["body"],
        )
    )
    opt = rec.get("optional_modules") or []
    if opt:
        rows = [["optional module"]] + [[m] for m in opt]
        t = Table(_wrap_cells(rows, styles), colWidths=[8 * cm])
        t.setStyle(_zebra_style(len(rows)))
        items.append(t)
        items.append(Spacer(1, 0.15 * cm))
    env = rec.get("env") or {}
    if env:
        rows = [["env var", "value"]] + [[k, str(v)] for k, v in env.items()]
        t = Table(_wrap_cells(rows, styles), colWidths=[6 * cm, 6 * cm])
        t.setStyle(_zebra_style(len(rows)))
        items.append(t)
        items.append(Spacer(1, 0.15 * cm))
    bundle = rec.get("bundle") or {}
    if bundle:
        rows = [["bundle field", "value"]] + [[k, str(v)] for k, v in bundle.items()]
        t = Table(_wrap_cells(rows, styles), colWidths=[6 * cm, 6 * cm])
        t.setStyle(_zebra_style(len(rows)))
        items.append(t)
        items.append(Spacer(1, 0.15 * cm))
    notes = (rec.get("notes") or "").strip()
    if notes:
        items.append(
            Paragraph(
                "<b>When to use:</b> "
                + _esc(" ".join(line.strip() for line in notes.splitlines() if line.strip())),
                styles["body"],
            )
        )
    items.append(Spacer(1, 0.3 * cm))
    return items


def _help_block(text: str, styles) -> list:
    # Wrap each line as a Code paragraph to preserve layout while still flowing.
    out: list = []
    for line in text.splitlines() or ["(no output)"]:
        out.append(Paragraph(_esc(line) or "&nbsp;", styles["code"]))
    return out


def _claims_flowables(claim_rows, claim_counts, styles) -> list:
    items: list = []
    if not claim_rows:
        items.append(
            Paragraph(
                "Claim matrix not parsed. Refer directly to "
                "<font face='Courier'>r_multiomics_factory/docs/NSCLC_CLAIM_VALIDATION_MATRIX.md</font>.",
                styles["body"],
            )
        )
        return items
    total = sum(claim_counts.values())
    items.append(
        Paragraph(
            f"Source: <font face='Courier'>r_multiomics_factory/docs/NSCLC_CLAIM_VALIDATION_MATRIX.md</font>. "
            f"<b>{total}</b> claims parsed: <b>{claim_counts.get('direct', 0)}</b> direct, "
            f"<b>{claim_counts.get('partial', 0)}</b> partial, "
            f"<b>{claim_counts.get('unsupported', 0)}</b> unsupported.",
            styles["body"],
        )
    )
    header = ["claim_id", "family", "module_family", "lane", "status", "evidence"]
    rows = [header] + claim_rows
    col_widths = [2.4 * cm, 2.8 * cm, 2.8 * cm, 2.0 * cm, 2.0 * cm, 4.5 * cm]
    t = Table(_wrap_cells(rows, styles), colWidths=col_widths, repeatRows=1)
    t.setStyle(_zebra_style(len(rows)))
    items.append(t)
    return items


def _tests_flowables(test_counts, styles) -> list:
    total = sum(test_counts.values())
    items: list = [
        Paragraph(
            f"Total collected on the named contract surface: <b>{total}</b> tests.",
            styles["body"],
        )
    ]
    rows = [["file", "tests", "protects against"]]
    for f in TEST_FILES:
        rows.append([f, str(test_counts.get(f, 0)), TEST_PURPOSE.get(f, "")])
    col_widths = [6 * cm, 1.5 * cm, 9 * cm]
    t = Table(_wrap_cells(rows, styles), colWidths=col_widths, repeatRows=1)
    t.setStyle(_zebra_style(len(rows)))
    items.append(t)
    return items


def _appendix_flowables(doc_index, styles) -> list:
    rows = [["path", "summary"]]
    for path, summary in doc_index:
        rows.append([path, summary])
    col_widths = [6.5 * cm, 10.0 * cm]
    t = Table(_wrap_cells(rows, styles), colWidths=col_widths, repeatRows=1)
    t.setStyle(_zebra_style(len(rows)))
    return [t]


# ---------------------------------------------------------------------------
# Main render
# ---------------------------------------------------------------------------


SECTION_NAMES = (
    "1. Executive summary",
    "2. Architecture",
    "3. Capability matrix",
    "4. Bundle v2.1 contract",
    "5. Workflow walkthroughs",
    "6. scfactory CLI reference",
    "7. NC2024 claim coverage",
    "8. Test coverage and quality gates",
    "9. Known limitations and roadmap",
    "Appendix: file index",
)


def render_pdf(
    out_path: Path,
    module_rows: list[ModuleRow],
    recipes: dict[str, dict],
    helps: dict[str, str],
    test_counts: dict[str, int],
    claim_rows: list[list[str]],
    claim_counts: dict[str, int],
    doc_index: list[tuple[str, str]],
    generated_at: str,
) -> None:
    styles = _make_styles()
    story: list = []

    # Cover
    story.append(Paragraph(DOC_TITLE, styles["title"]))
    story.append(Paragraph(DOC_SUBTITLE, styles["subtitle"]))
    story.append(
        Paragraph(
            f"<font face='Courier'>generated_at: {generated_at}</font>",
            styles["small"],
        )
    )
    story.append(Spacer(1, 0.5 * cm))
    story.append(
        Paragraph(
            "This document is auto-generated from on-disk project state (module catalog, "
            "recipes, README/docs index, and pytest --collect-only). Re-run "
            "<font face='Courier'>python scripts/generate_factories_report.py</font> to refresh.",
            styles["body"],
        )
    )
    story.append(PageBreak())

    # Static TOC. We deliberately avoid a two-pass build (which is fragile when
    # flowables are reused across builds in ReportLab) and just list the
    # sections in their fixed order.
    story.append(Paragraph("Table of contents", styles["h1"]))
    toc_rows = [["section"]] + [[name] for name in SECTION_NAMES]
    toc_table = Table(
        _wrap_cells(toc_rows, styles, "body"),
        colWidths=[PAGE_WIDTH - 2 * MARGIN],
    )
    toc_table.setStyle(_zebra_style(len(toc_rows)))
    story.append(toc_table)
    story.append(PageBreak())

    # 1. Executive summary
    story.append(Paragraph("1. Executive summary", styles["h1"]))
    story.append(
        Paragraph(
            "<b>singlecell_factory</b> is a Python/AnnData heavy-compute pipeline organized as "
            "33 concrete module files across the layers ingest, quality_control, latent_structure, "
            "markers, annotation, state_dynamics, biology, aggregation, reproduction, spatial, "
            "multimodal, genomic_optional, external_validation, covariates, and reporting. "
            "The catalog declares 29 active module specs (3 mandatory + 26 optional) with explicit "
            "dependencies, modality tags, and bridge_ready flags. The CLI entrypoint is "
            "<font face='Courier'>scfactory run</font>.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "<b>r_multiomics_factory</b> is the R-side downstream visualization workspace. It "
            "consumes the bundle v2.1 contract emitted by "
            "<font face='Courier'>scripts/export_singlecell_r_bundle.py</font>, reads it via "
            "<font face='Courier'>R_bundle/io_bundle.R</font>, and renders Seurat-based plots via "
            "the <font face='Courier'>R/*_module.R</font> files. Together: 97 collected tests on "
            "the named contract surface, 4 starter recipes, 3 modality extensions "
            "(protein/spatial/multimodal_obsm), and an NC2024 paper-faithful reproduction lane.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "<b>Who should use it:</b> pipeline users via "
            "<font face='Courier'>scfactory run</font>; agents via the documented bundle v2.1 "
            "contract; paper-faithful reproductions via "
            "<font face='Courier'>--recipe nc2024_paper</font>.",
            styles["body"],
        )
    )
    story.append(PageBreak())

    # 2. Architecture
    story.append(Paragraph("2. Architecture", styles["h1"]))
    story.extend(_architecture_flowables(module_rows, styles))
    story.append(PageBreak())

    # 3. Capability matrix (may span two pages — repeatRows=1 handles header)
    story.append(Paragraph("3. Capability matrix", styles["h1"]))
    story.append(
        Paragraph(
            f"Auto-extracted from "
            f"<font face='Courier'>workflow/modular/module_catalog.py</font> "
            f"({len(module_rows)} modules).",
            styles["body"],
        )
    )
    story.append(_capability_table(module_rows, styles))
    story.append(PageBreak())

    # 4. Bundle contract
    story.append(Paragraph("4. Bundle v2.1 contract", styles["h1"]))
    bundle_rows = [
        ["field / topic", "value"],
        ["schema_version (v2.1)", "singlecell_r_bundle_v2.1"],
        ["schema_version (back-compat)", "singlecell_r_bundle_v2"],
        ["accepted by ExportConfig", "v1, v2, v2.1"],
        ["top-level manifest", "schema_version, source, bundle, cell_alignment, expression, files, extensions"],
        ["KNOWN_EXTENSION_KEYS", "protein, spatial, multimodal_obsm"],
        ["protein extension", "CLR-normalized ADT (parquet); fields: normalization, isotype_controls"],
        ["spatial extension", "(x,y) coords + library metadata (parquet); image data is STRING-pointer only"],
        ["multimodal_obsm extension", "EXPERIMENTAL; one parquet per joint embedding (X_wnn, X_mofa)"],
        ["claim_guard literal", "not_for_de_or_new_quantitative_claims_without_full_object_validation"],
    ]
    t = Table(_wrap_cells(bundle_rows, styles), colWidths=[5.5 * cm, 11 * cm], repeatRows=1)
    t.setStyle(_zebra_style(len(bundle_rows)))
    story.append(t)
    story.append(Spacer(1, 0.3 * cm))
    story.append(
        Paragraph(
            "<b>CLAIM_GUARD explanation.</b> The literal string is surfaced in R as both a "
            "<font face='Courier'>message()</font> call and as "
            "<font face='Courier'>attr(expr_sparse, \"claim_guard\")</font>. It means bundle "
            "data is plotting-only and must not be used for new quantitative claims (DE, "
            "thresholding, etc.) without going back to the full AnnData and re-validating.",
            styles["body"],
        )
    )
    story.append(
        Paragraph(
            "<b>Atomicity.</b> Writers stage every file under a "
            "<font face='Courier'>TemporaryDirectory</font>, compute SHA256 per file (recorded "
            "in the manifest, cross-checked against on-disk byte size), and write the manifest "
            "last. <font face='Courier'>os.replace(temp_dir, final_output_dir)</font> performs "
            "the atomic swap; any pre-existing target is moved aside as a backup and restored "
            "on failure.",
            styles["body"],
        )
    )
    story.append(PageBreak())

    # 5. Recipes
    story.append(Paragraph("5. Workflow walkthroughs", styles["h1"]))
    for stem in ("nc2024_paper", "quick_explore", "cite_seq_full", "visium_neighborhoods"):
        story.extend(_recipe_flowables(stem, recipes.get(stem) or {}, styles))
    story.append(PageBreak())

    # 6. CLI reference
    story.append(Paragraph("6. scfactory CLI reference", styles["h1"]))
    for sub in ("_root", "run", "doctor", "report"):
        title = "scfactory" if sub == "_root" else f"scfactory {sub}"
        story.append(Paragraph(title, styles["h2"]))
        story.extend(_help_block(helps.get(sub, "(no output)"), styles))
        story.append(Spacer(1, 0.2 * cm))
    story.append(PageBreak())

    # 7. NC2024 claims
    story.append(Paragraph("7. NC2024 claim coverage", styles["h1"]))
    story.extend(_claims_flowables(claim_rows, claim_counts, styles))
    story.append(PageBreak())

    # 8. Tests
    story.append(Paragraph("8. Test coverage and quality gates", styles["h1"]))
    story.extend(_tests_flowables(test_counts, styles))
    story.append(PageBreak())

    # 9. Limitations
    story.append(Paragraph("9. Known limitations and roadmap", styles["h1"]))
    bullets = [
        "<b>multimodal_integration is EXPERIMENTAL.</b> Off by default, gated by "
        "<font face='Courier'>SC_MULTIMODAL_ENGINE</font>. Emits the bundle v2.1 "
        "<i>multimodal_obsm</i> extension only when explicitly enabled.",
        "<b>Tissue images are STRING-pointer only</b> in the spatial extension. Image bytes "
        "are never serialized into the bundle; readers must resolve paths themselves.",
        "<b>Two bundle-v1 schema tests are pinned to v1</b> (not removed) to preserve "
        "back-compat coverage.",
        "<b>No cell-cell communication / pathway scoring / RNA-velocity / trajectory parity</b> "
        "yet against R; listed as roadmap.",
        "<font face='Courier'>--include-multimodal-obsm</font> is opt-in only and is never "
        "auto-passed even on the multimodal modality.",
        "Pre-existing <i>repo-doc-sync</i> hook on singlecell_factory blocks commits until "
        "citation/doc references are staged.",
        "Roadmap: a dedicated \"How to add a new module\" extension guide is planned but is "
        "intentionally NOT included in this overview.",
    ]
    for b in bullets:
        story.append(Paragraph("- " + b, styles["body"]))
    story.append(PageBreak())

    # Appendix
    story.append(Paragraph("Appendix: file index", styles["h1"]))
    story.append(
        Paragraph(
            f"All top-level *.md and docs/*.md across both repos ({len(doc_index)} entries).",
            styles["body"],
        )
    )
    story.extend(_appendix_flowables(doc_index, styles))

    doc = _Doc(
        str(out_path),
        pagesize=A4,
        leftMargin=MARGIN,
        rightMargin=MARGIN,
        topMargin=MARGIN,
        bottomMargin=MARGIN,
        title=DOC_TITLE,
        author="generate_factories_report.py",
    )
    doc.build(story)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> int:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)

    module_rows = load_module_rows()
    recipes = load_recipes()
    helps = load_scfactory_help()
    test_counts = collect_test_counts()
    claim_rows, claim_counts, _ = load_claim_matrix()
    doc_index = list_doc_index()

    generated_at = (
        _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()
    )

    md_text = render_markdown(
        module_rows, recipes, helps, test_counts,
        claim_rows, claim_counts, doc_index, generated_at,
    )
    OUT_MD.write_text(md_text)

    render_pdf(
        OUT_PDF, module_rows, recipes, helps, test_counts,
        claim_rows, claim_counts, doc_index, generated_at,
    )

    print(f"wrote {OUT_MD} ({OUT_MD.stat().st_size} bytes)")
    print(f"wrote {OUT_PDF} ({OUT_PDF.stat().st_size} bytes)")
    if _WARNINGS:
        print(f"warnings: {len(_WARNINGS)}")
        for w in _WARNINGS:
            print(f"  - {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
