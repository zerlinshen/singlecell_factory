"""AST-driven __references__ insertion for workflow.modular.modules.

Usage:
    python scripts/dev/insert_module_references.py [--tier 1|2|3|all] [--dry-run]

Reads the canonical citation table embedded below (sourced verbatim from
docs/SCIENTIFIC_AUDIT_2026-05-15.md). For each target module:
  1. Parse with ast to locate the last top-level Import / ImportFrom.
  2. If the module already has a top-level __references__ assignment, SKIP
     (idempotent).
  3. Otherwise insert the __references__ literal immediately after the last
     top-level import, before any class/function defs.
  4. Write back.

Run-time invariants:
  - Pure file I/O on .py modules; no behavioural test execution.
  - Idempotent: re-running on a fully-patched repo is a no-op.
  - Order of inserted entries matches the audit doc's first-citation order
    per module.
"""
from __future__ import annotations

import argparse
import ast
import io
import json
import re
import sys
import textwrap
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MODULES_DIR = REPO_ROOT / "workflow" / "modular" / "modules"

# ---------------------------------------------------------------------------
# Canonical citation table (verbatim from docs/SCIENTIFIC_AUDIT_2026-05-15.md)
# Keys: module filename within MODULES_DIR.
# Values: dict mapping ref_key -> {title, authors, journal, year, doi, description}
# ---------------------------------------------------------------------------

TIER_1_REFS: dict[str, dict[str, dict[str, str]]] = {
    "cellranger.py": {
        "Zheng_10x_2017": {
            "title": "Massively parallel digital transcriptional profiling of single cells",
            "authors": "Zheng et al.",
            "journal": "Nature Communications",
            "year": "2017",
            "doi": "10.1038/ncomms14049",
            "description": "10x Chromium Single Cell 3' chemistry — describes the output format this module ingests.",
        },
        "CellRanger_software": {
            "title": "Cell Ranger Single Cell Software (v7+)",
            "authors": "10x Genomics",
            "journal": "Software documentation",
            "year": "2024",
            "doi": "https://support.10xgenomics.com/single-cell-gene-expression/software",
            "description": "Authoritative spec for the barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz layout consumed here.",
        },
    },
    "qc.py": {
        "scanpy": {
            "title": "SCANPY: large-scale single-cell gene expression data analysis",
            "authors": "Wolf, Angerer, Theis",
            "journal": "Genome Biology",
            "year": "2018",
            "doi": "10.1186/s13059-017-1382-0",
            "description": "scanpy.pp.calculate_qc_metrics / pp.filter_genes / pp.filter_cells semantics used here.",
        },
        "Luecken_Theis_2019": {
            "title": "Current best practices in single-cell RNA-seq analysis: a tutorial",
            "authors": "Luecken, Theis",
            "journal": "Molecular Systems Biology",
            "year": "2019",
            "doi": "10.15252/msb.20188746",
            "description": "Canonical QC threshold guidance (mito %, n_genes/cell, doublet detection).",
        },
    },
    "doublet_detection.py": {
        "Wolock_Scrublet_2019": {
            "title": "Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data",
            "authors": "Wolock, Lopez, Klein",
            "journal": "Cell Systems",
            "year": "2019",
            "doi": "10.1016/j.cels.2018.11.005",
            "description": "Reference implementation; supports whole-dataset and per-sample (grouped) doublet rate calibration.",
        },
    },
    "clustering.py": {
        "scanpy": {
            "title": "SCANPY: large-scale single-cell gene expression data analysis",
            "authors": "Wolf, Angerer, Theis",
            "journal": "Genome Biology",
            "year": "2018",
            "doi": "10.1186/s13059-017-1382-0",
            "description": "scanpy pp.pca / pp.neighbors / tl.umap / tl.leiden stack.",
        },
        "Halko_TruncatedSVD_2011": {
            "title": "Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions",
            "authors": "Halko, Martinsson, Tropp",
            "journal": "SIAM Review",
            "year": "2011",
            "doi": "10.1137/090771806",
            "description": "Truncated SVD method underlying sklearn.decomposition.TruncatedSVD used in CSS path.",
        },
        "McInnes_UMAP_2018": {
            "title": "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction",
            "authors": "McInnes, Healy, Melville",
            "journal": "arXiv preprint",
            "year": "2018",
            "doi": "arXiv:1802.03426",
            "description": "KNN-graph + UMAP layout used by scanpy.pp.neighbors / sc.tl.umap.",
        },
        "Traag_Leiden_2019": {
            "title": "From Louvain to Leiden: guaranteeing well-connected communities",
            "authors": "Traag, Waltman, van Eck",
            "journal": "Scientific Reports",
            "year": "2019",
            "doi": "10.1038/s41598-019-41695-z",
            "description": "Leiden community detection used as the clustering algorithm.",
        },
        "rapids_singlecell": {
            "title": "rapids-singlecell — GPU-accelerated scanpy",
            "authors": "scverse contributors",
            "journal": "Software (scverse)",
            "year": "2024",
            "doi": "https://github.com/scverse/rapids_singlecell",
            "description": "GPU acceleration path. Falls back to CPU scanpy when unavailable.",
        },
    },
    "batch_correction.py": {
        "Korsunsky_Harmony_2019": {
            "title": "Fast, sensitive and accurate integration of single-cell data with Harmony",
            "authors": "Korsunsky et al.",
            "journal": "Nature Methods",
            "year": "2019",
            "doi": "10.1038/s41592-019-0619-0",
            "description": "Default backend (rsc.pp.harmony_integrate GPU port, US-B3 parity ARI=1.000).",
        },
        "Polanski_BBKNN_2020": {
            "title": "BBKNN: fast batch alignment of single cell transcriptomes",
            "authors": "Polanski et al.",
            "journal": "Bioinformatics",
            "year": "2020",
            "doi": "10.1093/bioinformatics/btz625",
            "description": "BBKNN backend.",
        },
        "Johnson_ComBat_2007": {
            "title": "Adjusting batch effects in microarray expression data using empirical Bayes methods",
            "authors": "Johnson, Li, Rabinovic",
            "journal": "Biostatistics",
            "year": "2007",
            "doi": "10.1093/biostatistics/kxj037",
            "description": "ComBat backend (sc.pp.combat).",
        },
        "Hie_Scanorama_2019": {
            "title": "Efficient integration of heterogeneous single-cell transcriptomes using Scanorama",
            "authors": "Hie, Bryson, Berger",
            "journal": "Nature Biotechnology",
            "year": "2019",
            "doi": "10.1038/s41587-019-0113-3",
            "description": "Scanorama backend.",
        },
        "Lopez_scVI_2018": {
            "title": "Deep generative modeling for single-cell transcriptomics",
            "authors": "Lopez et al.",
            "journal": "Nature Methods",
            "year": "2018",
            "doi": "10.1038/s41592-018-0229-2",
            "description": "scVI backend via scvi-tools.",
        },
        "Haghverdi_MNN_2018": {
            "title": "Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors",
            "authors": "Haghverdi et al.",
            "journal": "Nature Biotechnology",
            "year": "2018",
            "doi": "10.1038/nbt.4091",
            "description": "MNN / fastMNN backends.",
        },
    },
    "annotation.py": {
        "Tirosh_marker_scoring_2016": {
            "title": "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq",
            "authors": "Tirosh et al.",
            "journal": "Science",
            "year": "2016",
            "doi": "10.1126/science.aad0501",
            "description": "Cluster-vote marker-mean scoring approach for cell type assignment.",
        },
        "Stuart_label_transfer_2019": {
            "title": "Comprehensive Integration of Single-Cell Data",
            "authors": "Stuart et al.",
            "journal": "Cell",
            "year": "2019",
            "doi": "10.1016/j.cell.2019.05.031",
            "description": "Reference-based label transfer principles; this module uses a simpler sklearn.NearestNeighbors KNN majority vote on a labeled reference.",
        },
    },
    "differential_expression.py": {
        "scanpy": {
            "title": "SCANPY: large-scale single-cell gene expression data analysis",
            "authors": "Wolf, Angerer, Theis",
            "journal": "Genome Biology",
            "year": "2018",
            "doi": "10.1186/s13059-017-1382-0",
            "description": "scanpy.tl.rank_genes_groups implementation (Wilcoxon/t-test/MAST/ROC).",
        },
        "Soneson_DE_benchmark_2018": {
            "title": "Bias, robustness and scalability in single-cell differential expression analysis",
            "authors": "Soneson, Robinson",
            "journal": "Nature Methods",
            "year": "2018",
            "doi": "10.1038/nmeth.4612",
            "description": "Benchmark showing Wilcoxon competitive with bespoke single-cell DE methods.",
        },
        "Benjamini_Hochberg_1995": {
            "title": "Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing",
            "authors": "Benjamini, Hochberg",
            "journal": "Journal of the Royal Statistical Society B",
            "year": "1995",
            "doi": "10.1111/j.2517-6161.1995.tb02031.x",
            "description": "BH FDR correction applied to per-gene p-values and per-substate marker scoring.",
        },
    },
}


# ---------------------------------------------------------------------------
# Insertion engine
# ---------------------------------------------------------------------------

def _has_module_level_references(src: str) -> bool:
    """True if the source already has a top-level __references__ assignment."""
    tree = ast.parse(src)
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and tgt.id == "__references__":
                    return True
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == "__references__":
            if node.value is not None:
                return True
    return False


def _find_last_top_import_lineno(src: str) -> int:
    """1-based line number of the last top-level Import or ImportFrom statement."""
    tree = ast.parse(src)
    last = 0
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            # node.end_lineno is the inclusive last line of the import statement
            last = max(last, getattr(node, "end_lineno", node.lineno))
    return last


def _format_references_literal(refs: dict[str, dict[str, str]]) -> str:
    """Render the __references__ dict as a Python literal block."""
    lines = ["", "", "__references__ = {"]
    for key, entry in refs.items():
        lines.append(f"    {json.dumps(key)}: {{")
        for k in ("title", "authors", "journal", "year", "doi", "description"):
            if k in entry:
                lines.append(f"        {json.dumps(k)}: {json.dumps(entry[k])},")
        lines.append("    },")
    lines.append("}")
    lines.append("")
    return "\n".join(lines)


def patch_module(module_path: Path, refs: dict[str, dict[str, str]], dry_run: bool = False) -> str:
    src = module_path.read_text(encoding="utf-8")
    if _has_module_level_references(src):
        return "skipped (already has __references__)"

    last_import_lineno = _find_last_top_import_lineno(src)
    if last_import_lineno == 0:
        return "error: no top-level imports found (cannot place __references__)"

    insertion_block = _format_references_literal(refs)

    src_lines = src.splitlines(keepends=True)
    # Insert AFTER the last import line (i.e., at index = last_import_lineno)
    new_src = "".join(src_lines[:last_import_lineno]) + insertion_block + "\n" + "".join(src_lines[last_import_lineno:])

    # Validate the patched source parses
    try:
        ast.parse(new_src)
    except SyntaxError as exc:
        return f"error: post-patch source fails to parse: {exc}"

    # Re-check that __references__ is now top-level
    if not _has_module_level_references(new_src):
        return "error: post-patch source does not expose __references__ at module level"

    if dry_run:
        return f"would patch (+{len(refs)} refs)"

    module_path.write_text(new_src, encoding="utf-8")
    return f"patched (+{len(refs)} refs)"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _apply_table(table: dict[str, dict[str, dict[str, str]]], dry_run: bool) -> int:
    rc = 0
    for filename, refs in table.items():
        path = MODULES_DIR / filename
        if not path.exists():
            print(f"  {filename}: MISSING ({path})", file=sys.stderr)
            rc = 1
            continue
        result = patch_module(path, refs, dry_run=dry_run)
        prefix = "DRY" if dry_run else "OK"
        print(f"  [{prefix}] {filename}: {result}")
        if result.startswith("error:"):
            rc = 1
    return rc


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tier", choices=["1", "2", "3", "all"], default="1")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    rc = 0
    if args.tier in ("1", "all"):
        print("=== Tier 1 (foundational) ===")
        rc |= _apply_table(TIER_1_REFS, args.dry_run)
    # Tier 2/3 tables registered in subsequent stories (US-002, US-003).
    if args.tier in ("2", "all"):
        try:
            from scripts.dev.tier2_refs import TIER_2_REFS  # type: ignore
            print("\n=== Tier 2 (downstream biology) ===")
            rc |= _apply_table(TIER_2_REFS, args.dry_run)
        except ImportError:
            print("\n=== Tier 2 (downstream biology) === SKIP (tier2_refs not yet defined; see US-002)")
    if args.tier in ("3", "all"):
        try:
            from scripts.dev.tier3_refs import TIER_3_REFS  # type: ignore
            print("\n=== Tier 3 (project-local utilities) ===")
            rc |= _apply_table(TIER_3_REFS, args.dry_run)
        except ImportError:
            print("\n=== Tier 3 (project-local utilities) === SKIP (tier3_refs not yet defined; see US-003)")

    return rc


if __name__ == "__main__":
    sys.exit(main())
