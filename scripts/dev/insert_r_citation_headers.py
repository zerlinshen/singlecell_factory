"""Insert `# Citations:` header block into each R module in r_multiomics_factory.

Idempotent: a file that already contains a `# Citations:` line in its top 20
lines is skipped.

Each file gets its canonical citations harvested from
docs/SCIENTIFIC_AUDIT_2026-05-15.md (R-side section). The header is placed
after any shebang / leading comment block, before the first `library(` call.
"""
from __future__ import annotations

import sys
from pathlib import Path

R_FACTORY = Path("/home/zerlinshen/r_multiomics_factory")

R_CITATIONS: dict[str, list[str]] = {
    "R/annotation_module.R": [
        "Tirosh et al. 2016. doi:10.1126/science.aad0501 (signature scoring methodology underlying AddModuleScore)",
        "Hao et al. 2021. doi:10.1016/j.cell.2021.04.048 (Seurat 4 / AddModuleScore implementation)",
    ],
    "R/atac_module.R": [
        "Cusanovich et al. 2018. doi:10.1016/j.cell.2018.06.052 (TF-IDF + LSI methodology for sparse peak matrices)",
        "Stuart et al. 2021. doi:10.1038/s41592-021-01282-5 (Signac multimodal chromatin analysis — methodological reference)",
    ],
    "R/batch_integration_module.R": [
        "Stuart et al. 2019. doi:10.1016/j.cell.2019.05.031 (Seurat CCA / anchor-based integration)",
        "Korsunsky et al. 2019. doi:10.1038/s41592-019-0619-0 (Harmony batch correction)",
        "Hao et al. 2021. doi:10.1016/j.cell.2021.04.048 (Seurat 4 / RPCA integration)",
    ],
    "R/composition_plots.R": [
        "Wickham 2016. ggplot2: Elegant Graphics for Data Analysis. Springer. ISBN 978-3-319-24277-4",
    ],
    "R/dim_plots.R": [
        "McInnes et al. 2018. arXiv:1802.03426 (UMAP)",
        "Wickham 2016. ggplot2 (visualization framework)",
        "Hao et al. 2021. doi:10.1016/j.cell.2021.04.048 (Seurat DimPlot conventions)",
    ],
    "R/expression_plots.R": [
        "Hao et al. 2021. doi:10.1016/j.cell.2021.04.048 (Seurat FeaturePlot / DotPlot / DoHeatmap conventions)",
        "Wickham 2016. ggplot2",
    ],
    "R/integration_module.R": [
        "Hao et al. 2021. doi:10.1016/j.cell.2021.04.048 (WNN multimodal integration)",
        "Argelaguet et al. 2020. doi:10.1186/s13059-020-02015-1 (MOFA multi-omics factor analysis)",
    ],
    "R/io_bridge.R": [
        "Zellkonverter / SeuratDisk h5ad-Seurat conversion conventions (Bioconductor + satijalab packages)",
        "Hao et al. 2021. doi:10.1016/j.cell.2021.04.048 (Seurat object model)",
    ],
    "R/marker_db_module.R": [
        "Wave 1 / P1A.S5 project module — see singlecell_factory/.omc/plans/multiomics-platform-evolution-consensus-plan.md",
        "Apache Arrow Project. doi:10.1145/3514221.3526057 (Arrow Parquet format)",
    ],
    "R/marker_module.R": [
        "Hao et al. 2021. doi:10.1016/j.cell.2021.04.048 (Seurat FindMarkers / FindAllMarkers)",
        "Finak et al. 2015. doi:10.1186/s13059-015-0844-5 (MAST DE for single-cell)",
        "Soneson & Robinson 2018. doi:10.1038/nmeth.4612 (Wilcoxon benchmark for scRNA DE)",
    ],
    "R/pipeline_steps.R": [
        "Project orchestrator — combines Seurat (Hao 2021, doi:10.1016/j.cell.2021.04.048) plotting modules into a single run_plot_suite call.",
    ],
    "R/preprocessing_module.R": [
        "Hao et al. 2021. doi:10.1016/j.cell.2021.04.048 (Seurat 4/5 preprocessing: NormalizeData / FindVariableFeatures / ScaleData / RunPCA / RunUMAP / FindClusters)",
    ],
    "R/protein_module.R": [
        "Stoeckius et al. 2017. doi:10.1038/nmeth.4380 (CITE-seq ADT)",
        "Mulè et al. 2022. doi:10.1038/s41467-022-29356-8 (DSB normalization for ADT)",
    ],
    "R/qc_plots.R": [
        "Hao et al. 2021. doi:10.1016/j.cell.2021.04.048 (Seurat QC violin / scatter conventions)",
        "Luecken & Theis 2019. doi:10.15252/msb.20188746 (QC threshold guidance)",
    ],
    "R/spatial_module.R": [
        "10x Genomics Visium platform — https://www.10xgenomics.com/products/spatial-gene-expression",
        "Chen et al. 2015. doi:10.1126/science.aaa6090 (MERFISH spatial transcriptomics)",
        "10x Genomics Xenium platform — https://www.10xgenomics.com/products/xenium-in-situ",
    ],
    "R/theme_config.R": [
        "Wickham 2016. ggplot2: Elegant Graphics for Data Analysis. Springer. ISBN 978-3-319-24277-4",
    ],
    "R/cli_utils.R": [
        "Project utility — see r_multiomics_factory/AGENTS.md for --project-root / --run-id contract.",
    ],
    "R_bundle/io_bundle.R": [
        "Project bundle schema — see singlecell_factory/contracts/bundle_schema.yaml (vendored copy at r_multiomics_factory/contracts/bundle_schema.yaml).",
        "Apache Arrow Project. doi:10.1145/3514221.3526057 (Arrow Parquet format used for v2/v2.1/v2.2 bundles).",
    ],
    "R_bundle/bundle_cache.R": [
        "Project utility — caches bundle SHA256 computation results to reduce redundant I/O.",
    ],
    "R_bundle/remote_bundle_manifest.R": [
        "Project utility — validates v1/v2/v2.1/v2.2 bundle manifests per singlecell_factory/contracts/bundle_schema.yaml.",
    ],
}


def _format_block(citations: list[str]) -> str:
    lines = ["# Citations:"]
    for c in citations:
        lines.append(f"#   - {c}")
    lines.append("")
    return "\n".join(lines) + "\n"


def patch_file(rel_path: str, citations: list[str]) -> str:
    path = R_FACTORY / rel_path
    if not path.exists():
        return f"MISSING ({path})"
    src = path.read_text(encoding="utf-8")
    top_20_lines = src.splitlines()[:20]
    if any("# Citations:" in line for line in top_20_lines):
        return "skipped (already has Citations block)"

    # Find insertion point: after any leading comment block (# ...), before first non-comment / non-blank line.
    lines = src.splitlines(keepends=True)
    insert_at = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped == "" or stripped.startswith("#"):
            insert_at = i + 1
            continue
        break

    block = _format_block(citations)
    new_src = "".join(lines[:insert_at]) + block + "".join(lines[insert_at:])
    path.write_text(new_src, encoding="utf-8")
    return f"patched (+{len(citations)} cites)"


def main() -> int:
    if not R_FACTORY.exists():
        print(f"ERROR: r-factory not found at {R_FACTORY}", file=sys.stderr)
        return 2

    rc = 0
    for rel, cites in R_CITATIONS.items():
        result = patch_file(rel, cites)
        print(f"  {rel:<40} {result}")
        if result.startswith("MISSING"):
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
