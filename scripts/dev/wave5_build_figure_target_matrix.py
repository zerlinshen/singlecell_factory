#!/usr/bin/env python3
"""Build the Wave-5 raw/public-data Trevino figure target matrix.

This is a planning/evidence utility, not a renderer.  It reads the local Cell
paper PDF text extraction plus the current public-resource manifests and writes a
machine-readable matrix that separates:

* panels that can be regenerated from public sequencing-derived resources,
* panels that are non-computational schematics/IHC images, and
* panels blocked by missing public models, variant catalogs, or external
  datasets.

Outputs are written under the project run directory so scientific artifacts do
not land in the factory repo.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class PanelTarget:
    figure_id: str
    panel_id: str
    panel_label: str
    panel_type: str
    required_public_inputs: str
    local_or_public_resources: str
    planned_reproduction_method: str
    initial_status: str
    status_reason: str
    quantitative_check: str


def _load_json(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {}


def _resource_status(
    resource_manifest: dict,
    author_inventory: dict,
    author_download_manifest: dict | None = None,
) -> dict[str, str]:
    status: dict[str, str] = {}
    for row in resource_manifest.get("geo_supplementary_files", []):
        name = Path(row.get("project_input_path", row.get("filename", ""))).name
        status[name] = row.get("status", "unknown")
    for row in author_inventory.get("rows", []):
        status[f"author:{row.get('filename')}"] = row.get("status", "unknown")
    for row in (author_download_manifest or {}).get("geo_supplementary_files", []):
        # Download status is stronger evidence than HEAD discovery.  Keep the
        # key namespace stable so downstream panel logic can simply query
        # ``author:<filename>`` regardless of whether evidence came from the
        # inventory or the download manifest.
        filename = row.get("filename")
        if filename:
            status[f"author:{filename}"] = row.get("status", "unknown")
    return status


def _main_panel_targets(status: dict[str, str], brain_asd_present: bool = False) -> list[PanelTarget]:
    """Curated main-figure matrix from the local paper caption text.

    The matrix intentionally starts conservatively.  A later renderer/validator
    must promote a panel to REPRODUCED only after it has generated outputs and
    recorded a metric.
    """

    def geo(name: str) -> str:
        return status.get(name, "missing")

    have_multiome = all(
        geo(name).startswith("present")
        for name in [
            "GSE162170_multiome_rna_counts.tsv.gz",
            "GSE162170_multiome_atac_counts.tsv.gz",
            "GSE162170_multiome_cell_metadata.txt.gz",
        ]
    )
    full_rna_pending = not geo("GSE162170_rna_counts.tsv.gz").startswith("present")
    full_atac_pending = not geo("GSE162170_atac_counts.tsv.gz").startswith("present")
    gene_activity_pending = not geo("GSE162170_atac_gene_activities.tsv.gz").startswith("present")
    chromvar_pending = not geo("GSE162170_atac_chromVAR_TF_activities.tsv.gz").startswith("present")
    author_fcm_status = status.get("author:FCM_Object.RDS", "missing")
    author_glial_rna_status = status.get("author:RNA_GlialPseudobulks.RDS", "missing")
    author_fcm_known = author_fcm_status in {"head_ok", "present"}
    author_fcm_ready = author_fcm_status.startswith("present") and author_glial_rna_status.startswith("present")
    author_bpnet_inputs_known = any(
        k.startswith("author:hft_c") and v in {"head_ok", "present"} for k, v in status.items()
    )

    pending = "DOWNLOAD_PENDING"
    multiome_status = "READY_FOR_LOADER" if have_multiome else pending
    full_singleome_status = pending if (full_rna_pending or full_atac_pending) else "READY_FOR_LOADER"
    ga_status = pending if gene_activity_pending else "READY_FOR_LOADER"
    chromvar_status = pending if chromvar_pending else "READY_FOR_LOADER"
    fcm_status = (
        "READY_FOR_LOADER"
        if author_fcm_ready
        else "AUTHOR_RESOURCE_DISCOVERED_NOT_DOWNLOADED"
        if author_fcm_known
        else "RESOURCE_SURVEY_PENDING"
    )
    fcm_reason = (
        "Author FCM and glial pseudobulk resources are downloaded; loader/render validation required."
        if author_fcm_ready
        else "Author FCM resources discovered but not yet downloaded/validated."
        if author_fcm_known
        else "Author FCM resources not yet surveyed."
    )
    velocity_ready = (
        not full_rna_pending
        and geo("GSE162170_rna_spliced_counts.tsv.gz").startswith("present")
        and geo("GSE162170_rna_unspliced_counts.tsv.gz").startswith("present")
    )
    velocity_status = "READY_FOR_LOADER" if velocity_ready else pending
    velocity_reason = (
        "RNA counts plus spliced/unspliced matrices are present; velocity loader/dependency validation required."
        if velocity_ready
        else "Spliced/unspliced matrices are pending and velocity dependency status still needs smoke."
    )
    bpnet_status = "METHOD_RESOURCE_GAP"
    if author_bpnet_inputs_known:
        bpnet_status = "PARTIAL_AUTHOR_INPUTS_DISCOVERED_MODEL_GAP"
    bpnet_resources = "author cluster bigWigs discovered; BPNet models/variant catalogs not found"
    bpnet_reason = "Author bigWigs are discoverable, but trained BPNet model weights and mutation catalogs are not yet in manifest."
    if brain_asd_present:
        bpnet_resources = "Brain_ASD code/peaks/GC negatives cloned; author cluster bigWigs discovered; trained model weights and scored mutation catalogs not found"
        bpnet_reason = "Public Brain_ASD training/scoring code is present, but finished BPNet model weights and Figure 7 mutation score tables are not in the cloned/public manifests."

    rows = [
        # Figure 1
        PanelTarget("1", "1A", "time/profiling/cell-type schematic", "schematic", "paper metadata only", "local PDF", "Do not regenerate as sequencing result; document as NON_COMPUTATIONAL.", "NON_COMPUTATIONAL", "Schematic, not a data-derived pipeline figure.", "presence in PDF text"),
        PanelTarget("1", "1B", "SOX9/CTIP2 IHC", "imaging", "source microscopy images", "local PDF raster only", "Cannot regenerate from scRNA/scATAC matrices; cite as non-computational source panel.", "NON_COMPUTATIONAL", "Raw microscopy images are not in discovered GEO/S3 matrices.", "classification evidence"),
        PanelTarget("1", "1C", "GFAP/KI67/PPP1R17 IHC", "imaging", "source microscopy images", "local PDF raster only", "Cannot regenerate from scRNA/scATAC matrices; cite as non-computational source panel.", "NON_COMPUTATIONAL", "Raw microscopy images are not in discovered GEO/S3 matrices.", "classification evidence"),
        PanelTarget("1", "1D", "RNA and ATAC UMAP by time", "embedding", "scRNA, scATAC and metadata matrices", "GEO GSE162170 RNA/ATAC resources", "Recompute or load embeddings, color by gestational age/time.", full_singleome_status, "Full scRNA/scATAC matrices are still being downloaded." if full_singleome_status == pending else "Inputs present; loader smoke required.", "cluster/time distribution and embedding existence"),
        PanelTarget("1", "1E", "SOX9/EOMES/NEUROD2/DLX2 multimodal overlays", "multiome marker overlay", "multiome RNA, gene activity, motif activity", "PCW21 multiome counts present; gene activity/chromVAR pending from GEO/S3", "Render RNA expression, ATAC gene activity, and TF motif activity overlays.", "PARTIAL_INPUTS_PENDING", "RNA/ATAC counts are present, but activity/motif inputs are incomplete until downloads finish.", "marker recovery and per-panel nonzero fraction"),
        PanelTarget("1", "1F", "cluster UMAP", "embedding", "multiome or full singleome cell metadata/counts", "PCW21 multiome counts and cluster names present", "Generate UMAP/clusters from public matrices and compare cluster names/counts.", multiome_status, "PCW21 multiome inputs are locally present." if have_multiome else "Local multiome inputs incomplete.", "ARI/NMI or cluster-count concordance"),
        PanelTarget("1", "1G", "marker-expression dotplot", "dotplot", "scRNA counts plus cluster metadata", "GEO RNA counts/metadata", "Compute dotplot for canonical cortical markers by cluster.", pending if full_rna_pending else "READY_FOR_LOADER", "Full scRNA counts are pending download." if full_rna_pending else "Inputs present; loader smoke required.", "marker-panel recall"),
        PanelTarget("1", "1H", "ATAC gene-activity marker dotplot", "dotplot", "scATAC gene activity plus cluster metadata", "GEO ATAC gene-activity matrix", "Compute marker gene activity dotplot by ATAC cluster.", ga_status, "ATAC gene-activity matrix pending download." if ga_status == pending else "Inputs present; loader smoke required.", "marker activity recall"),
        # Figure 2
        PanelTarget("2", "2A", "singleome integration schematic", "schematic", "paper metadata only", "local PDF", "Document design panel only.", "NON_COMPUTATIONAL", "Schematic, not a data-derived output.", "presence in PDF text"),
        PanelTarget("2", "2B", "matched scRNA/scATAC cluster UMAPs", "integration embedding", "full scRNA/scATAC matrices, metadata, CCA matching", "GEO matrices; author CCA_Matching.RDS discovered", "Recompute integration or load author matching object if needed, then compare cluster transfer.", full_singleome_status, "Full singleome inputs still pending." if full_singleome_status == pending else "Inputs present; integration smoke required.", "label-transfer agreement"),
        PanelTarget("2", "2C", "CRE-gene linkage heatmap", "peak-gene heatmap", "peak-gene links, ATAC accessibility, RNA expression", "local S2F link table; author PeakGeneLinks resources discovered; full matrices pending", "Rebuild pseudobulk heatmap of linked CRE-gene pairs.", "PARTIAL_INPUTS_PENDING", "Peak-gene evidence exists, but full matrices/downloads and loader validation are pending.", "top-link overlap and heatmap dimensions"),
        PanelTarget("2", "2D", "gene expression versus gene activity correlation", "correlation scatter", "scRNA expression, ATAC gene activity, link counts", "GEO RNA/ATAC gene activity; author GA_RNA_Correlations.RDS discovered", "Compute gene-level RNA/GA correlation and annotate TFs/GPCs.", "PARTIAL_INPUTS_PENDING", "Requires full RNA plus ATAC gene activity downloads.", "Spearman/Pearson concordance with author table"),
        PanelTarget("2", "2E", "GO enrichment of GPC genes", "enrichment", "GPC gene set and GO annotation", "author gene sets discovered; no GO DB recorded yet", "Run GO enrichment from public/installed annotation if available, otherwise document DB gap.", "METHOD_RESOURCE_GAP", "GO annotation database/tool has not yet been pinned in the run manifest.", "top-term overlap"),
        PanelTarget("2", "2F", "multiome generation schematic", "schematic", "paper metadata only", "local PDF", "Document design panel only.", "NON_COMPUTATIONAL", "Schematic, not a data-derived output.", "presence in PDF text"),
        PanelTarget("2", "2G", "multiome projections into singleome spaces", "projection", "multiome RNA/ATAC plus full singleome embeddings", "PCW21 multiome present; full singleome pending", "Project multiome RNA/ATAC into public singleome manifolds.", "PARTIAL_INPUTS_PENDING", "Multiome present but full singleome reference download/loader not finished.", "projection cluster agreement"),
        PanelTarget("2", "2H", "singleome versus multiome linkage overlap", "Venn/overlap", "singleome and multiome peak-gene links", "local S2F table; author PeakGeneLinks resources discovered", "Compute overlap once link tables are downloaded/parsed.", "PARTIAL_INPUTS_PENDING", "Requires complete/public link source mapping.", "Jaccard/overlap count"),
        PanelTarget("2", "2I", "predictive chromatin correspondence", "correlation scatter", "singleome and multiome linkage/correlation scores", "author GA_RNA_Correlations/PeakGeneLinks discovered", "Compare singleome vs multiome predictive chromatin statistics.", "PARTIAL_INPUTS_PENDING", "Requires full link/stat tables.", "correlation coefficient"),
        # Figure 3
        PanelTarget("3", "3A", "GluN RNA-velocity pseudotime", "trajectory/velocity", "scRNA counts, spliced/unspliced, velocity method", "GEO RNA spliced/unspliced", "Compute RNA velocity/pseudotime for GluN trajectory.", velocity_status, velocity_reason, "velocity graph and pseudotime marker monotonicity"),
        PanelTarget("3", "3B", "ATAC transferred pseudotime", "trajectory transfer", "scATAC embedding plus RNA pseudotime transfer", "GEO ATAC pending", "Transfer RNA pseudotime to ATAC cells by integration/nearest neighbors.", pending, "Depends on Figure 3A and full ATAC matrices.", "pseudotime transfer concordance"),
        PanelTarget("3", "3C", "GluN CRE-gene heatmap along pseudotime", "trajectory heatmap", "GluN pseudotime, links, RNA, ATAC", "GEO matrices and link tables", "Aggregate pseudobulks along GluN pseudotime and plot linked CRE-gene heatmap.", pending, "Depends on Figure 3A/B plus link tables.", "heatmap cluster sizes and known marker ordering"),
        PanelTarget("3", "3D", "GSEA for interaction clusters", "enrichment", "interaction cluster gene sets and pathway DB", "author gene_sets discovered", "Run gene-set enrichment for CRE-gene clusters.", "METHOD_RESOURCE_GAP", "Pathway database/version not yet pinned.", "top-term overlap"),
        PanelTarget("3", "3E", "TF motif enrichment in interaction clusters", "motif enrichment", "peak sets and motif annotations", "GEO chromVAR pending; author motif matrix discovered", "Run motif enrichment or use validated author motif matrix.", chromvar_status, "ChromVAR/motif inputs pending." if chromvar_status == pending else "Motif inputs present; loader smoke required.", "motif odds-ratio overlap"),
        PanelTarget("3", "3F", "dynamic TF expression and motif activity heatmaps", "dynamic heatmap", "TF expression, chromVAR activity, pseudotime bins", "GEO RNA/chromVAR; author TF_MotifExpressionCorrelation discovered", "Aggregate TF expression/activity along pseudotime.", "PARTIAL_INPUTS_PENDING", "Requires RNA/chromVAR downloads and pseudotime.", "dynamic TF panel recall"),
        PanelTarget("3", "3G", "motif correlation and synergy", "TF network", "motif activity clusters and synergy calculation", "author TF_MotifExpressionCorrelation discovered", "Recompute motif correlations; synergy may require original algorithm resources.", "METHOD_RESOURCE_GAP", "Synergy calculation method/resources require author-code audit.", "motif-cluster concordance"),
        PanelTarget("3", "3H", "motif activity versus gene expression correlations", "correlation heatmap", "TF expression and motif activity", "GEO RNA/chromVAR; author TF correlation object discovered", "Compute TF motif/gene-expression correlations.", "PARTIAL_INPUTS_PENDING", "Requires RNA/chromVAR matrices.", "correlation rank overlap"),
        PanelTarget("3", "3I", "gene-expression pseudotime versus motif synergy", "scatter", "pseudotime and motif synergy scores", "author code/resources incomplete", "Attempt after synergy audit; otherwise mark method-resource gap.", "METHOD_RESOURCE_GAP", "Depends on unavailable/unaudited synergy implementation.", "scatter correlation"),
        # Figure 4
        PanelTarget("4", "4A", "glial clustering/reprojection schematic and pseudobulks", "trajectory/reprojection", "glial RNA pseudobulks and FCM object", "author FCM/RNA_Glial resources discovered" if author_fcm_known else "resource pending", "Reconstruct or load public glial fuzzy-clustering state and pseudobulks.", fcm_status, fcm_reason, "pseudobulk count and cluster agreement"),
        PanelTarget("4", "4B", "module expression heatmap", "module heatmap", "glial pseudobulk module expression", "author FCM/RNA_Glial resources discovered", "Render module heatmap across pseudobulks.", fcm_status, fcm_reason, "module order concordance"),
        PanelTarget("4", "4C", "selected gene heatmap across glial pseudobulks", "gene heatmap", "glial pseudobulk expression", "author RNA_GlialPseudobulks discovered", "Render selected gene expression heatmap.", fcm_status, fcm_reason, "marker gene recovery"),
        PanelTarget("4", "4D", "module expression in UMAP embedding", "module embedding", "FCM embedding and module scores", "author FCM resources discovered", "Render module expression over fuzzy-clustering embedding.", fcm_status, fcm_reason, "module centroid concordance"),
        PanelTarget("4", "4E", "module centroids and overlap graph", "network overlay", "module centroids and Jaccard overlap", "author FCM resources discovered", "Compute/plot centroid graph thresholded by Jaccard.", fcm_status, fcm_reason, "edge-count/Jaccard summary"),
        PanelTarget("4", "4F", "ASCL1/HES4/OLIG1 module membership", "module membership", "gene modules and expression", "author FCM resources discovered", "Render membership/expression values.", fcm_status, fcm_reason, "target gene membership check"),
        PanelTarget("4", "4G", "EOMES/AQP4/MBP module membership", "module membership", "gene modules and expression", "author FCM resources discovered", "Render membership/expression values.", fcm_status, fcm_reason, "target gene membership check"),
        PanelTarget("4", "4H", "ASCL1/OLIG1 and EGFR membership", "module membership", "gene modules and expression", "author FCM resources discovered", "Render membership/expression values.", fcm_status, fcm_reason, "target gene membership check"),
        PanelTarget("4", "4I", "ASCL1/OLIG2/EGFR IHC", "imaging", "source microscopy images", "local PDF raster only", "Cannot regenerate from sequencing matrices.", "NON_COMPUTATIONAL", "Raw microscopy images not discovered as public data.", "classification evidence"),
        PanelTarget("4", "4J", "PDGFRA and SPARCL1 membership", "module membership", "gene modules and expression", "author FCM resources discovered", "Render membership/expression values.", fcm_status, fcm_reason, "target gene membership check"),
        PanelTarget("4", "4K", "SPARCL1/PDGFRA IHC", "imaging", "source microscopy images", "local PDF raster only", "Cannot regenerate from sequencing matrices.", "NON_COMPUTATIONAL", "Raw microscopy images not discovered as public data.", "classification evidence"),
        # Figure 5
        PanelTarget("5", "5A", "astrocyte-associated module membership", "module membership", "glial modules and expression", "author FCM resources discovered", "Render AQP4/TNC/ALDH2/APOE module membership.", fcm_status, fcm_reason, "target gene membership check"),
        PanelTarget("5", "5B", "motif enrichments in module-linked GREs", "motif enrichment", "module-linked GREs and motif annotations", "author motif/link resources discovered", "Compute motif enrichment for module 13 versus 14 GREs.", "PARTIAL_INPUTS_PENDING", "Requires module links plus motif matrix/annotations.", "motif enrichment overlap"),
        PanelTarget("5", "5C", "AQP4-positive glial reclustering", "clustering", "glial pseudobulks and FCM embedding", "author FCM resources discovered", "Recluster glial pseudobulks and identify A1-HES/A2-OLIG.", fcm_status, fcm_reason, "cluster count and marker specificity"),
        PanelTarget("5", "5D", "A1-HES versus A2-OLIG differential expression", "differential expression", "glial pseudobulk counts and DESeq2-equivalent method", "author RNA_GlialPseudobulks discovered", "Run DESeq2 if R package present or a documented equivalent fallback.", "METHOD_RESOURCE_GAP", "DESeq2/package availability and exact design need validation.", "top-DE gene overlap"),
        PanelTarget("5", "5E", "Bhaduri fetal scRNA UMAP by cortical area", "external reference projection", "Bhaduri et al. fetal scRNA reference", "not yet downloaded", "Download/locate Bhaduri reference or mark public-data gap.", "PUBLIC_DATA_GAP", "External reference dataset is not in current manifest.", "reference dataset presence and projection metrics"),
        PanelTarget("5", "5F", "module/DE gene expression in Bhaduri reference", "external reference validation", "Bhaduri reference and module gene sets", "not yet downloaded", "Score module/DE gene sets in Bhaduri reference.", "PUBLIC_DATA_GAP", "External reference dataset is not in current manifest.", "gene-set score concordance"),
        # Figure 6
        PanelTarget("6", "6A", "cell-cycle signature correlation with modules", "module correlation", "cell-cycle gene set and glial modules", "author gene set/FCM resources discovered", "Compute signature-module correlations across pseudobulks.", fcm_status, fcm_reason, "correlation matrix sanity"),
        PanelTarget("6", "6B", "ATAC projection schematic", "schematic", "paper metadata only", "local PDF", "Document design panel only.", "NON_COMPUTATIONAL", "Schematic, not a data-derived output.", "presence in PDF text"),
        PanelTarget("6", "6C", "ATAC pseudobulks projected into cycling modules", "projection", "ATAC pseudobulk GA/accessibility and FCM embedding", "author ATAC glial matrices discovered", "Project ATAC pseudobulks into FCM/cycling module space.", "PARTIAL_INPUTS_PENDING", "Requires author glial ATAC matrices or local recompute.", "projection neighborhood concordance"),
        PanelTarget("6", "6D", "branch-specific gene-activity heatmap", "gene-activity heatmap", "ATAC gene activity and branch labels", "GEO ATAC gene activities pending", "Compute branch-specific top active genes and GPC enrichment.", ga_status, "ATAC gene activity pending." if ga_status == pending else "Inputs present; branch labels still required.", "GPC enrichment and marker recall"),
        PanelTarget("6", "6E", "GPC motif and expression dynamics across branches", "dynamic heatmap", "GPC motif activity and RNA expression", "GEO chromVAR/RNA resources pending", "Aggregate motif activity/expression across branches.", "PARTIAL_INPUTS_PENDING", "Requires chromVAR, RNA, and branch labels.", "known GPC motif recovery"),
        PanelTarget("6", "6F", "GPC-only chromatin reprojection", "projection", "GPC-linked ATAC accessibility", "peak-gene links and ATAC counts pending", "Reproject branches using GPC-associated chromatin.", "PARTIAL_INPUTS_PENDING", "Requires link table and full ATAC matrix.", "branch separability metric"),
        PanelTarget("6", "6G", "multiome RNA projection into FCM embedding", "projection", "multiome RNA and FCM embedding", "multiome present; FCM discovered", "Project multiome RNA cells into glial FCM embedding.", "PARTIAL_INPUTS_PENDING", "Requires FCM resources plus multiome loader.", "projection cluster concordance"),
        # Figure 7
        PanelTarget("7", "7A", "mutation-prioritization schematic", "schematic", "paper metadata only", "local PDF", "Document design panel only.", "NON_COMPUTATIONAL", "Schematic, not a data-derived output.", "presence in PDF text"),
        PanelTarget("7", "7B", "cluster-specific BPNet enrichments on ATAC UMAP", "disease/BPNet", "trained BPNet models, variant catalogs, ATAC cluster bigWigs/peaks", bpnet_resources, "Attempt public resource discovery; do not approximate as reproduced without models/catalogs.", bpnet_status, bpnet_reason, "Fisher enrichment and UMAP overlay"),
        PanelTarget("7", "7C", "fetal-heart control enrichment", "disease/BPNet control", "fetal-heart BPNet/control enhancers and mutation catalogs", "not found", "Search public resources; otherwise mark method-resource gap.", "METHOD_RESOURCE_GAP", "Control model/resources are not in discovered GEO/S3 links.", "control OR/p-value"),
        PanelTarget("7", "7D", "SFARI nearest-gene enrichment", "disease enrichment", "prioritized mutation list and SFARI gene list", "author SFARI list present; mutation priorities not found", "Compute only if prioritized case/control mutations are public.", "METHOD_RESOURCE_GAP", "SFARI genes are present, but prioritized mutations are not yet public in manifest.", "case/control enrichment"),
        PanelTarget("7", "7E", "disrupted motif families", "variant motif disruption", "mutation list, BPNet scores, motif annotation", "Brain_ASD motif/scoring code present" if brain_asd_present else "not found", "Attempt with public mutation/model resources; otherwise mark method gap.", "METHOD_RESOURCE_GAP", "High-effect mutation and model-score table not found.", "motif family overlap"),
        PanelTarget("7", "7F", "NFIA mutation locus", "locus plot", "specific mutation, linked enhancer, BPNet prediction", "Brain_ASD code present but specific scored mutation table not found" if brain_asd_present else "not found", "Generate locus plot only if mutation/model/link resources are public.", "METHOD_RESOURCE_GAP", "Specific high-effect mutation/model evidence not found.", "locus evidence presence"),
        PanelTarget("7", "7G", "NPY mutation locus", "locus plot", "specific mutation, linked enhancer, BPNet prediction", "Brain_ASD code present but specific scored mutation table not found" if brain_asd_present else "not found", "Generate locus plot only if mutation/model/link resources are public.", "METHOD_RESOURCE_GAP", "Specific high-effect mutation/model evidence not found.", "locus evidence presence"),
    ]
    return rows


def _extract_supplemental_targets(pdf_text: str) -> list[PanelTarget]:
    """Create conservative supplemental panel rows from legend labels."""
    rows: list[PanelTarget] = []
    matches = list(re.finditer(r"Figure S(\d+)\.", pdf_text))
    for idx, match in enumerate(matches):
        fig = f"S{match.group(1)}"
        end = matches[idx + 1].start() if idx + 1 < len(matches) else len(pdf_text)
        caption = re.sub(r"\s+", " ", pdf_text[match.start() : end]).strip()
        panel_labels = sorted(set(re.findall(r"\(([A-Z])\)", caption)))
        if not panel_labels:
            panel_labels = ["all"]
        for panel in panel_labels:
            panel_caption = _panel_caption_fragment(caption, panel)
            panel_type, status, reason, inputs, method, metric = _classify_supp_panel(panel_caption)
            rows.append(
                PanelTarget(
                    figure_id=fig,
                    panel_id=f"{fig}{panel}",
                    panel_label=panel_caption[:220],
                    panel_type=panel_type,
                    required_public_inputs=inputs,
                    local_or_public_resources="local PDF legend; public matrices/resources to be resolved by main matrix where applicable",
                    planned_reproduction_method=method,
                    initial_status=status,
                    status_reason=reason,
                    quantitative_check=metric,
                )
            )
    return rows


def _panel_caption_fragment(caption: str, panel: str) -> str:
    if panel == "all":
        return caption[:500]
    pattern = rf"\({re.escape(panel)}\)\s*"
    m = re.search(pattern, caption)
    if not m:
        return caption[:500]
    next_m = re.search(r"\([A-Z]\)\s*", caption[m.end() :])
    end = m.end() + next_m.start() if next_m else min(len(caption), m.start() + 700)
    return caption[m.start() : end].strip()


def _classify_supp_panel(fragment: str) -> tuple[str, str, str, str, str, str]:
    low = fragment.lower()
    if "immunohistochem" in low or "image" in low or "scale bar" in low:
        return (
            "imaging",
            "NON_COMPUTATIONAL",
            "Supplemental microscopy/source image panel; not regenerated from sequencing matrices.",
            "source microscopy images",
            "Document/cite as non-computational unless raw microscopy files are found.",
            "classification evidence",
        )
    if "bpnet" in low or "mutation" in low or "sfari" in low or "conservation" in low:
        return (
            "disease/BPNet supplemental",
            "METHOD_RESOURCE_GAP",
            "Requires BPNet model outputs and/or mutation catalogs not yet present in public manifest.",
            "BPNet model resources, mutation catalogs, motif annotations",
            "Attempt public resource discovery; do not call reproduced without model/catalog evidence.",
            "enrichment or locus metric if inputs found",
        )
    if "external" in low or "bhaduri" in low or "nowakowski" in low or "adult" in low:
        return (
            "external reference validation",
            "PUBLIC_DATA_GAP",
            "Requires external reference dataset not yet in current run manifest.",
            "external public reference dataset plus current labels",
            "Download/locate reference and run projection/score validation.",
            "projection or gene-set concordance",
        )
    if "umap" in low or "heatmap" in low or "dotplot" in low or "expression" in low or "accessibility" in low or "chromvar" in low:
        return (
            "sequencing-derived supplemental",
            "PARTIAL_INPUTS_PENDING",
            "Likely reproducible from public matrices after downloads/loaders complete.",
            "GSE162170 public matrices and metadata",
            "Regenerate after public matrix loader and panel-specific method are validated.",
            "dimension/marker/cluster concordance",
        )
    return (
        "supplemental triage",
        "TRIAGE_PENDING",
        "Panel requires manual mapping after loader/resource inventory is complete.",
        "to be resolved",
        "Manual method mapping required before rendering.",
        "triage completeness",
    )


def _write_outputs(rows: list[PanelTarget], out_dir: Path, metadata: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {"metadata": metadata, "rows": [asdict(r) for r in rows]}
    (out_dir / "figure_target_matrix.json").write_text(json.dumps(payload, indent=2) + "\n")
    fieldnames = list(asdict(rows[0]).keys()) if rows else []
    with (out_dir / "figure_target_matrix.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    counts: dict[str, int] = {}
    for row in rows:
        counts[row.initial_status] = counts.get(row.initial_status, 0) + 1
    md = [
        "# Wave-5 Raw/Public Trevino Figure Target Matrix",
        "",
        f"Created UTC: {metadata['created_at_utc']}",
        f"Run dir: `{metadata['run_dir']}`",
        "",
        "## Initial status counts",
        "",
    ]
    for k in sorted(counts):
        md.append(f"- `{k}`: {counts[k]}")
    md.extend(
        [
            "",
            "## Main figure panel targets",
            "",
            "| Panel | Initial status | Type | Method | Reason |",
            "|---|---|---|---|---|",
        ]
    )
    for row in rows:
        if row.figure_id.startswith("S"):
            continue
        md.append(
            f"| {row.panel_id} | `{row.initial_status}` | {row.panel_type} | "
            f"{row.planned_reproduction_method} | {row.status_reason} |"
        )
    md.extend(
        [
            "",
            "## Supplemental coverage",
            "",
            "Supplemental figures S1-S8 are listed in the CSV/JSON at panel-label granularity when labels were recoverable from the local PDF text. "
            "They start as triage/resource statuses and must be promoted only after rendering plus metric evidence.",
            "",
        ]
    )
    (out_dir / "figure_target_matrix.md").write_text("\n".join(md) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--pdf-text", required=True, type=Path)
    ap.add_argument("--resource-manifest", type=Path)
    ap.add_argument("--author-inventory", type=Path)
    ap.add_argument("--author-download-manifest", type=Path)
    ap.add_argument("--brain-asd-inventory", type=Path)
    args = ap.parse_args()

    pdf_text = args.pdf_text.read_text(errors="replace")
    resource_manifest = _load_json(args.resource_manifest) if args.resource_manifest else {}
    author_inventory = _load_json(args.author_inventory) if args.author_inventory else {}
    author_download_manifest = _load_json(args.author_download_manifest) if args.author_download_manifest else {}
    status = _resource_status(resource_manifest, author_inventory, author_download_manifest)
    brain_asd_present = bool(args.brain_asd_inventory and args.brain_asd_inventory.exists())
    rows = _main_panel_targets(status, brain_asd_present=brain_asd_present)
    rows.extend(_extract_supplemental_targets(pdf_text))
    metadata = {
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "run_dir": str(args.run_dir),
        "pdf_text": str(args.pdf_text),
        "resource_manifest": str(args.resource_manifest) if args.resource_manifest else None,
        "author_inventory": str(args.author_inventory) if args.author_inventory else None,
        "author_download_manifest": str(args.author_download_manifest) if args.author_download_manifest else None,
        "brain_asd_inventory": str(args.brain_asd_inventory) if args.brain_asd_inventory else None,
        "main_panel_count": sum(1 for r in rows if not r.figure_id.startswith("S")),
        "supplemental_panel_count": sum(1 for r in rows if r.figure_id.startswith("S")),
        "status_contract": [
            "REPRODUCED is forbidden at this target-matrix stage.",
            "Promotion requires generated output path plus quantitative check evidence.",
        ],
    }
    _write_outputs(rows, args.run_dir / "python" / "figure_target_matrix", metadata)
    print(json.dumps({"rows": len(rows), **metadata}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
