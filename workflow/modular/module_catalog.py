from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ModuleSpec:
    """Stable module metadata used by orchestration, docs, and tests.

    The execution engine still imports module implementations lazily from
    pipeline.py; this catalog deliberately stays dependency-light so CLI/help
    paths can inspect module structure without importing scanpy/liana/scvelo.
    """

    name: str
    depends_on: tuple[str, ...] = ()
    layer: str = "analysis"
    modality: str = "singlecell_rna"
    owner: str = "singlecell_factory"
    bridge_ready: bool = False
    description: str = ""


MANDATORY_MODULES: tuple[str, ...] = (
    "cellranger",
    "qc",
    "ambient_correction",
    "doublet_detection",
)

DEFAULT_OPTIONAL_MODULES: tuple[str, ...] = (
    "clustering",
    "differential_expression",
    "annotation",
    "trajectory",
    "pseudo_velocity",
)


MODULE_SPECS: dict[str, ModuleSpec] = {
    "cellranger": ModuleSpec(
        name="cellranger",
        layer="ingest",
        description="Load prepared matrices or run Cell Ranger when requested.",
    ),
    "qc": ModuleSpec(
        name="qc",
        depends_on=("cellranger",),
        layer="quality_control",
        description="Cell/gene QC filtering and QC figures.",
    ),
    "ambient_correction": ModuleSpec(
        name="ambient_correction",
        depends_on=("qc",),
        layer="quality_control",
        description=(
            "Conditional DecontX ambient RNA correction (Yang 2020 Genome Biology). "
            "Per-sample trigger evaluation against 5 QC rules (top50, mt%, doublet rate, "
            "total/n_genes correlation, cross-sample housekeeping CV); self-skips when "
            "no trigger fires. R driver runs in r_multiomics conda env via subprocess. "
            "Policy: ops/policy/ambient_correction_policy.md Phase 2 (2026-05-20). "
            "Opt-out: SC_AMBIENT_TRIGGERS_DISABLE=1 or --ambient-disable-triggers."
        ),
    ),
    "doublet_detection": ModuleSpec(
        name="doublet_detection",
        depends_on=("qc", "ambient_correction"),
        layer="quality_control",
        description=(
            "Scrublet-style doublet scoring/removal. Runs on counts AFTER ambient "
            "correction (when triggered) so doublet calls use cleaned signal."
        ),
    ),
    "clustering": ModuleSpec(
        name="clustering",
        depends_on=("doublet_detection",),
        layer="latent_structure",
        bridge_ready=True,
        description="PCA/neighbors/UMAP/Leiden latent structure.",
    ),
    "cell_cycle": ModuleSpec(
        name="cell_cycle",
        depends_on=("clustering",),
        layer="covariates",
        description="Cell-cycle scoring and optional regression inputs.",
    ),
    "batch_correction": ModuleSpec(
        name="batch_correction",
        depends_on=("clustering",),
        layer="latent_structure",
        description="Batch integration, currently Harmony-first for NC2024.",
    ),
    "differential_expression": ModuleSpec(
        name="differential_expression",
        depends_on=("clustering",),
        layer="markers",
        bridge_ready=True,
        description="Cluster marker differential expression.",
    ),
    "annotation": ModuleSpec(
        name="annotation",
        depends_on=("clustering",),
        layer="annotation",
        bridge_ready=True,
        description="Cluster-voting cell type annotation.",
    ),
    "trajectory": ModuleSpec(
        name="trajectory",
        depends_on=("clustering",),
        layer="state_dynamics",
        description="Trajectory and pseudotime analysis.",
    ),
    "pseudo_velocity": ModuleSpec(
        name="pseudo_velocity",
        depends_on=("trajectory",),
        layer="state_dynamics",
        description="Trajectory-derived pseudo-velocity summaries.",
    ),
    "rna_velocity": ModuleSpec(
        name="rna_velocity",
        depends_on=("clustering",),
        layer="state_dynamics",
        modality="singlecell_rna_splicing",
        description="RNA velocity when spliced/unspliced modality is available.",
    ),
    "cnv_inference": ModuleSpec(
        name="cnv_inference",
        depends_on=("annotation",),
        layer="genomic_optional",
        modality="copy_number",
        description=(
            "CNV inference lane; runs after annotation so epithelial/tumor calls "
            "and normal-reference groups are auditable before CNV classification."
        ),
    ),
    "pathway_analysis": ModuleSpec(
        name="pathway_analysis",
        depends_on=("differential_expression",),
        layer="biology",
        description="Pathway enrichment from marker outputs.",
    ),
    "cell_communication": ModuleSpec(
        name="cell_communication",
        depends_on=("annotation",),
        layer="biology",
        description="Ligand-receptor communication summaries.",
    ),
    "gene_regulatory_network": ModuleSpec(
        name="gene_regulatory_network",
        depends_on=("clustering",),
        layer="biology",
        description="Gene regulatory network inference.",
    ),
    "validate_cbioportal": ModuleSpec(
        name="validate_cbioportal",
        depends_on=("differential_expression",),
        layer="external_validation",
        modality="external_cancer_genomics",
        description="cBioPortal validation of selected genes.",
    ),
    "immune_phenotyping": ModuleSpec(
        name="immune_phenotyping",
        depends_on=("annotation",),
        layer="biology",
        bridge_ready=True,
        description="Immune subtype scoring and visual summaries.",
    ),
    "tumor_microenvironment": ModuleSpec(
        name="tumor_microenvironment",
        depends_on=("annotation",),
        layer="biology",
        bridge_ready=True,
        description="TME subtype/signature summaries.",
    ),
    "gene_signature_scoring": ModuleSpec(
        name="gene_signature_scoring",
        depends_on=("clustering",),
        layer="biology",
        bridge_ready=True,
        description="Gene set/signature scoring.",
    ),
    "evolution": ModuleSpec(
        name="evolution",
        depends_on=("cnv_inference", "trajectory"),
        layer="genomic_optional",
        modality="copy_number_state_dynamics",
        description="Evolution summaries requiring CNV and trajectory evidence.",
    ),
    "pseudobulk_de": ModuleSpec(
        name="pseudobulk_de",
        depends_on=("differential_expression",),
        layer="markers",
        bridge_ready=True,
        description="Pseudobulk DE with explicit contrast contract.",
    ),
    "cell_fate": ModuleSpec(
        name="cell_fate",
        depends_on=("trajectory",),
        layer="state_dynamics",
        description="Cell fate summaries from trajectory outputs.",
    ),
    "composition": ModuleSpec(
        name="composition",
        depends_on=("annotation",),
        layer="reporting",
        bridge_ready=True,
        description="Cell composition tables and plots.",
    ),
    "metacell": ModuleSpec(
        name="metacell",
        depends_on=("clustering",),
        layer="aggregation",
        bridge_ready=True,
        description="Metacell aggregation for scale-safe summaries.",
    ),
    "paper_repro": ModuleSpec(
        name="paper_repro",
        depends_on=("clustering",),
        layer="reproduction",
        bridge_ready=True,
        description="Paper finding/parity reproduction checks.",
    ),
    "protein_adt": ModuleSpec(
        name="protein_adt",
        depends_on=("qc",),
        layer="ingest",
        modality="protein_adt",
        bridge_ready=True,
        description="CITE-seq / ADT QC + CLR normalization (bundle v2.1 protein extension).",
    ),
    "spatial_ingest": ModuleSpec(
        name="spatial_ingest",
        depends_on=("qc",),
        layer="spatial",
        modality="spatial_transcriptomics",
        bridge_ready=True,
        description=(
            "Spatial transcriptomics ingest: attach (x,y) coords + library "
            "metadata for Visium/Xenium/MERFISH (bundle v2.1 spatial extension)."
        ),
    ),
    "spatial_neighborhoods": ModuleSpec(
        name="spatial_neighborhoods",
        depends_on=("spatial_ingest",),
        layer="spatial",
        modality="spatial_transcriptomics",
        bridge_ready=True,
        description=(
            "Spatial neighborhood analytics (Moran's I, neighborhood "
            "enrichment, co-occurrence) via squidpy when available."
        ),
    ),
    "multimodal_integration": ModuleSpec(
        name="multimodal_integration",
        depends_on=("clustering",),
        layer="multimodal",
        modality="multimodal_joint_embedding",
        bridge_ready=True,
        description=(
            "[EXPERIMENTAL] Multimodal joint embedding via Seurat WNN (R "
            "subprocess) or MOFA (muon/mofapy2). Off by default; emits "
            "bundle v2.1 multimodal_obsm extension."
        ),
    ),
    "marker_db_loader": ModuleSpec(
        name="marker_db_loader",
        depends_on=(),
        layer="annotation",
        modality="singlecell_rna",
        bridge_ready=False,
        description=(
            "Load offline-vendored marker DBs (CellMarker2, PanglaoDB, CellTypist, scTypeDB) "
            "for (tissue, condition) routing. Writes adata.uns['marker_db_index']. "
            "No clustering dependency — consumed by context_aware_annotation."
        ),
    ),
    "context_aware_annotation": ModuleSpec(
        name="context_aware_annotation",
        depends_on=("clustering", "marker_db_loader"),
        layer="annotation",
        modality="singlecell_rna",
        bridge_ready=True,
        description=(
            "Context-aware cell type annotation driven by (tissue, condition) marker DBs. "
            "Writes adata.obs['context_aware_celltype'] and optionally "
            "adata.obs['context_aware_substate']. Opt-in; leaves legacy annotation.py untouched (AC-10)."
        ),
    ),
    "modality_registry": ModuleSpec(
        name="modality_registry",
        depends_on=(),
        layer="ingest",
        modality="agnostic",
        bridge_ready=False,
        description=(
            "Detect present modalities (RNA/ATAC/VDJ/Ribo/Hi-C/protein/spatial) by scanning "
            "obsm keys and obs columns. Writes adata.uns['modalities_present'] and "
            "runs/<run-id>/manifest/modalities.json. Wave 2 / P1B.S8."
        ),
    ),
    "cross_modality_qc": ModuleSpec(
        name="cross_modality_qc",
        depends_on=("modality_registry",),
        layer="quality_control",
        modality="agnostic",
        bridge_ready=False,
        description=(
            "Barcode overlap QC across modalities. Tri-state status: "
            "skipped_single_modality | ran | warn_divergence. Sparse-safe (no full-X densification). "
            "Wave 2 / P1B.S9."
        ),
    ),
    "atac_ingest": ModuleSpec(
        name="atac_ingest",
        depends_on=(),
        layer="ingest",
        modality="atac",
        bridge_ready=True,
        description=(
            "ATAC-seq sparse peak matrix ingest with TF-IDF + LSI embedding "
            "(Cusanovich 2018 methodology). Memory-safe: never densifies the peak matrix. "
            "Writes sparse adata.obsm['atac_peaks'], compatibility adata.obsm['X_atac'], "
            "and adata.uns['atac_peaks']. Wave 2 / P2.S11."
        ),
    ),
    "vdj_ingest": ModuleSpec(
        name="vdj_ingest",
        depends_on=(),
        layer="ingest",
        modality="vdj",
        bridge_ready=True,
        description=(
            "VDJ receptor ingest: parse cellranger-vdj filtered_contig_annotations.csv "
            "into per-cell clonotype_id + chain_pairing. Memory-safe (pandas-only, no AnnData "
            "duplication). Wave 2B / P2.S15."
        ),
    ),
    "vdj_metrics": ModuleSpec(
        name="vdj_metrics",
        depends_on=("vdj_ingest",),
        layer="metrics",
        modality="vdj",
        bridge_ready=True,
        description=(
            "Per-sample VDJ diversity (Shannon + Gini) and per-cell clonal_expansion class. "
            "Wave 2B / P2.S16."
        ),
    ),
    "atac_qc": ModuleSpec(
        name="atac_qc",
        depends_on=(),
        layer="quality_control",
        modality="atac",
        bridge_ready=False,
        description=(
            "Per-cell ATAC QC: TSS enrichment + FRiP from fragments.tsv.gz, ENCODE-style "
            "pass/warn/fail thresholds. Memory-safe (streams fragments). Wave 2B / P2.S12."
        ),
    ),
    "atac_lsi": ModuleSpec(
        name="atac_lsi",
        depends_on=("atac_ingest",),
        layer="embedding",
        modality="atac",
        bridge_ready=False,
        description=(
            "TF-IDF normalisation + truncated SVD (LSI) on adata.obsm['atac_peaks']. "
            "Writes adata.obsm['X_lsi'] (49 components, first dropped) and "
            "adata.uns['lsi_variance_explained']. Wave 5 / US-W5-1."
        ),
    ),
    "peak_to_gene": ModuleSpec(
        name="peak_to_gene",
        depends_on=("atac_ingest", "atac_lsi"),
        layer="annotation_prep",
        modality="atac",
        bridge_ready=False,
        description=(
            "Distance-based linkage of ATAC peaks to nearest gene TSS (Cicero-style). "
            "Writes adata.uns['peak_to_gene'] for joint RNA+ATAC analysis. Wave 2B / P2.S13."
        ),
    ),
    "hic_ingest": ModuleSpec(
        name="hic_ingest",
        depends_on=(),
        layer="ingest",
        modality="hic",
        bridge_ready=True,
        description=(
            "Hi-C / scHi-C contact map ingest. Reads .cool/.mcool (via cooler) or "
            "TSV contact-pair format. Stores contacts as scipy.sparse CSR. Wave 2B / P2.S18."
        ),
    ),
    "hic_tad": ModuleSpec(
        name="hic_tad",
        depends_on=("hic_ingest",),
        layer="annotation_prep",
        modality="hic",
        bridge_ready=True,
        description=(
            "TAD boundary detection via insulation score (Crane 2015) + A/B compartment "
            "scoring via first eigenvector of correlation matrix (Lieberman-Aiden 2009). "
            "Wave 2B / P2.S19."
        ),
    ),
    "ribo_ingest": ModuleSpec(
        name="ribo_ingest",
        depends_on=(),
        layer="ingest",
        modality="ribo",
        bridge_ready=True,
        description=(
            "Ribosome profiling (Ribo-seq) ingest. Computes per-gene per-sample "
            "translation efficiency (footprint / RNA + pseudocount). Memory-safe "
            "(sparse-aware RNA aggregation). Wave 2B / P2.S21."
        ),
    ),
}


def module_dependencies() -> dict[str, set[str]]:
    """Return dependency DAG in the legacy pipeline shape."""

    return {name: set(spec.depends_on) for name, spec in MODULE_SPECS.items()}


def optional_module_names() -> tuple[str, ...]:
    """Return every non-mandatory module in stable catalog order."""

    mandatory = set(MANDATORY_MODULES)
    return tuple(name for name in MODULE_SPECS if name not in mandatory)


def module_help_list() -> str:
    """Comma-separated optional module list for CLI help and docs."""

    return ", ".join(optional_module_names())


def modules_by_layer() -> dict[str, tuple[str, ...]]:
    """Group module names by architectural layer in catalog order."""

    grouped: dict[str, list[str]] = {}
    for name, spec in MODULE_SPECS.items():
        grouped.setdefault(spec.layer, []).append(name)
    return {layer: tuple(names) for layer, names in grouped.items()}


def bridge_ready_modules() -> tuple[str, ...]:
    """Modules whose outputs are intended for R/report bundle consumption."""

    return tuple(name for name, spec in MODULE_SPECS.items() if spec.bridge_ready)
