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
    "doublet_detection": ModuleSpec(
        name="doublet_detection",
        depends_on=("qc",),
        layer="quality_control",
        description="Scrublet-style doublet scoring/removal.",
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
        depends_on=("clustering",),
        layer="genomic_optional",
        modality="copy_number",
        description="CNV inference lane; keep scale-safe before full cohorts.",
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
