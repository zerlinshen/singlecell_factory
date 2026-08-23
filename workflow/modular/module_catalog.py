from __future__ import annotations

from dataclasses import dataclass

from .batch_risk import (
    BATCH_STRATEGY_INTEGRATE,
    BATCH_STRATEGY_SINGLE_BATCH,
)


@dataclass(frozen=True)
class ModuleSpec:
    """Stable module metadata used by orchestration, docs, and tests.

    The execution engine still imports module implementations lazily from
    pipeline.py; this catalog deliberately stays dependency-light so CLI/help
    paths can inspect module structure without importing scanpy/liana/scvelo.

    ``runs_after`` is an ORDERING-ONLY hint distinct from ``depends_on``. It
    declares that this module must be sequenced after the named modules WHEN
    BOTH are already in the requested set, but it must never force-pull a
    module into the run (no inclusion auto-pull) and never seed parallel-tier
    in-degree. ``depends_on`` carries hard dependency + auto-include semantics;
    ``runs_after`` only refines order among already-requested modules.
    """

    name: str
    depends_on: tuple[str, ...] = ()
    runs_after: tuple[str, ...] = ()
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
)

# Modality-aware defaults are catalog data, not launcher policy.  Every
# user-facing adapter must consume this mapping (through the accessor below)
# instead of copying module names into its own source file.
MODALITY_OPTIONAL_MODULES: dict[str, tuple[str, ...]] = {
    "rna_only": DEFAULT_OPTIONAL_MODULES,
    "cite_seq": DEFAULT_OPTIONAL_MODULES + ("protein_adt",),
    "spatial": DEFAULT_OPTIONAL_MODULES
    + ("spatial_ingest", "spatial_neighborhoods"),
    "multimodal": DEFAULT_OPTIONAL_MODULES
    + (
        "protein_adt",
        "spatial_ingest",
        "spatial_neighborhoods",
        "multimodal_integration",
    ),
}


@dataclass(frozen=True)
class AnalysisProfile:
    """A named, explicitly declared optional-module plan.

    Distinct from MODALITY_OPTIONAL_MODULES, which answers "what assay is
    this?". A profile answers "what design is this, and what has the operator
    declared about it?" — so it carries the batch declaration alongside the
    module list. ``batch_strategy`` values come from
    :mod:`workflow.modular.batch_risk`, never from a re-typed literal.
    """

    name: str
    optional_modules: tuple[str, ...]
    batch_strategy: str
    description: str


# Named designs. These exist because DEFAULT_OPTIONAL_MODULES contains no
# integration step: running it on multi-batch input is a legitimate choice only
# when it is a choice. Selecting a profile is that choice, recorded in the
# manifest. See batch_risk.py for the measured cost of getting this wrong.
ANALYSIS_PROFILES: dict[str, AnalysisProfile] = {
    "single_batch": AnalysisProfile(
        name="single_batch",
        optional_modules=DEFAULT_OPTIONAL_MODULES,
        batch_strategy=BATCH_STRATEGY_SINGLE_BATCH,
        description=(
            "Canonical defaults with an affirmative single-batch declaration. "
            "Fails at plan time if the input actually spans multiple batches."
        ),
    ),
    "multi_batch_harmony": AnalysisProfile(
        name="multi_batch_harmony",
        # batch_correction sits directly after clustering: it depends_on
        # clustering and OVERWRITES obs["leiden"], so every leiden consumer must
        # be sequenced after it (see the runs_after hints on
        # differential_expression and annotation).
        optional_modules=(
            "clustering",
            "batch_correction",
            "differential_expression",
            "annotation",
        ),
        batch_strategy=BATCH_STRATEGY_INTEGRATE,
        description=(
            "Canonical defaults plus Harmony integration for multi-batch input "
            "(batch key from --batch-key; method from --batch-method)."
        ),
    ),
}


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
            "Per-sample trigger evaluation against 4 QC rules (top50, mt%, "
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
        # Wave-3 option A: when cell_cycle is also requested, run it first so
        # regress_out can affect PCA/neighbors/UMAP/Leiden (Luecken & Theis).
        # Ordering-only — does not auto-include cell_cycle.
        runs_after=("cell_cycle",),
        layer="latent_structure",
        bridge_ready=True,
        description="PCA/neighbors/UMAP/Leiden latent structure.",
    ),
    "cell_cycle": ModuleSpec(
        name="cell_cycle",
        # Wave-3: no longer depends_on clustering. Hard dep is post-doublet
        # expression so scoring/regress can precede layout when both modules
        # are requested (see clustering.runs_after).
        depends_on=("doublet_detection",),
        layer="covariates",
        description=(
            "Cell-cycle scoring and optional regression. Prefer running before "
            "clustering when regress_cell_cycle=True so embeddings are cycle-corrected."
        ),
    ),
    "integration_select": ModuleSpec(
        name="integration_select",
        depends_on=("clustering",),
        layer="latent_structure",
        description=(
            "Per-run discovery integration-selection gate. Scores baseline/"
            "Harmony/scVI/shuffle candidates after clustering and SETS "
            "cfg.batch.method (routing Harmony to the working harmonypy-direct "
            "backend) before batch_correction. Fails loud on degraded "
            "candidates or a non-firing shuffle control. OPT-IN via "
            "--select-integration (runs an expensive scVI seed sweep; "
            "cache-bounded). Never writes inside the factory tree."
        ),
    ),
    "batch_correction": ModuleSpec(
        name="batch_correction",
        depends_on=("clustering",),
        runs_after=("integration_select",),
        layer="latent_structure",
        description=(
            "Batch integration, currently Harmony-first for NC2024. Ordering-only "
            "runs_after integration_select so a gate-chosen method/backend is "
            "applied before correction; does NOT auto-include integration_select."
        ),
    ),
    "differential_expression": ModuleSpec(
        name="differential_expression",
        depends_on=("clustering",),
        runs_after=("batch_correction",),
        layer="markers",
        bridge_ready=True,
        description=(
            "Cluster marker differential expression. Ordering-only runs_after "
            "batch_correction for the same reason annotation carries it: "
            "batch_correction declares provides obs.leiden and overwrites the "
            "labels, so markers computed before it describe pre-correction "
            "clusters that no longer exist in final_adata. Does NOT auto-include "
            "batch_correction."
        ),
    ),
    "annotation": ModuleSpec(
        name="annotation",
        depends_on=("clustering",),
        runs_after=("batch_correction",),
        layer="annotation",
        bridge_ready=True,
        description=(
            "Cluster-voting cell type annotation. Ordering-only runs_after "
            "batch_correction so cell_type binds to post-correction leiden "
            "(see leiden race fix); does NOT auto-include batch_correction."
        ),
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
            "Per-sample VDJ diversity (Shannon + Gini) and sample-local per-cell clonal_expansion class "
            "(keyed by biological sample when obs.sample is present, with _ALL_ fallback). Wave 2B / P2.S16."
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
            "Two-pass peak-to-gene linkage. (1) Distance-based assignment of "
            "ATAC peaks to gene TSSs within --peak-to-gene-window bp; writes "
            "adata.uns['peak_to_gene']. (2) Pearson + permutation-FDR "
            "scoring over the joint ATAC+RNA cell axis (Trevino 2021 / Ma "
            "2020 SHARE-seq; n_perms default 100, sparse-axis preserved); "
            "writes adata.uns['peak_to_gene_top1000'] and "
            "['peak_to_gene_linkages'] consumed by tf_network. NOT Cicero "
            "co-accessibility (see __references__ in module). Wave 2B / "
            "P2.S13 + Wave 5 / US-W5-6."
        ),
    ),
    "scatac_pseudobulk_da": ModuleSpec(
        name="scatac_pseudobulk_da",
        depends_on=("atac_ingest",),
        layer="analysis",
        modality="atac",
        bridge_ready=True,
        description=(
            "scATAC sample-level pseudobulk differential accessibility using DESeq2 "
            "primary and edgeR QL F-test cross-check."
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


def module_runs_after() -> dict[str, set[str]]:
    """Return ordering-only ``runs_after`` hints keyed by module name.

    Distinct from :func:`module_dependencies`: these edges only refine
    execution order among already-requested modules and never auto-include.
    """

    return {
        name: set(spec.runs_after)
        for name, spec in MODULE_SPECS.items()
        if spec.runs_after
    }


def optional_module_names() -> tuple[str, ...]:
    """Return every non-mandatory module in stable catalog order."""

    mandatory = set(MANDATORY_MODULES)
    return tuple(name for name in MODULE_SPECS if name not in mandatory)


def optional_modules_for_modality(modality: str) -> tuple[str, ...]:
    """Return the canonical optional-module defaults for ``modality``.

    Unknown modality names are rejected so adapters cannot silently relabel an
    unsupported input as RNA-only while presenting the result as auto-detected.
    """

    try:
        return MODALITY_OPTIONAL_MODULES[modality]
    except KeyError as exc:
        choices = ", ".join(MODALITY_OPTIONAL_MODULES)
        raise ValueError(
            f"unknown modality {modality!r}; expected one of: {choices}"
        ) from exc


def analysis_profile_names() -> tuple[str, ...]:
    """Return every named analysis profile in stable catalog order."""

    return tuple(ANALYSIS_PROFILES)


def analysis_profile(name: str) -> AnalysisProfile:
    """Return the named analysis profile.

    Unknown names are rejected so an adapter cannot advertise a profile that has
    no catalog entry, or silently fall back to the unintegrated defaults.
    """

    try:
        return ANALYSIS_PROFILES[name]
    except KeyError as exc:
        choices = ", ".join(ANALYSIS_PROFILES)
        raise ValueError(
            f"unknown analysis profile {name!r}; expected one of: {choices}"
        ) from exc


def optional_modules_for_profile(name: str) -> tuple[str, ...]:
    """Return the canonical optional-module list for a named profile."""

    return analysis_profile(name).optional_modules


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
