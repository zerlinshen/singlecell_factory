from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class CellRangerConfig:
    """Configuration for the mandatory Cell Ranger stage."""

    sample_root: Path
    outs_dir: Path
    fastq_dir: Path | None = None
    transcriptome_dir: Path | None = None
    sample_id: str = "lusc"
    localcores: int = 8
    localmem: int = 64
    force_run: bool = False
    run_if_missing: bool = True


@dataclass
class QCConfig:
    """Configuration for QC stage."""

    min_genes: int = 200
    max_genes: int = 7000
    min_counts: int = 500
    max_counts: int = 50000
    max_mito_pct: float = 20.0
    max_ribo_pct: float = 50.0
    min_cells: int = 3


@dataclass
class DoubletConfig:
    """Configuration for doublet detection stage."""

    expected_doublet_rate: float = 0.06
    remove_doublets: bool = True
    n_prin_comps: int = 30  # Scrublet internal PCA; paper uses 30


@dataclass
class ClusteringConfig:
    """Configuration for dimension reduction and clustering."""

    n_top_genes: int = 3000
    target_sum: float = 1e4
    n_pcs: int = 15  # paper: 15-PC Harmony space for Leiden clustering
    n_neighbors: int = 15
    leiden_resolution: float = 1.0  # paper: Leiden resolution=1.0
    random_state: int = 0
    scale_data: bool = False


@dataclass
class BatchConfig:
    """Configuration for batch correction."""

    batch_key: str = "sample"
    method: str = "harmony"  # harmony, bbknn, combat, scanorama, scvi, mnn, fastmnn
    # Harmony parameters (Korsunsky et al. 2019, Nat Methods). The defaults
    # are the upstream library defaults; expose them in the manifest so
    # cross-run drift in Harmony tuning is auditable.
    harmony_theta: float = 2.0
    harmony_sigma: float = 0.1
    # Plan F-Harmony / Principle 9 (2026-05-20): default bumped 10 -> 50.
    # Empirical: NC2024 877k tumor cohort with 75 batches does not converge at
    # 10 iters (G-C0 v2 launch run 2026-05-20 emitted "Harmony did not converge"
    # warning). Korsunsky 2019 recommends <= 100 max; 50 is the safe default
    # for large multi-batch cohorts and matches the harmonypy convergence
    # studies cited in the original paper.
    harmony_max_iter: int = 50
    # Plan F-Harmony / Principle 8: harmony backend selection. "auto" = GPU
    # via rapids-singlecell if available, else CPU harmonypy. "cpu" = explicit
    # scanpy_external.pp.harmony_integrate. "gpu" = explicit rsc (raises if
    # unavailable). G-C0 v3 contract requires the chosen backend(s) to
    # converge; non-convergence raises unless SC_ALLOW_HARMONY_NON_CONVERGENCE=1.
    harmony_backend: str = "auto"
    scvi_max_epochs: int = 200
    scvi_n_latent: int = 30
    scvi_early_stopping: bool = True


@dataclass
class CNVConfig:
    """Configuration for CNV inference."""

    reference_group: str | None = None  # e.g., "Fibroblast" for normal reference
    window_size: int = 100
    malignant_percentile: float = 75.0


@dataclass
class VelocityConfig:
    """Configuration for RNA velocity."""

    loom_path: Path | None = None
    bam_path: Path | None = None
    gtf_path: Path | None = None
    mode: str = "stochastic"  # stochastic, dynamical
    n_jobs: int = 4
    min_shared_counts: int = 20  # scVelo filter_and_normalize threshold
    n_pcs: int = 30  # PCA components for scVelo moments
    n_neighbors: int = 30  # neighbors for scVelo moments


@dataclass
class CbioPortalConfig:
    """Configuration for the optional validate_cbioportal module."""

    genes: list[str] = field(default_factory=list)
    study_id: str = "lusc_tcga_pan_can_atlas_2018"
    use_de_genes: bool = True
    top_n_de_genes: int = 20
    timeout: int = 30


@dataclass
class GeneSignatureConfig:
    """Configuration for gene signature scoring module."""

    signature_json: Path | None = None  # Optional user-provided signatures JSON
    use_builtin: bool = True  # Include built-in cancer signatures


@dataclass
class PaperReproConfig:
    """Configuration for the optional paper_repro module."""

    spec_json: Path | None = None
    strict: bool = False


@dataclass
class DEConfig:
    """Paper-aligned DE marker parameters."""

    marker_min_pct: float = 0.10       # paper: 0.30 — only paper-aligned launcher sets this
    logfc_threshold: float = 0.25      # paper: 0.0 — only paper-aligned launcher sets this
    marker_correction: str = "benjamini-hochberg"  # paper uses bonferroni; default stays BH


@dataclass
class PseudobulkDEConfig:
    """Paper-aligned pseudobulk DE thresholds."""

    padj_threshold: float = 0.05       # paper: median(padj) <= 0.05
    abs_logfc_threshold: float = 1.0   # paper: |median(logFC)| >= 1


@dataclass
class PseudobulkConfig:
    """Configuration for pseudobulk differential expression."""

    sample_col: str | None = None
    group_col: str = "cell_type"
    contrast_col: str | None = None
    contrast_a: str | None = None
    contrast_b: str | None = None
    contrast_json: Path | None = None
    exploratory_group_vs_rest: bool = False
    min_cells_per_sample: int = 3
    min_samples_per_condition: int = 2


@dataclass
class PipelineConfig:
    """Top-level modular workflow configuration."""

    project: str
    output_dir: Path
    cellranger: CellRangerConfig
    qc: QCConfig = field(default_factory=QCConfig)
    doublet: DoubletConfig = field(default_factory=DoubletConfig)
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    batch: BatchConfig = field(default_factory=BatchConfig)
    cnv: CNVConfig = field(default_factory=CNVConfig)
    velocity: VelocityConfig = field(default_factory=VelocityConfig)
    optional_modules: list[str] = field(
        default_factory=lambda: [
            "clustering",
            "differential_expression",
            "annotation",
            "trajectory",
            "pseudo_velocity",
        ]
    )
    markers: dict[str, list[str]] = field(default_factory=dict)
    cbioportal: CbioPortalConfig = field(default_factory=CbioPortalConfig)
    gene_signature: GeneSignatureConfig = field(default_factory=GeneSignatureConfig)
    paper_repro: PaperReproConfig = field(default_factory=PaperReproConfig)
    pseudobulk: PseudobulkConfig = field(default_factory=PseudobulkConfig)
    de_config: DEConfig = field(default_factory=DEConfig)
    pseudobulk_de: PseudobulkDEConfig = field(default_factory=PseudobulkDEConfig)
    regress_cell_cycle: bool = False
    trajectory_root_cluster: str | None = None
    checkpoint: bool = False
    resume_from: str | None = None
    parallel_workers: int = 1
    de_method: str = "wilcoxon"
    de_n_genes: int = 300
    de_pval_threshold: float = 0.05
    de_logfc_threshold: float = 0.25
    de_correction: str = "benjamini-hochberg"  # paper: bonferroni; default stays BH for backward compat
    de_min_pct: float = 0.10                   # paper: 0.30; only paper-aligned launcher overrides
    annotation_confidence_threshold: float = 0.1
    reference_adata: Path | None = None
    reference_label_key: str = "cell_type"
    reference_k: int = 15
    reference_min_confidence: float = 0.6
    reference_override_mode: str = "conservative"  # conservative, all
    random_state: int = 42  # global seed propagated to all stochastic modules via ctx.random_state
    gpu_mode: str = "auto"  # auto, off, force
    scale_mode: str = "standard"  # standard, large, massive — kept as preset bundle for backwards compat
    # Capability flags — set explicitly or expanded from scale_mode via scale_mode_to_capabilities()
    lazy_read: str = "auto"          # auto, true, false
    doublet_strategy: str = "auto"   # auto, grouped, whole, skip
    clustering_engine: str = "auto"  # auto, sparse_exact, css, gpu
    checkpoint_policy: str = "full"  # full
    # Wave 3 (US-W3-2). Policy for GPU clustering failure post-host-mutation.
    # raise (default): poison adata + raise ClusteringContractViolation.
    # restore-cpu: restore from adata.raw and route to CPU clustering (M2 only).
    # reload-checkpoint: reload adata from on-disk h5ad checkpoint (requires --checkpoint).
    # Overridable via SC_GPU_FAILURE_POLICY env var.
    gpu_failure_policy: str = "raise"
    # Cohort subset: obs_col=val1,val2 filter applied after loading (supports list for AND-chaining)
    cohort_subset: Optional[list[str]] = None
    annotation_strategy: str = "cluster_voting"  # cluster_voting, cell_argmax
    # Marker intelligence (P1A)
    tissue: str = "lung"
    condition: str = "NSCLC"
    validate_context: bool = False
    context_mismatch_threshold: float = 0.3
    context_min_cells: int = 20


# Maps scale_mode preset names to their capability flag bundles.
# "massive" = the canonical preset used by the NC2024 launch script.
_SCALE_MODE_PRESETS: dict[str, dict[str, str]] = {
    "standard": {
        "lazy_read": "auto",
        "doublet_strategy": "auto",
        "clustering_engine": "auto",
        "checkpoint_policy": "full",
    },
    "large": {
        "lazy_read": "auto",
        "doublet_strategy": "auto",
        "clustering_engine": "auto",
        "checkpoint_policy": "full",
    },
    # F-4 (Plan ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md, Principle 2):
    # `massive` preset previously routed clustering_engine -> "css". CSS is now
    # removed from the production science path; the preset is REMAPPED to
    # route clustering_engine -> "auto" instead. The operational flags
    # (lazy_read=true, doublet_strategy=grouped, checkpoint_policy=full) are
    # preserved — they are I/O / dispatch options, not scientific compromises.
    "massive": {
        "lazy_read": "true",
        "doublet_strategy": "grouped",
        "clustering_engine": "auto",
        "checkpoint_policy": "full",
    },
}


def scale_mode_to_capabilities(scale_mode: str) -> dict[str, str]:
    """Return the capability flag bundle for a given scale_mode preset name.

    Returns "standard" bundle if the preset is unknown, so old configs degrade safely.
    """
    return dict(_SCALE_MODE_PRESETS.get(scale_mode, _SCALE_MODE_PRESETS["standard"]))
