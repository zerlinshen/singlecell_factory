from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


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


@dataclass
class ClusteringConfig:
    """Configuration for dimension reduction and clustering."""

    n_top_genes: int = 3000
    target_sum: float = 1e4
    n_pcs: int = 40
    n_neighbors: int = 15
    leiden_resolution: float = 0.8
    random_state: int = 0
    scale_data: bool = False


@dataclass
class BatchConfig:
    """Configuration for batch correction."""

    batch_key: str = "sample"
    method: str = "harmony"  # harmony, bbknn, combat, scanorama, scvi, mnn, fastmnn
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
    regress_cell_cycle: bool = False
    trajectory_root_cluster: str | None = None
    checkpoint: bool = False
    resume_from: str | None = None
    parallel_workers: int = 1
    de_method: str = "wilcoxon"
    de_n_genes: int = 300
    de_pval_threshold: float = 0.05
    de_logfc_threshold: float = 0.25
    annotation_confidence_threshold: float = 0.1
    reference_adata: Path | None = None
    reference_label_key: str = "cell_type"
    reference_k: int = 15
    reference_min_confidence: float = 0.6
    reference_override_mode: str = "conservative"  # conservative, all
    gpu_mode: str = "auto"  # auto, off, force
    scale_mode: str = "standard"  # standard, large, massive — kept as preset bundle for backwards compat
    # Capability flags — set explicitly or expanded from scale_mode via scale_mode_to_capabilities()
    lazy_read: str = "auto"          # auto, true, false
    doublet_strategy: str = "auto"   # auto, grouped, whole, skip
    clustering_engine: str = "auto"  # auto, sparse_exact, css, gpu
    checkpoint_policy: str = "full"  # full, mandatory_only, metadata_only


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
    "massive": {
        "lazy_read": "true",
        "doublet_strategy": "grouped",
        "clustering_engine": "css",
        "checkpoint_policy": "mandatory_only",
    },
}


def scale_mode_to_capabilities(scale_mode: str) -> dict[str, str]:
    """Return the capability flag bundle for a given scale_mode preset name.

    Returns "standard" bundle if the preset is unknown, so old configs degrade safely.
    """
    return dict(_SCALE_MODE_PRESETS.get(scale_mode, _SCALE_MODE_PRESETS["standard"]))
