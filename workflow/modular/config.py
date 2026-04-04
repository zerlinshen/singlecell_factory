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
    method: str = "harmony"  # harmony, bbknn, combat


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
    regress_cell_cycle: bool = False
    trajectory_root_cluster: str | None = None
    checkpoint: bool = False
    resume_from: str | None = None
    parallel_workers: int = 1
