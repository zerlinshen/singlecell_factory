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
    # Cell-count-scaled expected doublet rate. 10x Chromium multiplet rate is
    # ~0.8% per 1000 cells loaded (10x v3 User Guide; scDblFinder, Germain 2021),
    # so a fixed 6% is wrong off ~6-7k cells. When True, per-sample (per-batch when
    # grouped) rate = min(0.08, 0.008 * n_obs/1000); expected_doublet_rate above
    # acts as the floor when scaling is disabled.
    scale_expected_doublet_rate: bool = True
    remove_doublets: bool = True
    n_prin_comps: int = 30  # Scrublet internal PCA; paper uses 30
    backend: str = "scrublet"  # scrublet, doubletfinder, scdblfinder, consensus
    consensus_logic: str = "or"  # or, and, rank
    consensus_pair: str = "scrublet_doubletfinder"
    r_conda_env: str = "r_multiomics"
    subprocess_timeout: int = 1800
    doubletfinder_pn: float = 0.25
    doubletfinder_pk: float = 0.09
    doubletfinder_pcs: int = 20
    scdblfinder_samples_col: str | None = None


@dataclass
class AmbientCorrectionConfig:
    """Configuration for conditional DecontX ambient RNA correction.

    Policy: ops/policy/ambient_correction_policy.md (Phase 2 decided 2026-05-20).
    Tool: DecontX (Yang et al. 2020, Genome Biology). R driver invoked from the
    r_multiomics conda env. Defaults below match policy factory thresholds.
    """

    # Activation switches
    disable_triggers: bool = False  # SC_AMBIENT_TRIGGERS_DISABLE=1 also forces skip
    dry_run_triggers_only: bool = False  # evaluate triggers only, no R call

    # Trigger thresholds (per ambient_correction_policy.md Phase 2).
    # NOTE: pct_counts_in_top_50 and pct_counts_mt follow scanpy convention
    # (PERCENT scale 0-100), not fraction. Defaults below are PERCENT.
    # T3 (doublet-excess) removed 2026-05-20 because the DAG places
    # ambient_correction BEFORE doublet_detection, leaving the doublet rate
    # unevaluable at trigger time.
    trigger_top50: float = 50.0              # T1: median pct_counts_in_top_50 > X% (scanpy percent)
    trigger_mt: float = 15.0                 # T2: median pct_counts_mt > X% (scanpy percent)
    trigger_count_correlation: float = 0.85  # T4: Spearman(total, n_genes) < X (fraction)
    trigger_cross_sample_cv: float = 0.50    # T5: cohort housekeeping CV > X (fraction)

    # R subprocess plumbing
    r_conda_env: str = "r_multiomics"
    decontx_max_iter: int = 200
    decontx_seed: int = 42
    subprocess_timeout: int = 1800           # 30 min per sample worst-case
    batch_obs_column: str | None = None      # if set, decontX --batch <col>
    keep_temp_on_failure: bool = True        # preserve temp dir if R fails


@dataclass
class ClusteringConfig:
    """Configuration for dimension reduction and clustering."""

    n_top_genes: int = 3000
    target_sum: float = 1e4
    n_pcs: int = 15  # paper: 15-PC Harmony space for Leiden clustering
    n_neighbors: int = 15
    leiden_resolution: float = 1.0  # paper: Leiden resolution=1.0
    leiden_resolution_sweep: tuple[float, ...] = ()
    random_state: int = 0
    scale_data: bool = False
    # HVG flavor (configurable; previously hardcoded "seurat"). Default stays
    # "seurat" because the suite's own hgmm species-mixing GT shows the
    # dispersion-based "seurat" flavor recovers the 2-population structure at
    # least as well as "seurat_v3" (ARI 0.293 vs ~0.27 at res=0.5; 0.173 vs
    # 0.127 at res=1.0 — real-run 2026-06-27:
    # governance/external_references/audits/round5_realrun_validation_singlecell_2026-06-27.md,
    # corroborating clustering_scoreboard.csv). "seurat_v3" (variance-stabilizing
    # HVG on RAW COUNTS; Stuart 2019, Heumos 2023 Nat Rev Genet) is exposed as an
    # option and is generally preferred for complex multi-batch tissue, but it is
    # NOT the GT-validated default on the local benchmark, so it is opt-in. When
    # "seurat_v3" is selected it MUST run on a raw-count layer; if none is
    # available at HVG time the call falls back to "seurat" and records a loud
    # hvg_flavor_fallback status (no silent log-data use).
    hvg_flavor: str = "seurat"
    # Leiden iterations. -1 = run to convergence (Traag 2019). The igraph default
    # (2 passes) may not converge and differs from the rapids GPU default,
    # breaking CPU/GPU parity; pinning this aligns both lanes.
    leiden_n_iterations: int = -1


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
    # Per-run discovery integration-selection gate (integration_select module,
    # plan per-run-discovery-integration-gate.md). OPT-IN only: the gate runs an
    # expensive scVI seed sweep, so it is never default-on. When enabled, the
    # gate requires baseline/harmony/scVI/shuffle candidates, fails loud on
    # degraded candidates or a non-firing shuffle control, and SETS `method`
    # (and, if it picks harmony, `harmony_backend="direct"`) before
    # batch_correction runs.
    select_integration: bool = False
    # Margin the candidate batch-mixing must beat baseline by (gate 1). Mirrors
    # select_integration.MARGIN_MIX; exposed so a run can tune the gate strictness.
    integration_margin_mix: float = 0.05
    # scVI seed sweep for the gate (>=3 seeds required for the band-aware tie).
    integration_scvi_seeds: tuple[int, ...] = (0, 1, 2)


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
    ambient: AmbientCorrectionConfig = field(default_factory=AmbientCorrectionConfig)
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
    # Hi-C / single-cell 3D genome optional vertical slice. These fields are
    # consumed by hic_ingest/hic_tad when those optional modules are selected.
    hic_contacts_path: Path | None = None
    hic_resolution_bp: int = 25_000
    hic_chromsizes_path: Path | None = None
    hic_tad_window_bins: int = 5
    hic_tad_boundary_k: float = 1.0
    allow_partial_run: bool = False


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
