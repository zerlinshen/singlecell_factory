from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import (
    BatchConfig,
    CbioPortalConfig,
    CellRangerConfig,
    ClusteringConfig,
    CNVConfig,
    DEConfig,
    DoubletConfig,
    GeneSignatureConfig,
    PaperReproConfig,
    PipelineConfig,
    PseudobulkConfig,
    PseudobulkDEConfig,
    QCConfig,
    VelocityConfig,
    scale_mode_to_capabilities,
)
from .module_catalog import DEFAULT_OPTIONAL_MODULES as DEFAULT_OPTIONAL_MODULE_NAMES
from .module_catalog import module_help_list
from .pipeline import MODULE_DEPENDENCIES, run_pipeline


DEFAULT_OPTIONAL_MODULES = ",".join(DEFAULT_OPTIONAL_MODULE_NAMES)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments for modular workflow."""
    parser = argparse.ArgumentParser(description="Modular single-cell workflow runner")
    parser.add_argument("--project", required=True)
    parser.add_argument("--sample-root", required=True, help="Dataset root, e.g. data/raw/lung_carcinoma_3k_count")
    parser.add_argument(
        "--outs-dir",
        default="",
        help="Path to Cell Ranger filtered_feature_bc_matrix. Defaults to <sample-root>/outs/filtered_feature_bc_matrix",
    )
    parser.add_argument("--output-dir", default="/home/zerlinshen/singlecell_factory/results")
    parser.add_argument(
        "--optional-modules",
        default=DEFAULT_OPTIONAL_MODULES,
        help=f"Comma-separated optional modules. Available: {module_help_list()}",
    )
    parser.add_argument("--markers-json", default="", help="Optional custom marker dictionary JSON file")

    # Cell Ranger
    parser.add_argument("--fastq-dir", default="")
    parser.add_argument(
        "--transcriptome-dir",
        default="",
        help=(
            "Cell Ranger reference directory. Used by cellranger count and for RNA velocity "
            "GTF auto-discovery (<ref>/genes/genes.gtf[.gz])."
        ),
    )
    parser.add_argument("--sample-id", default="lusc")
    parser.add_argument("--localcores", type=int, default=8)
    parser.add_argument("--localmem", type=int, default=64)
    parser.add_argument("--force-cellranger", action="store_true")
    parser.add_argument("--no-run-cellranger-if-missing", action="store_true")

    # QC
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-genes", type=int, default=7000)
    parser.add_argument("--min-counts", type=int, default=500)
    parser.add_argument("--max-counts", type=int, default=50000)
    parser.add_argument("--max-mito-pct", type=float, default=20.0)
    parser.add_argument("--max-ribo-pct", type=float, default=50.0)
    parser.add_argument("--min-cells", type=int, default=3)

    # Doublet detection
    parser.add_argument("--expected-doublet-rate", type=float, default=0.06)
    parser.add_argument("--no-remove-doublets", action="store_true", help="Keep doublets (mark but don't remove)")

    # Clustering
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=40)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--leiden-resolution", type=float, default=0.8)
    parser.add_argument("--scale-data", action="store_true", help="Apply sc.pp.scale() before PCA")
    parser.add_argument(
        "--gpu-mode",
        default="auto",
        choices=["auto", "off", "force"],
        help="GPU policy: auto (try GPU, fall back to CPU), off (CPU only), force (require GPU)",
    )
    parser.add_argument(
        "--scale-mode",
        default="standard",
        choices=["standard", "large", "massive"],
        help=(
            "Dataset-size execution profile: standard (default), large (safer defaults for ~100k+ cells), "
            "massive (minimal-memory first-pass for several-hundred-thousand to million-cell runs). "
            "Acts as a preset bundle that expands to --lazy-read / --doublet-strategy / "
            "--clustering-engine / --checkpoint-policy. Explicit capability flags override the preset."
        ),
    )
    parser.add_argument(
        "--lazy-read",
        default="",
        choices=["", "auto", "true", "false"],
        help="Control lazy zarr loading: auto (>5 GB triggers lazy), true, false. Overrides scale-mode preset.",
    )
    parser.add_argument(
        "--doublet-strategy",
        default="",
        choices=["", "auto", "grouped", "whole", "skip"],
        help=(
            "Doublet detection strategy: auto (grouped when n_obs>=100k and sample column present), "
            "grouped, whole, skip. Overrides scale-mode preset."
        ),
    )
    parser.add_argument(
        "--clustering-engine",
        default="",
        choices=["", "auto", "sparse_exact", "css", "gpu"],
        help="Clustering engine: auto, sparse_exact, css, gpu. Overrides scale-mode preset.",
    )
    parser.add_argument(
        "--checkpoint-policy",
        default="",
        choices=["", "full", "mandatory_only", "metadata_only"],
        help="Checkpoint policy: full (save all), mandatory_only (skip early modules), metadata_only. Overrides scale-mode preset.",
    )

    # Cohort subset
    parser.add_argument(
        "--cohort-subset",
        action="append",
        default=None,
        dest="cohort_subset",
        metavar="OBS_COL=VAL1,VAL2",
        help=(
            "Filter cells to a subset: <obs_col>=<value1>,<value2>. "
            "Multiple --cohort-subset flags are AND-combined. "
            "E.g. --cohort-subset disease=lung_adenocarcinoma,lung_squamous_cell_carcinoma"
        ),
    )

    # Differential expression
    parser.add_argument(
        "--de-method",
        default="wilcoxon",
        choices=["wilcoxon", "t-test", "t-test_overestim_var", "logreg"],
        help="Method for scanpy.tl.rank_genes_groups (default: wilcoxon)",
    )
    parser.add_argument(
        "--de-n-genes",
        type=int,
        default=300,
        help="Maximum number of genes ranked per cluster in DE (default: 300)",
    )
    parser.add_argument(
        "--de-correction",
        default="benjamini-hochberg",
        choices=["benjamini-hochberg", "bonferroni"],
        help="Multiple-testing correction for marker DE (default: benjamini-hochberg; paper uses bonferroni)",
    )
    parser.add_argument(
        "--de-min-pct",
        type=float,
        default=0.10,
        help="Min fraction of cells expressing a gene for DE marker (default: 0.10; paper: 0.30)",
    )

    # Cell cycle
    parser.add_argument("--regress-cell-cycle", action="store_true", help="Regress out cell cycle effects")

    # Batch correction
    parser.add_argument("--batch-key", default="sample", help="Column in obs for batch labels")
    parser.add_argument(
        "--batch-method",
        default="harmony",
        choices=["harmony", "bbknn", "combat", "scanorama", "scvi", "mnn", "fastmnn"],
    )
    parser.add_argument(
        "--scvi-max-epochs",
        type=int,
        default=200,
        help="Max training epochs for scVI when --batch-method scvi (default: 200)",
    )
    parser.add_argument(
        "--scvi-n-latent",
        type=int,
        default=30,
        help="Latent dimension for scVI when --batch-method scvi (default: 30)",
    )
    parser.add_argument(
        "--no-scvi-early-stopping",
        action="store_true",
        help="Disable early stopping for scVI when --batch-method scvi.",
    )

    # Trajectory
    parser.add_argument("--trajectory-root-cluster", default=None, help="Leiden cluster ID for DPT root")

    # CNV inference
    parser.add_argument("--cnv-reference-group", default=None, help="Cell type to use as normal reference for CNV")
    parser.add_argument("--cnv-window-size", type=int, default=100)

    # RNA velocity
    parser.add_argument("--velocity-loom", default="", help="Path to loom file with spliced/unspliced counts")
    parser.add_argument("--velocity-bam", default="", help="Path to possorted_genome_bam.bam for spliced/unspliced extraction")
    parser.add_argument(
        "--velocity-gtf",
        default="",
        help=(
            "Path to genes.gtf(.gz) for spliced/unspliced extraction. "
            "Optional if --transcriptome-dir is set (auto-discovery enabled)."
        ),
    )
    parser.add_argument("--velocity-mode", default="stochastic", choices=["stochastic", "dynamical"])
    parser.add_argument("--velocity-n-jobs", type=int, default=4, help="Parallel workers for BAM extraction and scVelo dynamics (default: 4)")
    parser.add_argument("--velocity-min-shared-counts", type=int, default=20, help="Minimum shared counts for scVelo gene filtering (default: 20)")
    parser.add_argument("--velocity-n-pcs", type=int, default=30, help="PCA components for scVelo moments (default: 30)")
    parser.add_argument("--velocity-n-neighbors", type=int, default=30, help="Neighbors for scVelo moments (default: 30)")

    # Gene signature scoring
    parser.add_argument(
        "--signature-json", default="",
        help="JSON file with custom gene signatures: {name: [gene1, gene2, ...]}",
    )
    parser.add_argument(
        "--paper-spec-json",
        default="",
        help=(
            "JSON spec for paper-driven reproduction tracking "
            "(paper metadata + figure parity checks)."
        ),
    )
    parser.add_argument(
        "--paper-repro-strict",
        action="store_true",
        help="Fail the run if paper_repro reports missing provenance or failed figure checks.",
    )

    # cBioPortal validation
    parser.add_argument(
        "--cbioportal-genes",
        default="",
        help="Comma-separated gene symbols for cBioPortal validation, e.g. ELF3,TP53",
    )
    parser.add_argument(
        "--cbioportal-study",
        default="lusc_tcga_pan_can_atlas_2018",
        help="cBioPortal study ID (default: lusc_tcga_pan_can_atlas_2018)",
    )
    parser.add_argument(
        "--no-cbioportal-de-genes",
        action="store_true",
        help="Do not automatically pull top DE genes into cBioPortal validation",
    )
    parser.add_argument(
        "--cbioportal-top-n",
        type=int,
        default=20,
        help="Number of top DE genes to include in cBioPortal validation (default: 20)",
    )

    # Pseudobulk DE
    parser.add_argument(
        "--pseudobulk-sample-col",
        default="",
        help="obs column identifying biological samples for pseudobulk DE. Defaults to --batch-key, then sample/batch/donor/patient.",
    )
    parser.add_argument(
        "--pseudobulk-group-col",
        default="cell_type",
        help="obs column used for stratified pseudobulk contrasts or exploratory group-vs-rest output (default: cell_type).",
    )
    parser.add_argument(
        "--pseudobulk-contrast-col",
        default="",
        help="obs column containing the confirmatory pseudobulk condition labels, e.g. condition or disease.",
    )
    parser.add_argument(
        "--pseudobulk-contrast-a",
        default="",
        help="First condition label for confirmatory pseudobulk DE.",
    )
    parser.add_argument(
        "--pseudobulk-contrast-b",
        default="",
        help="Second condition label for confirmatory pseudobulk DE.",
    )
    parser.add_argument(
        "--pseudobulk-contrast-json",
        default="",
        help="Optional JSON list of pseudobulk contrasts with contrast_col/contrast_a/contrast_b/name fields.",
    )
    parser.add_argument(
        "--pseudobulk-exploratory-group-vs-rest",
        action="store_true",
        help="Allow exploratory group-vs-rest pseudobulk when no explicit contrast is provided.",
    )
    parser.add_argument("--pseudobulk-min-cells-per-sample", type=int, default=3)
    parser.add_argument("--pseudobulk-min-samples-per-condition", type=int, default=2)

    # Configurable thresholds
    parser.add_argument("--de-pval-threshold", type=float, default=0.05,
                        help="Adjusted p-value threshold for significant DE genes (default: 0.05)")
    parser.add_argument("--de-logfc-threshold", type=float, default=0.25,
                        help="Minimum absolute log fold change for DE genes (default: 0.25)")
    parser.add_argument("--annotation-confidence-threshold", type=float, default=0.1,
                        help="Minimum score for confident cell type assignment (default: 0.1)")
    parser.add_argument(
        "--annotation-strategy",
        default="cluster_voting",
        choices=["cluster_voting", "cell_argmax"],
        help="Annotation strategy: cluster_voting (aggregate per cluster then score) or cell_argmax (per-cell score then majority vote). Default: cluster_voting",
    )
    parser.add_argument(
        "--reference-adata",
        default="",
        help=(
            "Optional reference h5ad for annotation label transfer. "
            "When provided, annotation will run marker scoring + KNN reference mapping."
        ),
    )
    parser.add_argument(
        "--reference-label-key",
        default="cell_type",
        help="Column in reference obs used for label transfer (default: cell_type)",
    )
    parser.add_argument(
        "--reference-k",
        type=int,
        default=15,
        help="K neighbors for reference label transfer (default: 15)",
    )
    parser.add_argument(
        "--reference-min-confidence",
        type=float,
        default=0.6,
        help="Min reference confidence to override marker label (default: 0.6)",
    )
    parser.add_argument(
        "--reference-override-mode",
        default="conservative",
        choices=["conservative", "all"],
        help=(
            "How reference mapping overrides marker labels: "
            "'conservative' (Unknown/low-confidence only) or 'all'."
        ),
    )

    # Checkpointing & resume
    parser.add_argument(
        "--checkpoint", action="store_true",
        help="Save checkpoints after each module for crash recovery",
    )
    parser.add_argument(
        "--resume-from",
        default=None,
        help="Resume pipeline from this module using saved checkpoints",
    )

    # Parallel execution
    parser.add_argument(
        "--parallel-workers",
        type=int,
        default=1,
        help="Number of parallel workers for independent modules (default: 1 = sequential)",
    )
    return parser.parse_args()


def _validate_modules(optional_modules_str: str) -> list[str]:
    """Parse and validate optional module names."""
    modules = [m.strip() for m in optional_modules_str.split(",") if m.strip()]
    unknown = set(modules) - MODULE_DEPENDENCIES.keys()
    if unknown:
        raise SystemExit(
            f"Error: unknown module(s): {sorted(unknown)}. "
            f"Available: {sorted(MODULE_DEPENDENCIES.keys())}"
        )
    return modules


def _load_markers(markers_json: str) -> dict[str, list[str]]:
    if not markers_json:
        return {}
    payload = json.loads(Path(markers_json).read_text(encoding="utf-8"))
    return {str(k): [str(g) for g in v] for k, v in payload.items()}


def _apply_scale_mode(args: argparse.Namespace) -> argparse.Namespace:
    """Expand scale_mode preset into capability flags, then apply numeric tuning.

    Explicit capability flags (non-empty) always win over the preset bundle.
    This preserves the contract that --scale-mode massive produces identical
    behavior to the NC2024 launch script while allowing per-flag overrides.
    """
    # Expand the preset bundle first; explicit flags override below.
    preset = scale_mode_to_capabilities(args.scale_mode)
    if not args.lazy_read:
        args.lazy_read = preset["lazy_read"]
    if not args.doublet_strategy:
        args.doublet_strategy = preset["doublet_strategy"]
    if not args.clustering_engine:
        args.clustering_engine = preset["clustering_engine"]
    if not args.checkpoint_policy:
        args.checkpoint_policy = preset["checkpoint_policy"]

    if args.scale_mode == "standard":
        return args

    if args.scale_mode == "large":
        if args.optional_modules == DEFAULT_OPTIONAL_MODULES:
            args.optional_modules = "clustering,annotation,differential_expression"
        if args.n_top_genes == 3000:
            args.n_top_genes = 2000
        if args.n_pcs == 40:
            args.n_pcs = 30
        if args.n_neighbors == 15:
            args.n_neighbors = 12
        if args.leiden_resolution == 0.8:
            args.leiden_resolution = 0.6
        if args.de_n_genes == 300:
            args.de_n_genes = 200
        return args

    # massive
    if args.optional_modules == DEFAULT_OPTIONAL_MODULES:
        args.optional_modules = "clustering"
    if args.n_top_genes == 3000:
        args.n_top_genes = 1000
    if args.n_pcs == 40:
        args.n_pcs = 20
    if args.n_neighbors == 15:
        args.n_neighbors = 10
    if args.leiden_resolution == 0.8:
        args.leiden_resolution = 0.4
    if args.de_n_genes == 300:
        args.de_n_genes = 100
    return args


def _validate_args(args: argparse.Namespace) -> None:
    """Validate CLI arguments before pipeline execution."""
    if args.min_genes >= args.max_genes:
        raise SystemExit(f"Error: --min-genes ({args.min_genes}) must be < --max-genes ({args.max_genes})")
    if args.min_counts >= args.max_counts:
        raise SystemExit(f"Error: --min-counts ({args.min_counts}) must be < --max-counts ({args.max_counts})")
    if not (0 <= args.max_mito_pct <= 100):
        raise SystemExit(f"Error: --max-mito-pct must be in [0, 100], got {args.max_mito_pct}")
    if not (0 <= args.max_ribo_pct <= 100):
        raise SystemExit(f"Error: --max-ribo-pct must be in [0, 100], got {args.max_ribo_pct}")
    if args.parallel_workers < 1:
        raise SystemExit(f"Error: --parallel-workers must be >= 1, got {args.parallel_workers}")
    if args.n_pcs < 1:
        raise SystemExit(f"Error: --n-pcs must be >= 1, got {args.n_pcs}")
    if args.n_neighbors < 1:
        raise SystemExit(f"Error: --n-neighbors must be >= 1, got {args.n_neighbors}")
    if args.leiden_resolution <= 0:
        raise SystemExit(f"Error: --leiden-resolution must be > 0, got {args.leiden_resolution}")
    if args.de_n_genes < 1:
        raise SystemExit(f"Error: --de-n-genes must be >= 1, got {args.de_n_genes}")
    if args.scvi_max_epochs < 1:
        raise SystemExit(f"Error: --scvi-max-epochs must be >= 1, got {args.scvi_max_epochs}")
    if args.scvi_n_latent < 1:
        raise SystemExit(f"Error: --scvi-n-latent must be >= 1, got {args.scvi_n_latent}")
    if args.reference_k < 1:
        raise SystemExit(f"Error: --reference-k must be >= 1, got {args.reference_k}")
    if not (0 <= args.reference_min_confidence <= 1):
        raise SystemExit(
            "Error: --reference-min-confidence must be in [0, 1], "
            f"got {args.reference_min_confidence}"
        )


def main() -> None:
    """CLI entrypoint."""
    import faulthandler
    import logging
    from ._run_ledger import RunLedger
    from ._shutdown import run_shutdown_cleanup

    faulthandler.enable()

    args = _apply_scale_mode(parse_args())
    _validate_args(args)
    sample_root = Path(args.sample_root)
    outs_dir = (
        Path(args.outs_dir)
        if args.outs_dir
        else sample_root / "outs" / "filtered_feature_bc_matrix"
    )

    cfg = PipelineConfig(
        project=args.project,
        output_dir=Path(args.output_dir),
        cellranger=CellRangerConfig(
            sample_root=sample_root,
            outs_dir=outs_dir,
            fastq_dir=Path(args.fastq_dir) if args.fastq_dir else None,
            transcriptome_dir=Path(args.transcriptome_dir) if args.transcriptome_dir else None,
            sample_id=args.sample_id,
            localcores=args.localcores,
            localmem=args.localmem,
            force_run=args.force_cellranger,
            run_if_missing=not args.no_run_cellranger_if_missing,
        ),
        qc=QCConfig(
            min_genes=args.min_genes,
            max_genes=args.max_genes,
            min_counts=args.min_counts,
            max_counts=args.max_counts,
            max_mito_pct=args.max_mito_pct,
            max_ribo_pct=args.max_ribo_pct,
            min_cells=args.min_cells,
        ),
        doublet=DoubletConfig(
            expected_doublet_rate=args.expected_doublet_rate,
            remove_doublets=not args.no_remove_doublets,
        ),
        clustering=ClusteringConfig(
            n_top_genes=args.n_top_genes,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            leiden_resolution=args.leiden_resolution,
            scale_data=args.scale_data,
        ),
        batch=BatchConfig(
            batch_key=args.batch_key,
            method=args.batch_method,
            scvi_max_epochs=args.scvi_max_epochs,
            scvi_n_latent=args.scvi_n_latent,
            scvi_early_stopping=not args.no_scvi_early_stopping,
        ),
        cnv=CNVConfig(
            reference_group=args.cnv_reference_group,
            window_size=args.cnv_window_size,
        ),
        velocity=VelocityConfig(
            loom_path=Path(args.velocity_loom) if args.velocity_loom else None,
            bam_path=Path(args.velocity_bam) if args.velocity_bam else None,
            gtf_path=Path(args.velocity_gtf) if args.velocity_gtf else None,
            mode=args.velocity_mode,
            n_jobs=args.velocity_n_jobs,
            min_shared_counts=args.velocity_min_shared_counts,
            n_pcs=args.velocity_n_pcs,
            n_neighbors=args.velocity_n_neighbors,
        ),
        optional_modules=_validate_modules(args.optional_modules),
        markers=_load_markers(args.markers_json),
        gene_signature=GeneSignatureConfig(
            signature_json=Path(args.signature_json) if args.signature_json else None,
        ),
        paper_repro=PaperReproConfig(
            spec_json=Path(args.paper_spec_json) if args.paper_spec_json else None,
            strict=args.paper_repro_strict,
        ),
        pseudobulk=PseudobulkConfig(
            sample_col=args.pseudobulk_sample_col or None,
            group_col=args.pseudobulk_group_col,
            contrast_col=args.pseudobulk_contrast_col or None,
            contrast_a=args.pseudobulk_contrast_a or None,
            contrast_b=args.pseudobulk_contrast_b or None,
            contrast_json=Path(args.pseudobulk_contrast_json) if args.pseudobulk_contrast_json else None,
            exploratory_group_vs_rest=args.pseudobulk_exploratory_group_vs_rest,
            min_cells_per_sample=args.pseudobulk_min_cells_per_sample,
            min_samples_per_condition=args.pseudobulk_min_samples_per_condition,
        ),
        cbioportal=CbioPortalConfig(
            genes=[g.strip() for g in args.cbioportal_genes.split(",") if g.strip()],
            study_id=args.cbioportal_study,
            use_de_genes=not args.no_cbioportal_de_genes,
            top_n_de_genes=args.cbioportal_top_n,
        ),
        de_config=DEConfig(
            marker_min_pct=args.de_min_pct,
            logfc_threshold=args.de_logfc_threshold,
            marker_correction=args.de_correction,
        ),
        pseudobulk_de=PseudobulkDEConfig(),
        regress_cell_cycle=args.regress_cell_cycle,
        trajectory_root_cluster=args.trajectory_root_cluster,
        checkpoint=args.checkpoint,
        resume_from=args.resume_from,
        parallel_workers=args.parallel_workers,
        de_method=args.de_method,
        de_n_genes=args.de_n_genes,
        de_pval_threshold=args.de_pval_threshold,
        de_logfc_threshold=args.de_logfc_threshold,
        de_correction=args.de_correction,
        de_min_pct=args.de_min_pct,
        annotation_confidence_threshold=args.annotation_confidence_threshold,
        reference_adata=Path(args.reference_adata) if args.reference_adata else None,
        reference_label_key=args.reference_label_key,
        reference_k=args.reference_k,
        reference_min_confidence=args.reference_min_confidence,
        reference_override_mode=args.reference_override_mode,
        gpu_mode=args.gpu_mode,
        scale_mode=args.scale_mode,
        lazy_read=args.lazy_read,
        doublet_strategy=args.doublet_strategy,
        clustering_engine=args.clustering_engine,
        checkpoint_policy=args.checkpoint_policy,
        cohort_subset=args.cohort_subset,
        annotation_strategy=args.annotation_strategy,
    )
    ledger = None
    try:
        _ledger_ctx = type("_LedgerCtx", (), {"cfg": cfg})()
        ledger = RunLedger(_ledger_ctx, args.project, Path(args.output_dir))
        ledger.record_start()
    except Exception as _ledger_exc:
        logging.getLogger(__name__).warning("RunLedger.record_start failed: %s", _ledger_exc)
        ledger = None

    import os as _os
    _watchdog_thread = None
    if _os.environ.get("SC_MEM_WATCHDOG", "").lower() == "on":
        try:
            from ._mem_watchdog import start as _watchdog_start
            _watchdog_ctx = type("_WatchdogCtx", (), {"metadata": {}})()
            _watchdog_thread = _watchdog_start(_watchdog_ctx)
        except Exception as _wd_exc:
            logging.getLogger(__name__).warning("MemoryWatchdog start failed: %s", _wd_exc)

    try:
        manifest = run_pipeline(cfg, ledger=ledger)
        print(manifest)
    finally:
        run_shutdown_cleanup()


if __name__ == "__main__":
    main()
