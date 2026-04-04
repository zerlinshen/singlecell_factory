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
    DoubletConfig,
    GeneSignatureConfig,
    PipelineConfig,
    QCConfig,
    VelocityConfig,
)
from .pipeline import run_pipeline


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
        default="clustering,differential_expression,annotation,trajectory,pseudo_velocity",
        help=(
            "Comma-separated optional modules. Available: clustering, cell_cycle, "
            "batch_correction, differential_expression, annotation, trajectory, "
            "pseudo_velocity, rna_velocity, cnv_inference, pathway_analysis, "
            "cell_communication, gene_regulatory_network, validate_cbioportal, "
            "immune_phenotyping, tumor_microenvironment, gene_signature_scoring"
        ),
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

    # Cell cycle
    parser.add_argument("--regress-cell-cycle", action="store_true", help="Regress out cell cycle effects")

    # Batch correction
    parser.add_argument("--batch-key", default="sample", help="Column in obs for batch labels")
    parser.add_argument("--batch-method", default="harmony", choices=["harmony", "bbknn", "combat", "scanorama"])

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


def _load_markers(markers_json: str) -> dict[str, list[str]]:
    if not markers_json:
        return {}
    payload = json.loads(Path(markers_json).read_text(encoding="utf-8"))
    return {str(k): [str(g) for g in v] for k, v in payload.items()}


def main() -> None:
    """CLI entrypoint."""
    args = parse_args()
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
        optional_modules=[m.strip() for m in args.optional_modules.split(",") if m.strip()],
        markers=_load_markers(args.markers_json),
        gene_signature=GeneSignatureConfig(
            signature_json=Path(args.signature_json) if args.signature_json else None,
        ),
        cbioportal=CbioPortalConfig(
            genes=[g.strip() for g in args.cbioportal_genes.split(",") if g.strip()],
            study_id=args.cbioportal_study,
            use_de_genes=not args.no_cbioportal_de_genes,
            top_n_de_genes=args.cbioportal_top_n,
        ),
        regress_cell_cycle=args.regress_cell_cycle,
        trajectory_root_cluster=args.trajectory_root_cluster,
        checkpoint=args.checkpoint,
        resume_from=args.resume_from,
        parallel_workers=args.parallel_workers,
    )
    manifest = run_pipeline(cfg)
    print(manifest)


if __name__ == "__main__":
    main()
