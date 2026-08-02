from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import (
    AmbientCorrectionConfig,
    BatchConfig,
    CbioPortalConfig,
    CellRangerConfig,
    ClusteringConfig,
    CompositionConfig,
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
from .legacy_output import (
    LEGACY_DEFAULT_OUTPUT_DIR,
    legacy_output_warning_message,
    record_legacy_output_access,
)
from .pipeline import MODULE_DEPENDENCIES, run_pipeline


DEFAULT_OPTIONAL_MODULES = ",".join(DEFAULT_OPTIONAL_MODULE_NAMES)


# ---------------------------------------------------------------------------
# Canonical scientific profile — the single source of truth for every surface
# that can start an analysis.
#
# Defined ahead of parse_args() on purpose: the argparse defaults below read
# from this dict rather than repeating the literals, and the config.py dataclass
# defaults are pinned to it by tests/test_paper_param_alignment.py. Before
# 2026-07-28 all three were independent copies and two of them disagreed
# (dataclass n_pcs=15/resolution=1.0 vs CLI 40/0.8), so identical science
# launched programmatically and via the CLI silently diverged.
#
# Departures from canonical are only sanctioned through a NAMED profile in
# _SCIENTIFIC_PROFILE_OVERRIDES, and cost --acknowledge-scientific-non-equivalence.
# ---------------------------------------------------------------------------
_CANONICAL_SCIENTIFIC_PARAMETERS = {
    "optional_modules": DEFAULT_OPTIONAL_MODULES,
    "n_top_genes": 3000,
    # 40 PCs, not the NC2024 paper's 15. Heterogeneous multi-batch tumour tissue
    # keeps biologically meaningful variance well past the first ~15 components,
    # and rare populations are exactly what a low truncation discards; the
    # benchmark guidance is to err high because over-inclusion is much cheaper
    # than truncation (Luecken & Theis 2019, Mol Syst Biol 15:e8746; Heumos 2023,
    # Nat Rev Genet 24:550-572).
    "n_pcs": 40,
    "n_neighbors": 15,
    # 0.8 is the long-standing canonical CLI value and was NOT retuned here; the
    # 2026-07-28 change aligned the dataclass to it. Sweep per cohort with
    # --leiden-resolution-sweep before overriding.
    "leiden_resolution": 0.8,
    "de_n_genes": 300,
    "doublet_strategy": "auto",
    "clustering_engine": "auto",
}

# Named, explicitly non-equivalent departures from canonical. Keys double as the
# --scientific-profile choices, so a profile can never be advertised without an
# implementation (or implemented without being reachable).
_SCIENTIFIC_PROFILE_OVERRIDES = {
    "canonical": {},
    "legacy-large": {
        "optional_modules": "clustering,annotation,differential_expression",
        "n_top_genes": 2000,
        "n_pcs": 30,
        "n_neighbors": 12,
        "leiden_resolution": 0.6,
        "de_n_genes": 200,
    },
    "legacy-massive": {
        "optional_modules": "clustering",
        "n_top_genes": 1000,
        "n_pcs": 20,
        "n_neighbors": 10,
        "leiden_resolution": 0.4,
        "de_n_genes": 100,
        "doublet_strategy": "grouped",
    },
    # NC2024 reproduction: Leiden at resolution 1.0 on a 15-PC Harmony embedding
    # (docs/PUBLICATION_READY.md, AUDIT_2026-04-26_v2.md). These were the
    # ClusteringConfig dataclass defaults until 2026-07-28, where they shadowed
    # every programmatic caller instead of being opted into. Kept as a named
    # profile so the reproduction stays selectable, auditable in the resolved
    # diff, and self-documenting. The paper's DE settings (bonferroni,
    # min_pct=0.30) stay launcher-level and are deliberately NOT bundled here —
    # this profile changes clustering geometry only.
    "paper-15pc": {
        "n_pcs": 15,
        "leiden_resolution": 1.0,
    },
}


def parse_args() -> argparse.Namespace:
    """Parse command line arguments for modular workflow."""
    parser = argparse.ArgumentParser(description="Modular single-cell workflow runner")
    parser.add_argument("--project", required=True)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--sample-root",
        help="Dataset root, e.g. data/raw/lung_carcinoma_3k_count",
    )
    input_group.add_argument(
        "--input-h5ad",
        help=(
            "Exact prepared AnnData file to ingest. The filename is arbitrary; "
            "the canonical loader reads this path directly instead of searching "
            "for <sample-root>/prepared_input.h5ad."
        ),
    )
    parser.add_argument(
        "--outs-dir",
        default="",
        help="Path to Cell Ranger filtered_feature_bc_matrix. Defaults to <sample-root>/outs/filtered_feature_bc_matrix",
    )
    parser.add_argument("--output-dir", default=str(LEGACY_DEFAULT_OUTPUT_DIR))
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
    parser.add_argument(
        "--doublet-backend",
        default="scrublet",
        choices=["scrublet", "doubletfinder", "scdblfinder", "consensus"],
        help="Doublet backend. Default stays scrublet; use consensus/second-opinion lanes for benchmarking.",
    )
    parser.add_argument(
        "--doublet-consensus-logic",
        default="or",
        choices=["or", "and", "rank"],
        help="Consensus merge logic when --doublet-backend consensus is used.",
    )
    parser.add_argument(
        "--doublet-consensus-pair",
        default="scrublet_doubletfinder",
        choices=[
            "scrublet_doubletfinder",
            "scrublet_scdblfinder",
            "scrublet_doubletfinder_scdblfinder",
        ],
        help="Backends to combine for consensus doublet calling.",
    )
    parser.add_argument("--doubletfinder-pn", type=float, default=0.25)
    parser.add_argument("--doubletfinder-pk", type=float, default=0.09)
    parser.add_argument("--doubletfinder-pcs", type=int, default=20)
    parser.add_argument(
        "--scdblfinder-samples-col",
        default="",
        help="Optional obs column passed to scDblFinder as sample labels. Empty = no sample grouping.",
    )
    parser.add_argument("--doublet-subprocess-timeout", type=int, default=1800)

    # Ambient RNA correction (DecontX, conditional per-sample triggers)
    # Policy: ops/policy/ambient_correction_policy.md (Phase 2 decided 2026-05-20)
    parser.add_argument(
        "--ambient-disable-triggers",
        action="store_true",
        help=(
            "Skip ambient correction entirely (paper-faithful mode where the source "
            "paper did not do ambient RNA correction). Equivalent to SC_AMBIENT_TRIGGERS_DISABLE=1."
        ),
    )
    parser.add_argument(
        "--ambient-dry-run-triggers-only",
        action="store_true",
        help="Evaluate ambient triggers and record decision but skip the R DecontX call (introspection mode).",
    )
    parser.add_argument(
        "--ambient-trigger-top50",
        type=float,
        default=50.0,
        help="T1 trigger: median pct_counts_in_top_50 (scanpy percent 0-100) > this fires DecontX (default 50.0).",
    )
    parser.add_argument(
        "--ambient-trigger-mt",
        type=float,
        default=15.0,
        help="T2 trigger: median pct_counts_mt (scanpy percent 0-100) > this fires DecontX (default 15.0).",
    )
    # T3 (doublet-excess) removed 2026-05-20: DAG ordering puts ambient_correction
    # before doublet_detection, so predicted_doublet is never present when triggers
    # are evaluated. Doublet-driven ambient re-trigger is deferred to a potential
    # future Phase 3. See ops/policy/ambient_correction_policy.md Phase 2 notes.
    parser.add_argument(
        "--ambient-trigger-count-corr",
        type=float,
        default=0.85,
        help="T4 trigger: Spearman(total_counts, n_genes_by_counts) < this fires DecontX (default 0.85).",
    )
    parser.add_argument(
        "--ambient-trigger-cohort-cv",
        type=float,
        default=0.50,
        help="T5 trigger: cohort housekeeping-gene CV > this fires DecontX on all samples (default 0.50).",
    )
    parser.add_argument(
        "--ambient-decontx-max-iter",
        type=int,
        default=200,
        help="DecontX max EM iterations (default 200; Yang 2020 Genome Biology).",
    )
    parser.add_argument(
        "--ambient-batch-column",
        default="",
        help="obs column to pass to DecontX as --batch (per-sample contamination modelling). Empty = no batch.",
    )

    # Clustering. Defaults are read from the canonical profile rather than
    # re-typed, so the flag and the profile cannot drift apart.
    parser.add_argument(
        "--n-top-genes", type=int, default=_CANONICAL_SCIENTIFIC_PARAMETERS["n_top_genes"]
    )
    parser.add_argument("--n-pcs", type=int, default=_CANONICAL_SCIENTIFIC_PARAMETERS["n_pcs"])
    parser.add_argument(
        "--n-neighbors", type=int, default=_CANONICAL_SCIENTIFIC_PARAMETERS["n_neighbors"]
    )
    parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=_CANONICAL_SCIENTIFIC_PARAMETERS["leiden_resolution"],
    )
    parser.add_argument(
        "--leiden-resolution-sweep",
        default="",
        help=(
            "Comma-separated diagnostic Leiden resolutions to audit, e.g. 0.5,0.8,1.0. "
            "Writes clustering/leiden_resolution_sweep.csv without changing --leiden-resolution."
        ),
    )
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
            "Resource-only dataset-size strategy. Expands only resource execution "
            "settings such as --lazy-read / --checkpoint-policy and never "
            "changes HVGs, PCs, neighbors, resolution, DE limits, or module selection."
        ),
    )
    parser.add_argument(
        "--scientific-profile",
        default="canonical",
        choices=list(_SCIENTIFIC_PROFILE_OVERRIDES),
        help=(
            "Scientific parameter profile, separate from --scale-mode. Non-canonical "
            "profiles change the analysis and require "
            "--acknowledge-scientific-non-equivalence. 'paper-15pc' reproduces the "
            "NC2024 clustering geometry (15-PC Harmony space, Leiden resolution 1.0)."
        ),
    )
    parser.add_argument(
        "--acknowledge-scientific-non-equivalence",
        action="store_true",
        help=(
            "Acknowledge that resolved changes to HVGs, PCs, neighbors, Leiden "
            "resolution, DE limits, or module selection are not scientifically equivalent."
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
            "grouped, whole, skip. Scientific choice independent of --scale-mode."
        ),
    )
    parser.add_argument(
        "--clustering-engine",
        default="",
        choices=["", "auto", "sparse_exact", "css", "gpu"],
        help=(
            "Clustering engine: auto, sparse_exact, css, gpu. Scientific choice "
            "independent of --scale-mode."
        ),
    )
    parser.add_argument(
        "--checkpoint-policy",
        default="",
        choices=["", "full"],
        help="Checkpoint policy: full (save all checkpoints). Overrides scale-mode preset.",
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

    parser.add_argument(
        "--multimodal-engine",
        default=None,
        choices=["off", "wnn", "mofa"],
        help="Engine for multimodal_integration module (default off). Set SC_MULTIMODAL_ENGINE to override; module must also be in --optional-modules.",
    )
    parser.add_argument(
        "--second-obsm-key",
        default=None,
        help=(
            "Second-modality obsm key for multimodal_integration (default: protein_clr). "
            "When set, passing the key is treated as an explicit opt-in: missing key raises "
            "ModuleContractError for wnn/mofa engines. Set SC_MULTIMODAL_SECOND_OBSM to override."
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
        default=_CANONICAL_SCIENTIFIC_PARAMETERS["de_n_genes"],
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
        "--harmony-backend",
        default="auto",
        choices=["auto", "cpu", "gpu", "direct"],
        help=(
            "Harmony backend selection (Plan F-Harmony, G-C0 v3). 'auto' uses "
            "the proven harmonypy direct path on this environment. 'direct' "
            "calls harmonypy.run_harmony. 'cpu' forces harmonypy via "
            "scanpy_external. 'gpu' forces rsc (raises if unavailable). "
            "Non-convergence raises unless SC_ALLOW_HARMONY_NON_CONVERGENCE=1 "
            "opt-in."
        ),
    )
    parser.add_argument(
        "--harmony-max-iter",
        type=int,
        default=50,
        help=(
            "Maximum Harmony iterations (Korsunsky 2019). Default bumped to "
            "50 in 2026-05-20 — empirical: 10 does not converge on 75-batch "
            "tumor cohorts."
        ),
    )
    parser.add_argument(
        "--harmony-theta",
        type=float,
        default=2.0,
        help="Harmony theta diversity parameter (Korsunsky 2019).",
    )
    parser.add_argument(
        "--harmony-sigma",
        type=float,
        default=0.1,
        help="Harmony sigma kernel width parameter (Korsunsky 2019).",
    )
    parser.add_argument(
        "--select-integration",
        action="store_true",
        help=(
            "Opt-in: run the per-run discovery integration-selection gate "
            "(integration_select module) after clustering and before "
            "batch_correction. It requires baseline/Harmony/scVI/shuffle "
            "candidates, fails loud on degraded candidates or a non-firing "
            "shuffle control, and SETS cfg.batch.method (routing a Harmony "
            "pick to the working harmonypy-direct backend). Runs an expensive "
            "scVI seed sweep (cache-bounded); NOT default-on. Adds "
            "integration_select to the requested modules."
        ),
    )
    parser.add_argument(
        "--integration-margin-mix",
        type=float,
        default=0.05,
        help=(
            "Margin the candidate batch-mixing must beat baseline by in the "
            "integration_select gate (default: 0.05)."
        ),
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
    parser.add_argument(
        "--trajectory-root-justification",
        default=None,
        help=(
            "Biological rationale for the DPT root cluster. Both this value and "
            "--trajectory-root-cluster are required for claim-capable pseudotime."
        ),
    )

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

    # Hi-C / single-cell 3D genome
    parser.add_argument(
        "--hic-contacts-path",
        default="",
        help=(
            "Optional .cool/.mcool or tab-separated Hi-C contact-pair file for "
            "hic_ingest. TSV columns: chrom1,pos1,chrom2,pos2,count."
        ),
    )
    parser.add_argument(
        "--hic-resolution-bp",
        type=int,
        default=25_000,
        help="Hi-C bin size in bp for TSV ingest or metadata reporting (default: 25000).",
    )
    parser.add_argument(
        "--hic-chromsizes-path",
        default="",
        help="Optional chrom.sizes file for TSV Hi-C ingest; inferred from contacts when omitted.",
    )
    parser.add_argument(
        "--hic-tad-window-bins",
        type=int,
        default=5,
        help="Window radius in bins for insulation-score TAD boundary calling (default: 5).",
    )
    parser.add_argument(
        "--hic-tad-boundary-k",
        type=float,
        default=1.0,
        help="Boundary threshold in SD units below mean insulation score (default: 1.0).",
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

    # Replicate-aware composition analysis
    parser.add_argument(
        "--composition-sample-col",
        default="",
        help=(
            "obs column identifying independent biological samples. When omitted, "
            "the composition module may resolve --batch-key or a standard sample "
            "column for descriptive proportions only."
        ),
    )
    parser.add_argument(
        "--composition-condition-col",
        default="",
        help=(
            "Explicit obs condition/covariate column for scCODA inference. Empty "
            "means descriptive sample-level proportions only."
        ),
    )
    parser.add_argument(
        "--composition-contrast-a",
        default="",
        help="First condition label for a two-level composition contrast.",
    )
    parser.add_argument(
        "--composition-contrast-b",
        default="",
        help="Second condition label for a two-level composition contrast.",
    )
    parser.add_argument(
        "--composition-covariates",
        default="",
        help="Comma-separated sample-level covariate columns added to the scCODA formula.",
    )
    parser.add_argument(
        "--composition-min-samples-per-condition",
        type=int,
        default=2,
        help="Minimum independent biological samples in every modeled condition (default: 2).",
    )

    # Pseudobulk DE
    parser.add_argument(
        "--pseudobulk-sample-col",
        default="",
        help=(
            "obs column identifying biological samples for pseudobulk DE. "
            "Required explicitly for confirmatory contrasts; exploratory mode "
            "may fall back to --batch-key, then sample/batch/donor/patient."
        ),
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

    # Marker intelligence (P1A)
    parser.add_argument(
        "--tissue",
        default="lung",
        help="Tissue type for marker DB routing (default: lung)",
    )
    parser.add_argument(
        "--condition",
        default="NSCLC",
        help="Disease/condition for marker DB routing (default: NSCLC)",
    )
    parser.add_argument(
        "--validate-context",
        action="store_true",
        help="Compute per-cluster context validation scores and write context_validation.json (AC-4)",
    )
    parser.add_argument(
        "--context-mismatch-threshold",
        type=float,
        default=0.3,
        help="Score threshold below which a cluster triggers a context mismatch warning (default: 0.3)",
    )
    parser.add_argument(
        "--context-min-cells",
        type=int,
        default=20,
        help="Minimum cells per cluster for context validation (default: 20)",
    )

    # Reproducibility
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Global random seed propagated to all stochastic modules (default: 42)",
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

    # Project-root layout (additive; legacy --output-dir still works when absent)
    parser.add_argument(
        "--project-root",
        default=None,
        metavar="PATH",
        help=(
            "Project directory root. When set, all artifacts are written to "
            "<project-root>/runs/<run-id>/python/. "
            "Preferred over legacy --output-dir."
        ),
    )
    parser.add_argument(
        "--run-id",
        default=None,
        metavar="STR",
        help=(
            "Run identifier (format: YYYY-MM-DDTHHMMZ-<7hex>). "
            "Auto-generated from UTC timestamp + factory SHA when absent."
        ),
    )
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help=(
            "Allow pipeline to run when factory git tree has uncommitted changes. "
            "The diff SHA256 is recorded in manifest.json for forensics."
        ),
    )
    parser.add_argument(
        "--allow-partial-run",
        action="store_true",
        help=(
            "Return exit status 0 even when a requested optional module failed. "
            "Failure remains recorded in both manifests; intended only for "
            "explicit recovery workflows."
        ),
    )
    return parser.parse_args()


def _load_pipeline_manifest(result: object) -> tuple[dict, str]:
    """Normalize the preserved run_pipeline Path contract for CLI consumption."""
    if isinstance(result, dict):
        return dict(result), ""
    if isinstance(result, (str, Path)):
        path = Path(result)
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError(f"pipeline manifest must contain a JSON object: {path}")
        return payload, str(path)
    raise TypeError(
        f"run_pipeline returned unsupported manifest type: {type(result).__name__}"
    )


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


def _resolve_optional_modules(args: argparse.Namespace) -> list[str]:
    """Validate the requested module list and fold in opt-in flag modules.

    ``--select-integration`` adds ``integration_select`` to the requested set
    (idempotent; never default-on). Dependency/ordering resolution is handled by
    the pipeline (integration_select depends_on clustering; batch_correction
    runs_after integration_select).
    """
    modules = _validate_modules(args.optional_modules)
    if getattr(args, "select_integration", False) and "integration_select" not in modules:
        modules.append("integration_select")
    contrast_json = getattr(args, "pseudobulk_contrast_json", "")
    contrast_a = getattr(args, "pseudobulk_contrast_a", "")
    contrast_b = getattr(args, "pseudobulk_contrast_b", "")
    contrast_requested = bool(contrast_json or contrast_a or contrast_b)
    if contrast_requested and not (contrast_json or (contrast_a and contrast_b)):
        raise SystemExit(
            "Error: a pseudobulk condition contrast requires both "
            "--pseudobulk-contrast-a and --pseudobulk-contrast-b, or a "
            "--pseudobulk-contrast-json contract."
        )
    if contrast_requested and "pseudobulk_de" not in modules:
        modules.append("pseudobulk_de")

    composition_condition = getattr(args, "composition_condition_col", "")
    composition_a = getattr(args, "composition_contrast_a", "")
    composition_b = getattr(args, "composition_contrast_b", "")
    if bool(composition_a) != bool(composition_b):
        raise SystemExit(
            "Error: a composition contrast requires both "
            "--composition-contrast-a and --composition-contrast-b."
        )
    if (composition_a or composition_b) and not composition_condition:
        raise SystemExit(
            "Error: a composition contrast requires --composition-condition-col; "
            "sample identifiers are not condition labels."
        )
    composition_requested = any((
        getattr(args, "composition_sample_col", ""),
        composition_condition,
        composition_a,
        composition_b,
        getattr(args, "composition_covariates", ""),
    ))
    if composition_requested and "composition" not in modules:
        modules.append("composition")
    return modules


def _load_markers(markers_json: str) -> dict[str, list[str]]:
    if not markers_json:
        return {}
    payload = json.loads(Path(markers_json).read_text(encoding="utf-8"))
    return {str(k): [str(g) for g in v] for k, v in payload.items()}


def _parse_resolution_sweep(value: str) -> tuple[float, ...]:
    if not value:
        return ()
    parsed: list[float] = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            resolution = float(token)
        except ValueError as exc:
            raise SystemExit(f"Error: invalid --leiden-resolution-sweep token: {token!r}") from exc
        if resolution <= 0:
            raise SystemExit(
                "Error: --leiden-resolution-sweep values must be > 0, "
                f"got {resolution}"
            )
        parsed.append(resolution)
    return tuple(dict.fromkeys(parsed))


def _apply_scale_mode(args: argparse.Namespace) -> argparse.Namespace:
    """Expand a resource-only scale preset into operational resource flags.

    Explicit resource flags (non-empty) always win over the preset bundle.
    Scientific parameters are intentionally untouched.
    """
    # F-4 (Plan ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md):
    # `--scale-mode massive` previously routed clustering_engine -> "css"
    # (a hard scientific compromise per Principle 2). The preset is now remapped
    # to clustering_engine -> "auto". The subsequent resource/science split
    # also moved grouped doublet selection behind the explicit
    # scientific_profile=legacy-massive acknowledgement path.
    if args.scale_mode == "massive":
        import warnings as _warnings
        _warnings.warn(
            "--scale-mode=massive: semantics changed as of 2026-05-19. The "
            "clustering_engine slot of this preset was remapped from 'css' "
            "(removed from production per Plan F-1, Principle 2) to 'auto'. "
            "Resource flags (lazy_read=true, checkpoint_policy=full) are preserved. "
            "Grouped doublet handling and CSS are scientific/method choices and are "
            "not selected by scale mode. Use an acknowledged scientific profile for "
            "non-equivalent method changes. See "
            "~/.omc/plans/nc-cell-clustering-final-strategy-plan.md",
            DeprecationWarning,
            stacklevel=2,
        )
    # Expand the preset bundle first; explicit flags override below.
    preset = scale_mode_to_capabilities(args.scale_mode)
    if not args.lazy_read:
        args.lazy_read = preset["lazy_read"]
    if not args.checkpoint_policy:
        args.checkpoint_policy = preset["checkpoint_policy"]

    return args


def _apply_scientific_profile(args: argparse.Namespace) -> argparse.Namespace:
    """Resolve scientific parameters and require explicit non-equivalence consent."""
    profile = getattr(args, "scientific_profile", "canonical")
    for field_name, canonical_value in _CANONICAL_SCIENTIFIC_PARAMETERS.items():
        if getattr(args, field_name) == "":
            setattr(args, field_name, canonical_value)
    overrides = _SCIENTIFIC_PROFILE_OVERRIDES[profile]
    for field_name, resolved_value in overrides.items():
        if getattr(args, field_name) == _CANONICAL_SCIENTIFIC_PARAMETERS[field_name]:
            setattr(args, field_name, resolved_value)

    resolved_diff = {}
    for field_name, canonical_value in _CANONICAL_SCIENTIFIC_PARAMETERS.items():
        resolved_value = getattr(args, field_name)
        if resolved_value != canonical_value:
            resolved_diff[field_name] = {
                "canonical": canonical_value,
                "resolved": resolved_value,
            }

    acknowledged = bool(
        getattr(args, "acknowledge_scientific_non_equivalence", False)
    )
    if resolved_diff and not acknowledged:
        changed = ", ".join(sorted(resolved_diff))
        raise SystemExit(
            "Error: resolved scientific settings differ from the canonical profile "
            f"({changed}). Re-run with --acknowledge-scientific-non-equivalence "
            "to record that the analyses are not scientifically equivalent."
        )

    args.resolved_scientific_parameter_diff = resolved_diff
    args.scientific_non_equivalence_acknowledged = acknowledged
    return args


def _validate_args(args: argparse.Namespace) -> None:
    """Validate CLI arguments before pipeline execution."""
    if args.input_h5ad:
        input_h5ad = Path(args.input_h5ad).expanduser()
        if input_h5ad.suffix.lower() != ".h5ad":
            raise SystemExit(
                f"Error: --input-h5ad must name a .h5ad file, got {input_h5ad}"
            )
        if not input_h5ad.is_file():
            raise SystemExit(f"Error: --input-h5ad file not found: {input_h5ad}")
        if args.outs_dir:
            raise SystemExit(
                "Error: --outs-dir cannot be combined with --input-h5ad; "
                "the exact AnnData file is the ingest source."
            )
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
    if args.hic_resolution_bp <= 0:
        raise SystemExit(f"Error: --hic-resolution-bp must be > 0, got {args.hic_resolution_bp}")
    if args.hic_tad_window_bins < 1:
        raise SystemExit(f"Error: --hic-tad-window-bins must be >= 1, got {args.hic_tad_window_bins}")
    if args.hic_tad_boundary_k < 0:
        raise SystemExit(f"Error: --hic-tad-boundary-k must be >= 0, got {args.hic_tad_boundary_k}")
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
    import os
    import socket
    import sys
    from ._run_ledger import RunLedger
    from ._shutdown import run_shutdown_cleanup
    from .project_paths import resolve_run_dir, python_dir
    from .manifest_writer import factory_git_state, write_manifest

    faulthandler.enable()

    args = _apply_scientific_profile(_apply_scale_mode(parse_args()))

    # Propagate --multimodal-engine to env var the module reads at runtime.
    if args.multimodal_engine is not None and "SC_MULTIMODAL_ENGINE" not in os.environ:
        os.environ["SC_MULTIMODAL_ENGINE"] = args.multimodal_engine
    # Propagate --second-obsm-key: env var wins if already set; CLI flag is explicit opt-in.
    if args.second_obsm_key is not None and "SC_MULTIMODAL_SECOND_OBSM" not in os.environ:
        os.environ["SC_MULTIMODAL_SECOND_OBSM"] = args.second_obsm_key
    _validate_args(args)

    # GOV-2 cutover semantics:
    # When `SC_REQUIRE_PROJECT_ROOT` is UNSET: a missing `--project-root` emits
    # a contracted warning, records external JSONL access telemetry, and falls
    # back to legacy `results/`. When
    # `SC_REQUIRE_PROJECT_ROOT=1`: a missing `--project-root` is a hard error
    # (exit code 2). Warning and hard-error are mutually exclusive (no
    # double-fire). Retirement is controlled by the explicit deprecation
    # contract and its measured access window.
    if args.project_root is None:
        if os.environ.get("SC_REQUIRE_PROJECT_ROOT") == "1":
            print(
                "ERROR: --project-root is required (SC_REQUIRE_PROJECT_ROOT=1 is set). "
                "Pass --project-root or unset the env var.",
                file=sys.stderr,
            )
            sys.exit(2)

    input_h5ad = (
        Path(args.input_h5ad).expanduser().resolve()
        if args.input_h5ad
        else None
    )
    sample_root = (
        Path(args.sample_root).expanduser()
        if args.sample_root
        else input_h5ad.parent
    )
    outs_dir = (
        Path(args.outs_dir)
        if args.outs_dir
        else sample_root / "outs" / "filtered_feature_bc_matrix"
    )

    # Resolve effective output directory from --project-root or legacy --output-dir.
    _FACTORY_ROOT = Path(__file__).resolve().parent.parent.parent
    _run_dir = None
    _run_id = None
    if args.project_root is not None:
        # Dirty-tree gate: abort unless --allow-dirty is set.
        _py_state = factory_git_state(_FACTORY_ROOT)
        if _py_state["dirty"] and not args.allow_dirty:
            print(
                "ERROR: singlecell_factory has uncommitted changes. "
                "Commit or stash them, or pass --allow-dirty to record the diff.",
                file=sys.stderr,
            )
            raise SystemExit(1)
        _run_dir = resolve_run_dir(
            Path(args.project_root),
            run_id=args.run_id,
            factory_sha=_py_state["sha"] or None,
        )
        _run_id = _run_dir.name
        effective_output_dir = python_dir(_run_dir)
    else:
        import warnings
        try:
            telemetry_path = record_legacy_output_access(
                output_dir=Path(args.output_dir),
                project=args.project,
            )
        except (OSError, ValueError) as exc:
            print(
                "ERROR: legacy output telemetry could not be recorded; the "
                "deprecation contract blocks this compatibility launch: "
                f"{exc}",
                file=sys.stderr,
            )
            raise SystemExit(2) from exc
        warning_message = legacy_output_warning_message(telemetry_path)
        warnings.warn(
            warning_message,
            DeprecationWarning,
            stacklevel=2,
        )
        print(f"DeprecationWarning: {warning_message}", file=sys.stderr)
        effective_output_dir = Path(args.output_dir)

    cfg = PipelineConfig(
        project=args.project,
        output_dir=effective_output_dir,
        cellranger=CellRangerConfig(
            sample_root=sample_root,
            outs_dir=outs_dir,
            input_h5ad=input_h5ad,
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
            backend=args.doublet_backend,
            consensus_logic=args.doublet_consensus_logic,
            consensus_pair=args.doublet_consensus_pair,
            subprocess_timeout=args.doublet_subprocess_timeout,
            doubletfinder_pn=args.doubletfinder_pn,
            doubletfinder_pk=args.doubletfinder_pk,
            doubletfinder_pcs=args.doubletfinder_pcs,
            scdblfinder_samples_col=args.scdblfinder_samples_col or None,
        ),
        ambient=AmbientCorrectionConfig(
            disable_triggers=args.ambient_disable_triggers,
            dry_run_triggers_only=args.ambient_dry_run_triggers_only,
            trigger_top50=args.ambient_trigger_top50,
            trigger_mt=args.ambient_trigger_mt,
            trigger_count_correlation=args.ambient_trigger_count_corr,
            trigger_cross_sample_cv=args.ambient_trigger_cohort_cv,
            decontx_max_iter=args.ambient_decontx_max_iter,
            batch_obs_column=args.ambient_batch_column or None,
        ),
        clustering=ClusteringConfig(
            n_top_genes=args.n_top_genes,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            leiden_resolution=args.leiden_resolution,
            leiden_resolution_sweep=_parse_resolution_sweep(args.leiden_resolution_sweep),
            scale_data=args.scale_data,
        ),
        batch=BatchConfig(
            batch_key=args.batch_key,
            method=args.batch_method,
            harmony_theta=args.harmony_theta,
            harmony_sigma=args.harmony_sigma,
            harmony_max_iter=args.harmony_max_iter,
            harmony_backend=args.harmony_backend,
            scvi_max_epochs=args.scvi_max_epochs,
            scvi_n_latent=args.scvi_n_latent,
            scvi_early_stopping=not args.no_scvi_early_stopping,
            select_integration=args.select_integration,
            integration_margin_mix=args.integration_margin_mix,
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
        optional_modules=_resolve_optional_modules(args),
        markers=_load_markers(args.markers_json),
        gene_signature=GeneSignatureConfig(
            signature_json=Path(args.signature_json) if args.signature_json else None,
        ),
        paper_repro=PaperReproConfig(
            spec_json=Path(args.paper_spec_json) if args.paper_spec_json else None,
            strict=args.paper_repro_strict,
        ),
        composition=CompositionConfig(
            sample_col=args.composition_sample_col or None,
            condition_col=args.composition_condition_col or None,
            contrast_a=args.composition_contrast_a or None,
            contrast_b=args.composition_contrast_b or None,
            covariates=tuple(
                value.strip()
                for value in args.composition_covariates.split(",")
                if value.strip()
            ),
            min_samples_per_condition=args.composition_min_samples_per_condition,
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
        trajectory_root_justification=args.trajectory_root_justification,
        checkpoint=args.checkpoint,
        resume_from=args.resume_from,
        allow_partial_run=args.allow_partial_run,
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
        scientific_profile=args.scientific_profile,
        scientific_non_equivalence_acknowledged=(
            args.scientific_non_equivalence_acknowledged
        ),
        resolved_scientific_parameter_diff=args.resolved_scientific_parameter_diff,
        lazy_read=args.lazy_read,
        doublet_strategy=args.doublet_strategy,
        clustering_engine=args.clustering_engine,
        checkpoint_policy=args.checkpoint_policy,
        cohort_subset=args.cohort_subset,
        annotation_strategy=args.annotation_strategy,
        random_state=args.random_state,
        tissue=args.tissue,
        condition=args.condition,
        validate_context=args.validate_context,
        context_mismatch_threshold=args.context_mismatch_threshold,
        context_min_cells=args.context_min_cells,
        hic_contacts_path=Path(args.hic_contacts_path) if args.hic_contacts_path else None,
        hic_resolution_bp=args.hic_resolution_bp,
        hic_chromsizes_path=Path(args.hic_chromsizes_path) if args.hic_chromsizes_path else None,
        hic_tad_window_bins=args.hic_tad_window_bins,
        hic_tad_boundary_k=args.hic_tad_boundary_k,
    )
    ledger = None
    try:
        _ledger_ctx = type("_LedgerCtx", (), {"cfg": cfg})()
        ledger_root = _run_dir if _run_dir is not None else Path(args.output_dir)
        ledger = RunLedger(_ledger_ctx, args.project, ledger_root)
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
        manifest_result = run_pipeline(cfg, ledger=ledger)
        print(manifest_result)
        manifest, producer_manifest = _load_pipeline_manifest(manifest_result)
        if _run_dir is not None:
            try:
                modules_run = list(manifest.get("modules_run", []))
                bundle_sha256 = manifest.get("bundle_sha256", "")
                write_manifest(
                    _run_dir,
                    project_id=args.project,
                    run_id=_run_id,
                    modules_run=modules_run,
                    bundle_sha256=bundle_sha256 or None,
                    requested_modules=list(manifest.get("requested_modules", modules_run)),
                    planned_modules=list(manifest.get("planned_modules", modules_run)),
                    executed_modules=list(manifest.get("executed_modules", modules_run)),
                    completed_modules=list(manifest.get("completed_modules", modules_run)),
                    skipped_modules=list(manifest.get("skipped_modules", [])),
                    failed_modules=list(manifest.get("failed_modules", [])),
                    overall_status=str(manifest.get("overall_status", "")),
                    producer_manifest=producer_manifest,
                    factory_python_state=manifest.get("factory_python"),
                    produced_on=socket.gethostname(),
                    factory_python_path=_FACTORY_ROOT,
                    extra={"allow_partial_run": bool(args.allow_partial_run)},
                )
            except Exception as _mf_exc:
                logging.getLogger(__name__).warning("manifest write failed: %s", _mf_exc)
                raise RuntimeError("project-root manifest write failed") from _mf_exc
        if manifest.get("failed_modules") and not args.allow_partial_run:
            print(
                "ERROR: requested module failure(s): "
                + ", ".join(str(name) for name in manifest["failed_modules"]),
                file=sys.stderr,
            )
            raise SystemExit(1)
    finally:
        run_shutdown_cleanup()


if __name__ == "__main__":
    main()
