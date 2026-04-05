# singlecell_factory v5.0 — Modular scRNA-seq Pipeline

A comprehensive, production-ready single-cell RNA-seq analysis framework with **mandatory QC + 21 optional analysis modules + automatic dependency resolution + GPU acceleration + categorized output**.

Designed for 10X Genomics datasets. Tested on lung squamous cell carcinoma (LUSC) 3K cells.

## Features

- **3 mandatory modules** (cellranger, QC, doublet detection) ensure data quality baseline
- **21 optional analysis modules** covering the full scRNA-seq workflow
- Automatic topological dependency resolution — just list what you want, dependencies are auto-included
- **GPU acceleration** — rapids-singlecell backend for PCA/UMAP/neighbors (auto-detected)
- **Categorized output** — each module's figures and tables in its own subfolder
- **Checkpoint & resume** — zarr-accelerated checkpoints with h5ad fallback
- **Parallel execution** — thread-safe parallel tiers with cost-aware scheduling
- **Module runtime telemetry** — per-module wall-time automatically stored in `run_manifest.json`
- Multi-backend support: each module auto-detects the best available tool
- Validated against cBioPortal mutation data
- Engineering principle: **accuracy and reproducibility first, performance second**

## Architecture

### Mandatory Modules (always run)

| Module | Function |
|---|---|
| `cellranger` | Load `filtered_feature_bc_matrix`, unified data entry |
| `qc` | Cell/gene filtering (mito/ribo/hemoglobin/n_genes/UMI), QC visualizations |
| `doublet_detection` | Scrublet doublet detection and removal |

### Optional Modules (auto-dependency resolution)

| Module | Function | Depends on |
|---|---|---|
| `clustering` | Normalization, HVG, PCA, UMAP, Leiden clustering | doublet_detection |
| `cell_cycle` | Cell cycle scoring (S/G2M), optional regression | clustering |
| `batch_correction` | Multi-sample batch correction (Harmony/BBKNN/Combat/Scanorama) | clustering |
| `differential_expression` | Cluster marker genes (`wilcoxon` default, configurable), significance filtering | clustering |
| `annotation` | Marker-based cell type annotation with confidence scores | clustering |
| `trajectory` | PAGA trajectory graph + DPT pseudotime + gene expression dynamics | clustering |
| `pseudo_velocity` | Pseudo-RNA velocity with arrow/stream plots (vectorized) | trajectory |
| `rna_velocity` | Real RNA velocity (scVelo: stochastic/dynamical) | clustering |
| `cnv_inference` | Expression-based CNV inference (infercnvpy / vectorized sliding window) | clustering |
| `pathway_analysis` | Gene set enrichment (gseapy/decoupler/built-in Hallmark, BH FDR) | differential_expression |
| `cell_communication` | Ligand-receptor interactions (LIANA / manual L-R scoring) | annotation |
| `gene_regulatory_network` | TF activity inference (decoupler+DoRothEA / manual TF-target) | clustering |
| `validate_cbioportal` | Cross-reference DE genes with cBioPortal mutation data | differential_expression |
| `immune_phenotyping` | 15 immune subtypes + exhaustion/cytotoxicity/activation scores | annotation |
| `tumor_microenvironment` | TME scoring (CYT/TIS/IFN-gamma/ESTIMATE) + checkpoint profiling | annotation |
| `gene_signature_scoring` | 10 built-in cancer signatures + custom JSON signatures | clustering |
| `evolution` | CNV-based clonal clustering, phylogenetic dendrogram, pseudotime-ordered evolution | cnv_inference + trajectory |
| `pseudobulk_de` | Pseudobulk DE (pydeseq2 / Mann-Whitney fallback) — statistically proper multi-sample comparison | differential_expression |
| `cell_fate` | Probabilistic cell fate mapping (CellRank / diffusion-based fallback) | trajectory |
| `composition` | Differential cell type composition analysis (pertpy/scCODA / chi-squared fallback) | annotation |
| `metacell` | Metacell aggregation (SEACells / MiniBatchKMeans fallback) — noise reduction for large datasets | clustering |

### Module Dependency DAG

```
cellranger -> qc -> doublet_detection -> clustering -+-> differential_expression -+-> pathway_analysis
                                                      |                            +-> validate_cbioportal
                                                      |                            +-> pseudobulk_de
                                                      +-> annotation -+-> cell_communication
                                                      |               +-> immune_phenotyping
                                                      |               +-> tumor_microenvironment
                                                      |               +-> composition
                                                      +-> trajectory -+-> pseudo_velocity
                                                      |               +-> cell_fate
                                                      +-> cnv_inference --+--> evolution
                                                      |                        (requires both trajectory + cnv_inference)
                                                      +-> cell_cycle
                                                      +-> batch_correction
                                                      +-> rna_velocity
                                                      +-> gene_regulatory_network
                                                      +-> gene_signature_scoring
                                                      +-> metacell
```

## Project Structure

```
singlecell_factory/
├── workflow/modular/
│   ├── cli.py              # CLI entry point
│   ├── config.py           # Dataclass configuration
│   ├── context.py          # Runtime context (per-module output dirs, checkpointing)
│   ├── pipeline.py         # Pipeline orchestration + dependency resolution
│   ├── perf_baseline.py    # Performance baselines
│   └── modules/            # 24 analysis modules
├── data/raw/               # Input datasets
├── results/                # Pipeline output (each run = timestamped folder with analysis subfolders)
├── tests/                  # Test suite
├── reports/                # Analysis reports
├── environment.yml         # Conda environment
└── pyproject.toml          # Python project config
```

## Installation

```bash
cd /home/zerlinshen/singlecell_factory
conda env create -f environment.yml
conda activate sc10x
export MPLCONFIGDIR=$PWD/.mplconfig
export NUMBA_CACHE_DIR=/tmp/numba_cache
```

`MPLCONFIGDIR` and `NUMBA_CACHE_DIR` are strongly recommended for stable `scanpy/scvelo` startup in some Conda environments.

## Dependency Requirements

- Core: `python>=3.10`, `scanpy`, `anndata`, `numpy`, `scipy`, `pandas`, `matplotlib`, `scikit-learn`
- Mandatory-module runtime: `scrublet`
- Optional backends (auto-detected at runtime):
  - `rapids-singlecell`, `cupy` (GPU clustering)
  - `harmonypy` / `bbknn` / `scanorama` (batch correction)
  - `infercnvpy`, `pybiomart` (CNV)
  - `gseapy`, `decoupler` (pathway / TF activity)
  - `liana` (cell-cell communication)
  - `cellrank` (cell fate)
  - `pydeseq2` (pseudobulk DE)
  - `SEACells` (metacell)
  - `scvelo`, `pysam` (RNA velocity)

If a backend is missing, the corresponding module falls back to an implemented alternative whenever possible.

## Data Preparation (LUSC)

The LUSC 3K dataset should be placed at:

```
data/raw/lung_carcinoma_3k_count/outs/filtered_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

## Quick Start

### Full Analysis (recommended)

```bash
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,cell_cycle,differential_expression,annotation,\
trajectory,pseudo_velocity,cnv_inference,pathway_analysis,cell_communication,\
gene_regulatory_network,immune_phenotyping,tumor_microenvironment,\
gene_signature_scoring
```

### With Checkpointing (crash recovery)

```bash
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression,annotation \
  --checkpoint
```

Resume from a failed module:

```bash
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression,annotation \
  --checkpoint --resume-from annotation
```

### Parallel Execution

```bash
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression,annotation,immune_phenotyping \
  --parallel-workers 4
```

### Immuno-Oncology Deep Analysis

```bash
python -m workflow.modular.cli \
  --project LUSC_immuno \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression,annotation,\
immune_phenotyping,tumor_microenvironment,gene_signature_scoring,pathway_analysis
```

### RNA Velocity from BAM (no loom file)

```bash
python -m workflow.modular.cli \
  --project LUSC_velocity \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,rna_velocity \
  --velocity-bam data/raw/lung_carcinoma_3k_count/outs/possorted_genome_bam.bam \
  --transcriptome-dir /path/to/refdata-gex-GRCh38-2020-A
```

If `--transcriptome-dir` is set, `genes.gtf(.gz)` is auto-discovered from:
`<transcriptome-dir>/genes/genes.gtf(.gz)` (or `<transcriptome-dir>/genes.gtf(.gz)`).
You can still pass `--velocity-gtf` explicitly to override auto-discovery.

## Output Structure

Each run creates an independent timestamped folder under `results/`. Subfolders are named by analysis type:

```
results/LUSC_3k_Analysis_{timestamp}/
├── final_adata.h5ad                          # Final AnnData object
├── run_manifest.json                         # Run metadata
├── module_status.csv                         # Module execution status
│
├── cellranger/                               # Data loading
├── qc/                                       # Quality control
│   ├── qc_violin_pre_filter.png
│   ├── qc_violin_post_filter.png
│   ├── qc_scatter_pre_filter.png
│   └── qc_scatter_post_filter.png
├── doublet_detection/
│   └── doublet_scores.png
├── clustering/
│   ├── pca_variance_explained.png
│   └── umap_leiden.png
├── annotation/
│   ├── cell_type_annotation.csv
│   ├── cluster_majority_cell_type.csv
│   ├── umap_cell_type.png
│   ├── umap_annotation_confidence.png
│   └── cell_type_composition.png
├── cell_cycle/
│   ├── cell_cycle_scores.csv
│   └── cell_cycle_umap.png
├── cnv_inference/
│   ├── cnv_scores.csv
│   ├── cnv_classification.json
│   ├── cnv_score_umap.png
│   └── cnv_heatmap.png
├── differential_expression/
│   ├── marker_genes.csv
│   ├── marker_genes_all.csv
│   ├── marker_top5_by_cluster.csv
│   ├── de_dotplot_top5.png
│   ├── de_heatmap_top5.png
│   └── de_volcano.png
├── gene_regulatory_network/
│   ├── tf_activity_per_cluster.csv
│   ├── tf_top_per_cluster.json
│   └── tf_activity_heatmap.png
├── gene_signature_scoring/
│   ├── gene_signature_scores.csv
│   ├── gene_signature_per_cluster.csv
│   ├── signature_heatmap.png
│   ├── signature_umap.png
│   └── signature_correlation.png
├── trajectory/
│   ├── dpt_pseudotime.csv
│   ├── pseudotime_per_cluster.csv
│   ├── pseudotime_top_genes.csv
│   ├── pseudotime_dpt_umap.png
│   ├── paga_trajectory.png
│   ├── pseudotime_gene_heatmap.png
│   └── pseudotime_violin_per_cluster.png
├── cell_communication/
│   ├── cell_communication_lr.csv
│   └── cell_communication_heatmap.png
├── immune_phenotyping/
│   ├── immune_phenotyping.csv
│   ├── immune_subtype_summary.csv
│   ├── umap_immune_subtype.png
│   ├── immune_exhaustion_umap.png
│   ├── immune_cytotoxicity_umap.png
│   ├── immune_subtype_composition.png
│   └── immune_signature_heatmap.png
├── tumor_microenvironment/
│   ├── tme_scores_per_cell.csv
│   ├── tme_scores_per_cluster.csv
│   ├── checkpoint_expression.csv
│   ├── tme_cyt_umap.png
│   ├── tme_tis_umap.png
│   ├── tme_signature_heatmap.png
│   ├── checkpoint_dotplot.png
│   └── tme_immune_stromal_bar.png
├── pathway_analysis/
│   ├── pathway_enrichment.csv
│   └── pathway_enrichment_bar.png
├── pseudo_velocity/
│   ├── pseudo_velocity_speed.csv
│   ├── pseudo_velocity_per_cluster.csv
│   ├── pseudo_velocity_arrows.png
│   ├── pseudo_velocity_stream.png
│   ├── pseudo_velocity_speed_umap.png
│   └── pseudo_velocity_speed_boxplot.png
├── rna_velocity/
│   ├── velocity_confidence.csv
│   ├── velocity_top_genes.csv
│   ├── velocity_stream_umap.png
│   ├── velocity_grid_umap.png
│   ├── velocity_length_distribution.png
│   ├── velocity_latent_time_umap.png         # dynamical mode only
│   └── velocity_phase_portraits.png          # dynamical mode only
└── evolution/
    ├── evolution_clone_assignment.csv
    ├── evolution_clone_stats.csv
    ├── evolution_clone_markers.csv
    ├── evolution_clone_umap.png
    ├── evolution_clone_composition.png
    ├── evolution_phylo_dendrogram.png
    ├── evolution_timeline.png
    └── evolution_cnv_by_clone.png
```

Each run is fully independent. Multiple runs accumulate under `results/`:

```
results/
├── LUSC_3k_Analysis_20260404_062118/    # Run 1
├── LUSC_immuno_20260405_091500/         # Run 2
└── LUSC_tme_20260405_140000/            # Run 3
```

## CLI Parameters

### General

| Parameter | Description |
|---|---|
| `--project` | Run name (used in output directory naming) |
| `--sample-root` | Dataset root directory |
| `--output-dir` | Output root (default: `results/`) |
| `--optional-modules` | Comma-separated module list (dependencies auto-included) |
| `--checkpoint` | Save checkpoints after each module for crash recovery |
| `--resume-from MODULE` | Resume from a specific module using saved checkpoints |
| `--parallel-workers N` | Number of parallel workers (default: 1 = sequential) |
| `--de-method` | DE method for `scanpy.tl.rank_genes_groups` (default: `wilcoxon`) |
| `--de-n-genes` | Max genes ranked per cluster in DE (default: 300) |
| `--de-pval-threshold` | Adjusted p-value cutoff for DE significance (default: 0.05) |
| `--de-logfc-threshold` | Minimum absolute log fold change for DE (default: 0.25) |
| `--annotation-confidence-threshold` | Minimum score for confident cell type assignment (default: 0.1) |

### QC

| Parameter | Default | Description |
|---|---|---|
| `--min-genes` | 200 | Minimum genes per cell |
| `--max-genes` | 7000 | Maximum genes per cell |
| `--min-counts` | 500 | Minimum UMI counts |
| `--max-counts` | 50000 | Maximum UMI counts |
| `--max-mito-pct` | 20 | Maximum mitochondrial % |
| `--max-ribo-pct` | 50 | Maximum ribosomal % |

### Clustering

| Parameter | Default | Description |
|---|---|---|
| `--n-top-genes` | 3000 | Number of HVGs |
| `--n-pcs` | 40 | PCA dimensions |
| `--n-neighbors` | 15 | k-NN neighbors |
| `--leiden-resolution` | 0.8 | Leiden clustering resolution |

### Batch Correction

| Parameter | Default | Description |
|---|---|---|
| `--batch-key` | sample | Batch column in adata.obs |
| `--batch-method` | harmony | Method: harmony/bbknn/combat/scanorama |

### RNA Velocity

| Parameter | Default | Description |
|---|---|---|
| `--velocity-loom` | - | Loom file with spliced/unspliced counts |
| `--velocity-bam` | - | Path to possorted_genome_bam.bam (Cell Ranger BAM output) |
| `--velocity-gtf` | - | Path to genes.gtf(.gz); optional if `--transcriptome-dir` is set (auto-discovery) |
| `--velocity-mode` | stochastic | scVelo mode: stochastic/dynamical |
| `--velocity-n-jobs` | 4 | Parallel workers for BAM extraction and scVelo dynamics |
| `--velocity-min-shared-counts` | 20 | Minimum shared counts for scVelo gene filtering |
| `--velocity-n-pcs` | 30 | PCA components for scVelo moments |
| `--velocity-n-neighbors` | 30 | Neighbors for scVelo moments |

If no loom file is provided, the module can extract spliced/unspliced counts directly from Cell Ranger BAM output using pysam (parallelized by chromosome). `genes.gtf(.gz)` is taken from `--velocity-gtf` or auto-discovered from `--transcriptome-dir`. The module runs scVelo on an internal adata copy and transfers only cell-level results back, preserving the shared gene index and enabling parallel execution with other modules.
RNA velocity BAM classification uses strict exon/intron rules only (no lossy fast-path).
Policy update (April 5, 2026): all lossy velocity shortcuts were removed from pipeline code.

Velocity BAM extraction cache:
- Default cache dir: `/tmp/singlecell_factory_velocity_cache`
- Override: `SCF_VELOCITY_CACHE_DIR=/path/to/cache`
- Cache key includes BAM/GTF identity + barcode/gene ordering + velocity extraction mode; repeated runs with the same inputs can skip BAM parsing.

On hosts with >=16 logical CPU cores, if `--velocity-n-jobs` is left at default (`4`), BAM extraction workers are auto-bumped to `8` and recorded in `run_manifest.json -> metadata.velocity_extract_n_jobs`.

**Dynamical mode** additionally outputs latent time UMAP and phase portraits for top velocity genes.

### Other

| Parameter | Description |
|---|---|
| `--regress-cell-cycle` | Regress out cell cycle effects |
| `--trajectory-root-cluster` | Leiden cluster ID for DPT root |
| `--cnv-reference-group` | Normal reference cell type for CNV (e.g., Fibroblast) |
| `--signature-json` | Custom gene signature JSON file |
| `--cbioportal-genes` | Comma-separated gene list for cBioPortal validation |

## Multi-Backend Support

| Module | Primary | Fallback |
|---|---|---|
| `pathway_analysis` | gseapy (MSigDB) | decoupler (PROGENy) -> built-in Hallmark (BH FDR) |
| `cell_communication` | LIANA (multi-method consensus) | Manual L-R scoring (18 TME pairs) |
| `gene_regulatory_network` | decoupler + DoRothEA | Manual TF-target scoring (12 key TFs) |
| `cnv_inference` | infercnvpy | Vectorized sliding-window smoothing |
| `batch_correction` | Harmony | BBKNN / Combat / Scanorama |

## Performance Optimizations (v5.0)

### GPU Acceleration
- **rapids-singlecell**: Auto-detected GPU backend for PCA, neighbors, UMAP, and Leiden clustering. Falls back to CPU scanpy when unavailable. Expected **10-50x speedup** on datasets >10K cells with NVIDIA GPU.

### Memory & I/O
- **Memory guard**: Parallel worker count is constrained by estimated AnnData copy size + available RAM to avoid OOM in branch execution.
- **Raw matrix policy**: Normalized/log1p matrix is stored in `adata.raw` before downstream analyses to improve biological interpretability for marker/score modules.
- **Zarr checkpoints**: 3-5x faster checkpoint I/O via `adata.write_zarr()` with automatic h5ad fallback for incompatible key names.
- **Async figure I/O**: Background disk writes with error-resilient flush.

### Computation
- **Clustering**: PCA on HVGs only, optional scaling via `--scale-data`, CPU/GPU branches now both produce UMAP outputs.
- **QC**: Optimized `sc.pp.calculate_qc_metrics` (`log1p=False`, `percent_top=None`).
- **Differential Expression**: Default switched to `wilcoxon` with configurable `--de-method` and `--de-n-genes`.
- **Trajectory**: Pseudotime-gene correlation refactored to sparse-friendly computation, avoiding full-matrix densification.
- **Cell communication**: Vectorized `np.where` scoring — eliminates O(n^2) Python loops.
- **CNV inference**: `scipy.ndimage.uniform_filter1d` vectorized sliding window.
- **RNA velocity**: BAM extraction parallelized by chromosome (~4-5x speedup).

### Pipeline Orchestration
- **Cost-aware tier scheduling**: Heavy modules (rna_velocity, cnv_inference) start first in parallel tiers.
- **Thread-safe status tracking**: `threading.Lock` on module status writes.
- **Copy-on-write branching**: Parallel modules get isolated AnnData copies with proper merge-back of obs/obsm/uns/metadata/module_dirs.
- **Runtime telemetry**: `metadata.module_runtime_sec` + `metadata.pipeline_wall_seconds` are persisted in each `run_manifest.json`.
- **Checkpoint/resume**: `--checkpoint` + `--resume-from` for crash recovery.

### Performance Benchmark Results (LUSC 3K Dataset)

Benchmark script output (`output/performance_benchmark_20260405_optfix/benchmark_compare.json`):

| Metric | Baseline | Optimized | Change |
|---|---|---|---|
| Wall-clock time | 10.04 s | 9.16 s | **-8.79%** |
| Peak RSS memory | 799 MB | 564 MB | **-29.45%** |
| Clusters found | 14 | 15 | ARI = 0.54 |

Latest end-to-end real-data run on **April 5, 2026** (`results/LUSC_FULL_STRICT_ALLMODULES_20260405_173510`):
- Pipeline wall time: **34.513 s**
- Cells: `2588 -> 2382` after QC -> `2377` after doublet removal
- Clusters: `17`
- Module status: **24/24 ok** (3 mandatory + 21 optional)
- Heaviest modules by wall-time: `evolution (12.206s)`, `differential_expression (11.548s)`, `rna_velocity (8.873s)`, `clustering (8.776s)`
- Detailed module usage + result/performance report (LUSC example): `reports/LUSC_模块使用说明与结果性能报告_20260405.md`

RNA velocity bottleneck benchmark on real LUSC data (strict mode, same inputs, **April 5, 2026**):
- Cold run: `results/LUSC_VEL_STRICT_ONLY_20260405_173025`
- Warm run (same cache dir): `results/LUSC_VEL_STRICT_ONLY_WARM_20260405_173113`

| Metric | Cold | Warm | Change |
|---|---:|---:|---:|
| Pipeline wall time | 37.143 s | 15.116 s | **-59.3%** |
| `rna_velocity` module time | 26.422 s | 4.418 s | **-83.3%** |
| Velocity layer loading (`metadata.velocity_extract_seconds`) | 22.006 s | 0.090 s | **-99.6%** |
| Mean velocity confidence | 0.86117 | 0.86117 | identical |

Output reproducibility checks (cold vs warm): identical file hashes for `rna_velocity/velocity_confidence.csv`, `rna_velocity/velocity_top_genes.csv`, and `clustering/umap_leiden.png`.

RNA velocity extraction parallel scaling on the same dataset (strict mode, cold cache):
- `n_jobs=4`: `38.323 s`
- `n_jobs=8`: `26.148 s` (**31.8% faster**)

### Best Practice Configuration Recommendations

1. **Large datasets (>50k cells)**: start with `--parallel-workers 2-4`, then scale up only if RAM headroom is sufficient.
2. **DE for publication**: keep `--de-method wilcoxon`; use `t-test_overestim_var` only when speed is the top priority.
3. **Single-sample runs**: `pseudobulk_de` is expected to skip (requires >=2 biological samples).
4. **Batch correction**: only enable when a real batch column exists; otherwise keep skipped to avoid unnecessary recomputation.
5. **Runtime tuning loop**: inspect `run_manifest.json -> metadata.module_runtime_sec` and optimize the top 3 slowest modules first.
6. **RNA velocity repeated runs**: set `SCF_VELOCITY_CACHE_DIR` to a fast local SSD path and reuse it across reruns.
7. **RNA velocity policy**: keep strict exon/intron classification for all runs; optimize with cache and parallelism instead of lossy shortcuts.
8. **High-core hosts**: default velocity extraction auto-uses 8 workers on >=16 cores; set `--velocity-n-jobs` explicitly to override.

## Quality Verification Checklist

1. Check `module_status.csv` — all modules should show "ok"
2. Review `qc/qc_violin_*.png` — verify filter thresholds are reasonable
3. Review `clustering/umap_leiden.png` — clusters should be well-separated
4. Check `annotation/umap_cell_type.png` — cell types should be biologically coherent
5. Review `immune_phenotyping/immune_signature_heatmap.png` — CD8_exhausted should correlate with exhaustion score
6. Check `tumor_microenvironment/checkpoint_dotplot.png` — checkpoint expression patterns

## Methodology & References

Every analysis module uses publicly recognized, peer-reviewed methods. Below is the complete methodology audit with citations.

---

### 01. Cell Ranger Data Loading (`cellranger`)

| Item | Detail |
|---|---|
| **Method** | 10X Genomics Cell Ranger count matrix loading |
| **Implementation** | `scanpy.read_10x_mtx()` |
| **Input format** | Market Matrix (MTX) sparse format from Cell Ranger |
| **Reference** | Zheng et al., *Nature Communications*, 2017. DOI: [10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049) |

---

### 02. Quality Control (`qc`)

| Item | Detail |
|---|---|
| **Method** | Threshold-based cell and gene filtering |
| **Metrics** | n_genes_by_counts, total_counts, pct_counts_mt, pct_counts_ribo, pct_counts_hb |
| **Implementation** | `scanpy.pp.calculate_qc_metrics()`, `scanpy.pp.filter_genes()` |
| **Gene classes** | Mitochondrial (MT-), Ribosomal (RPS/RPL), Hemoglobin (HBA/HBB) |
| **Defaults** | 200-7000 genes, 500-50000 UMI, <=20% mito, <=50% ribo |
| **Reference** | Luecken & Theis, *Molecular Systems Biology*, 2019. DOI: [10.15252/msb.20188746](https://doi.org/10.15252/msb.20188746) |

---

### 03. Doublet Detection (`doublet_detection`)

| Item | Detail |
|---|---|
| **Method** | Scrublet — synthetic doublet simulation |
| **Algorithm** | Simulates doublets by averaging random cell pairs, scores real cells by k-NN similarity to synthetic doublets |
| **Implementation** | `scrublet.Scrublet.scrub_doublets(min_counts=2, min_cells=3, n_prin_comps=30)` |
| **Fallback** | For tiny/degenerate datasets, pipeline falls back to all-singlets (records `doublet_method=fallback_all_singlets`) instead of aborting |
| **Parameters** | expected_doublet_rate=0.06, min_gene_variability_pctl=85 |
| **Reference** | **Wolock et al., *Cell Systems*, 2019.** DOI: [10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005) |

---

### 04. Clustering (`clustering`)

| Item | Detail |
|---|---|
| **Normalization** | CPM (counts-per-million) + log1p, via `scanpy.pp.normalize_total()` |
| **HVG selection** | Seurat flavor HVG selection (`scanpy.pp.highly_variable_genes`) |
| **PCA** | Truncated SVD via ARPACK solver on HVGs only (memory-optimized) |
| **Neighbor graph** | k-nearest neighbors (k=15 default) in PCA space |
| **Dimensionality reduction** | UMAP (Uniform Manifold Approximation and Projection) |
| **Clustering** | Leiden community detection algorithm (igraph, undirected, resolution=0.8) |
| **References** | **McInnes et al., *JOSS*, 2018.** DOI: [10.21105/joss.00861](https://doi.org/10.21105/joss.00861) (UMAP); **Traag et al., *Scientific Reports*, 2019.** DOI: [10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z) (Leiden); **Stuart et al., *Cell*, 2019.** DOI: [10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031) (Seurat v3 HVG) |

---

### 05. Cell Type Annotation (`annotation`)

| Item | Detail |
|---|---|
| **Method** | Gene set scoring with maximum-score assignment |
| **Implementation** | `scanpy.tl.score_genes()` per cell type marker panel |
| **Confidence metric** | max_score - second_max_score (vectorized via `np.partition`) |
| **Marker panels** | 10 cell types: Tumor epithelial, T cell, NK cell, B cell, Myeloid/Macro, Fibroblast, Endothelial, Plasma cell, Mast cell, Dendritic cell |
| **Unknown threshold** | Cells with confidence < 0.1 labeled "Unknown" |
| **Reference** | **Tirosh et al., *Science*, 2016.** DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |

---

### 06. Cell Cycle Scoring (`cell_cycle`)

| Item | Detail |
|---|---|
| **Method** | S-phase and G2M-phase gene set scoring |
| **Gene sets** | 47 S-phase genes + 49 G2M-phase genes (Tirosh/Regev lab) |
| **Implementation** | `scanpy.tl.score_genes_cell_cycle()` |
| **Optional** | `scanpy.pp.regress_out(["S_score", "G2M_score"])` to remove cell cycle effects |
| **Reference** | **Tirosh et al., *Science*, 2016.** DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |

---

### 07. CNV Inference (`cnv_inference`)

| Item | Detail |
|---|---|
| **Method (primary)** | Expression-based sliding-window CNV inference |
| **Algorithm** | Center gene expression by reference mean, smooth per chromosome using `scipy.ndimage.uniform_filter1d`, compute per-cell variance as CNV score |
| **Gene positions** | Auto-fetched from Ensembl BioMart via `pybiomart` if not pre-annotated |
| **Classification** | Percentile thresholding (default 75th) for malignant vs normal |
| **Method (fallback)** | `infercnvpy.tl.infercnv()` + `infercnvpy.tl.cnv_score()` |
| **References** | **Patel et al., *Science*, 2014.** DOI: [10.1126/science.1254257](https://doi.org/10.1126/science.1254257); **Tirosh et al., *Science*, 2016.** DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |

---

### 08. Differential Expression (`differential_expression`)

| Item | Detail |
|---|---|
| **Statistical test** | Wilcoxon rank-sum by default (configurable via `--de-method`) |
| **Implementation** | `scanpy.tl.rank_genes_groups(method=<de_method>, use_raw=False, n_genes=<de_n_genes>, pts=True)` |
| **Multiple testing** | Benjamini-Hochberg FDR correction (scanpy default) |
| **Significance filter** | adjusted p-value < 0.05 AND |log2FC| > 0.25 (configurable via `--de-pval-threshold`, `--de-logfc-threshold`) |
| **Output** | Ranked genes per cluster (`--de-n-genes`, default 300) with scores/p-values/fold changes |
| **References** | **Wilcoxon, *Biometrics Bulletin*, 1945.** DOI: [10.2307/3001968](https://doi.org/10.2307/3001968); **Benjamini & Hochberg, *JRSS-B*, 1995.** DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) |

---

### 09. Gene Regulatory Network (`gene_regulatory_network`)

| Item | Detail |
|---|---|
| **Method (primary)** | decoupler + DoRothEA regulon database |
| **Algorithm** | Univariate Linear Model (ULM) — infers TF activity from target gene expression |
| **Regulon confidence** | Levels A (high), B (medium), C (curated) from DoRothEA |
| **Implementation** | `decoupler.get_dorothea(organism="human", levels=["A","B","C"])`, `decoupler.run_ulm()` |
| **Fallback** | Gene set scoring for 12 curated TF-target sets (TP53, MYC, STAT1, NFKB1, HIF1A, SOX2, FOXP3, E2F1, NOTCH1, SMAD3, JUN, AP-1) |
| **References** | **Garcia-Alonso et al., *Genome Research*, 2019.** DOI: [10.1101/gr.240663.118](https://doi.org/10.1101/gr.240663.118) (DoRothEA); **Badia-i-Mompel et al., *Bioinformatics Advances*, 2022.** DOI: [10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016) (decoupler) |

---

### 10. Gene Signature Scoring (`gene_signature_scoring`)

| Item | Detail |
|---|---|
| **Method** | Scanpy gene set scoring (`sc.tl.score_genes`) |
| **Built-in signatures** | 10 cancer hallmark panels (see table below) |
| **Custom input** | User-provided JSON: `{"name": ["GENE1", "GENE2", ...]}` |
| **Minimum genes** | Requires >= 2 genes present per signature |

**Built-in Cancer Signatures:**

| Signature | Genes | Reference |
|---|---|---|
| Proliferation | MKI67, TOP2A, PCNA, MCM2, MCM6, CDK1, CCNB1, CCNB2 | Standard oncology panel |
| Apoptosis resistance | BCL2, BCL2L1, MCL1, BIRC5, XIAP, CFLAR | BCL2 family |
| Angiogenesis | VEGFA, VEGFB, FLT1, KDR, PECAM1, ANGPT2, NRP1 | VEGF/FLT pathway |
| Invasion/metastasis | MMP2, MMP9, MMP14, SNAI1, TWIST1, VIM, CDH2 | MMP/EMT markers |
| EMT mesenchymal | VIM, CDH2, FN1, SNAI2, ZEB1, ZEB2, TWIST1, MMP2 | Tan et al., *EMBO Mol Med*, 2014. DOI: [10.15252/emmm.201404208](https://doi.org/10.15252/emmm.201404208) |
| EMT epithelial | CDH1, EPCAM, KRT8, KRT18, KRT19, CLDN4, OCLN | Tan et al., *EMBO Mol Med*, 2014. DOI: [10.15252/emmm.201404208](https://doi.org/10.15252/emmm.201404208) |
| Stemness | POU5F1, NANOG, SOX2, KLF4, MYC, LIN28A, SALL4, BMI1 | Malta et al., *Cell*, 2018. DOI: [10.1016/j.cell.2018.03.034](https://doi.org/10.1016/j.cell.2018.03.034) |
| Hypoxia | VEGFA, SLC2A1, HK2, LDHA, PGK1, CA9, BNIP3, ENO1 | Buffa et al., *Br J Cancer*, 2010. DOI: [10.1038/sj.bjc.6605450](https://doi.org/10.1038/sj.bjc.6605450) |
| DNA damage response | BRCA1, BRCA2, ATM, ATR, RAD51, CHEK1, CHEK2, TP53 | DDR pathway |
| Glycolysis | HK2, PFKP, PKM, LDHA, ENO1, GAPDH, TPI1, ALDOA | Warburg effect |

---

### 11. Trajectory / Pseudotime (`trajectory`)

| Item | Detail |
|---|---|
| **Methods** | PAGA trajectory graph + Diffusion Pseudotime (DPT) + gene expression dynamics |
| **PAGA** | Partition-based graph abstraction — discovers trajectory topology between clusters |
| **DPT** | Diffusion map embedding + geodesic distance from root cell for temporal ordering |
| **Gene trends** | Top 30 genes correlated with pseudotime, smoothed heatmap visualization |
| **Implementation** | `scanpy.tl.paga()`, `scanpy.tl.diffmap()`, `scanpy.tl.dpt()` |
| **References** | **Wolf et al., *Genome Biology*, 2019.** DOI: [10.1186/s13059-019-1663-x](https://doi.org/10.1186/s13059-019-1663-x) (PAGA); **Haghverdi et al., *Nature Methods*, 2016.** DOI: [10.1038/nmeth.3971](https://doi.org/10.1038/nmeth.3971) (DPT) |

---

### 12. Cell Communication (`cell_communication`)

| Item | Detail |
|---|---|
| **Method (primary)** | LIANA multi-method consensus scoring |
| **Algorithm** | Aggregates 4+ L-R scoring methods (CellPhoneDB, NATMI, SingleCellSignalR, Connectome) into consensus rank |
| **Implementation** | `liana.mt.rank_aggregate(groupby="cell_type", resource_name="consensus")` |
| **Fallback** | Manual L-R scoring: pre-computed mean expression matrix, outer-product scoring for 18 curated TME L-R pairs (PD-L1/PD-1, VEGFA/KDR, TGFB1/TGFBR2, etc.) |
| **Reference** | **Dimitrov et al., *Nature Communications*, 2022.** DOI: [10.1038/s41467-022-30755-0](https://doi.org/10.1038/s41467-022-30755-0) |

---

### 13. Immune Phenotyping (`immune_phenotyping`)

| Item | Detail |
|---|---|
| **Method** | Gene set scoring for 15 immune subtypes + 3 functional signatures |
| **Implementation** | `scanpy.tl.score_genes()` per subtype panel |
| **Assignment** | Argmax scoring restricted to cells annotated as immune parent types |

**15 Immune Subtypes:**

| Subtype | Key Markers | Source |
|---|---|---|
| CD4 naive | CCR7, LEF1, SELL, TCF7, IL7R | Zheng et al., *Cell*, 2017 |
| CD4 memory | IL7R, S100A4, ANXA1, CCR6, AQP3 | Zheng et al., *Cell*, 2017 |
| Treg | FOXP3, IL2RA, CTLA4, IKZF2, TNFRSF18 | Zheng et al., *Cell*, 2017 |
| Th1 / Th2 / Th17 | TBX21 / GATA3 / RORC + cytokines | Zhang et al., *Nature*, 2018 |
| CD8 effector | CD8A, GZMB, PRF1, NKG7, GNLY, IFNG | Zheng et al., *Cell*, 2017 |
| CD8 memory | CD8A, IL7R, GPR183, SELL, TCF7 | Zhang et al., *Nature*, 2018 |
| CD8 exhausted | CD8A, PDCD1, LAG3, HAVCR2, TIGIT, TOX, CTLA4 | Zheng et al., *Cell*, 2017 |
| NK cytotoxic | NKG7, GNLY, KLRD1, KLRF1, NCAM1, FCGR3A | Tirosh et al., *Science*, 2016 |
| Macro M1 | CD68, NOS2, IL1B, TNF, CXCL10, CD80 | Zhang et al., *Nature*, 2018 |
| Macro M2 | CD68, CD163, MRC1, MSR1, TGFB1, IL10 | Zhang et al., *Nature*, 2018 |
| cDC1 / cDC2 / pDC | CLEC9A / CD1C / LILRA4 + subtype markers | Zhang et al., *Nature*, 2018 |

**Functional Scores:**
- **Exhaustion**: LAG3, HAVCR2, TIGIT, PDCD1, CTLA4, TOX, TOX2, ENTPD1
- **Cytotoxicity**: GZMA, GZMB, GZMK, GZMH, PRF1, NKG7, GNLY, FASLG
- **Activation**: CD69, CD38, HLA-DRA, ICOS, TNFRSF9, TNFRSF4

| Reference | DOI |
|---|---|
| **Zheng et al., *Cell*, 2017** — "Landscape of Infiltrating T Cells in Liver Cancer" | [10.1016/j.cell.2017.05.035](https://doi.org/10.1016/j.cell.2017.05.035) |
| **Zhang et al., *Nature*, 2018** — "Lineage tracking reveals dynamic relationships of T cells in colorectal cancer" | [10.1038/s41586-018-0694-x](https://doi.org/10.1038/s41586-018-0694-x) |
| **Tirosh et al., *Science*, 2016** — "Dissecting the multicellular ecosystem of metastatic melanoma" | [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |

---

### 14. Tumor Microenvironment (`tumor_microenvironment`)

| Signature | Genes | Method | Reference |
|---|---|---|---|
| **CYT** (Cytolytic Activity) | GZMA, PRF1 | Geometric mean: sqrt(GZMA x PRF1) | Rooney et al., *Cell*, 2015. DOI: [10.1016/j.cell.2014.12.033](https://doi.org/10.1016/j.cell.2014.12.033) |
| **TIS** (T-cell Inflamed, 18-gene) | CD27, CD274, CD276, CD8A, CMKLR1, CXCL9, CXCR6, HLA-DQA1, HLA-DRB1, HLA-E, IDO1, LAG3, NKG7, PDCD1LG2, PSMB10, STAT1, TIGIT, TGFB1 | Gene set score | Ayers et al., *JCI*, 2017. DOI: [10.1172/JCI91190](https://doi.org/10.1172/JCI91190) |
| **IFN-gamma** (10-gene) | IFNG, STAT1, CCR5, CXCL9, CXCL10, CXCL11, IDO1, PRF1, GZMA, HLA-DRA | Gene set score | Ayers et al., *JCI*, 2017. DOI: [10.1172/JCI91190](https://doi.org/10.1172/JCI91190) |
| **Immunosuppression** | IDO1, CD274, PDCD1LG2, HAVCR2, LAG3, TGFB1, IL10, VEGFA, ARG1 | Gene set score | Composite |
| **Immune ESTIMATE** | CD2, CD3D, CD3E, CD3G, CD8A, CD19, CD79A, MS4A1, LCK, GZMB, NKG7, PRF1, GNLY | Gene set score | Yoshihara et al., *Nat Commun*, 2013. DOI: [10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612) |
| **Stromal ESTIMATE** | DCN, LUM, COL1A1, COL3A1, COL5A1, FAP, ACTA2, VIM, FN1, SPARC, POSTN, THY1 | Gene set score | Yoshihara et al., *Nat Commun*, 2013. DOI: [10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612) |

**Checkpoint Molecules Profiled:** PD-1 (PDCD1), PD-L1 (CD274), PD-L2 (PDCD1LG2), CTLA-4, LAG-3, TIM-3 (HAVCR2), TIGIT, VISTA (VSIR), IDO1, B7-H3 (CD276)

---

### 15. Pathway Analysis (`pathway_analysis`)

| Backend | Method | Implementation | Reference |
|---|---|---|---|
| **gseapy** (primary) | Over-Representation Analysis with MSigDB Hallmark gene sets | `gseapy.enrich(gene_sets="MSigDB_Hallmark_2020")` | Subramanian et al., *PNAS*, 2005. DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102); Liberzon et al., *Cell Systems*, 2015. DOI: [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004) |
| **decoupler + PROGENy** | Multivariate Linear Model (MLM) pathway activity scoring | `decoupler.run_mlm()` with PROGENy top 300 footprint genes | Schubert et al., *Nature Communications*, 2018. DOI: [10.1038/s41467-017-02391-6](https://doi.org/10.1038/s41467-017-02391-6) |
| **Fallback** | Hypergeometric overlap test with built-in Hallmark gene sets | `scipy.stats.hypergeom.sf()` + Benjamini-Hochberg FDR | Standard hypergeometric test |

---

### 16. Pseudo-Velocity (`pseudo_velocity`)

| Item | Detail |
|---|---|
| **Method** | Local pseudotime gradient estimation in UMAP space with full visualization |
| **Algorithm** | For each cell, compute k-NN (k=15), calculate velocity = mean(spatial_direction x pseudotime_gradient) across neighbors. Produces quiver arrows, IDW-interpolated stream plots, and per-cluster speed statistics |
| **Implementation** | Vectorized NumPy broadcasting + `sklearn.neighbors.NearestNeighbors` + `scipy.spatial.cKDTree` for grid interpolation |
| **Visualizations** | Arrow field, stream plot, speed UMAP overlay, per-cluster speed boxplot |
| **Reference** | Inspired by La Manno et al., *Nature*, 2018. DOI: [10.1038/s41586-018-0414-6](https://doi.org/10.1038/s41586-018-0414-6) (RNA velocity concept adapted to pseudotime gradients) |

---

### 17. Batch Correction (`batch_correction`)

| Method | Algorithm | Reference |
|---|---|---|
| **Harmony** (default) | Iterative PCA-based soft clustering with batch diversity maximization | Korsunsky et al., *Nature Methods*, 2019. DOI: [10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0) |
| **BBKNN** | Batch-balanced k-nearest neighbors graph construction | Polanski et al., *Bioinformatics*, 2020. DOI: [10.1093/bioinformatics/btz625](https://doi.org/10.1093/bioinformatics/btz625) |
| **ComBat** | Empirical Bayes linear model batch adjustment | Johnson et al., *Biostatistics*, 2007. DOI: [10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037) |
| **Scanorama** | Mutual nearest neighbor panoramic stitching | Hie et al., *Nature Biotechnology*, 2019. DOI: [10.1038/s41587-019-0113-3](https://doi.org/10.1038/s41587-019-0113-3) |

---

### 18. RNA Velocity (`rna_velocity`)

| Item | Detail |
|---|---|
| **Method** | scVelo — RNA velocity from splicing kinetics |
| **Stochastic mode** | Moment-based estimation of splicing/degradation rates |
| **Dynamical mode** | Full transcriptome kinetics model with latent time inference |
| **Implementation** | `scvelo.tl.velocity()`, `scvelo.tl.velocity_graph()`, `scvelo.tl.velocity_confidence()` |
| **Architecture** | Copy-on-write execution (`adata.copy()`), then transfer only cell-level outputs back to main adata (non-mutating, parallel-safe) |
| **Requirements** | Spliced/unspliced count layers (from loom file, velocyto, or BAM extraction) |
| **BAM extraction** | If no loom file is provided, spliced/unspliced counts are extracted directly from Cell Ranger BAM output (`possorted_genome_bam.bam`) using pysam + GTF-based exon/intron classification |
| **Outputs** | `velocity_stream_umap.png`, `velocity_grid_umap.png`, `velocity_length_distribution.png`, `velocity_confidence.csv`, `velocity_top_genes.csv` (+ `velocity_latent_time_umap.png`, `velocity_phase_portraits.png` in dynamical mode) |
| **Reference** | **Bergen et al., *Nature Biotechnology*, 2020.** DOI: [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3) |

---

### 19. cBioPortal Validation (`validate_cbioportal`)

| Item | Detail |
|---|---|
| **Method** | Cross-reference DE genes with TCGA large-cohort mutation data |
| **API** | cBioPortal REST API v2 (public, no authentication required) |
| **Endpoints** | Gene lookup, sample counts, mutation frequency profiling |
| **Default study** | LUSC TCGA Pan-Cancer Atlas 2018 |
| **Reference** | **Cerami et al., *Cancer Discovery*, 2012.** DOI: [10.1158/2159-8290.CD-12-0095](https://doi.org/10.1158/2159-8290.CD-12-0095); **Gao et al., *Science Signaling*, 2013.** DOI: [10.1126/scisignal.2004088](https://doi.org/10.1126/scisignal.2004088) |

---

### 20. Tumor Evolution (`evolution`)

| Item | Detail |
|---|---|
| **Method** | CNV-based clonal clustering + pseudotime-ordered evolution + phylogenetic reconstruction |
| **Clone identification** | Hierarchical clustering (Ward's method) on per-cell CNV profiles using correlation distance |
| **Phylogenetic tree** | Clone centroid dendrogram from `scipy.cluster.hierarchy` |
| **Evolution timeline** | Pseudotime distribution per clone reveals temporal ordering of clonal expansion |
| **Clone markers** | Wilcoxon rank-sum test for clone-specific DE genes with BH FDR correction |
| **Implementation** | `scipy.cluster.hierarchy.linkage()`, `scipy.spatial.distance.pdist()`, `scanpy.tl.rank_genes_groups()` |
| **References** | **Patel et al., *Science*, 2014.** DOI: [10.1126/science.1254257](https://doi.org/10.1126/science.1254257); **Gao et al., *Nature Biotechnology*, 2021.** DOI: [10.1038/s41587-020-00795-2](https://doi.org/10.1038/s41587-020-00795-2) (CopyKAT); **Tirosh et al., *Science*, 2016.** DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |

## Automated Reference Management

The pipeline includes an automated mechanism to keep the `Complete Citation List` up-to-date:

1. **Module-Level Citations**: New modules should follow [module_template.py](file:///home/zerlinshen/singlecell_factory/workflow/modular/modules/module_template.py) and define a `__references__` dictionary.
2. **Auto Discovery**: [update_references.py](file:///home/zerlinshen/singlecell_factory/scripts/update_references.py) parses `__references__` and also scans module source for DOI patterns, then enriches metadata via Crossref API when available.
3. **Pre-commit Hook**: [pre-commit](file:///home/zerlinshen/singlecell_factory/.githooks/pre-commit) enforces sync; if references changed, commit is blocked until README is staged.
4. **Enable Hook**: run `git config core.hooksPath .githooks`.

## Validation Protocol (LUSC Real Run)

Use [validate_optimizations.py](file:///home/zerlinshen/singlecell_factory/scripts/validate_optimizations.py) to validate acceleration without sacrificing reliability.

- **Baseline set**: LUSC 3K (`data/raw/lung_carcinoma_3k_count/outs/filtered_feature_bc_matrix`).
- **Control experiment**: baseline vs optimized with same hardware, same seed, repeated runs (`--repeats`).
- **Core metrics**: ARI, pairwise Precision/Recall/F1, marker-gene F1, runtime, data-quality audit (`data_quality_audit.csv`).
- **Stats**: 95% CI + significance test (`runtime_p_value`, `accuracy_p_value`) with threshold `p < 0.05`.
- **Tolerance**: max metric loss threshold configurable (`--loss-threshold-pct`, default `0.5`).
- **Robustness tests**: extreme dropout, low cell count, NaN injection sanitization.
- **Rollback policy**: if thresholds fail, report emits `ROLLBACK_REQUIRED` and automated rollback steps.

Example:

```bash
python scripts/validate_optimizations.py --repeats 5 --warmup-runs 1 --seed 42
```

Outputs are written to `reports/validation/`:
- `optimization_validation_report.json`
- `optimization_validation_report.md`
- `perf_runs.csv`
- `accuracy_runs.csv`
- `data_quality_audit.csv`
- `robustness_results.csv`
- `runtime_boxplot.png`
- `metric_distribution.png`
- `marker_f1_distribution.png`

## Standardized Module Docs

Use [MODULE_TECH_DOC_TEMPLATE.md](file:///home/zerlinshen/singlecell_factory/docs/MODULE_TECH_DOC_TEMPLATE.md) for every new module technical document. The template standardizes:
- Functional description
- Implementation details and complexity
- Benchmark and validation evidence
- Citation and best-practice references

---

### Complete Citation List

| # | Reference | DOI | Used by |
|---|---|---|---|
| 1 | Zheng et al., *Nature Communications*, 2017 | [10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049) | `cellranger` |
| 2 | Luecken & Theis, *Molecular Systems Biology*, 2019 | [10.15252/msb.20188746](https://doi.org/10.15252/msb.20188746) | `qc` |
| 3 | Wolock et al., *Cell Systems*, 2019 | [10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005) | `doublet_detection` |
| 4 | Stuart et al., *Cell*, 2019 | [10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031) | `clustering` (Seurat HVG) |
| 5 | McInnes et al., *JOSS*, 2018 | [10.21105/joss.00861](https://doi.org/10.21105/joss.00861) | `clustering` (UMAP) |
| 6 | Traag et al., *Scientific Reports*, 2019 | [10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z) | `clustering` (Leiden) |
| 7 | Tirosh et al., *Science*, 2016 | [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) | `cell_cycle`, `annotation`, `immune_phenotyping`, `cnv_inference` |
| 8 | Korsunsky et al., *Nature Methods*, 2019 | [10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0) | `batch_correction` (Harmony) |
| 9 | Polanski et al., *Bioinformatics*, 2020 | [10.1093/bioinformatics/btz625](https://doi.org/10.1093/bioinformatics/btz625) | `batch_correction` (BBKNN) |
| 10 | Johnson et al., *Biostatistics*, 2007 | [10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037) | `batch_correction` (ComBat) |
| 11 | Hie et al., *Nature Biotechnology*, 2019 | [10.1038/s41587-019-0113-3](https://doi.org/10.1038/s41587-019-0113-3) | `batch_correction` (Scanorama) |
| 12 | Wilcoxon, *Biometrics Bulletin*, 1945 | [10.2307/3001968](https://doi.org/10.2307/3001968) | `differential_expression` |
| 13 | Benjamini & Hochberg, *JRSS-B*, 1995 | [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) | `differential_expression` (FDR) |
| 14 | Haghverdi et al., *Nature Methods*, 2016 | [10.1038/nmeth.3971](https://doi.org/10.1038/nmeth.3971) | `trajectory` (DPT) |
| 15 | Wolf et al., *Genome Biology*, 2019 | [10.1186/s13059-019-1663-x](https://doi.org/10.1186/s13059-019-1663-x) | `trajectory` (PAGA) |
| 16 | La Manno et al., *Nature*, 2018 | [10.1038/s41586-018-0414-6](https://doi.org/10.1038/s41586-018-0414-6) | `pseudo_velocity` |
| 17 | Bergen et al., *Nature Biotechnology*, 2020 | [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3) | `rna_velocity` (scVelo) |
| 18 | Patel et al., *Science*, 2014 | [10.1126/science.1254257](https://doi.org/10.1126/science.1254257) | `cnv_inference`, `evolution` |
| 19 | Subramanian et al., *PNAS*, 2005 | [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102) | `pathway_analysis` (GSEA) |
| 20 | Liberzon et al., *Cell Systems*, 2015 | [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004) | `pathway_analysis` (MSigDB) |
| 21 | Schubert et al., *Nature Communications*, 2018 | [10.1038/s41467-017-02391-6](https://doi.org/10.1038/s41467-017-02391-6) | `pathway_analysis` (PROGENy) |
| 22 | Dimitrov et al., *Nature Communications*, 2022 | [10.1038/s41467-022-30755-0](https://doi.org/10.1038/s41467-022-30755-0) | `cell_communication` (LIANA) |
| 23 | Garcia-Alonso et al., *Genome Research*, 2019 | [10.1101/gr.240663.118](https://doi.org/10.1101/gr.240663.118) | `gene_regulatory_network` (DoRothEA) |
| 24 | Badia-i-Mompel et al., *Bioinformatics Advances*, 2022 | [10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016) | `gene_regulatory_network` (decoupler) |
| 25 | Cerami et al., *Cancer Discovery*, 2012 | [10.1158/2159-8290.CD-12-0095](https://doi.org/10.1158/2159-8290.CD-12-0095) | `validate_cbioportal` |
| 26 | Gao et al., *Science Signaling*, 2013 | [10.1126/scisignal.2004088](https://doi.org/10.1126/scisignal.2004088) | `validate_cbioportal` |
| 27 | Zheng et al., *Cell*, 2017 | [10.1016/j.cell.2017.05.035](https://doi.org/10.1016/j.cell.2017.05.035) | `immune_phenotyping` |
| 28 | Zhang et al., *Nature*, 2018 | [10.1038/s41586-018-0694-x](https://doi.org/10.1038/s41586-018-0694-x) | `immune_phenotyping` |
| 29 | Rooney et al., *Cell*, 2015 | [10.1016/j.cell.2014.12.033](https://doi.org/10.1016/j.cell.2014.12.033) | `tumor_microenvironment` (CYT) |
| 30 | Ayers et al., *JCI*, 2017 | [10.1172/JCI91190](https://doi.org/10.1172/JCI91190) | `tumor_microenvironment` (TIS, IFN-gamma) |
| 31 | Yoshihara et al., *Nature Communications*, 2013 | [10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612) | `tumor_microenvironment` (ESTIMATE) |
| 32 | Malta et al., *Cell*, 2018 | [10.1016/j.cell.2018.03.034](https://doi.org/10.1016/j.cell.2018.03.034) | `gene_signature_scoring` (stemness) |
| 33 | Tan et al., *EMBO Mol Med*, 2014 | [10.15252/emmm.201404208](https://doi.org/10.15252/emmm.201404208) | `gene_signature_scoring` (EMT) |
| 34 | Buffa et al., *Br J Cancer*, 2010 | [10.1038/sj.bjc.6605450](https://doi.org/10.1038/sj.bjc.6605450) | `gene_signature_scoring` (hypoxia) |
| 35 | Gao et al., *Nature Biotechnology*, 2021 | [10.1038/s41587-020-00795-2](https://doi.org/10.1038/s41587-020-00795-2) | `evolution` (CopyKAT clonal analysis) |
---

## Results

Each pipeline run creates an independent folder under `results/`. Example verified run:

```
results/LUSC_REVIEW_OPT_FIX_20260405_20260405_142200/
├── cellranger/                   (data loading)
├── qc/                           (4 QC plots)
├── doublet_detection/            (1 doublet score plot)
├── clustering/                   (2 plots: PCA elbow, UMAP)
├── annotation/                   (3 plots + 2 CSV)
├── cell_cycle/                   (1 plot + 1 CSV)
├── cnv_inference/                (2 plots + 1 CSV + 1 JSON)
├── differential_expression/      (3 plots + 3 CSV)
├── gene_regulatory_network/      (1 plot + 1 CSV + 1 JSON)
├── gene_signature_scoring/       (3 plots + 2 CSV)
├── trajectory/                   (4 plots + 3 CSV)
├── cell_communication/           (1 plot + 1 CSV)
├── immune_phenotyping/           (5 plots + 2 CSV)
├── tumor_microenvironment/       (5 plots + 3 CSV)
├── pathway_analysis/             (1 plot + 1 CSV)
├── pseudo_velocity/              (4 plots + 2 CSV)
├── final_adata.h5ad
├── evolution/                    (4 plots + 3 CSV)
├── module_status.csv             (22/22 modules: all OK)
└── run_manifest.json
```

## Legacy Workflows

Old entry points are still available: `workflow/standard.py`, `workflow/velocity.py`, `workflow/benchmark.py`.

Recommended: use `python -m workflow.modular.cli` for all new analysis.
