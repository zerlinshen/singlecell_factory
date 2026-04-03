# Single-Cell RNA-seq Analysis Protocol

A step-by-step guide for beginners using the **singlecell_factory** pipeline.

---

## 1. What Is This Pipeline?

### 1.1 Background: Single-Cell RNA Sequencing

Traditional RNA sequencing ("bulk RNA-seq") measures gene expression averaged across millions of cells. You get one number per gene — but tumors, immune systems, and developing tissues contain dozens of different cell types mixed together. You lose all that complexity.

**Single-cell RNA sequencing (scRNA-seq)** measures gene expression in *individual cells*. A typical experiment captures 1,000-100,000 cells, each with expression levels for 20,000-30,000 genes. This lets you:

- Discover cell types you didn't know were there
- See how cancer cells differ from immune cells within the same tumor
- Track how cells change over time (pseudotime / trajectory analysis)
- Identify which genes are turned on/off in each cell type

### 1.2 What This Pipeline Does

This pipeline takes raw 10X Genomics scRNA-seq data and runs a complete analysis workflow:

```
Raw count matrix (cells x genes)
    |
    v
Quality control --> Doublet removal --> Clustering --> Cell type annotation
    |                                       |
    v                                       v
Differential expression         Trajectory / Pseudotime
    |                                       |
    v                                       v
Pathway analysis              Pseudo-velocity (flow arrows)
    |
    v
Immune phenotyping --> TME scoring --> Tumor evolution
```

The pipeline produces **77 output files** (figures + tables) organized into 18 analysis folders.

### 1.3 Who Is This For?

- Bioinformatics students learning scRNA-seq analysis
- Biologists who want to analyze their own 10X data
- Researchers who need a reproducible, citable analysis pipeline

### 1.4 Example Dataset

Throughout this protocol, we use **lung squamous cell carcinoma (LUSC)** data from 10X Genomics — approximately 3,000 cells from a human lung tumor sample. By the end, we identified:

- **2,377 cells** passing quality filters (from 2,588 raw)
- **17 clusters** of cells
- **10 cell types** (T cells, macrophages, B cells, tumor epithelial, etc.)
- **3 evolutionary clones** with distinct CNV profiles
- **6,551 differentially expressed genes** across clusters

---

## 2. Prerequisites

### 2.1 Software

```bash
# Create conda environment
conda env create -f environment.yml
conda activate sc10x
```

**Key dependencies** (installed automatically):
| Package | Purpose |
|---|---|
| scanpy | Core scRNA-seq analysis |
| anndata | Data structure for single-cell data |
| scrublet | Doublet detection |
| igraph + leidenalg | Community detection (clustering) |
| infercnvpy | Copy number variation inference |
| pybiomart | Gene position annotation |
| matplotlib | Visualization |

### 2.2 Input Data

You need a **Cell Ranger output directory** containing the filtered count matrix:

```
data/raw/your_dataset/outs/filtered_feature_bc_matrix/
    barcodes.tsv.gz    # Cell barcodes (one per cell)
    features.tsv.gz    # Gene names and IDs
    matrix.mtx.gz      # Sparse count matrix (cells x genes)
```

This is the standard output from 10X Genomics Cell Ranger `count` pipeline. If you have FASTQ files instead, run Cell Ranger first.

### 2.3 Hardware

| Dataset size | RAM needed | Time estimate |
|---|---|---|
| 1,000-5,000 cells | 8 GB | 5-15 minutes |
| 5,000-20,000 cells | 16 GB | 15-45 minutes |
| 20,000-100,000 cells | 32-64 GB | 1-4 hours |

---

## 3. Pipeline Overview

### 3.1 Module Dependency Graph

Modules are connected by dependencies. When you request a module, all its upstream dependencies are automatically included:

```
cellranger --> qc --> doublet_detection --> clustering
                                              |
              +-------------------------------+-------------------------------+
              |               |               |               |               |
              v               v               v               v               v
        annotation    diff_expression    trajectory      cnv_inference    cell_cycle
              |               |               |               |
              v               v               v               v
      immune_pheno     pathway_analysis  pseudo_velocity   evolution
      tumor_microenv   validate_cbio
      cell_communic
```

### 3.2 Mandatory vs Optional

| Type | Modules | Why mandatory? |
|---|---|---|
| **Mandatory** (always run) | cellranger, qc, doublet_detection | Every analysis needs clean, de-duplicated data |
| **Optional** (you choose) | All 18 others | Pick what's relevant to your biological question |

### 3.3 How to Run

```bash
# Minimal run (just QC + clustering)
python -m workflow.modular.cli \
  --project my_analysis \
  --sample-root data/raw/my_dataset \
  --optional-modules clustering

# Full analysis (all 18 optional modules)
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,cell_cycle,differential_expression,annotation,\
trajectory,pseudo_velocity,cnv_inference,pathway_analysis,cell_communication,\
gene_regulatory_network,compare_10x,immune_phenotyping,tumor_microenvironment,\
gene_signature_scoring,evolution
```

---

## 4. Step-by-Step Module Guide

### Module 1: Cell Ranger Data Loading (`cellranger`)

**What it does:** Loads the 10X Cell Ranger count matrix into memory as an AnnData object — the standard data structure for scRNA-seq in Python.

**Method:** `scanpy.read_10x_mtx()` reads the sparse matrix files (barcodes, features, matrix).

**What you get:**
- An AnnData object with shape (n_cells x n_genes)
- In our LUSC example: 2,588 cells x 38,606 genes

**How to interpret:** This is just data loading — nothing to interpret yet. If this fails, check that your `filtered_feature_bc_matrix/` directory contains all three `.gz` files.

---

### Module 2: Quality Control (`qc`)

**What it does:** Removes low-quality cells and uninformative genes. Bad cells include:
- Empty droplets (too few genes detected)
- Dying cells (high mitochondrial gene %, because mRNA leaks from damaged cells)
- Multiplets/debris (abnormally high gene counts)

**Method:** Threshold-based filtering on QC metrics computed by scanpy.
- Reference: Luecken & Theis, *Molecular Systems Biology*, 2019. DOI: [10.15252/msb.20188746](https://doi.org/10.15252/msb.20188746)

**Key parameters and when to change them:**

| Parameter | Default | When to increase | When to decrease |
|---|---|---|---|
| `--min-genes` | 200 | Very noisy data | High-quality data with low gene detection |
| `--max-genes` | 7000 | Cell types with naturally high complexity | Suspected doublets |
| `--max-mito-pct` | 20% | Metabolically active tissues (heart, muscle) | Strict filtering |
| `--max-ribo-pct` | 50% | Immune cells (high ribo) | Non-immune tissues |

**Output files:**
| File | What to look for |
|---|---|
| `qc_violin_pre_filter.png` | Distribution of genes/cell, UMI/cell, mito%, ribo% BEFORE filtering |
| `qc_violin_post_filter.png` | Same metrics AFTER filtering — distributions should be cleaner |
| `qc_scatter_pre_filter.png` | Genes vs UMI counts, colored by mito% |
| `qc_scatter_post_filter.png` | Same after filtering |

**How to interpret:**
- **Good:** Post-filter violin plots show tight, unimodal distributions. Mito% should be low (< 10% for most cells).
- **Bad:** Bimodal distributions in post-filter plots suggest your thresholds are too loose. Very few cells remaining suggests thresholds are too strict.

**LUSC example:** 2,588 raw cells -> 2,382 after QC (206 removed, 8.0%).

---

### Module 3: Doublet Detection (`doublet_detection`)

**What it does:** Identifies "doublets" — droplets that accidentally captured two cells instead of one. Doublets appear as hybrid cell types and can create false clusters.

**Method:** Scrublet (Wolock et al., *Cell Systems*, 2019. DOI: [10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005))

Scrublet works by:
1. Simulating artificial doublets by averaging random pairs of cells
2. Building a shared k-NN graph of real + simulated cells
3. Scoring each real cell by how similar it is to simulated doublets
4. Automatically finding a threshold to call doublets

**Output:** `doublet_scores.png` — histogram showing the distribution of doublet scores. There should be a clear bimodal distribution with a small peak of doublets on the right.

**How to interpret:**
- **Good:** Clear separation between singlet peak (left) and doublet peak (right). Low doublet rate (< 10%).
- **Bad:** No clear threshold — the two peaks overlap heavily. Consider adjusting `--expected-doublet-rate`.

**LUSC example:** 5 doublets detected (0.2%), final cell count: 2,377.

---

### Module 4: Clustering (`clustering`)

**What it does:** This is the core analysis step that groups cells with similar expression profiles into clusters. It involves:

1. **Normalization** — Makes cell-to-cell comparisons fair (CPM + log transform)
2. **Highly Variable Gene (HVG) selection** — Finds genes that vary most across cells (3,000 by default)
3. **PCA** — Reduces 20,000+ genes to 40 principal components
4. **Neighbor graph** — Builds a k-nearest-neighbor graph (k=15) in PCA space
5. **UMAP** — Projects the high-dimensional data to 2D for visualization
6. **Leiden clustering** — Finds communities (clusters) in the neighbor graph

**Methods and citations:**
- UMAP: McInnes et al., *JOSS*, 2018. DOI: [10.21105/joss.00861](https://doi.org/10.21105/joss.00861)
- Leiden: Traag et al., *Scientific Reports*, 2019. DOI: [10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)
- Seurat HVG: Stuart et al., *Cell*, 2019. DOI: [10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031)

**Key parameters:**

| Parameter | Default | Effect of increasing | Effect of decreasing |
|---|---|---|---|
| `--leiden-resolution` | 0.8 | More clusters (finer subtypes) | Fewer clusters (broader groups) |
| `--n-top-genes` | 3000 | More genes in analysis | Fewer, more variable genes |
| `--n-pcs` | 40 | Captures more variance | Faster, less noise |
| `--n-neighbors` | 15 | Smoother clusters | More local structure preserved |

**Output files:**
| File | What it shows |
|---|---|
| `pca_variance_explained.png` | Elbow plot — how many PCs to keep. Look for the "elbow" where the curve flattens |
| `umap_leiden.png` | UMAP colored by cluster. This is the most important overview figure |

**How to interpret:**
- **Good UMAP:** Distinct, well-separated clusters. Each cluster is a potential cell type.
- **Bad UMAP:** One big blob (resolution too low), or 50+ tiny clusters (resolution too high), or fragmented clouds (poor normalization).

**LUSC example:** 17 clusters found at resolution 0.8. PCA: 40 PCs capture 90% of variance.

---

### Module 5: Cell Type Annotation (`annotation`)

**What it does:** Assigns biological names (T cell, Macrophage, Tumor cell, etc.) to each cluster by scoring cells against known marker gene panels.

**Method:** Scanpy gene set scoring (`sc.tl.score_genes`) with marker panels for 10 cell types. Each cell gets a score for each cell type; the highest score wins. Confidence = max_score - second_best_score.

**Reference:** Tirosh et al., *Science*, 2016. DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)

**Built-in cell type markers:**

| Cell type | Key markers |
|---|---|
| Tumor epithelial | EPCAM, KRT7, KRT8, KRT18, KRT19, MUC1 |
| T cell | CD3D, CD3E, CD4, CD8A, IL7R |
| NK cell | NKG7, GNLY, KLRD1 |
| B cell | MS4A1, CD79A, CD19 |
| Myeloid/Macrophage | CD68, CD14, LYZ |
| Fibroblast | DCN, LUM, COL1A1 |
| Endothelial | PECAM1, VWF |
| Plasma cell | MZB1, JCHAIN |
| Mast cell | TPSAB1, TPSB2 |
| Dendritic cell | FCER1A, CD1C |

**Output files:**
| File | What it shows |
|---|---|
| `umap_cell_type.png` | UMAP colored by assigned cell type |
| `umap_annotation_confidence.png` | UMAP colored by confidence score (brighter = more confident) |
| `cell_type_composition.png` | Stacked bar chart: what fraction of each cluster is each cell type |
| `cell_type_annotation.csv` | Per-cell annotation with confidence scores |

**How to interpret:**
- **Good:** Each cluster is dominated by one cell type. Confidence scores are high (> 0.2). Few "Unknown" cells (< 10%).
- **Bad:** Many "Unknown" cells, or one cell type assigned to clusters that look very different on UMAP. This means your marker panels don't match your tissue type — consider custom markers via `--markers-json`.

**LUSC example:** 10 cell types identified, 4.4% Unknown. T cells dominate (5 clusters), with macrophages, B cells, tumor epithelial, and smaller populations.

---

### Module 6: Cell Cycle Scoring (`cell_cycle`)

**What it does:** Assigns each cell a cell cycle phase (G1, S, or G2M) and continuous S/G2M scores. This helps you determine whether clusters are driven by biology or just by proliferation state.

**Method:** Gene set scoring using 47 S-phase genes and 49 G2M-phase genes from the Tirosh/Regev lab.
- Reference: Tirosh et al., *Science*, 2016. DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)

**Output:** `cell_cycle_umap.png` — UMAP colored by cell cycle phase.

**How to interpret:**
- If a cluster is entirely S/G2M cells, it may be a "proliferating" subpopulation rather than a distinct cell type.
- Use `--regress-cell-cycle` if cell cycle effects dominate your UMAP.

**LUSC example:** G1: 1,683 cells (71%), S: 280 (12%), G2M: 414 (17%).

---

### Module 7: CNV Inference (`cnv_inference`)

**What it does:** Infers copy number variations (gains/losses of chromosomal regions) from gene expression patterns. This is critical for distinguishing **malignant** cells from **normal** cells in tumor samples.

**Method:** Sliding-window smoothing of gene expression along chromosomal positions. Gene positions are auto-fetched from Ensembl BioMart. Per-cell CNV score = variance of the smoothed signal.
- Reference: Patel et al., *Science*, 2014. DOI: [10.1126/science.1254257](https://doi.org/10.1126/science.1254257)

**Output files:**
| File | What it shows |
|---|---|
| `cnv_heatmap.png` | Rows = cells, columns = genes ordered by chromosome. Red/blue = gains/losses |
| `cnv_score_umap.png` | UMAP colored by CNV score. Bright = likely malignant |
| `cnv_scores.csv` | Per-cell CNV score |
| `cnv_classification.json` | Malignant vs normal cell counts |

**How to interpret:**
- **cnv_heatmap.png:** Look for horizontal bands of red (gains) or blue (losses) in specific chromosomal regions. Malignant cells often show clear chromosome-arm level events (e.g., 3q gain in LUSC).
- **cnv_score_umap.png:** Malignant cells should cluster together with high CNV scores.

**LUSC example:** 594 cells classified as malignant (25%). Mean CNV score: 0.0025.

---

### Module 8: 10X Validation (`compare_10x`)

**What it does:** Compares your pipeline's results against Cell Ranger's built-in analysis to quantify agreement.

**Metrics:**
- **ARI** (Adjusted Rand Index) — measures cluster assignment agreement (0 = random, 1 = perfect). DOI: [10.1007/BF01908075](https://doi.org/10.1007/BF01908075)
- **NMI** (Normalized Mutual Information) — measures information shared between clusterings (0 = none, 1 = identical)

**How to interpret:**
- ARI > 0.5: Good agreement
- ARI 0.3-0.5: Moderate (different resolution or method choices)
- ARI < 0.3: Poor — investigate why

**LUSC example:** ARI = 0.42, NMI = 0.59. Moderate agreement — expected because we use Leiden (not Cell Ranger's graph-based clustering) at different resolution.

---

### Module 9: Differential Expression (`differential_expression`)

**What it does:** Finds "marker genes" — genes that are significantly higher or lower in one cluster compared to all others. These markers define what makes each cluster biologically distinct.

**Method:** Wilcoxon rank-sum test (non-parametric, no distribution assumptions) with Benjamini-Hochberg FDR correction.
- References: Wilcoxon, *Biometrics Bulletin*, 1945. DOI: [10.2307/3001968](https://doi.org/10.2307/3001968); Benjamini & Hochberg, *JRSS-B*, 1995. DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

**Significance filters:** adjusted p-value < 0.05 AND |log2 fold change| > 0.25.

**Output files:**
| File | What it shows |
|---|---|
| `marker_genes.csv` | All significant markers (filtered) |
| `marker_top5_by_cluster.csv` | Top 5 markers per cluster — your "cheat sheet" for naming clusters |
| `de_dotplot_top5.png` | Dot plot: size = fraction of cells expressing, color = mean expression |
| `de_heatmap_top5.png` | Heatmap of top markers across clusters |
| `de_volcano.png` | Volcano plot: log2FC vs -log10(p-value) |

**How to interpret:**
- **de_dotplot_top5.png** is the most informative. Each row is a gene, each column is a cluster. Big, dark dots = strong, consistent markers.
- Compare your top markers to known biology. For example, if cluster 5 shows high EPCAM, KRT7, KRT18 — that's an epithelial/tumor cluster.

**LUSC example:** 6,551 significant marker genes found across 17 clusters.

---

### Module 10: Gene Regulatory Network (`gene_regulatory_network`)

**What it does:** Infers which transcription factors (TFs) are active in each cluster. TFs are master regulators that control downstream gene programs.

**Method (primary):** decoupler + DoRothEA. DoRothEA provides a curated database of TF-target gene relationships (confidence levels A/B/C). decoupler uses Univariate Linear Models to infer TF activity from target gene expression.
- References: Garcia-Alonso et al., *Genome Research*, 2019. DOI: [10.1101/gr.240663.118](https://doi.org/10.1101/gr.240663.118); Badia-i-Mompel et al., *Bioinformatics Advances*, 2022. DOI: [10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016)

**Output files:**
| File | What it shows |
|---|---|
| `tf_activity_heatmap.png` | TF activity per cluster — which regulators drive each cluster? |
| `tf_activity_per_cluster.csv` | Numeric TF activity scores |

**How to interpret:** High TF activity in a cluster suggests that TF's program is active there. Example: high STAT1 in immune clusters (interferon signaling), high MYC in proliferating tumor cells.

---

### Module 11: Gene Signature Scoring (`gene_signature_scoring`)

**What it does:** Scores each cell against 10 cancer hallmark gene signatures. This reveals which biological processes are active where.

**Built-in signatures:**

| Signature | What it measures | Key reference |
|---|---|---|
| Proliferation | Active cell division (MKI67, TOP2A) | Standard oncology panel |
| Apoptosis resistance | Resistance to cell death (BCL2 family) | BCL2 family biology |
| Angiogenesis | Blood vessel formation (VEGF pathway) | VEGF/FLT pathway |
| EMT mesenchymal | Epithelial-to-mesenchymal transition | Tan et al., *EMBO Mol Med*, 2014. DOI: [10.15252/emmm.201404208](https://doi.org/10.15252/emmm.201404208) |
| Stemness | Cancer stem cell features | Malta et al., *Cell*, 2018. DOI: [10.1016/j.cell.2018.03.034](https://doi.org/10.1016/j.cell.2018.03.034) |
| Hypoxia | Low oxygen response | Buffa et al., *Br J Cancer*, 2010. DOI: [10.1038/sj.bjc.6605450](https://doi.org/10.1038/sj.bjc.6605450) |
| Glycolysis | Warburg effect (tumor metabolism) | Warburg effect biology |
| DDR | DNA damage response | DDR pathway |

**Output files:**
| File | What it shows |
|---|---|
| `signature_heatmap.png` | Signature scores per cluster |
| `signature_umap.png` | Top 4 most variable signatures on UMAP |
| `signature_correlation.png` | Which signatures co-occur? (Pearson correlation matrix) |

**How to interpret:**
- Tumor clusters should score high on proliferation, glycolysis, EMT
- Immune clusters should score low on tumor-associated signatures
- High EMT + low epithelial = mesenchymal phenotype (more invasive)

---

### Module 12: Trajectory / Pseudotime (`trajectory`)

**What it does:** Orders cells along a continuous "pseudotime" axis that represents biological progression — for example, from stem-like to differentiated, or from naive to exhausted T cells. Also builds a PAGA graph showing how clusters connect.

**Methods:**
- PAGA: Wolf et al., *Genome Biology*, 2019. DOI: [10.1186/s13059-019-1663-x](https://doi.org/10.1186/s13059-019-1663-x)
- DPT: Haghverdi et al., *Nature Methods*, 2016. DOI: [10.1038/nmeth.3971](https://doi.org/10.1038/nmeth.3971)

**Output files:**
| File | What it shows |
|---|---|
| `pseudotime_dpt_umap.png` | UMAP colored by pseudotime (yellow = early, purple = late) + diffusion map |
| `paga_trajectory.png` | PAGA graph: nodes = clusters, edges = connectivity strength |
| `pseudotime_gene_heatmap.png` | Top 30 genes correlated with pseudotime, ordered by time |
| `pseudotime_violin_per_cluster.png` | Pseudotime distribution per cluster |
| `pseudotime_per_cluster.csv` | Per-cluster pseudotime statistics |
| `pseudotime_top_genes.csv` | Genes most correlated with pseudotime |

**How to interpret:**
- **pseudotime_gene_heatmap.png** is the key figure. It shows how gene expression changes along the trajectory. Look for waves of activation/repression — these represent biological programs turning on/off.
- **paga_trajectory.png** shows which clusters are connected. Thick edges = strong transitions. This reveals the "path" cells take.
- **pseudotime_violin_per_cluster.png** shows where each cluster falls on the timeline. Clusters with early pseudotime are "source" populations; late pseudotime are "destination" populations.

**LUSC example:** Cluster 0 has the lowest mean pseudotime (0.02) — likely the starting population. Cluster 12 has the highest (0.85) — the most differentiated.

---

### Module 13: Cell Communication (`cell_communication`)

**What it does:** Identifies ligand-receptor interactions between cell types — which cell types are "talking" to which, and through what signaling pathways.

**Method (primary):** LIANA multi-method consensus (Dimitrov et al., *Nature Communications*, 2022. DOI: [10.1038/s41467-022-30755-0](https://doi.org/10.1038/s41467-022-30755-0))

**Fallback:** Manual scoring of 18 curated TME ligand-receptor pairs (PD-L1/PD-1, VEGFA/KDR, TGFB1/TGFBR2, etc.)

**Output files:**
| File | What it shows |
|---|---|
| `cell_communication_lr.csv` | All scored L-R interactions with source/target cell types |
| `cell_communication_heatmap.png` | Top L-R interactions as a heatmap |

**How to interpret:**
- Look for **immune checkpoint** interactions: PD-L1 (tumor) -> PD-1 (T cell) indicates immune evasion
- **VEGFA -> KDR** between tumor and endothelial cells indicates angiogenesis signaling
- High TGFB1 interactions suggest immunosuppressive microenvironment

**LUSC example:** 50 significant L-R interaction pairs identified.

---

### Module 14: Immune Phenotyping (`immune_phenotyping`)

**What it does:** Assigns fine-grained immune subtypes to immune cells (15 subtypes) and computes functional scores for exhaustion, cytotoxicity, and activation.

**References:**
- Zheng et al., *Cell*, 2017. DOI: [10.1016/j.cell.2017.05.035](https://doi.org/10.1016/j.cell.2017.05.035)
- Zhang et al., *Nature*, 2018. DOI: [10.1038/s41586-018-0694-x](https://doi.org/10.1038/s41586-018-0694-x)

**15 immune subtypes identified:**

| Category | Subtypes |
|---|---|
| CD4 T cells | CD4 naive, CD4 memory, Treg, Th1, Th2, Th17 |
| CD8 T cells | CD8 effector, CD8 memory, CD8 exhausted |
| NK cells | NK cytotoxic |
| Macrophages | M1 (pro-inflammatory), M2 (anti-inflammatory) |
| Dendritic cells | cDC1, cDC2, pDC |

**Functional scores:**
- **Exhaustion** (LAG3, PDCD1, TIGIT, TOX, ...) — are T cells "tired"?
- **Cytotoxicity** (GZMB, PRF1, NKG7, ...) — are T cells actively killing?
- **Activation** (CD69, CD38, ICOS, ...) — are T cells recently activated?

**Output files:**
| File | What it shows |
|---|---|
| `umap_immune_subtype.png` | UMAP colored by immune subtype |
| `immune_subtype_composition.png` | Stacked bar: immune subtype proportions per cluster |
| `immune_signature_heatmap.png` | Exhaustion/cytotoxicity/activation scores per subtype |
| `immune_exhaustion_umap.png` | UMAP colored by exhaustion score |
| `immune_cytotoxicity_umap.png` | UMAP colored by cytotoxicity score |

**How to interpret:**
- High CD8 exhausted + high exhaustion score = T cells are being suppressed by the tumor
- High M2/low M1 ratio = immunosuppressive macrophage polarization
- High Treg proportion = active immune suppression

**LUSC example:** CD4 memory (461 cells), CD8 effector (303), M2 macrophages (298), CD8 exhausted (121), Tregs (84).

---

### Module 15: Tumor Microenvironment (`tumor_microenvironment`)

**What it does:** Computes established immunotherapy-relevant scores that predict response to immune checkpoint inhibitors.

**Published signatures:**

| Score | What it predicts | Reference |
|---|---|---|
| **CYT** (Cytolytic Activity) | Active immune killing in the tumor | Rooney et al., *Cell*, 2015. DOI: [10.1016/j.cell.2014.12.033](https://doi.org/10.1016/j.cell.2014.12.033) |
| **TIS** (T-cell Inflamed, 18-gene) | Response to anti-PD-1 therapy | Ayers et al., *JCI*, 2017. DOI: [10.1172/JCI91190](https://doi.org/10.1172/JCI91190) |
| **IFN-gamma** (10-gene) | Interferon gamma signaling | Ayers et al., *JCI*, 2017 |
| **Immune ESTIMATE** | Immune cell infiltration level | Yoshihara et al., *Nat Commun*, 2013. DOI: [10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612) |
| **Stromal ESTIMATE** | Stromal/fibroblast content | Yoshihara et al., 2013 |

**Checkpoint molecules profiled:** PD-1, PD-L1, CTLA-4, LAG-3, TIM-3, TIGIT, VISTA, IDO1, B7-H3

**Output files:**
| File | What it shows |
|---|---|
| `tme_signature_heatmap.png` | TME scores per cluster |
| `checkpoint_dotplot.png` | Checkpoint molecule expression per cell type (size = % expressing, color = mean level) |
| `tme_cyt_umap.png` | CYT score on UMAP |
| `tme_immune_stromal_bar.png` | Immune vs stromal ESTIMATE per cluster |

**How to interpret:**
- **checkpoint_dotplot.png** is clinically actionable. If PD-L1 is high on tumor cells and PD-1 is high on T cells, the patient may benefit from anti-PD-1/PD-L1 therapy.
- High TIS score = "inflamed" tumor (better immunotherapy response)
- High stromal + low immune ESTIMATE = "cold" tumor (poor immunotherapy candidate)

**LUSC example:** Mean CYT = 0.15. 6 TME signatures scored across all clusters.

---

### Module 16: Pathway Analysis (`pathway_analysis`)

**What it does:** Identifies which biological pathways are enriched in the differentially expressed genes.

**Method (primary):** Over-representation analysis with MSigDB Hallmark gene sets via gseapy.
- References: Subramanian et al., *PNAS*, 2005. DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102); Liberzon et al., *Cell Systems*, 2015. DOI: [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)

**Output files:**
| File | What it shows |
|---|---|
| `pathway_enrichment.csv` | All enriched pathways with p-values and overlapping genes |
| `pathway_enrichment_bar.png` | Bar chart of top enriched pathways |

**How to interpret:** Look for pathways consistent with your tissue type. In LUSC: EMT, interferon-gamma response, angiogenesis, TNF-alpha signaling, MYC targets, and hypoxia are all expected in a lung tumor.

---

### Module 17: Pseudo-Velocity (`pseudo_velocity`)

**What it does:** Computes velocity vectors showing the *direction* cells are moving in gene expression space, based on pseudotime gradients. This produces flow arrows and stream plots on UMAP.

**Method:** For each cell, calculate the average direction to k-nearest neighbors weighted by their pseudotime difference. Speed = magnitude of the velocity vector.

**Output files:**
| File | What it shows |
|---|---|
| `pseudo_velocity_arrows.png` | Quiver plot: arrows show direction of cell state transitions |
| `pseudo_velocity_stream.png` | Streamline plot: continuous flow field interpolated on a grid |
| `pseudo_velocity_speed_umap.png` | UMAP colored by velocity speed (fast-changing vs stable cells) |
| `pseudo_velocity_speed_boxplot.png` | Speed distribution per cluster |

**How to interpret:**
- **Arrows** should flow from early (stem/naive) to late (differentiated/exhausted) populations
- **High speed** cells are actively transitioning between states
- **Low speed** cells are in stable terminal states
- Arrows pointing in multiple directions suggest a branching decision point

---

### Module 18: Tumor Evolution (`evolution`)

**What it does:** Reconstructs the clonal architecture of the tumor by clustering cells based on their CNV profiles, then orders clones along pseudotime to reveal the evolutionary trajectory.

**Method:** Hierarchical clustering (Ward's method) on CNV profiles to identify clones. Clone-specific marker genes via Wilcoxon test. Phylogenetic dendrogram from clone centroid distances.
- Reference: Gao et al., *Nature Biotechnology*, 2021. DOI: [10.1038/s41587-020-00795-2](https://doi.org/10.1038/s41587-020-00795-2)

**Output files:**
| File | What it shows |
|---|---|
| `evolution_clone_umap.png` | UMAP colored by clone assignment |
| `evolution_phylo_dendrogram.png` | Phylogenetic tree of clones based on CNV similarity |
| `evolution_timeline.png` | Pseudotime distribution per clone (violin + density) — reveals temporal ordering |
| `evolution_clone_composition.png` | Clone proportions per Leiden cluster |
| `evolution_cnv_by_clone.png` | CNV score distribution per clone |
| `evolution_clone_stats.csv` | Per-clone summary: size, CNV score, pseudotime, dominant cell type |
| `evolution_clone_markers.csv` | Clone-specific differentially expressed genes |

**How to interpret:**
- **evolution_phylo_dendrogram.png** is the key figure. It shows how clones are related — like a family tree of cancer. Closely branched clones share more CNV events.
- **evolution_timeline.png** reveals temporal ordering. If Clone_1 has early pseudotime and Clone_3 has late pseudotime, Clone_3 likely evolved *from* Clone_1.
- **evolution_clone_markers.csv** identifies genes that distinguish clones — potential drivers of clonal evolution.

**LUSC example:** 3 clones identified. Clone_1 (54% of cells, early pseudotime, T cell dominant), Clone_2 (20%, late pseudotime, Myeloid dominant), Clone_3 (26%, early pseudotime, T cell dominant).

---

## 5. Running the Pipeline

### 5.1 Common Commands

```bash
# Activate environment
conda activate sc10x

# Full analysis (recommended for first run)
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,cell_cycle,differential_expression,annotation,\
trajectory,pseudo_velocity,cnv_inference,pathway_analysis,cell_communication,\
gene_regulatory_network,compare_10x,immune_phenotyping,tumor_microenvironment,\
gene_signature_scoring,evolution

# Quick exploratory run (just clustering + DE)
python -m workflow.modular.cli \
  --project quick_look \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression

# Immuno-oncology focused
python -m workflow.modular.cli \
  --project immuno_deep \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,annotation,differential_expression,\
immune_phenotyping,tumor_microenvironment,cell_communication,pathway_analysis

# With crash recovery
python -m workflow.modular.cli \
  --project my_run \
  --sample-root data/raw/my_dataset \
  --optional-modules clustering,annotation \
  --checkpoint

# Resume after a crash
python -m workflow.modular.cli \
  --project my_run \
  --sample-root data/raw/my_dataset \
  --optional-modules clustering,annotation \
  --checkpoint --resume-from annotation
```

### 5.2 Adjusting Parameters

```bash
# Finer clustering (more clusters)
python -m workflow.modular.cli ... --leiden-resolution 1.5

# Stricter QC
python -m workflow.modular.cli ... --max-mito-pct 10 --min-genes 500

# Custom cell type markers
python -m workflow.modular.cli ... --markers-json my_markers.json

# Use normal fibroblasts as CNV reference
python -m workflow.modular.cli ... --cnv-reference-group Fibroblast
```

---

## 6. Interpreting Results

### 6.1 First Things to Check

After a run completes, follow this order:

| Step | File | What to verify |
|---|---|---|
| 1 | `module_status.csv` | All modules show "ok" |
| 2 | `qc/qc_violin_post_filter.png` | Distributions are clean, no bimodality |
| 3 | `clustering/umap_leiden.png` | Clusters are well-separated |
| 4 | `annotation/umap_cell_type.png` | Cell types make biological sense |
| 5 | `differential_expression/de_dotplot_top5.png` | Marker genes match expected biology |
| 6 | `trajectory/paga_trajectory.png` | Cluster connectivity is biologically plausible |
| 7 | `evolution/evolution_phylo_dendrogram.png` | Clone relationships make sense with CNV data |

### 6.2 Red Flags

| Symptom | Likely cause | Fix |
|---|---|---|
| > 50% cells removed by QC | Thresholds too strict, or bad sample quality | Relax thresholds; check sample prep |
| One giant cluster + many tiny ones | Leiden resolution wrong | Try `--leiden-resolution 0.5` or `1.2` |
| > 20% "Unknown" cell types | Marker panels don't match tissue | Use `--markers-json` with tissue-specific markers |
| All cells labeled same type | Too few marker genes present, or normalization issue | Check `--n-top-genes`, try `--scale-data` |
| CNV inference fails | No gene position data | Auto-fetched from Ensembl; check internet connection |
| UMAP is a single blob | Insufficient biological variation, or batch effects | Try `--batch-method harmony` if multi-sample |

### 6.3 Output File Types

| Extension | What it is | How to open |
|---|---|---|
| `.png` | Figure / plot | Any image viewer |
| `.csv` | Table (comma-separated) | Excel, R, Python (pandas) |
| `.json` | Structured metadata | Text editor, Python (json) |
| `.h5ad` | AnnData object (all data + results) | Python: `import anndata; adata = anndata.read_h5ad("final_adata.h5ad")` |

---

## 7. FAQ / Troubleshooting

**Q: How long does the full pipeline take?**
A: For 3,000 cells with all 18 modules: approximately 5-10 minutes. Larger datasets scale roughly linearly.

**Q: Can I add my own gene signatures?**
A: Yes. Create a JSON file: `{"my_signature": ["GENE1", "GENE2", "GENE3"]}` and pass `--signature-json my_sigs.json`.

**Q: What if I don't have a loom file for RNA velocity?**
A: The `rna_velocity` module requires spliced/unspliced layers from velocyto. Without a loom file, it will skip gracefully. Use `pseudo_velocity` instead — it works from pseudotime alone.

**Q: Can I run on mouse data?**
A: The pipeline is designed for human data (human gene names, human marker panels, Ensembl human gene positions). For mouse, you would need to customize marker panels and change the BioMart organism.

**Q: How do I compare two conditions (e.g., treated vs untreated)?**
A: Load both samples with a shared batch key column. Use `--batch-method harmony --batch-key condition` to integrate, then look for condition-specific clusters and DE genes.

**Q: What does "checkpoint" do?**
A: `--checkpoint` saves the AnnData object after each module completes. If the pipeline crashes at module 15, you can `--resume-from` module 15 instead of re-running everything from scratch.

**Q: How do I cite this pipeline?**
A: Cite the methods used by each module. The README contains a complete citation list with 33 DOIs. At minimum, cite scanpy (Wolf et al., Genome Biology, 2018), plus the specific methods you used (e.g., Leiden, Scrublet, DPT, etc.).

---

## 8. Complete Output Inventory

A full pipeline run produces this structure:

```
results/{project}_{timestamp}/
├── final_adata.h5ad              # Complete AnnData (load in Python for further analysis)
├── run_manifest.json             # Run parameters and metadata
├── module_status.csv             # Pass/fail status for each module
│
├── cellranger/                   # [M1] Data loading
├── qc/                           # [M2] 4 QC plots
├── doublet_detection/            # [M3] Doublet score histogram
├── clustering/                   # [M4] PCA elbow + UMAP
├── annotation/                   # [M5] Cell type UMAP + composition + CSVs
├── cell_cycle/                   # [M6] Cell cycle UMAP + scores CSV
├── cnv_inference/                # [M7] CNV heatmap + UMAP + classification
├── compare_10x/                  # [M8] ARI/NMI metrics
├── differential_expression/      # [M9] Volcano + dotplot + heatmap + CSVs
├── gene_regulatory_network/      # [M10] TF activity heatmap + CSVs
├── gene_signature_scoring/       # [M11] Signature heatmap + UMAP + correlation
├── trajectory/                   # [M12] PAGA + pseudotime heatmap + violin + CSVs
├── cell_communication/           # [M13] L-R heatmap + CSV
├── immune_phenotyping/           # [M14] Immune UMAP + composition + signatures
├── tumor_microenvironment/       # [M15] TME heatmap + checkpoint dotplot + CSVs
├── pathway_analysis/             # [M16] Enrichment bar + CSV
├── pseudo_velocity/              # [M17] Arrows + stream + speed UMAP + boxplot
└── evolution/                    # [M18] Clone UMAP + dendrogram + timeline + CSVs
```

**Total: 77 files across 18 module folders.**
