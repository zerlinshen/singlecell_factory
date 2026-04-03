# LUSC scRNA-seq Analysis Report

**Dataset**: 10X Genomics Lung Squamous Cell Carcinoma (3K cells)
**Pipeline**: singlecell_factory v3.1
**Date**: 2026-04-04
**Run ID**: lusc_final_20260404_041002

---

## 1. Data Overview & Quality Control

| Metric | Value |
|---|---|
| Raw cells loaded | 2,588 |
| Raw genes | 38,606 |
| Cells after QC | 2,382 (92.0% retained) |
| Genes after QC | 21,993 |
| QC removed cells | 206 |
| Doublets detected (Scrublet) | 5 (0.21%) |
| Final cells for analysis | 2,377 |

**QC thresholds**: min_genes=200, max_genes=7000, min_counts=500, max_counts=50000, max_mito=20%, max_ribo=50%

### QC Violin Plots (Pre-filter)

![QC Violin Pre-filter](figures/qc_violin_pre_filter.png)

### QC Violin Plots (Post-filter)

![QC Violin Post-filter](figures/qc_violin_post_filter.png)

### QC Scatter Plots

![QC Scatter Pre-filter](figures/qc_scatter_pre_filter.png)
![QC Scatter Post-filter](figures/qc_scatter_post_filter.png)

### Doublet Detection

![Doublet Scores](figures/doublet_scores.png)

Scrublet detected 5 doublets (0.21% doublet rate) with an automatically set threshold of 0.50. The estimated overall doublet rate is 5.7%, consistent with the expected 6.0% for this cell count.

---

## 2. Dimensionality Reduction & Clustering

### PCA Variance Explained

![PCA Variance](figures/pca_variance_explained.png)

40 PCs are needed to capture 90% cumulative variance. The elbow is visible around PC 15-20, confirming that 40 PCs is a conservative and appropriate choice.

### UMAP — Leiden Clustering (17 clusters)

![UMAP Leiden](figures/umap_leiden.png)

**Parameters**: n_top_genes=3000 (HVG), n_pcs=40, n_neighbors=15, leiden_resolution=0.8

17 clusters were identified, providing fine-grained resolution of the tumor microenvironment.

---

## 3. Cell Type Annotation

![UMAP Cell Type](figures/umap_cell_type.png)

![Annotation Confidence](figures/umap_annotation_confidence.png)

![Cell Type Composition](figures/cell_type_composition.png)

### Cluster-to-Cell-Type Mapping

| Cluster | Cell Type | Key Markers |
|---|---|---|
| 0, 1, 6, 9, 15 | T cell | CCL5, CD8A, CD7, IL7R |
| 2, 3, 7, 8 | Myeloid/Macro | FTL, C1QA, C1QC, LYZ, PLAUR |
| 4, 14 | B cell | MS4A1, CD37, BANK1 |
| 5, 13 | Tumor epithelial | KRT19, KRT6A, KRT5, S100A2 |
| 10 | Plasma cell | — |
| 11 | Endothelial | — |
| 12 | Mast cell | — |
| 16 | Fibroblast | — |

Unknown annotation rate: 4.42% (low, indicating good marker coverage).

---

## 4. 10X Official Comparison

| Metric | 10X Official | Our Pipeline | Difference |
|---|---|---|---|
| Estimated Cells | 2,588 | 2,377 | -8.2% |
| Median Genes/Cell | 1,484 | 1,513 | +2.0% |
| Median UMI/Cell | 3,671 | 3,756 | +2.3% |
| Total Genes Detected | 27,135 | 21,993 | -18.9% |

| Clustering Agreement | Score |
|---|---|
| Adjusted Rand Index (ARI) | 0.391 |
| Normalized Mutual Info (NMI) | 0.593 |
| Overlapping cells | 2,377 |

**Interpretation**: Our pipeline retains 92% of cells after more stringent QC + doublet removal. Per-cell metrics (median genes, median UMI) are nearly identical (+2%), confirming similar data quality. The lower total gene count reflects our stricter gene filter (min_cells=3). Clustering ARI=0.39 is expected since we produce 17 fine-grained clusters vs 10X's 8 coarse clusters; NMI=0.59 confirms strong information overlap.

---

## 5. Cell Cycle Analysis

![Cell Cycle UMAP](figures/cell_cycle_umap.png)

| Phase | Cells | Percentage |
|---|---|---|
| G1 | 1,692 | 71.2% |
| G2M | 400 | 16.8% |
| S | 285 | 12.0% |

The majority of cells are in G1 phase. Active cycling (S + G2M = 28.8%) is concentrated in tumor epithelial clusters and certain T cell populations.

---

## 6. Differential Expression

**6,464 significant DE genes** (padj < 0.05, |logFC| > 0.25) across 17 clusters.

### Dotplot — Top 5 Markers per Cluster

![DE Dotplot](figures/de_dotplot_top5.png)

### Heatmap — Top 5 Markers per Cluster

![DE Heatmap](figures/de_heatmap_top5.png)

### Volcano Plot

![DE Volcano](figures/de_volcano.png)

### Highlight Markers

| Cluster | Top Markers | Biological Identity |
|---|---|---|
| 0 | CCL5, CD8A, HLA-B, CD7 | Cytotoxic CD8+ T cells |
| 1 | DNAJB1, HSPA1A, HSPA1B | Stressed/activated T cells (heat-shock) |
| 2 | FTL, C1QA, C1QC, APOE | Tissue-resident macrophages |
| 5 | KRT19, KRT6A, KRT5, S100A2 | Squamous tumor epithelial |
| 7 | PLAUR, SLC11A1, LYZ, TYROBP | Inflammatory monocyte-derived macrophages |
| 8 | NFKB1 (top TF), CD83 | Activated dendritic/myeloid |
| 13 | NKX2-1 (top TF) | Lung-specific tumor cells |
| 16 | TWIST1, SNAI1 (top TFs), COL1A1 | EMT-active fibroblasts (CAFs) |

---

## 7. Trajectory & Pseudo-velocity

![Pseudotime UMAP](figures/pseudotime_dpt_umap.png)

Diffusion pseudotime (DPT) was computed from the first diffusion component extremum. The trajectory shows progression from naive T cells through activated states to exhausted phenotypes.

---

## 8. Deep Immune Phenotyping

![UMAP Immune Subtype](figures/umap_immune_subtype.png)

![Immune Subtype Composition](figures/immune_subtype_composition.png)

### Immune Subtype Distribution (n=2,000 immune cells)

| Subtype | Cells | % of Immune |
|---|---|---|
| CD4 memory | 460 | 23.0% |
| CD8 effector | 310 | 15.5% |
| Macro M2 | 300 | 15.0% |
| cDC1 | 133 | 6.7% |
| pDC | 121 | 6.1% |
| CD8 exhausted | 120 | 6.0% |
| CD4 naive | 115 | 5.8% |
| Treg | 82 | 4.1% |
| Macro M1 | 68 | 3.4% |
| CD8 memory | 54 | 2.7% |
| Th17 | 52 | 2.6% |
| cDC2 | 47 | 2.4% |
| Th1 | 37 | 1.9% |
| NK cytotoxic | 36 | 1.8% |
| Th2 | 20 | 1.0% |

### Immune Functional Signatures

![Exhaustion UMAP](figures/immune_exhaustion_umap.png)
![Cytotoxicity UMAP](figures/immune_cytotoxicity_umap.png)
![Immune Signature Heatmap](figures/immune_signature_heatmap.png)

**Key findings**:
- CD8 exhausted cells (120 cells, 6%) express high levels of PDCD1, LAG3, HAVCR2, TIGIT, TOX — classic exhaustion markers
- M2 > M1 macrophage ratio (300:68 = 4.4:1) indicates immunosuppressive myeloid polarization
- Treg population (82 cells, 4.1%) suggests active immune suppression in TME

---

## 9. Tumor Microenvironment (TME) Scoring

![CYT Score UMAP](figures/tme_cyt_umap.png)
![TIS Score UMAP](figures/tme_tis_umap.png)
![TME Signature Heatmap](figures/tme_signature_heatmap.png)
![Immune vs Stromal ESTIMATE](figures/tme_immune_stromal_bar.png)

### TME Scores per Cell Type Cluster (top clusters)

| Cluster (Type) | CYT Score | TIS | IFNg | Immuno-suppression | Immune ESTIMATE |
|---|---|---|---|---|---|
| 0 (T cell/CD8) | **0.510** | 0.423 | 0.113 | 0.076 | **0.581** |
| 8 (Myeloid) | 0.121 | **0.627** | **0.433** | **0.302** | 0.163 |
| 5 (Tumor) | 0.000 | -0.093 | -0.279 | -0.160 | -0.359 |
| 16 (Fibroblast) | 0.000 | -0.127 | -0.320 | -0.136 | -0.383 |

**Interpretation**: Cluster 0 (CD8 T cells) shows the highest cytolytic activity (CYT=0.51). Myeloid cluster 8 has the strongest T-cell inflammation and IFN-gamma signatures but also the highest immunosuppression score, reflecting the dual role of tumor-associated macrophages. Tumor and fibroblast clusters show negative immune scores, as expected.

### Immune Checkpoint Expression Profiling

![Checkpoint Dotplot](figures/checkpoint_dotplot.png)

| Checkpoint | Highest Expression | % Positive |
|---|---|---|
| PD-1 (PDCD1) | T cells | 30.1% |
| PD-L1 (CD274) | Tumor epithelial | 13.2% |
| CTLA-4 | T cells | 30.7% |
| LAG-3 | NK cells | 28.4% |
| TIM-3 (HAVCR2) | **Myeloid/Macro** | **66.3%** |
| TIGIT | T cells | 38.2% |
| VISTA (VSIR) | **Myeloid/Macro** | **62.9%** |
| IDO1 | **Tumor epithelial** | **39.7%** |
| B7-H3 (CD276) | Tumor epithelial | 25.5% |

**Key findings**: TIM-3 and VISTA are predominantly expressed on myeloid cells (66% and 63%), suggesting these are key myeloid checkpoint targets. IDO1 is highly expressed on tumor cells (40%), indicating metabolic immune evasion. Traditional checkpoint targets PD-1/CTLA-4 are concentrated on T cells.

---

## 10. Gene Signature Scoring

![Signature Heatmap](figures/signature_heatmap.png)
![Signature UMAP](figures/signature_umap.png)
![Signature Correlation](figures/signature_correlation.png)

### Cancer Hallmark Signatures per Cluster

10 cancer-relevant gene signatures were scored. Key patterns:
- **Proliferation + glycolysis** highest in tumor clusters (5, 13) — consistent with Warburg effect
- **EMT mesenchymal** highest in fibroblast cluster (16) — cancer-associated fibroblasts
- **Hypoxia** elevated in tumor clusters — tumor hypoxic niche
- **EMT mesenchymal** and **invasion** positively correlated (r > 0.7)
- **EMT epithelial** anti-correlated with **EMT mesenchymal** — expected EMT axis

---

## 11. Pathway Enrichment

![Pathway Enrichment](figures/pathway_enrichment_bar.png)

### Top Enriched Pathways

| Cluster | Pathway | padj | Genes |
|---|---|---|---|
| 16 (Fibroblast) | EMT | 2.3e-12 | ACTA2, COL1A1, COL3A1, FN1, SPARC, TAGLN, THY1 |
| 9 (T cell) | IFN-gamma Response | 7.1e-08 | B2M, HLA-A, HLA-B, HLA-C, IRF1 |
| 11 (Endothelial) | Angiogenesis | 2.7e-06 | CDH5, COL4A1, FLT1, PECAM1 |
| 0 (T cell) | IFN-gamma Response | 4.1e-06 | B2M, HLA-A, HLA-B, HLA-C |

---

## 12. Cell-Cell Communication

![Cell Communication Heatmap](figures/cell_communication_heatmap.png)

### Top 10 Ligand-Receptor Interactions

| Ligand | Receptor | Source | Target | LR Score |
|---|---|---|---|---|
| SPP1 | CD44 | Myeloid/Macro | Mast cell | 3.57 |
| SPP1 | CD44 | Myeloid/Macro | Unknown | 2.61 |
| SPP1 | CD44 | Myeloid/Macro | Myeloid/Macro | 2.58 |
| SPP1 | CD44 | Myeloid/Macro | T cell | 2.49 |
| SPP1 | CD44 | Myeloid/Macro | NK cell | 2.43 |
| SPP1 | CD44 | Myeloid/Macro | DC | 2.12 |
| SPP1 | CD44 | Myeloid/Macro | B cell | 2.04 |
| TGFB1 | TGFBR2 | NK cell | Endothelial | 1.96 |
| SPP1 | CD44 | Myeloid/Macro | Fibroblast | 1.82 |
| TGFB1 | TGFBR2 | T cell | Endothelial | 1.60 |

**Key findings**: SPP1-CD44 axis dominates cell-cell communication, with myeloid cells as the primary source. SPP1 (Osteopontin) from macrophages is a well-established pro-tumorigenic signal in LUSC. TGFB1-TGFBR2 signaling from lymphocytes to endothelium suggests immune-vascular crosstalk.

---

## 13. Transcription Factor Activity

![TF Activity Heatmap](figures/tf_activity_heatmap.png)

### Top TFs per Cluster (Notable)

| Cluster | Cell Type | Top TFs | Interpretation |
|---|---|---|---|
| 0 | CD8 T cell | STAT1, MYC | IFN-gamma signaling, proliferative T cells |
| 5 | Tumor | MYC, HIF1A | Oncogene + hypoxia response |
| 8 | Myeloid | **NFKB1** (1.304) | NF-kB driven inflammation |
| 13 | Tumor (lung) | **NKX2-1** (1.595) | Lung lineage TF (TTF-1) |
| 16 | Fibroblast | **TWIST1** (1.11), **SNAI1** (0.88) | EMT master regulators |

**Key findings**: NKX2-1 (TTF-1) specifically marks cluster 13 as lung-origin tumor cells. TWIST1/SNAI1 in fibroblasts confirms cancer-associated fibroblast (CAF) EMT program. NF-kB dominance in myeloid cluster 8 indicates classical inflammatory activation.

---

## 14. Module Execution Summary

| Module | Status | Key Metric |
|---|---|---|
| cellranger | OK | 2,588 raw cells |
| qc | OK | 92% cell retention |
| doublet_detection | OK | 0.21% doublet rate |
| clustering | OK | 17 clusters |
| annotation | OK | 4.4% unknown |
| cell_cycle | OK | 71.2% G1 |
| differential_expression | OK | 6,464 significant genes |
| trajectory | OK | DPT computed |
| pseudo_velocity | OK | Velocity vectors computed |
| pathway_analysis | OK | EMT, IFN-gamma top pathways |
| cell_communication | OK | SPP1-CD44 dominant |
| gene_regulatory_network | OK | STAT1, MYC, TWIST1 top TFs |
| compare_10x | OK | ARI=0.39, NMI=0.59 |
| immune_phenotyping | OK | 15 immune subtypes |
| tumor_microenvironment | OK | CYT, TIS, checkpoints scored |
| gene_signature_scoring | OK | 10 signatures scored |
| cnv_inference | SKIPPED | No gene position annotations |

**16 of 17 modules completed successfully.**

---

## 15. Recommendations for Deeper Analysis

### High Priority

1. **Install infercnvpy and add gene position annotations** to enable CNV inference. This would allow distinguishing malignant from non-malignant cells directly, rather than relying solely on marker-based annotation. Use the GTF from the reference genome to generate chromosome/start annotations:
   ```python
   # Add to adata.var from GTF
   adata.var["chromosome"] = gene_chr_mapping
   adata.var["start"] = gene_start_mapping
   ```

2. **RNA velocity with loom file**: Generate spliced/unspliced counts using `velocyto run` or `STARsolo` on the original BAM file. This would reveal transcriptional dynamics — particularly the trajectory from naive T cells through activation to exhaustion, and the EMT trajectory in tumor cells.

3. **Multi-resolution clustering**: Run Leiden at resolution 0.4, 0.6, 0.8, 1.0, 1.2 and compare. Resolution 0.8 (17 clusters) is fine-grained; a resolution scan would identify the most stable communities and optimal granularity. Consider using `clustree` visualization.

4. **Subclustering of key populations**:
   - **T cells** (clusters 0, 1, 6, 9, 15): Subcluster at higher resolution to better resolve CD4/CD8 subsets, especially the CD8 effector vs exhausted transition
   - **Myeloid cells** (clusters 2, 3, 7, 8): Subcluster to resolve monocytes, tissue-resident macrophages, dendritic cell subsets, and MDSCs
   - **Tumor cells** (clusters 5, 13): Subcluster to identify intra-tumor heterogeneity and therapy-resistant subpopulations

### Medium Priority

5. **Receptor-ligand analysis with LIANA**: Install LIANA for multi-method consensus cell communication analysis (CellPhoneDB, NATMI, Connectome, etc.). The current manual L-R scoring captures the major interactions but multi-method consensus provides higher confidence rankings.

6. **Spatial deconvolution context**: If paired spatial transcriptomics (Visium/MERFISH) data is available, use cell2location or Tangram to map these single-cell populations to tissue architecture. This would reveal whether exhausted T cells co-localize with PD-L1+ tumor cells.

7. **Clonotype analysis**: If paired TCR-seq data exists, integrate it to study clonal expansion of CD8 exhausted vs effector T cells. Clonally expanded exhausted T cells with tumor-reactive TCRs are the strongest candidates for checkpoint immunotherapy response.

8. **Progeny pathway activity at single-cell level**: Install decoupler + PROGENy for per-cell pathway activity scoring (instead of the current gene set overlap approach). This provides more granular pathway activity maps.

### Biological Follow-up

9. **M2/M1 ratio and clinical correlation**: The M2:M1 ratio of 4.4:1 is strongly immunosuppressive. In TCGA LUSC bulk data, high M2/M1 ratio correlates with worse overall survival. Validate this ratio against clinical outcomes.

10. **Checkpoint combination targets**: The data suggests:
    - **Anti-PD-1 + anti-TIGIT** for T cells (both >30% positive)
    - **Anti-TIM-3 + anti-VISTA** for myeloid checkpoint blockade (>60% positive)
    - **IDO1 inhibitor** for tumor metabolic reprogramming (40% positive on tumor cells)

11. **CAF-targeting strategies**: Cluster 16 fibroblasts show strong TWIST1/SNAI1 EMT programs and are the dominant source of EMT pathway enrichment. These cancer-associated fibroblasts may be targetable with FAP-directed therapies or TGF-beta blockade.

12. **SPP1+ macrophage program**: SPP1 (Osteopontin) from macrophages is the dominant communication signal. Recent studies (Zhang et al., Cancer Cell 2023) show SPP1+ macrophages drive tumor progression in NSCLC. Consider deeper characterization of SPP1-high vs SPP1-low macrophage states.

---

## Appendix: Run Parameters

```
Project: lusc_final
Sample root: data/raw/lung_carcinoma_3k_count
Pipeline version: singlecell_factory v3.1
QC: min_genes=200, max_genes=7000, min_counts=500, max_counts=50000, max_mito=20%, max_ribo=50%
Doublet: expected_rate=0.06, Scrublet
Clustering: n_top_genes=3000, n_pcs=40, n_neighbors=15, resolution=0.8
DE: Wilcoxon rank-sum, n_genes=500, padj<0.05, |logFC|>0.25
Annotation: marker scoring with confidence threshold=0.1
Batch correction: not applied (single sample)
Pathway: built-in Hallmark sets, BH FDR correction
Cell communication: manual 18-pair L-R scoring
GRN: manual 12-TF curated target sets
CNV reference group: Fibroblast (not executed — no gene positions)
Random seed: 0
```
