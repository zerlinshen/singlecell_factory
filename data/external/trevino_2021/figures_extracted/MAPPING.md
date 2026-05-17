# Trevino 2021 Figure Mapping (Stage C0)

**Source PDF:** `paper.pdf` (SHA `19dda596...456be`)
**Extraction:** `pdfimages -png` (717 embedded images, `img-NNN.png`) + `pdftoppm -r 200 -png` (41 page rasters, `page-NN.png`)
**Authored:** 2026-05-17 by ralph autonomous execution against plan v5.0.4 Stage C0
**Status:** **DRAFT — requires fresh-subagent review per AC-V5-TREVINO-PDF-2** before Stage C5 PDF assembly.

## Visual identification of figure pages (whole-page rasters preferred for comparison-PDF embedding)

Page rasters are more reliable than embedded-image extracts because:
- Embedded-image extraction (`pdfimages`) often splits composite figure panels into many sub-images (717 images for a 41-page paper) and loses text labels/legends/captions.
- Whole-page rasters (`pdftoppm`) preserve the complete published figure as the reader sees it.

| Cell paper Fig | Page raster | Plan v5.0.4 figure ID(s) | Visual content (verified by ralph 2026-05-17) | Status |
|---|---|---|---|---|
| Figure 1 | `page-04.png` | F-1B, F-1C, F-1D | A=schematic; B=IHC SOX9/CTIP2; C=IHC GFAP/Ki67; D=UMAP × Sample.Age; E=marker dot plot; F=scRNA-seq UMAP × cluster (38,500 cells); G=gene activity scores dot plot; H=ATAC dot plot | OK |
| Figure 2 | `page-05.png` | F-2A, F-2BD, F-5 (peak-gene linkage) | A=CCA schematic; B=scRNA+scATAC integration UMAP; C=normalized ATAC vs RNA heatmap (4,679 cells); D=GA-RNA correlation vs # linked enhancers; E=GO enrichment; F=multiome schematic; G=multiome RNA+ATAC UMAP; H=peak/gene cluster sizes; I=GPC singleome vs multiome scatter | OK |
| Figure 3 | `page-07.png` | F-3 (NOT-RUN — chromVAR TF motif) | A=pseudotime UMAP (mRNA velocity); B=mapped pseudotime; C=ATAC+RNA pseudotime heatmaps; D=gene/peak/GO enrichment; E=TF motif enrichment (peaks); F=TF motif activity heatmap (pseudotime); G=motif cluster heatmap; H=cluster enrichment; I=mean motif energy scatter | Reference only; F-3 is NOT-RUN placeholder |
| Figure 4 | `page-09.png` | F-4A, F-4D, F-4EG | A=schematic; B=fuzzy C-means clustering heatmap; C=heatmap; D=module gene UMAPs; E=connectivity dendrogram; F-K=UMAP+IHC panels (ASCL1, EOMES, HES1, AQP4, OLIG1, MBP, EGFR, PDGFRA, SPARCL1) | OK |
| Figure 5 | `page-11.png` | (supplementary — not in v5.0.4 reproducibility audit) | A=astrocyte module UMAPs; B=motif enrichment scatter; C=Louvain UMAP; D=differential expression; E=heterogeneity UMAP; F=fuzzy clusters | Reference only |
| Figure 6 | `page-12.png` | F-6 (NOT-RUN — TF-gene network) | A=cell-cycle signature UMAP; B=ATAC pseudobulk fuzzy clustering; C=projection chromatin → cycling; D=early-branch GPCs heatmap; E=motif enrichment; F=GPC dynamics fuzzy clusters; G=projection scRNA into ATAC | Reference only; F-6 is NOT-RUN placeholder |
| Figure 7 | `page-13.png` | F-7 (NOT-RUN — disease GWAS overlay) | A=model architecture; B=cluster ASD de novo enrichment UMAP; C=enrichment bar; D=ASD-mutated genes; E=motif family frequency; F=NFIX-NFIA/GluN4; G=MYF6-CTCF/Early RG | Reference only; F-7 is NOT-RUN placeholder |

## Plan v5.0.4 figure-ID → Trevino source mapping (binding for Stream C5 PDF assembly)

| Plan figure ID | Trevino source | Notes |
|---|---|---|
| **F-1B** UMAP × cell type | `page-04.png` panel F (scRNA-seq UMAP) | best match for our v4.2 UMAP × cell_type render |
| **F-1C** Marker gene overlay | `page-04.png` panels G+H (gene activity + ATAC dot plots) | per-marker dot/UMAP comparison |
| **F-1D** Cell-type comp × Sample.Age | `page-04.png` panel D (UMAP × age, our render flips to bar chart) | conceptual match (cell-type-by-age compositional info) |
| **F-2A** ATAC LSI + lineage tree | `page-05.png` panel B or panel G (multiome ATAC UMAP) | partial — Trevino includes a lineage tree we don't reproduce |
| **F-2BD** ATAC peak landscape | `page-05.png` panel C (ATAC/RNA heatmap) | our render aggregates per cell-type |
| **F-4A** Trajectory inference | `page-09.png` panel A + `page-07.png` panel A (pseudotime UMAP) | Trevino's main trajectory figure is Fig 4; pseudotime UMAP is Fig 3A |
| **F-4D** UMAP × DPT pseudotime | `page-07.png` panel A (pseudotime overlay UMAP) | direct match |
| **F-4EG** Branch-specific gene dynamics | `page-09.png` panels D, F-K (module gene UMAPs) + `page-12.png` panel F (GPC fuzzy clusters) | branch dynamics scattered across Cell Fig 4 + Fig 6 |
| **F-5** (supp) Peak-gene linkage | `page-05.png` panel D (GA-RNA correlation scatter) | direct match |
| **F-3** (NOT-RUN) chromVAR TF motif | `page-07.png` (Fig 3 reference) | placeholder page in our PDF |
| **F-6** (NOT-RUN) TF-gene network | `page-12.png` (Fig 6 reference) | placeholder page in our PDF |
| **F-7** (NOT-RUN) Disease GWAS | `page-13.png` (Fig 7 reference) | placeholder page in our PDF |

**RNA-only Stream A2 reference sources:**

| Plan figure ID | Trevino source | Notes |
|---|---|---|
| **F-RNA-1B** RNA-only UMAP × Leiden | `page-04.png` panel F (single-cell scRNA UMAP) — Trevino's was scRNA-seq alone before multiome integration | direct conceptual match (RNA-only seurat_clusters res=0.3 baseline) |
| **F-RNA-1C** RNA-only marker overlay | `page-04.png` panel G (gene activity dot plot) | direct match |
| **F-RNA-4A** RNA-only trajectory | `page-07.png` panel A (mRNA velocity pseudotime UMAP) | direct match if computable |

## Renaming policy

For Stream C5 PDF assembly, the comparison-PDF builder MUST reference these source files using their existing page-NN.png filenames. **No renaming is performed** in this iteration — the mapping above is the binding lookup table.

If a future iteration requires symbolic renames (e.g., `trevino_2021_fig1.png`), it must:
1. Use symlinks (preserve SHA chain back to MANIFEST.json)
2. Record the rename in this MAPPING.md under a "Renames" section with new ledger entry

## Review pending

This MAPPING.md is authored by ralph autonomous execution. Per AC-V5-TREVINO-PDF-2, a fresh-subagent (`code-reviewer`) MUST review before Stream C5 PDF assembly. Review file expected at `.omc/research/wave5/trevino_figure_mapping_review_2026-05-17.md`.

**Known risks (R-PDF-EXTRACTION-1 mitigations):**
- Page rasterization at 200 DPI; for final PDF embed may upscale to 600 DPI or re-rasterize at higher resolution
- Figure-panel-letter ambiguity in plan v5.0.4 (plan calls 'F-1B' for "UMAP × cell type" which corresponds to Trevino Fig 1 panel **F**, not panel B — see "Plan v5.0.4 figure-ID → Trevino source mapping" table above for binding interpretation)
