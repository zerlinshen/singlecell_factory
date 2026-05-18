# Multi-Omics Module Rationale

Scope: this note documents the current Phase 2 multi-omics contract between
`singlecell_factory`, `r_multiomics_factory` (R-native analysis), and
`plotting_factory` (introduced 2026-05-18; `r/modality/*_plots.R` owns the plot
halves extracted from each R modality module in Phase 2.1). ATAC is the active
v2.2 vertical slice. VDJ, Ribo-seq, and Hi-C remain reserved bundle slots until
their exporters have full Python-to-R round-trip tests.

## ATAC ingest

Rationale: store the peak-count matrix as sparse `adata.obsm["atac_peaks"]`,
peak coordinates in `adata.uns["atac_peaks"]` / `adata.uns["atac_var"]`, and a
compatibility LSI embedding for downstream plotting. Peak matrices are sparse by
nature; keeping them sparse preserves memory safety on large cohorts.

Literature support:

- Cusanovich et al. 2018, Cell, doi:10.1016/j.cell.2018.06.052: sparse
  single-cell chromatin accessibility matrices and TF-IDF/LSI-style reduction.
- Stuart et al. 2021, Nature Methods, doi:10.1038/s41592-021-01282-5: Signac
  conventions for scATAC storage, QC, and multimodal chromatin analysis.
- Zhang et al. 2024, Nature Methods, doi:10.1038/s41592-023-02139-9: SnapATAC2
  as the Python-side ATAC analysis reference.

## ATAC LSI

Rationale: apply TF-IDF to the sparse peak matrix, run TruncatedSVD, and drop the
first component by default because it is commonly depth-correlated. Component
count is bounded by the number of peaks, preserving the 50-component production
default when legal while allowing small smoke fixtures.

Literature support:

- Cusanovich et al. 2018, Cell, doi:10.1016/j.cell.2018.06.052.
- Stuart et al. 2021, Nature Methods, doi:10.1038/s41592-021-01282-5.
- Granja et al. 2021, Nature Genetics, doi:10.1038/s41588-021-00790-6: ArchR
  iterative LSI and first-component dropping convention.

## Peak-to-gene links

Rationale: link ATAC peaks to nearby gene TSS windows using RNA expression and
sparse peak accessibility. This is an interpretable first-pass bridge, not a
replacement for full co-accessibility modeling.

Literature support:

- Pliner et al. 2018, Molecular Cell, doi:10.1016/j.molcel.2018.06.044:
  Cicero peak co-accessibility and gene activity framing.
- Trevino et al. 2021, Cell, doi:10.1016/j.cell.2021.07.039: developmental
  enhancer-to-gene linking from single-cell multiome data.
- Ma et al. 2020, Cell, doi:10.1016/j.cell.2020.09.056: SHARE-seq joint
  accessibility/expression evidence for cis-regulatory links.

## R bundle ATAC extension

Rationale: the compact R bundle exports only plotting/reporting-scale ATAC data:
LSI coordinates and peak metadata. It intentionally does not serialize the full
sparse peak-count matrix into the R plotting bundle. The full AnnData remains
the authority for quantitative ATAC claims.

Supporting contracts:

- `scripts/export_singlecell_r_bundle.py::maybe_export_atac(...)` writes
  `extensions/atac/lsi.parquet` with a `cell` column and
  `extensions/atac/peaks.parquet` with `peak_id`, `chrom`, `start`, `end`.
- `r_multiomics_factory/R_bundle/io_bundle.R` attaches
  `bundle$extensions$atac$data`.
- `r_multiomics_factory/R/atac_module.R` fails soft for absent/reserved slots.

## VDJ reserved slot

Rationale: VDJ ingest and repertoire metrics exist in Python, but the bundle
slot stays `reserved` until exporter, R loader, and round-trip tests are
complete. The intended semantics are per-cell clonotype annotations, chain
pairing, clonal expansion, and per-sample diversity metrics.

Literature support:

- Scirpy, doi:10.1093/bioinformatics/btaa611: AnnData-native VDJ analysis
  conventions.
- Bagaev et al. 2015, PLoS Computational Biology, doi:10.1371/journal.pcbi.1004503:
  VDJtools repertoire metrics.
- Shannon 1948, doi:10.1002/j.1538-7305.1948.tb01338.x: entropy metric used for
  repertoire diversity.

## Ribo-seq reserved slot

Rationale: Ribo-seq ingest can compute per-gene, per-sample translation
efficiency from footprint and RNA counts, but the bundle slot stays `reserved`
until exporter and R reader tests are present.

Literature support:

- Ingolia et al. 2009, Science, doi:10.1126/science.1168978: ribosome
  profiling and footprint count semantics.
- Ingolia 2014, Nature Reviews Genetics, doi:10.1038/nrg3950: ribosome
  profiling interpretation and translation efficiency cautions.

## Hi-C reserved slot

Rationale: Hi-C ingest stores sparse binned contact matrices, and Hi-C TAD
scoring derives insulation boundaries plus first-eigenvector compartments. The
bundle slot remains `reserved` because exporting sparse contact maps and
compartment tables requires a dedicated R contract and scale-aware tests.

Literature support:

- Lieberman-Aiden et al. 2009, Science, doi:10.1126/science.1181369: original
  Hi-C contact-map framing.
- Crane et al. 2015, Nature, doi:10.1038/nature14450: insulation score for TAD
  boundary detection.
- Nora et al. 2012, Nature, doi:10.1038/nature11049: TAD organization.
- Abdennur and Mirny 2020, Bioinformatics, doi:10.1093/bioinformatics/btz540:
  Cooler storage conventions for genomic contact matrices.

## Densify and scale policy

Rationale: intentional dense conversions must be bounded and explicit. Large
cell-by-feature sparse matrices route through `plan_densify(...)` or carry a
`# densify-allowed:` marker for small local operations. This keeps ATAC, Hi-C,
and differential-expression paths from silently converting large sparse
matrices into dense arrays.
