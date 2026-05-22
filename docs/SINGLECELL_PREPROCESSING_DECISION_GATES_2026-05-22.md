# Single-Cell Preprocessing Decision Gates — 2026-05-22

This note records the current pipeline facts behind QC, normalization, clustering
resolution, batch correction, annotation, and CNV ordering. It is a guardrail for
future optimization rounds: do not choose a new default from one benchmark alone.
Any new default needs either a second dataset with a different shape or a
data-shape conditional that activates it.

## Current Facts

- QC thresholds are fixed, configurable CLI values: `min_genes=200`,
  `max_genes=7000`, `min_counts=500`, `max_counts=50000`,
  `max_mito_pct=20`, `max_ribo_pct=50`, `min_cells=3`. The QC module now writes
  `qc_threshold_audit.json` with thresholds, metric quantiles, and per-filter
  failure counts so threshold choices can be compared across datasets.
- Normalization is Scanpy total-count normalization plus log transform:
  `sc.pp.normalize_total(target_sum=1e4)` followed by `sc.pp.log1p`. It is not
  Seurat `SCTransform`, and the Python AnnData pipeline does not create a Seurat
  `SCT` assay. If an R/Seurat lane uses SCT in the future, PCA/neighbors/UMAP and
  clustering in that lane must explicitly use the SCT assay, and the Python lane
  must compare against it as a separate method, not silently replace defaults.
- Leiden resolution is still controlled by `--leiden-resolution`; use
  `--leiden-resolution-sweep 0.5,0.8,1.0` to write a diagnostic
  `leiden_resolution_sweep.csv` without changing final labels. The table records
  cluster count and small-cluster burden per resolution for visual/manual review.
  If batch correction reruns Leiden, the batch module writes its own
  corrected-representation sweep under `batch_correction/`.
- Batch correction supports Harmony by default and scVI as an explicit backend.
  The batch module already writes before/after UMAPs and now also writes
  `batch_mixing_metrics.json`. A good correction increases neighborhood batch
  entropy and lowers same-batch neighbor fraction while preserving true cell-type
  separation. If cells of the same biological type remain split by sample after
  correction, that is residual batch effect; if different cell types remain
  separate, that is expected biology.
- Annotation uses cluster-voting marker scoring plus optional KNN reference label
  transfer. It already writes `cluster_score_matrix.csv`,
  `cell_type_annotation.csv`, and cluster-majority labels. It now also writes
  epithelial marker support files for EPICAM/KRT8/KRT18 when present.
- CNV now depends on annotation in the module DAG. CNV writes
  `cnv_annotation_qc.json`, including cell-type counts and reference-group
  usability, before producing CNV classifications.

## Required Optimization Shape

1. Run at least two data shapes before changing defaults: for example a small
   curated 10x/reference dataset and a larger multi-sample tumor cohort.
2. Compare QC threshold retention, mitochondrial/ribosomal tails, cluster marker
   coherence, batch-mixing metrics, epithelial marker support, and CNV reference
   validity together.
3. Treat SCT, scVI, Harmony, and alternative resolution choices as method lanes.
   Promote a lane only with cross-dataset evidence or a clear shape-conditional
   rule.
4. Final publication figures should be regenerated through the project figure
   quality policy; diagnostic PNGs prove module behavior but are not manuscript
   panels by themselves.
