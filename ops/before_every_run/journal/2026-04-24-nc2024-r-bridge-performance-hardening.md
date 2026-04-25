# 2026-04-24 NC2024 R bridge performance hardening

## Task / objective

Harden the remote `singlecell_factory` to R/reporting bridge so large NC2024 outputs can be plotted downstream without copying or densifying the full `final_adata.h5ad`, while preserving scientific fidelity.

## Execution mode

- `execution_mode = benchmark/design-hardening`
- This was not a new full-cohort pipeline run.
- This did not rerun prepare.
- This did not create a new `NC2024_NSCLC_FULL_COHORT...` result directory.

## What changed

- Added sparse-safe compact bundle exporter:
  - `/home/zerlinshen/singlecell_factory/scripts/export_singlecell_r_bundle.py`
- Added targeted exporter tests:
  - `/home/zerlinshen/singlecell_factory/tests/test_singlecell_r_bundle_export.py`
- Added R-side bundle manifest validation helper:
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/R/remote_bundle_manifest.R`
- Updated R bundle plotting scripts to validate required files/schema before loading:
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/plot_remote_bundle.R`
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/plot_remote_bundle_large.R`
- Updated the remote pull-and-plot wrapper so it exports and transfers a compact bundle instead of copying `final_adata.h5ad`:
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/pull_and_plot_remote_result.sh`
- Updated R bridge README to recommend the compact bundle path for large cohorts:
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/README.md`

## Verification performed

- Remote Python unit tests:
  - command: `/home/zerlinshen/conda/bin/conda run -n sc_gpu pytest -q -o addopts="" tests/test_singlecell_r_bundle_export.py`
  - first-pass result: `3 passed in 0.04s`
  - post-architect-hardening result: `4 passed in 0.04s`
- Local R parser/smoke test:
  - `Rscript` parsed the updated helper and plotting scripts.
  - A tiny manifest-backed bundle successfully generated:
    - `umap_by_group.png`
    - `umap_by_cluster.png`
    - `cell_fraction.png`
    - `marker_dot.png`
    - a PDF
- Full NC2024 fresh-run compact bundle export:
  - source:
    - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/final_adata.h5ad`
  - output:
    - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_bundle`
  - first-pass elapsed wall time: `0:58.43`
  - first-pass max RSS: `13741976 KB` (~13.1 GiB)
  - refreshed post-contract-hardening elapsed wall time: `0:59.37`
  - refreshed post-contract-hardening max RSS: `13740932 KB` (~13.1 GiB)
  - exit status: `0`

## Resulting bundle

- Bundle path:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_bundle`
- Bundle size:
  - `126M`
- Manifest schema:
  - `singlecell_r_bundle_v1`
- Expression contract:
  - source slot: `X`
  - value scale: `source_X_as_stored`
  - export dtype: `float32`
  - intended use: `plotting_and_visual_summary_only`
  - claim guard: `not_for_de_or_new_quantitative_claims_without_full_object_validation`
- Source shape:
  - `810218 x 30374`
- Cells exported:
  - `810218`
- Exported obs columns:
  - `cell_type`
  - `leiden`
  - `sample`
  - `patient`
  - `condition`
  - `total_counts`
  - `n_genes_by_counts`
  - `pct_counts_mt`
  - `doublet_score`
  - `predicted_doublet`
  - `immune_subtype`
- Exported marker genes:
  - `CD3E`
  - `LYZ`
  - `MS4A1`
  - `NKG7`
  - `ELF3`
  - `EPCAM`
  - `KRT8`
  - `KRT18`
  - `PTPRC`
- Bundle files:
  - `obs.csv.gz` (`21823347` bytes, `810218 x 11`)
  - `X_umap.csv.gz` (`12104428` bytes, `810218 x 2`)
  - `X_pca.csv.gz` (`82683770` bytes, `810218 x 20`)
  - `marker_expr.csv.gz` (`15018952` bytes, `810218 x 9`)
  - `bundle_manifest.json`
  - `bundle_manifest.tsv`
  - `README.md`

## What succeeded

- The R bridge now has a compact, manifest-backed handoff artifact for the fresh canonical NC2024 full cohort.
- The bridge no longer requires pulling or converting the whole `.h5ad` for routine R plotting.
- The exporter copies only metadata, embeddings, and selected marker expression.
- The exporter records provenance, schema, selected columns, marker list, generated files, sizes, and SHA-256 hashes.
- The R large plotting path now fails early for malformed bundles, missing manifest, wrong schema, missing file contract keys, missing expression semantics, or cell ID/order mismatches.
- The R large plotting path no longer has the `labs(y = cluster)` runtime bug.
- The Seurat-based compact plotting path now refuses non-raw-count-like compact bundles and points operators to `plot_remote_bundle_large.R` instead.

## What remains risky

- The remote server does not appear to expose `Rscript`; full R plotting remains a Mac/local bridge responsibility unless R is installed remotely.
- The full compact bundle export peaked around 13 GiB RSS, which is safe relative to the full dense matrix problem but still not "zero memory"; future exporter work can stream embeddings/obs more aggressively if needed.
- This lane improves bridge reliability and performance; it does not by itself upgrade any biological claim verdict.

## Artifact classification

- `canonical`:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_bundle`
- `canonical`:
  - `/home/zerlinshen/singlecell_factory/scripts/export_singlecell_r_bundle.py`
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/R/remote_bundle_manifest.R`
- `evidence-only`:
  - local tiny R smoke-test temporary outputs under `/tmp/nc2024_r_bundle_plot.*`
- `superseded`:
  - direct full `.h5ad` copy as default handoff strategy for large cohorts
- `failed exploratory`:
  - none

## Next operator reminder

- For large cohort R plots, prefer:
  - `scripts/export_singlecell_r_bundle.py`
  - then `bridges/local_r_pipeline_macbook/scripts/plot_remote_bundle_large.R`
- Do not use direct `.h5ad` -> Seurat conversion as the default path for NC2024-scale objects.
- If more R-side analyses require additional genes, add them to `--markers` and regenerate the compact bundle rather than loading all genes into R.
- Treat compact bundle marker values as `source_X_as_stored` visual summaries, not as DE-ready counts or a basis for new quantitative expression claims.
