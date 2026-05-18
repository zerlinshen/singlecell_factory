# Wave-5 Trevino public RNA matrix viability run — 2026-05-18

Run id: `2026-05-17T2004Z-13c2c88`
Project root: `/home/zerlinshen/projects/wave5-trevino/`
Pipeline output: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458`
Autopilot phase: `ralph`

## Objective

Continue the Trevino Cell all-figure reproduction mainline from the earliest public computable input currently available. Public FASTQ/SRA raw reads remain not found in GEO/BioProject inspection, so this run starts from GEO public processed scRNA count matrices while preserving that boundary explicitly.

## Preflight

- Read `ops/before_every_run/LATEST.md` and newest Wave-5 journal before launching.
- Disk before pipeline: about `219G` available on `/dev/nvme0n1p5`.
- Scientific outputs written under `/home/zerlinshen/projects/wave5-trevino/`, not inside the factory tree.
- Factory tree was dirty; pipeline launched with `--allow-dirty` so provenance is recorded by the run manifest.

## What succeeded

### Public RNA prepared input

- Prepared input: `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_public_rna_prepared/prepared_input.h5ad`
- Shape before QC: `57868 x 33355`
- Layers mode: `none` (spliced/unspliced/ambiguous deferred to separate velocity loader audit)
- Gene symbol mapping for annotation: `{'mapped_to_symbol_var_names': 32366, 'fallback_to_gene_id': 989, 'unmapped_gene_ids': 976, 'duplicate_symbols_fell_back': 13}`
- Brain marker config: `/home/zerlinshen/projects/wave5-trevino/configs/wave5_brain_development_markers.json`

### Pipeline viability run

- Command log: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/pipeline_public_rna_viability_2026-05-17T2004Z-13c2c88.log`
- `module_status.csv`: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/module_status.csv`
- `run_manifest.json`: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/run_manifest.json`
- `final_adata.h5ad`: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/final_adata.h5ad`
- Modules all `ok`: `True`
- Final object after QC/doublet filtering: `55653 x 25519`
- Leiden clusters: `24`
- Wall time / memory evidence: `['\tElapsed (wall clock) time (h:mm:ss or m:ss): 3:05.96', '\tMaximum resident set size (kbytes): 5438124', '\tExit status: 0', 'PIPELINE_EXIT=0', 'finished_at=2026-05-18T04:08:04+08:00']`

## Annotation sanity snapshot

Cell-type counts from the brain marker panel:
- Excitatory neuron: 36503
- Radial glia / neural progenitor: 6973
- Inhibitory interneuron: 4969
- Cycling progenitor: 3020
- Intermediate progenitor: 1990
- Oligodendrocyte lineage: 827
- Microglia: 817
- Pericyte / smooth muscle: 386
- Endothelial: 168


Annotation confidence summary: `{'count': 55653.0, 'mean': 0.42570947460839476, 'std': 0.385695935275982, 'min': 0.0, '25%': 0.15585762103398637, '50%': 0.337431697845459, '75%': 0.5792499423086792, 'max': 3.5936075587272645}`

## What remains risky / next

- This is not FASTQ-level reproduction; it is public matrix-level reproduction because public FASTQs are still not located.
- Velocity-related layers were deliberately not loaded into the first pipeline object; layer gene-order mismatch is a separate loader/audit task.
- Next operator should use this run as the first successful public RNA pipeline viability baseline, then continue figure-specific reproduction gates from `pipeline_viability_summary.json` and the figure target matrix.

## Artifact classification

- canonical for current public RNA viability lane:
  - `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_public_rna_prepared/prepared_input.h5ad`
  - `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/final_adata.h5ad`
  - `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/run_manifest.json`
  - `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/module_status.csv`
  - `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/pipeline_viability/pipeline_viability_summary.json`
- evidence-only:
  - `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/pipeline_public_rna_viability_2026-05-17T2004Z-13c2c88.log`
  - smoke inputs under `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_public_rna_prepared_smoke*`
- failed exploratory:
  - earlier all-layer smoke attempts that exposed spliced/unspliced gene-order mismatch
- superseded:
  - counts-only full prepared input without gene-symbol mapping; replaced by annotation-ready build at the same path

## RNA velocity layer follow-up — 2026-05-17T20:35:41.091747+00:00

- Layer gene audit: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/rna_velocity_layer_audit/rna_layer_gene_audit.json`. Counts has `33355` genes; spliced/unspliced/ambiguous have `32648` genes; common intersection `32648`; counts-only genes dropped for velocity object `707`.
- Full velocity-ready AnnData: `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_public_rna_velocity_prepared/prepared_input.h5ad` with shape `57868 x 32648` and layers `['ambiguous', 'spliced', 'unspliced']`.
- scVelo smoke PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/rna_velocity_smoke/scvelo_smoke_retry.json`; subset `3000 x 1930` produced `velocity` layer and `velocity_graph` shape `[3000, 3000]`.
- First scVelo attempt failed only because installed scvelo 0.3.4 does not accept `n_top_genes` in `filter_and_normalize`; retry used manual top expressed genes without installing new dependencies.

## Figure 3A reproduction follow-up — 2026-05-17T20:47:44.166570+00:00

- Figure 3A public-matrix velocity outputs written under `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3A_velocity`.
- Validation artifact: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3a_public_matrix_velocity_reproduction.json` with `all_checks_pass=True` and status `REPRODUCED_PUBLIC_MATRIX_SUBSET`.
- Computed subset: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3A_velocity/figure3a_velocity_computed_subset.h5ad` shape `[6000, 2432]`; velocity graph `[6000, 6000]`; velocity pseudotime finite cells `6000`.
- Quantitative marker direction checks: [{"marker": "SOX2", "expected": "negative", "rho": -0.11768205604254171, "pass": true}, {"marker": "PAX6", "expected": "negative", "rho": -0.06033879966210802, "pass": true}, {"marker": "NEUROD1", "expected": "positive", "rho": 0.17167049594404055, "pass": true}, {"marker": "NEUROD6", "expected": "positive", "rho": 0.24301496370478703, "pass": true}, {"marker": "SLC17A7", "expected": "positive", "rho": 0.15297624508934532, "pass": true}, {"marker": "SATB2", "expected": "positive", "rho": 0.16986815690920778, "pass": true}]
- Note: velocity stream PDF uses a PNG-embedded fallback because matplotlib PDF rejected non-finite scVelo stream path vertices; PNG/SVG rendered natively and all output hashes are recorded in the summary JSON.

## Figure 3B ATAC transferred pseudotime follow-up — 2026-05-17T20:55:42.269821+00:00

- Output dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3B_atac_pseudotime`.
- Summary: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3B_atac_pseudotime/figure3b_atac_pseudotime_summary.json`.
- Validation artifact: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3b_public_matrix_atac_pseudotime_transfer.json`; generated-output checks pass `True`.
- Status is `GENERATED_PUBLIC_MATRIX_SUBSET_WEAK_MARKER_CONCORDANCE`: ATAC gene-activity transfer produced figures and finite pseudotime for 12,000 cells, but PAX6 had positive rather than expected negative Spearman rho (`+0.053877`).
- Passed marker directions: SOX2 negative; NEUROD1, NEUROD6, SLC17A7, SATB2 positive.
- Boundary: public processed-matrix transfer only; no fragments/FASTQ-level scATAC reprocessing and no author model parity is claimed.
- Artifact classification: Figure 3B output dir and validation JSON are evidence-only for current public-matrix reproduction; usable as input for Figure 3C, but not a full biological-concordance PASS.
- Next operator: begin Figure 3C from 3A/3B outputs plus Table S3B/FCM/GA RNA correlations; keep weak PAX6 caveat visible.

## Figure 3C CRE-gene heatmap follow-up — 2026-05-17T21:05:57.233741+00:00

- Output dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3C_cre_gene_heatmap`.
- Summary: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3C_cre_gene_heatmap/figure3c_cre_gene_heatmap_summary.json`.
- Validation artifact: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3c_public_matrix_cre_gene_heatmap.json`; generated-output checks pass `True`.
- Coverage: `12927/13989` usable links; derived cluster counts `{'k1': 2512, 'k2': 949, 'k3': 2666, 'k4': 3286, 'k5': 3514}`.
- Caveat: `Paper caption/methods state k=5, but local Table S3B extraction has 10 Link cluster labels.` The generated panel derives k=5 from public expression profiles and does not claim exact author 363-pseudobulk or original cluster-label parity.
- Artifact classification: evidence-only for public-matrix Figure 3C; suitable input for 3D/3E enrichment attempts.

## Figure 3D gene-set enrichment follow-up — 2026-05-17T21:08:46.960387+00:00

- Output dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3D_geneset_enrichment`.
- Summary: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3D_geneset_enrichment/figure3d_geneset_enrichment_summary.json`.
- Validation: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3d_public_matrix_local_geneset_enrichment.json`.
- Loaded `16` local gene sets; FDR<=0.05 rows `3`.
- Caveat: not exact topGO/GO parity; evidence-only local gene-set panel.

## Figure 3E motif enrichment follow-up — 2026-05-17T21:11:40.302894+00:00

- Output dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3E_motif_enrichment`.
- Summary: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3E_motif_enrichment/figure3e_motif_enrichment_summary.json`.
- Validation: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3e_public_matrix_derived_k5_motif_enrichment.json`.
- Motif matrix `657930 x 452`; Bonferroni<=0.05 positive rows `160`.
- Caveat: not exact LOLA/JASPAR/topGO parity; derived-k5 public-matrix motif panel.


## Figure 3F TF/motif heatmap follow-up — 2026-05-17T21:19:46.757259+00:00

- Output dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3F_tf_motif_heatmap`.
- Summary: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3F_tf_motif_heatmap/figure3f_tf_motif_heatmap_summary.json`.
- Validation: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3f_public_matrix_tf_motif_heatmap.json`; all checks pass `True`.
- RDS export: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3F_tf_motif_heatmap/tf_motif_expression_correlation.tsv` from `TF_MotifExpressionCorrelation.RDS`; export log records `rows=66426 cols=5`.
- Selected `31` TF/motif pairs from `1668` candidates after filters; `31` unique genes and `22` unique motifs.
- Binned cells: RNA `6000` and ATAC `12000` over `60` bins; finite expression/motif matrices validated.
- Artifact classification: evidence-only public-matrix Figure 3F panel; suitable for current honest reproduction package, not exact author pseudobulk/motif-cluster parity.
- Caveat: TF/motif heatmap from public RNA/chromVAR matrices binned by Figure 3A/3B pseudotime and TF_MotifExpressionCorrelation-derived pairs. Not exact author 363-pseudobulk, motif-cluster, or original 31/24 object parity.
- Next operator: continue Figure 3G/H/I resource audit and run code-quality checks for all new Figure 3 scripts.


## Figure 3G/3H/3I motif-correlation/synergy follow-up — 2026-05-17T21:26:43.574498+00:00

- Figure 3H output dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3H_tf_motif_gene_correlation`.
- Figure 3H summary: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3H_tf_motif_gene_correlation/figure3h_tf_motif_gene_correlation_summary.json`.
- Figure 3H validation: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3h_public_resource_tf_motif_gene_correlation.json`; all checks pass `True`.
- Figure 3H generated a motif-level TF/gene correlation heatmap from `TF_MotifExpressionCorrelation` and Figure 3F selected pairs: observed pairs `158/682`, strong |rho| pairs `136`.
- Figure 3G/3I audit dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3G_3I_motif_synergy_audit`.
- Figure 3G/3I validation: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3g_3i_public_resource_motif_synergy_gap.json`; statuses `{'3G': 'METHOD_RESOURCE_GAP_SYNERGY_NOT_REPRODUCED', '3I': 'METHOD_RESOURCE_GAP_MEAN_SYNERGY_NOT_REPRODUCED'}`.
- Audit root cause: basic input files are present, but chromVAR and motifmatchr are unavailable in the current R environment and no local/public 24 motif-cluster/synergy mapping file was found.
- Artifact classification: 3H is evidence-only public/author-resource correlation panel; 3G/3I are honest method-resource gaps, not failed fabricated figures.
- Next operator: proceed to code-quality checks and then remaining target-matrix panels; do not claim 3G/3I until chromVAR/motif-cluster resources are installed/pinned.


## Figure 3 script quality gate — 2026-05-17T21:27:46.971589+00:00

- Python syntax checks passed for public RNA loader and Figure 3A/3B/3C/3D/3E/3F/3H/gap-audit scripts.
- R parse passed for `scripts/dev/wave5_figure3_motif_enrichment.R`.
- `git diff --check` passed.
- Trevino loader regression passed: `8 passed in 0.78s` for `tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py`.
- Current Figure 3 lane is evidence-rich but not fully reproduced at the synergy level: 3G/3I remain method-resource gaps, not generated figures.


## Figure 1D/1F/1G/1H public-matrix follow-up — 2026-05-17T21:33:36.955085+00:00

- Output dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_1_public_panels`.
- Summary: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_1_public_panels/figure1_public_matrix_summary.json`.
- Validation: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure1_public_matrix_1d_1f_1g_1h.json`; all checks pass `True`.
- Generated outputs: 1D RNA/ATAC UMAP by age; 1F RNA annotation/ATAC cluster UMAP; 1G RNA marker dotplot; 1H ATAC gene-activity marker dotplot.
- Metrics: RNA cells `55653`, RNA markers `28/28`; ATAC UMAP subset cells `12000`; ATAC marker dotplot cells `31304`, ATAC markers `28/28`.
- Artifact classification: evidence-only public-matrix Figure 1 panel set; not exact FASTQ/fragments/author-embedding parity.
- Next operator: run quality checks, then continue Figure 1E or Figure 2B.


## Figure 1 script quality gate — 2026-05-17T21:34:24.636432+00:00

- Python syntax checks passed for Figure 1 and current Figure 3 public scripts.
- `git diff --check` passed.
- Trevino loader regression passed: `8 passed in 0.77s` for `tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py`.


## Figure 1E multimodal marker overlay follow-up — 2026-05-17T21:37:33.411243+00:00

- Output dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_1E_multimodal_marker_overlays`.
- Summary: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_1E_multimodal_marker_overlays/figure1e_multimodal_marker_overlay_summary.json`.
- Validation: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure1e_public_matrix_multimodal_marker_overlays.json`; all checks pass `True`.
- Generated RNA expression, ATAC gene activity, and chromVAR motif activity overlays for SOX9/EOMES/NEUROD2/DLX2.
- Caveat: DLX2 has no exact chromVAR row locally; MA0882.1_DLX6 is recorded as a DLX-family proxy, not exact DLX2 motif parity.
- Artifact classification: evidence-only public-matrix Figure 1E panel.

## 2026-05-18 05:44 CST — Figure 2B author CCA matching
- Generated Figure 2B public processed-matrix reconstruction from author `CCA_Matching.RDS`.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2b_public_author_cca_matching.json`.
- Output directory: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2B_cca_matching/`.
- Validation: status `REPRODUCED_AUTHOR_CCA_PUBLIC_MATRIX_FIGURE2B`; all checks pass; CCA rows 89,172; RNA cells 55,653; ATAC cells 12,000; coverage 1.000/1.000.
- Boundary: author matching table is used, but embeddings are current public RNA pipeline + bounded public ATAC gene-activity embedding, not raw FASTQ/fragments or exact author peak-LSI UMAP parity.
- Next: Figure 2C/2D CRE-gene linkage heatmap/correlation evidence.

## 2026-05-18 05:46 CST — Figure 2B quality gate
- Verification: `python -m py_compile ...` passed; `git diff --check` passed; `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> 8 passed in 0.79s.

## 2026-05-18 05:56 CST — Figure 2C/2D public author-resource evidence
- Generated Figure 2C proxy + Figure 2D author-resource panel.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2cd_public_author_cre_gene_resources.json`.
- Output directory: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2C_2D_cre_gene/`.
- Status 2C: `GENERATED_AUTHOR_RESOURCE_GENE_ACTIVITY_RNA_HEATMAP_NOT_CRE_ACCESSIBILITY_PAIR_PARITY`; Status 2D: `REPRODUCED_AUTHOR_RESOURCE_GA_RNA_CORRELATION_WITH_LINK_COUNT_DISCREPANCY`.
- Checks: generated-output checks pass; exact 2C paper link count does not match public RDS (64,030 vs 64,878) and exact CRE-accessibility pair heatmap is not claimed because peak-name row metadata is missing from the local accessibility pseudobulk RDS.
- Next: Figure 2E GO/gene-set enrichment / 2F-I multiome panels.

## 2026-05-18 05:58 CST — Figure 2C/2D quality gate
- Verification: `python -m py_compile ...` passed; `git diff --check` passed; `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> 8 passed in 0.79s.

## 2026-05-18 06:05 CST — Figure 2E public/local GPC enrichment
- Generated Figure 2E GPC enrichment from the public author GA/RNA + linked-CRE table.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2e_public_local_gpc_enrichment.json`.
- Output directory: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2E_gpc_enrichment/`.
- Status: `GENERATED_LOCAL_GPC_GENESET_ENRICHMENT_NOT_TOPGO_PARITY`; count-aligned GPCs 185 vs paper 185; published-rule candidates 225; gene sets tested 8.
- Top term: `GO_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY` overlap 46/185, BH q=2.1e-10.
- Boundary: not exact topGO v2.36.0 parity and not exact Table S2 variable-gene universe parity.

## 2026-05-18 06:05 CST — Figure 2B-2E quality gate
- Verification: `python -m py_compile ...` passed; `git diff --check` passed; `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> 8 passed in 0.78s.

## 2026-05-18 06:20 CST — Figure 2F/2G/2H/2I public multiome evidence
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2fghi_public_multiome_evidence.json`.
- Output directory: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2F_2I_multiome/`.
- Statuses: `{'2F': 'REPRODUCED_PUBLIC_AUTHOR_MULTIOME_QC_SUMMARY', '2G': 'GENERATED_PUBLIC_AUTHOR_RNA_AND_BOUNDED_ATAC_MULTIOME_PROJECTION', '2H': 'METHOD_RESOURCE_GAP_MULTIOME_LINKAGE_VENN_NOT_REPRODUCED', '2I': 'GENERATED_PUBLIC_MULTIOME_GPC_CORRELATION_CORRESPONDENCE'}`.
- Metrics: 2G multiome cells 8981; ATAC marker genes 70; 2I GPC n=185, Spearman rho=0.290.
- Boundary: 2H exact Venn not reproduced due missing separate multiome significant link table; 2G ATAC projection is marker-NN approximation.
- Downloaded resource checksums: `/home/zerlinshen/projects/wave5-trevino/inputs/brainchromatin_s3/processed/downloaded_key_resources.sha256`.

## 2026-05-18 06:20 CST — Figure 2F-I quality gate
- Verification: `python -m py_compile ...` passed; `git diff --check` passed; `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> 8 passed in 0.78s.


## Figure 4/5 public author FCM evidence — 2026-05-17T22:35:49.077762+00:00

- Evidence JSON: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure4_5_public_author_fcm_evidence.json`.
- Output directory: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_4_5_fcm_public_panels`.
- Statuses: `{'4A': 'GENERATED_AUTHOR_FCM_PUBLIC_PSEUDOBULK_OVERVIEW', '4B': 'REPRODUCED_AUTHOR_FCM_MODULE_HEATMAP_FROM_CENTERS', '4C': 'REPRODUCED_AUTHOR_FCM_SELECTED_GENE_HEATMAP', '4D': 'GENERATED_AUTHOR_FCM_MODULE_SCORE_UMAPS', '4E': 'REPRODUCED_AUTHOR_FCM_CENTROID_JACCARD_NETWORK', '4F': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4G': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4H': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4I': 'IMAGE_PANEL_NOT_COMPUTATIONAL_REPRODUCTION_LOCAL_PDF_ONLY', '4J': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4K': 'IMAGE_PANEL_NOT_COMPUTATIONAL_REPRODUCTION_LOCAL_PDF_ONLY', '5A': 'REPRODUCED_AUTHOR_FCM_ASTRO_GENE_MEMBERSHIP_EXPRESSION', '5B': 'GENERATED_AUTHOR_FCM_LINKED_PEAK_MOTIF_ENRICHMENT_M13_VS_M14', '5C': 'GENERATED_AUTHOR_FCM_AQP4_POSITIVE_RECLUSTERING_LOCAL_RULE', '5D': 'GENERATED_PUBLIC_FCM_PSEUDOBULK_DE_NOT_DESEQ2_PARITY', '5E': 'METHOD_RESOURCE_GAP_BHADURI_PRIMARY10X_NOT_STAGED', '5F': 'METHOD_RESOURCE_GAP_BHADURI_PRIMARY10X_NOT_STAGED'}`.
- Checks: `{'pseudobulks_eq_1267': True, 'modules_eq_14': True, 'jaccard_edges_gt_0_2_eq_32': True, 'motif_matrix_motifs_ge_400': True, 'module13_and_14_have_linked_peaks': True, 'a1_and_a2_labels_nonempty': True, 'selected_fig4_5_genes_found_full_ge_18': True, 'all_output_hashes_recorded': True}`.
- Metrics: pseudobulks `1267`, modules `14`, Jaccard links >0.2 `32`, motif matrix `657930` peaks x `452` motifs, A1-HES pseudobulks `356`, A2-OLIG pseudobulks `105`.
- Motif highlights: ASCL1/NHLH1 trend positive for m13 vs m14; SOX21 trends negative/m14-side, but local Fisher/BH enrichment is not strongly significant. Exact paper differential motif statistic parity is not claimed.
- Boundary: 4I/4K are immunohistochemistry image panels, not computational matrix outputs; 5D uses a local Welch-test pseudo-bulk proxy, not DESeq2/Table S5 parity; 5E/5F require UCSC Cell Browser `organoidreportcard/primary10X` Bhaduri data and are recorded as not staged yet.
- Next: continue Figure 6 chromatin-state/GPC branch resource audit and attempt public-resource panels.

## Figure 4/5 quality gate — 2026-05-17T22:35:49.077762+00:00

- `SC_REQUIRE_PROJECT_ROOT=1 python scripts/dev/wave5_render_public_fcm_figure4_5.py --project-root /home/zerlinshen/projects/wave5-trevino --run-id 20260517T1436Z-13c2c88`: passed.
- `python -m py_compile scripts/dev/wave5_render_public_fcm_figure4_5.py scripts/dev/wave5_render_public_multiome_figure2fghi.py scripts/dev/wave5_render_public_gpc_enrichment_figure2e.py`: passed.
- `git diff --check`: passed.
- PNG integrity check: 14 Figure 4/5 PNG outputs verified.
- `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py`: `8 passed in 0.77s`.
- Disk free `/home/zerlinshen/projects`: `207.1 GB`.


## Figure 6 public FCM/GPC branch evidence — 2026-05-17T22:42:40.390186+00:00

- Evidence JSON: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure6_public_fcm_gpc_branch_evidence.json`.
- Output directory: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_6_fcm_gpc_branch`.
- Statuses: `{'6A': 'REPRODUCED_PUBLIC_FCM_CELL_CYCLE_MODULE_CORRELATION', '6B': 'GENERATED_PUBLIC_FCM_ATAC_PROJECTION_SCHEMATIC', '6C': 'REPRODUCED_AUTHOR_FCM_ATAC_BRANCH_PROJECTION_FROM_GA_ANCHORS', '6D': 'GENERATED_PUBLIC_FCM_BRANCH_GENE_ACTIVITY_HEATMAP_LOCAL_SPECIFICITY', '6E': 'GENERATED_PUBLIC_FCM_GPC_TF_MOTIF_GENE_BRANCH_HEATMAP', '6F': 'REPRODUCED_AUTHOR_FCM_GPC_ANCHOR_REPROJECTION_BRANCH_ARROWS', '6G': 'GENERATED_BOUNDED_MULTIOME_RNA_TO_FCM_NN_PROJECTION_COLORED_BY_ATAC_CLUSTER'}`.
- Checks: `{'pseudobulks_eq_1267': True, 'branches_a_b_c_nonempty': True, 'cell_cycle_genes_found_ge_50': True, 'gpc_ids_in_gene_activity_ge_100': True, 'figure6d_top_genes_eq_50': True, 'motif_tfs_found_ge_8': True, 'multiome_cells_projected_eq_8981': True, 'multiome_projection_genes_ge_300': True, 'all_output_hashes_recorded': True}`.
- Metrics: branch A/B/C `61/27/126`, cell-cycle genes found `152`, GPC IDs in gene activity `138`, multiome cells projected `8981`, projection genes `500`, ATAC clusters `14`.
- Boundary: Figure 6G is bounded nearest-neighbor multiome RNA→FCM projection colored by public multiome ATAC clusters, not exact author uwot projection parity; Figure 6D/6E use local specificity/branch means over public branch labels.
- Next: continue Figure 7 disease/BPNet/ASD mutation resource audit.

## Figure 6 quality gate — 2026-05-17T22:42:40.390186+00:00

- `SC_REQUIRE_PROJECT_ROOT=1 python scripts/dev/wave5_render_public_fcm_figure6.py --project-root /home/zerlinshen/projects/wave5-trevino --run-id 20260517T1436Z-13c2c88`: passed.
- `python -m py_compile scripts/dev/wave5_render_public_fcm_figure6.py scripts/dev/wave5_render_public_fcm_figure4_5.py scripts/dev/wave5_render_public_multiome_figure2fghi.py`: passed.
- `git diff --check`: passed.
- PNG integrity check: 7 Figure 6 PNG outputs verified.
- `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py`: `8 passed in 0.85s`.
- Disk free `/home/zerlinshen/projects`: `207.1 GB`.

## 2026-05-18 07:06 CST — Figure 7 public ASD/BPNet evidence
- Generated/audited Figure 7 panels from official Cell Supplementary Table S6 and public ATAC subset.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure7_public_asd_bpnet_evidence.json`.
- Output directory: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_7_asd_bpnet_public_panels`.
- Statuses: `{'7A': 'GENERATED_FROM_METHODS_AND_PUBLIC_TABLE_S6_COUNTS', '7B': 'GENERATED_S6D_CLUSTER_HIGH_EFFECT_SUMMARY_ON_PUBLIC_ATAC_UMAP_NOT_EXACT_AUTHOR_FISHER_DENOMINATORS', '7C': 'GENERATED_S6_COUNT_AUDIT_EXACT_OR_1_909_DENOMINATORS_NOT_PUBLIC_IN_S6', '7D': 'GENERATED_SFARI_AUDIT_CAPTION_BENCHMARK_PLUS_RECOMPUTED_COUNTS_NOT_EXACT_24_17_FROM_LOCAL_TABLES', '7E': 'REPRODUCED_FROM_PUBLIC_TABLE_S6E_MOTIF_OVERLAP_EXCESS', '7F': 'GENERATED_NFIA_S6_LOCUS_SUMMARY_NOT_EXACT_BPNET_LOGO_TRACK_PARITY', '7G': 'GENERATED_NPY_S6_LOCUS_SUMMARY_NOT_EXACT_BPNET_LOGO_TRACK_PARITY'}`.
- Metrics: S6A/S6B mutations 42931/42049; S6D rows 800; S6D join missing 0; caption-like high-effect case/control 261/232; motif overlaps case/control 106/83.
- Boundary: exact 7B/7C denominator universe and 7F/7G trained BPNet/ref-alt/DeepLift tracks are not available in the staged public resources; 261 vs caption 262 case-count mismatch is retained as an honest resource/header ambiguity.

## 2026-05-18 07:06 CST — Figure 7 quality gate
- Verification: render command passed; py_compile passed; `git diff --check` passed; 7 PNGs verified; Trevino loader tests -> 8 passed in 0.78s.

## 2026-05-18 07:13 CST — Final quality gate / completion handoff
- Completion report: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/TREVINO_ALL_FIGURES_COMPLETION_REPORT.md`.
- Final quality gate: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/all_figures_final_quality_gate_review.json`; code review `APPROVE/CLEAR`; scientific decision `conditional`.
- Evidence audit: 17 JSON files; 53 panel statuses; missing paths 0; hash mismatches 0; unexpected false checks 0.
- Autopilot state marked complete; future work should start from optional raw FASTQ/fragments/BPNet-model acquisition if exact parity is required.
