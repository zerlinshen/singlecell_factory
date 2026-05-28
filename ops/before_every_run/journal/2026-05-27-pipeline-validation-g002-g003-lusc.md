# Pipeline Validation G002/G003 LUSC Marker/Mixing Run - 2026-05-27

- Context: next ultragoal round for integration biology and multiomics validation.
- Evidence root: `/home/zerlinshen/projects/pipeline-validation-20260527/`.
- Hidden Codex goal slot could not be recreated because the previous completed goal still occupies the thread slot; repo-native ultragoal artifacts were used.

## Completed

- G001 evidence and claim ledger created:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/claim_ledger.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/evidence_index.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/G001_claim_ledger_report.md`
- G002/G003/G004 LUSC real-data validation run completed on:
  - `/home/zerlinshen/projects/lusc-integration-gate-20260527/prepared/lusc_squamous_dataset_axis.h5ad`
  - 92,430 cells, 17,764 genes, 9 dataset batches, 24 `cell_type_major` labels.
- Validation script:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/scripts/run_lusc_g002_g003_marker_mixing.py`
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/marker_retention.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/embedding_label_purity.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/cell_type_mixing.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/rare_population_preservation.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g002_g003_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/G002_G003_lusc_marker_mixing_report.md`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/G004_lusc_rare_population_report.md`

## Result

- Marker retention: LUSC real-data pass.
  - 23/23 canonical/edge marker panels passed after adding rare/edge panels.
  - Marker gene coverage was 100% after mapping Ensembl IDs via `var.feature_name`.
  - Minimum marker AUROC among evaluated panels: 0.743973 (`Club`).
- Cell-type-specific mixing: LUSC partial pass with flags.
  - Harmony improved within-label same-batch neighbor fraction by >0.05 for 18/24 labels.
  - Strong improvements included T cell CD8, plasma cell, NK cell, monocyte, and macrophage.
  - Biological purity flags: `Alveolar cell type 1`, `DC mature`, and `cDC2` lost >0.10 same-label neighbor purity after Harmony.
- Rare/edge preservation: LUSC partial pass.
  - 13/19 rare/edge labels passed.
  - Review labels: `Alveolar cell type 1`, `DC mature`, `T cell regulatory`, `cDC1`, `cDC2`, `other`.

## Cautions

- Do not call G002/G003/G004 globally passed yet. Trevino validation remains pending.
- Do not ignore purity drops in small/edge labels; carry them into G004 rare/edge-population preservation.
- The first script attempt produced 0 marker coverage because the h5ad uses Ensembl IDs as `var_names`; fixed by using `var.feature_name`.
- Shuffle control collapses label purity and remains the structure-destroying negative control, not an integration candidate.

## Next

1. Repeat G002/G003/G004 on Trevino before upgrading `partial_lusc_*` to passed.
2. Run a focused LUSC sensitivity check for AT1, mature DC, cDC1, cDC2, and Treg if Trevino repeats the same pattern.
3. Continue to annotation/DE sanity only after the rare/edge-population caveats are resolved or explicitly carried forward.

## Trevino Ralph Follow-Up - 2026-05-27 23:43 CST

- Objective: complete the pending Trevino branch for G002/G003/G004 using real public RNA data before moving to downstream annotation/DE sanity.
- Input:
  - `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/final_adata.h5ad`
  - 55,653 cells, 25,519 genes, 8 samples, 9 `cell_type` labels.
- Script:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/scripts/run_trevino_g002_g003_g004_marker_mixing.py`
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_marker_retention.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_embedding_label_purity.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_cell_type_mixing.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_rare_population_preservation.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/trevino_g002_g003_g004_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/trevino_g002_g003_g004/G002_G003_G004_trevino_marker_mixing_report.md`

### Result

- Marker retention: Trevino pass, 9/9 canonical brain marker panels passed.
- Cell-type-specific mixing: Trevino partial pass with flags.
  - Harmony improved same-sample neighbor fraction by >0.05 for 8/9 labels.
  - Same-label purity drops >0.10 occurred for `Inhibitory interneuron` and `Intermediate progenitor`.
- Rare/edge preservation: Trevino partial pass.
  - 5/7 rare/edge labels passed.
  - Review labels: `Inhibitory interneuron`, `Intermediate progenitor`.

### Validation

- `python3 -m py_compile` passed for the Trevino validation script.
- `python3 -m json.tool` parsed the summary JSON.
- All four Trevino TSV outputs passed per-file column-count integrity checks.

### Carry-Forward

- G002 can now be treated as passed across LUSC and Trevino.
- G003/G004 are validation-complete but not clean scientific passes; keep the LUSC and Trevino purity/rare-population flags in the final review rather than hiding them.
- Next Ralph goal: G005 downstream annotation and sample-level DE sanity on LUSC.

## G005 LUSC Annotation and DE Sanity - 2026-05-27 23:47 CST

- Objective: validate downstream annotation/DE sanity using real LUSC count data and sample-level metadata.
- Input:
  - `/home/zerlinshen/projects/lusc-integration-gate-20260527/prepared/lusc_squamous_dataset_axis.h5ad`
  - 92,430 cells, 17,764 genes, 87 samples, 41 donors, 9 datasets.
- Script:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/scripts/run_lusc_g005_annotation_de_sanity.py`
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_annotation_summary.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_pseudobulk_sample_metadata.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_de_origin_tumor_primary_vs_normal_adjacent.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_de_tumor_stage_advanced_vs_early.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_sentinel_gene_checks.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_dataset_confounding.tsv`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/g005_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g005_annotation_de_sanity/G005_lusc_annotation_de_sanity_report.md`

### Result

- Annotation support: 23/24 labels passed marker/sample/dataset support; `other` remains review-only.
- Annotation confidence is unavailable in this object and must not be claimed.
- Origin tumor/normal sentinel direction sanity: 10/11 genes pass direction; `AQP5` mismatches.
- Pseudobulk DE is sample-level and computable, but both tested contrasts have dataset-confounding flags:
  - `origin:tumor_primary_vs_normal_adjacent` top50 DE dataset eta2 >0.50 fraction = 0.82.
  - `tumor_stage:advanced_vs_early` dataset/group Cramer's V = 0.893208 and top50 DE dataset eta2 >0.50 fraction = 1.0.

### Validation

- `python3 -m py_compile` passed for the G005 script.
- `python3 -m json.tool` parsed `g005_summary.json`.
- All six G005 TSV outputs passed column-count integrity checks.

### Carry-Forward

- G005 status: `conditional_pass_annotation_de_sanity_with_dataset_confounding_flags`.
- Do not use the current sample-level DE tables as final biological claims without a stronger model/design that handles dataset effects.
- Next Ralph goal: G006 real 3D genome / Hi-C factory bridge validation.

## G006 Real 3D Contact / Hi-C Factory Bridge - 2026-05-27 23:55 CST

- Objective: validate the 3D genome factory bridge with a real NG2025/LUSC contact artifact rather than synthetic contacts.
- Source:
  - `/home/zerlinshen/projects/ng2025-3d-genome/data/open/hic_lung/TCGA_HiChIP_hic/LUSC_H3K27ac.allValidPairs.hic`
  - Extracted chr21 at 100 kb through `hicstraw`.
- Script:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/scripts/run_g006_hicstraw_factory_bridge.py`
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/LUSC_H3K27ac_chr21_100000bp.contacts.tsv.gz`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/factory_run/hic_ingest/hic_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/factory_run/hic_tad/tad_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/g006_real_hic_factory_bridge.h5ad`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/singlecell_r_bundle_v22_hic_real/bundle_manifest.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/g006_summary.json`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g006_hic_factory_bridge/G006_hic_factory_bridge_report.md`

### Result

- Real extraction: 49,616 hicstraw records, total contact count 1,766,659.
- Factory modules: `hic_ingest` and `hic_tad` status `ok`.
- Bundle: schema `singlecell_r_bundle_v2.2`; Hi-C extension `active`; 467 bins; 98,860 sparse matrix contacts.
- Reference comparison:
  - Bundle contacts equal matrix nnz.
  - Ingest total contacts equal hicstraw total contacts.
- Regression check:
  - `python3 -m pytest tests/test_wave2b_hic_smoke.py tests/test_singlecell_r_bundle_export.py::test_export_bundle_v22_hic_extension -q --no-cov`
  - Result: `13 passed in 0.54s`.

### R-Side Load and Boundaries

- R-side load passed using `/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript`.
- stdout confirmed `STATUS=active`, `CONTACTS_NROW=98860`, and `COMPARTMENT_STATUS=low_information`.
- The real source is H3K27ac HiChIP, not unbiased Hi-C. This validates factory I/O and bundle extension wiring, but not final unbiased TAD/compartment biology.
- Compartment status is `low_information` for chr21; do not promote compartment biology claims from this run.

### Carry-Forward

- G006 status: `passed_with_hichip_scientific_boundary`.
- Next Ralph goal: G007 cross-factory handoff and documentation sync.

## G007 Cross-Factory Doc Sync - 2026-05-28 00:02 CST

- Objective: make the human-facing and agent-facing docs agree with the current G001-G006 claim ledger.
- Outputs:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/G007_cross_factory_doc_sync_report.md`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g007_doc_sync_summary.json`
- Updated `singlecell_factory` surfaces:
  - `README.md`
  - `AGENTS.md`
  - `AI_AGENT_PROTOCOL.md`
  - `CODEX.md`
  - `CLAUDE.md`
  - `PROTOCOL.md`
  - `ops/before_every_run/LATEST.md`
  - `ops/before_every_run/journal/2026-05-27-pipeline-validation-g002-g003-lusc.md`
- Updated `r_multiomics_factory` surfaces:
  - `README.md`
  - `AGENTS.md`
  - `AI_AGENT_PROTOCOL.md`
  - `CODEX_PROFILE.md`
  - `CLAUDE.md`

### Validation

- `repo-doc-sync --strict` for `singlecell_factory`: NO DRIFT.
- `repo-doc-sync --strict` for `r_multiomics_factory`: NO DRIFT.
- `git diff --check` passed for touched docs in both repos.
- `python3 -m py_compile` passed for all three validation scripts.
- `g007_doc_sync_summary.json` parsed successfully.

### Carry-Forward

- G007 status: `passed`.
- Docs now preserve G002 pass, G003/G004 purity/rare-population flags, G005 dataset-confounded DE, and G006 H3K27ac HiChIP low-information 3D-genome boundary.
- Next Ralph goal: G008 final independent code and science review.

## G008 Post-Review Reconciliation - 2026-05-28 00:25 CST

- Objective: fix the review blockers before final synthesis instead of declaring
  the validation round complete on stale evidence.
- Fixed script gates:
  - G006 no longer treats missing R-side validation as a technical pass; decision
    now requires `hic_ingest`, `hic_tad`, bundle schema, active Hi-C extension,
    contact-count parity, ingest total-count parity, and R-side load.
  - G005 no longer falls back from `layers["counts"]` to `adata.X` for pseudobulk
    DE.
  - G002/G003/G004 LUSC and Trevino now emit explicit marker label-shuffle and
    embedding-shuffle negative-control metrics.
- Real-data rerun evidence:
  - LUSC marker retention: 23/23 panels pass; 23/23 marker negative controls pass.
    Embedding negative control fired on 24/24 labels. Rare/edge result remains
    13/19 pass with review labels retained.
  - Trevino marker retention: 9/9 panels pass; 9/9 marker negative controls pass.
    Embedding negative control fired on 9/9 labels. Rare/edge result remains
    5/7 pass with review labels retained.
  - G005 raw-count gate found 3,070 fractional-count cells across 16 samples
    (`Guo_Zhang_2018`, `Maynard_Bivona_2020`). Those samples were excluded rather
    than rounded. DE ran on 71 raw-count-compatible samples and remains
    dataset-confounded for both tested contrasts.
  - G006 all 7 technical gates pass, including R-side load, but the scientific
    boundary remains H3K27ac HiChIP chr21 and `low_information` compartments.
- Updated artifacts:
  - `claim_ledger.tsv`
  - `g007_doc_sync_summary.json`
  - `G007_cross_factory_doc_sync_report.md`
  - current-state docs in both `singlecell_factory` and `r_multiomics_factory`.
- Carry-forward:
  - Final verdict cannot be "pipeline scientifically final-ready"; it should be
    partial/conditional with explicit blockers for rare-label purity, DE
    confounding/count compatibility, and unbiased 3D-genome biology.

## G008 Final Review - 2026-05-28 00:40 CST

- Final artifacts:
  - `/home/zerlinshen/projects/pipeline-validation-20260527/G008_final_code_science_review_report.md`
  - `/home/zerlinshen/projects/pipeline-validation-20260527/g008_final_review_summary.json`
- Independent review result:
  - Code-review lane: `APPROVE`.
  - Architect/science lane: `WATCH`.
- Registered verdict: `PASS_SUPPORTED_NOT_FINAL`.
- Interpretation:
  - The pipeline is supported for bounded tested claims.
  - It is not global final-claim-ready because G003/G004 retain rare/purity flags,
    G005 DE remains raw-count-subset and dataset-confounded, and G006 is a
    H3K27ac HiChIP low-information technical bridge rather than unbiased
    3D-genome biology.
