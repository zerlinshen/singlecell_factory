# Latest NC2024 legacy artifact cleanup update — 2026-05-18

- Off-mainline cleanup completed after the module-flow validation was recorded as `evidence-only`.
- Deleted evidence-only run: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88`.
- Deleted factory legacy/smoke residue: `/home/zerlinshen/singlecell_factory/results/test_20260516_*`, `/home/zerlinshen/singlecell_factory/results/ops`, and `/home/zerlinshen/singlecell_factory/output`.
- Preserved canonical NC2024 source of truth: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.
- Preserved module-flow validation report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-module-flow-validation/REPORT.md`.
- Preserved lightweight historical validation evidence: `/home/zerlinshen/singlecell_factory/results/small_real_validate_20260426`.
- Reclaimed bytes: `1614159948` (~1.50 GiB).
- Cleanup record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-legacy-artifact-cleanup/`.
- New journal: `before-every-run/journal/2026-05-18-nc2024-legacy-artifact-cleanup.md`.

---

# Latest NC2024 module flow validation update — 2026-05-18

- Controlled project-root validation run completed under `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T1013Z-13c2c88`.
- Verdict: `MODULE_FLOW_VALIDATION_PASS_WITH_TEST_CONTRACT_DRIFT`.
- Module status: 8 `ok`, 1 expected `skipped` (`batch_correction`, one `sample` batch only). Final AnnData `5281 x 19504`.
- Flow details: marker DB loaded 5 DBs / 198 entries; GPU clustering completed; context-aware annotation had 0 validation mismatches; DE used expected CPU fallback because current `rapids_singlecell.tl` lacks `rank_genes_groups`; doublet detection used fallback all-singlets after RAPIDS scrublet dtype rejection.
- Project governance validator passed for the new run: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-module-flow-validation/`.
- Module/process tests: `verify_module_references` OK (`45/45 modules verified; 2 project-local`), DAG/project-root/governance tests `29 passed`, minimal modular pipeline test `1 passed`.
- Known test-contract drift: broad targeted pytest including all `tests/test_modular.py` produced `78 passed`, `5 failed`; failures are stale DAG-size/GPU-fallback/copy-object expectations versus the current expanded DAG and M2 GPU failure policy.
- Artifact classification: new run is `evidence-only`; current NC2024 source of truth remains `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`.
- New journal: `before-every-run/journal/2026-05-18-nc2024-module-flow-validation.md`.

---

# Latest paper reproduction ladder policy update — 2026-05-18

- Reproduction workflow is now faithful-first: clone/stage upstream repo, pin commit/data/license, run raw-data reproduction when possible, reproduce both data objects and figures, classify parity/resource gaps, then perform module-gap analysis and context optimization.
- Governance record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-paper-reproduction-ladder-policy/`.
- Project policies for Cell/Trevino and NC2024 now include `upstream_repository`, `raw_data_reproduction`, `data_object_reproduction`, `figure_reproduction`, `module_gap_decisions`, and `context_optimization_decisions`.
- Skills were updated both locally and on `ubuntu-tail` so future agents use this ladder before adding/updating reusable modules.
- No pipeline run was launched and no data was deleted in this policy round.

---

# Latest NC2024 project retention and Cell status update — 2026-05-18

- NC2024/cancer project retention rule adopted: keep only the latest validated project run result; older bulky project run outputs can be deleted once inventory/provenance and replacement source of truth are recorded.
- Applied to `/home/zerlinshen/projects/nc-reproduction`: kept `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/`; deleted 3 older run directories totaling `142437966218` bytes.
- Cleanup record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-project-latest-run-retention/`.
- Post-cleanup project governance validation: `pass`, severity `{}`, 1 run inspected.
- Cell/Trevino reproduction status checked: quality gate decision `conditional`; 17 evidence JSON files, 31 direct evidence panels, no unexpected false checks, no missing referenced paths, no hash mismatches.
- Cell status record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-cell-reproduction-status/REPORT.md`.

---

# Latest NC2024 real run after structure-change update — 2026-05-18

- Legacy factory NC2024 `results/` cleanup completed after user approval: 48 entries, `30667532081` bytes removed; `results/` is now `27M` and has zero NC2024/nc2024 leftovers.
- True NC2024 pipeline run completed under project root: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/`.
- Run verdict: `REAL_RUN_PASS_PROJECT_ARCHITECTURE_CONFIRMED`.
- Report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-real-run-validation/REPORT.md`.
- Module status: 8 `ok`, 1 expected `skipped` (`batch_correction`, one batch only). Final AnnData `5281 x 19504`.
- Recovery note: old `launch_wave3_nc2024.sh` generated invalid non-hex run ids under the new architecture; launcher was patched to use the current factory git short SHA.
- Rule going forward: NC2024 scientific outputs live under `/home/zerlinshen/projects/nc-reproduction/runs/`; factory keeps only code/contracts/validators/docs/governance records.

---

# Latest NC2024 architecture validation update — 2026-05-18

- New architecture validation using `/home/zerlinshen/projects/nc-reproduction` completed.
- Verdict: `PROJECT_ARCHITECTURE_PASS_WITH_FACTORY_LEGACY_DEBT`.
- Report: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-nc2024-architecture-validation/REPORT.md`.
- Project validator: `pass`, severity `{'info': 2, 'warning': 1}`, 3 runs inspected.
- No governance files were written under the project root; validation outputs are factory control-plane only.
- Legacy debt: 48 NC2024 result paths remain under `singlecell_factory/results/`, total `30667532081` bytes. Do not create new NC2024 scientific outputs there; perform a separate retention cleanup before deleting legacy results.

---

# Latest before/after-run note — Wave-5 Trevino public Cell reproduction

## 2026-05-18 05:56 CST — Figure 2C/2D public author-resource evidence
- Mainline: Cell paper all-figure reproduction from public earliest-computable inputs; still not FASTQ-level.
- Completed this checkpoint: Figure 2C proxy heatmap and Figure 2D GA/RNA correlation scatter from author/public resources.
- Status 2C: `GENERATED_AUTHOR_RESOURCE_GENE_ACTIVITY_RNA_HEATMAP_NOT_CRE_ACCESSIBILITY_PAIR_PARITY`.
- Status 2D: `REPRODUCED_AUTHOR_RESOURCE_GA_RNA_CORRELATION_WITH_LINK_COUNT_DISCREPANCY`.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2cd_public_author_cre_gene_resources.json`.
- Outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2C_2D_cre_gene/`.
- Metrics: significant links 64,030 vs paper expected 64,878 (delta -848); heatmap genes 400; pseudobulks 1267; GA/RNA table rows 19,173; unique genes 2,739.
- Boundary: Figure 2D uses author GA_RNA_Correlations.RDS joined to public author PeakGeneLinks_Significant.RDS. Figure 2C is a proxy heatmap using author glial KNN pseudobulk ATAC gene activity and RNA matrices for genes with significant linked CREs. Exact Figure 2C CRE-accessibility pair rows are not claimed because the locally available accessibility pseudobulk RDS lacks peak names and therefore cannot be safely joined to peak.name without additional author row-order metadata.
- Next immediate action: Figure 2E GO/gene-set enrichment, then 2F-I multiome panels.

# Latest before/after-run note — Wave-5 Trevino public Cell reproduction

## 2026-05-18 05:44 CST — Figure 2B author CCA matching generated/validated
- Mainline: Cell paper all-figure reproduction from public earliest-computable inputs; still not FASTQ-level because public FASTQ/SRA raw reads are not available in the current evidence set.
- Completed this checkpoint: Figure 2B matched RNA/ATAC cluster UMAPs from author public `CCA_Matching.RDS`.
- Status: `REPRODUCED_AUTHOR_CCA_PUBLIC_MATRIX_FIGURE2B`; all validation checks pass.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2b_public_author_cca_matching.json`.
- Outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2B_cca_matching/`.
- Metrics: CCA rows 89,172; RNA cells 55,653; ATAC cells 12,000; RNA/ATAC match-cluster coverage 1.000/1.000; RNA cell-type proxy agreement 0.808.
- Truth boundary: Author public CCA_Matching.RDS was used for RNA↔ATAC nearest-neighbor links. RNA UMAP is the current singlecell_factory public RNA pipeline embedding; ATAC UMAP is the bounded public gene-activity embedding produced for Figure 3B. This is not FASTQ/fragments-level reprocessing and not exact author peak-LSI UMAP parity.
- Verification: `python -m py_compile ...` passed; `git diff --check` passed; `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> 8 passed in 0.79s.
- Next immediate action: Figure 2C/2D CRE-gene linkage heatmap/correlation evidence using public/author peak-gene resources.

# Latest Wave-5 Figure 1E quality gate update (2026-05-18)

- Quality gate PASS after Figure 1E addition: Python syntax checks, `git diff --check`, and Trevino loader regression tests passed.
- Test evidence: `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> `8 passed in 0.79s`.
- Next immediate action: Figure 2B public RNA/ATAC integration evidence.

---

# Latest Wave-5 Figure 1E multimodal overlay update (2026-05-18)

- Figure 1E public-matrix multimodal marker overlays are available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_1E_multimodal_marker_overlays`.
- Validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure1e_public_matrix_multimodal_marker_overlays.json`; all checks pass `True`.
- Metrics: RNA cells `55653`, ATAC overlay cells `12000`, markers `SOX9, EOMES, NEUROD2, DLX2`.
- Honest caveat: RNA expression, ATAC gene activity, and chromVAR motif activity were rendered from public processed matrices. ATAC panels use the 12k public ATAC subset embedding. DLX2 uses a DLX6 motif row as a DLX-family proxy because no exact DLX2 chromVAR row was found.
- Next immediate action: Figure 2B public RNA/ATAC integration evidence.

---

# Latest Wave-5 Figure 1 quality gate update (2026-05-18)

- Quality gate PASS after Figure 1D/1F/1G/1H addition: Python syntax checks, `git diff --check`, and Trevino loader regression tests passed.
- Test evidence: `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> `8 passed in 0.77s`.
- Next immediate action: Figure 1E multimodal marker overlays or Figure 2B public integration evidence.

---

# Latest Wave-5 Figure 1 public-matrix update (2026-05-18)

- Figure 1D/1F/1G/1H public-matrix outputs are available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_1_public_panels`.
- Validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure1_public_matrix_1d_1f_1g_1h.json`; all checks pass `True`.
- Metrics: RNA `55653` cells, RNA markers `28/28`; ATAC UMAP subset `12000` cells; ATAC marker dotplot `31304` cells, ATAC markers `28/28`.
- Honest caveat: Public processed RNA and ATAC gene-activity matrices were used. ATAC UMAP panels use the existing 12k Figure 3B subset embedding; ATAC marker dotplot uses all public 31,304 ATAC metadata-aligned cells. This is not FASTQ/fragments-level or exact author embedding parity.
- Next immediate action: Figure 1E multimodal overlays or Figure 2B public integration evidence.

---

# Latest Wave-5 Figure 3 quality gate update (2026-05-18)

- Quality gate PASS after Figure 3F/3H/gap-audit additions: Python syntax checks, R parse, `git diff --check`, and Trevino loader regression tests all passed.
- Test evidence: `pytest --no-cov -q tests/test_trevino_loader.py tests/test_trevino_synthetic_fixture.py` -> `8 passed in 0.78s`.
- Current Figure 3 status: 3A/3C/3D/3E/3F/3H generated/validated; 3B generated with weak PAX6 concordance caveat; 3G/3I honestly marked motif-synergy method-resource gaps.
- Next immediate action: continue remaining target-matrix panels outside Figure 3, prioritizing main Figure 1/2 READY_FOR_LOADER panels.

---

# Latest Wave-5 Figure 3G/3H/3I motif-correlation/synergy update (2026-05-18)

- Figure 3H motif/gene correlation heatmap outputs are available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3H_tf_motif_gene_correlation`.
- Figure 3H validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3h_public_resource_tf_motif_gene_correlation.json`; metrics: `158/682` observed gene-motif correlations, `136` strong |rho| pairs.
- Figure 3G/3I synergy audit: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3g_3i_public_resource_motif_synergy_gap.json`; statuses `{'3G': 'METHOD_RESOURCE_GAP_SYNERGY_NOT_REPRODUCED', '3I': 'METHOD_RESOURCE_GAP_MEAN_SYNERGY_NOT_REPRODUCED'}`.
- Honest caveat: 3H is author-resource motif-level correlation evidence, not exact 24 motif-cluster parity; 3G/3I are not reproduced because chromVAR/motifmatchr and the motif-cluster synergy mapping/workflow are absent locally.
- Next immediate action: run quality checks for new Figure 3 scripts, then proceed to remaining main-figure target-matrix panels.

---

# Latest Wave-5 Figure 3F TF/motif heatmap update (2026-05-18)

- Figure 3F public-matrix TF expression / motif activity heatmap outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3F_tf_motif_heatmap`.
- Validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3f_public_matrix_tf_motif_heatmap.json`; all checks pass `True`.
- Metrics: candidate pairs after filters `1668`, selected pairs `31`, unique genes `31`, unique motifs `22`; RNA cells binned `6000`, ATAC cells binned `12000`.
- Honest caveat: status `GENERATED_PUBLIC_MATRIX_TF_MOTIF_HEATMAP_DERIVED_PAIRS`; TF/motif heatmap from public RNA/chromVAR matrices binned by Figure 3A/3B pseudotime and TF_MotifExpressionCorrelation-derived pairs. Not exact author 363-pseudobulk, motif-cluster, or original 31/24 object parity.
- Next immediate action: Figure 3G/H/I resource audit or remaining Figure 3 panels, then code-quality checks for newly added scripts.

---

# Latest Wave-5 Figure 3E motif enrichment update (2026-05-18)

- Figure 3E public-matrix derived-k5 motif enrichment outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3E_motif_enrichment`.
- Validation PASS for generated outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3e_public_matrix_derived_k5_motif_enrichment.json`.
- Motif matrix readiness: `657930 x 452`; cluster count `5`; all link rows mapped to peak IDs.
- Enrichment: `1021` positive rows in `all_consensus_peaks`; Bonferroni<=0.05 positive rows `160`.
- Honest caveat: status `GENERATED_PUBLIC_MATRIX_DERIVED_K5_MOTIF_ENRICHMENT`; not exact LOLA/JASPAR/topGO environment parity and inherits Figure 3C derived-k5 caveat.
- Next immediate action: Figure 3F TF expression/motif activity heatmap resource audit/build.

---

# Latest Wave-5 Figure 3D gene-set enrichment update (2026-05-18)

- Figure 3D public-matrix local gene-set enrichment outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3D_geneset_enrichment`.
- Validation PASS for generated outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3d_public_matrix_local_geneset_enrichment.json`.
- Loaded `16` local author gene-set files/terms; tested `35` enrichment rows; FDR<=0.05 rows `3`.
- Honest caveat: status `GENERATED_PUBLIC_MATRIX_LOCAL_GENESET_ENRICHMENT_NOT_TOPGO_PARITY`; this is not exact topGO v2.36.0 / GO database-version parity.
- Next immediate action: Figure 3E motif enrichment rework on Figure 3C derived k5 clusters.

---

# Latest Wave-5 Figure 3C CRE-gene heatmap update (2026-05-18)

- Figure 3C public-matrix derived-k5 CRE-gene heatmap outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3C_cre_gene_heatmap`.
- Validation PASS for generated outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3c_public_matrix_cre_gene_heatmap.json`; PNG/SVG/PDF, link-row TSV, matrix NPZ, summary JSON, and log are present/nonempty.
- Coverage: `12927/13989` Table S3B links usable (`0.924`); peak-coordinate mapping complete; finite ATAC accessibility and RNA expression matrices.
- Honest caveat: the paper declares k=5, but local Table S3B has `10` cluster labels; this run derives k=5 from public RNA expression profiles and does not claim exact author 363-pseudobulk/CCA/fragment-level parity.
- Next immediate action: Figure 3D/3E enrichment panels using `figure3c_link_row_summary.tsv`, author motif resources, and local gene-set files.

---

# Latest Wave-5 Figure 3B ATAC pseudotime transfer update (2026-05-18)

- Figure 3B public-matrix ATAC transferred pseudotime outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3B_atac_pseudotime`.
- Generated-output validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3b_public_matrix_atac_pseudotime_transfer.json`; all expected PNG/SVG/PDF/H5AD/JSON/log files exist and are nonempty.
- Honest biological status: `GENERATED_PUBLIC_MATRIX_SUBSET_WEAK_MARKER_CONCORDANCE`, because PAX6 direction failed (`rho=+0.053877` where pre-registered expectation was negative). SOX2 and NEUROD1/NEUROD6/SLC17A7/SATB2 directions passed.
- Computed ATAC subset: `12000 / 31304` cells; transferred from Figure 3A RNA velocity subset with `6000` finite pseudotime reference cells; `25` shared transfer genes; finite transferred pseudotime cells `12000`.
- Boundary: public processed ATAC gene activity + RNA velocity transfer; not fragments/FASTQ-level scATAC and not author integration-model parity.
- Next immediate action: start Figure 3C linkage/co-accessibility evidence using Figure 3A/3B plus Table S3B/FCM/link resources.

---

# Latest Wave-5 Figure 3A velocity reproduction update (2026-05-18)

- Figure 3A public-matrix RNA velocity/pseudotime outputs are now available at `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_3A_velocity`.
- Validation PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure3a_public_matrix_velocity_reproduction.json`; status `REPRODUCED_PUBLIC_MATRIX_SUBSET`; all output PNG/SVG/PDF files exist and nonempty.
- Computed subset `[6000, 2432]` from `55653` overlapping QC-passed public RNA cells; velocity graph `[6000, 6000]`; velocity pseudotime finite cells `6000`.
- Marker direction check passed: SOX2/PAX6 decrease with pseudotime, NEUROD1/NEUROD6/SLC17A7/SATB2 increase with pseudotime.
- Boundary: this is public processed-matrix, 6000-cell stratified-subset reproduction, not FASTQ-level or full-cell velocity graph.

---

# Latest Wave-5 RNA velocity input update (2026-05-18)

- Full velocity-ready public RNA input: `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_public_rna_velocity_prepared/prepared_input.h5ad`; shape `57868 x 32648`; layers `['ambiguous', 'spliced', 'unspliced']`.
- Layer audit: counts has `33355` genes, velocity layers have `32648` genes, intersection `32648`; `707` counts-only genes are intentionally dropped only for the velocity object.
- scVelo smoke PASS: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/rna_velocity_smoke/scvelo_smoke_retry.json`; 3k-cell subset produced `velocity` layer and `velocity_graph` `[3000, 3000]` with scvelo `0.3.4`.
- Keep the counts-only pipeline object as the clustering/annotation baseline; use this layer-aligned object for Figure 3A/RNA velocity reproduction.

---

# Latest Wave-5 public RNA viability update (2026-05-18)

- Current public-matrix pipeline run: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/`.
- Prepared public RNA input: `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_public_rna_prepared/prepared_input.h5ad`; shape before QC `57868 x 33355`; gene-symbol var_names mapped for marker annotation (`{'mapped_to_symbol_var_names': 32366, 'fallback_to_gene_id': 989, 'unmapped_gene_ids': 976, 'duplicate_symbols_fell_back': 13}`).
- Pipeline output: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458` with `final_adata.h5ad`, `run_manifest.json`, and `module_status.csv` present.
- Module status: `cellranger -> qc -> doublet_detection -> clustering -> annotation -> differential_expression -> composition` all `ok`.
- Final object after QC/doublet filtering: `55653 x 25519`, `24` Leiden clusters; top annotations: Excitatory neuron=36503, Radial glia / neural progenitor=6973, Inhibitory interneuron=4969, Cycling progenitor=3020, Intermediate progenitor=1990.
- Truth boundary: no public FASTQ/SRA raw-read source found yet; this lane is GEO public processed scRNA matrix-level reproduction, with velocity layers deferred to a separate loader audit.

---


# Latest Wave-5.1 Final Completion Update (2026-05-17)

- Wave-5.1 Trevino Cell reproduction evidence package is complete with honest MV3 deferral.
- Final audit: `.omx/reports/wave5-v51-final-completion-audit-20260517.json`.
- v5.1 ledger: `ops/run_ledger/wave5_trevino_20260517T1005Z-7b539a5.v5.1.json` validates against `wave5_v5_1.schema.json`; v4.2/v5.0 schemas reject it.
- v5.1 comparison PDF: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1005Z-7b539a5/python/figures/v5/comparison_3fa006c876fae372_v5.1.pdf`, SHA `3fa006c876fae37289661bc5875766b76bff25634bd1346e0276c0158b239d6d`.
- MV3 F-3/F-6/F-7 status: `HONEST-DEFERRED` at `python/figures/v5/mv3/mv3_disposition.json`; no new dependencies or external GWAS/motif/genome resources were installed or fabricated.
- Final targeted gate: `75 passed, 11 deselected in 1.23s`; contract parity OK; MV2 ledger status is `PASS-CONTRACT-ONLY` (not full analytical R parity).

---

# Latest Wave-5.1 HV1 Update (2026-05-17)

- HV1 real Layer-3 metrics completed for run `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1005Z-7b539a5/`.
- Outputs: `python/figures/v5/quantitative_metrics/*.json` (`15` JSON files) and `python/figures/v5/reproduction_rate.json`.
- Real-metric figure count: `10` main/RNA-only figures with `n_scored_metrics > 0`; AC-V51-HV1-1 PASS (threshold `>=8`).
- Reproduction rate: v5.1 `0.4031615505603515` vs v5.0.4 baseline `0.338`, delta `+0.0651615505603515`; AC-V51-HV1-3 PASS.
- GEO SHA verification: all local GSE162170 files OK; manifest at `data/external/trevino_2021/geo/MANIFEST.json`.
- ATAC peak Spearman remains honest `NOT-COMPUTED`: the public ATAC counts file lacks cell-barcode headers in this local copy, so F-5 is not used to inflate the headline rate.
- Post-HV1 audit: `.omx/reports/wave5-v51-cell-mainline-gap-audit-20260517-post-hv1.json`; remaining mainline gaps are MV3 disposition, comparison_v5.1 PDF, and v5.1 ledger.

---
# Latest Wave-5.1 Trevino Cell Reproduction Update (2026-05-17)

- Current v5.1 evidence run: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1005Z-7b539a5/`.
- HV2 WNN pre-registered sweep completed: `288/288` OK rows, selected `resolution=0.5`, `n_neighbors=20`, `rna_weighted_2x`, `mutual`; ARI `0.5202309143816428`; selected config at `python/figures/v5/rna_only/selected_config.json`.
- MV1a processed-count replay completed: `python/final_adata.h5ad`, `python/run_manifest.json`, and `python/module_status.csv` exist under the project root; parity vs v4.2 PASS at `python/parity/parity_vs_v4_2.json`.
- R bundle exported under the project run at `python/bundle/` with `singlecell_r_bundle_v2.2` and `r_factory_sha_at_export=cba6162`.
- Operational caution: this adopted run id is legacy/no-hyphen (`20260517T1005Z-7b539a5`) and does not satisfy newer `workflow.modular.project_paths` run-id regex; direct project-run paths were used for bundle export.
- New journal: `ops/before_every_run/journal/2026-05-17-wave5-v51-hv2-mv1a.md`.

---

# BeforeEveryRun Latest Summary

## Current Canonical Execution State

- Canonical prepared input:
  `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- Retained fresh stage-1 evidence run:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
- Current clean full-cohort rerun source of truth:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652`
- Current clean rerun launch log:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO.launch.log`
- Current clean rerun local package:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase5_clean_full_cohort_rerun_20260424`
- Current human-facing workspace entrance:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/00_HUMAN_START_HERE.md`
- Current agent-facing run ledger:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/agent_runs/README.md`

## Current Methodology / Governance State

- As of the `2026-04-25` methodology audit, the next paper-aligned optimization
  lane is sparse-exact rather than CSS approximation:
  `/home/zerlinshen/singlecell_factory/ops/nc2024_methodology_audit/AUDIT_2026-04-25.md`
- That audit records paper-faithful Scrublet, Harmony 15-PC clustering, Leiden
  resolution `1.0`, paper-aligned DE parameters, and tumor-vs-background /
  healthy cohort splitting through `--cohort-subset`.
- Paper-aligned launcher:
  `/home/zerlinshen/singlecell_factory/scripts/run_nc2024_paper_aligned_20260425.sh`
- Sparse-exact exploratory launcher:
  `/home/zerlinshen/singlecell_factory/scripts/run_nc2024_full_cohort_sparse_exact_20260425.sh`
- Observed `2026-04-25` `NC2024_NSCLC_FULL_COHORT_SPARSE_EXACT_REAL_AUTO_*`
  result directories are not canonical successful runs unless a later audit
  finds `final_adata.h5ad`, `run_manifest.json`, and `module_status.csv`.
  Current read-only inspection saw only early mandatory outputs/checkpoints and
  a zero-byte sparse-exact launch log.
- Both Mac-led SSH orchestration and direct remote operation are valid. The Mac
  remains the review/organization/report-packaging surface; the remote repo
  remains the compute, remote-R, and run-truth surface.
- As of the `2026-04-30` architecture contract round, module hierarchy metadata
  lives in `/home/zerlinshen/singlecell_factory/workflow/modular/module_catalog.py`.
  `workflow/modular/pipeline.py` derives the legacy dependency DAG from that
  catalog, and CLI optional-module help derives from the same source.
- No full-cohort rerun was launched in the architecture contract round.
  Controller-validation smoke now lives at:
  `/home/zerlinshen/singlecell_factory/scripts/validate_nc2024_architecture_contract.py`
  It checks current NC2024 v2 run dirs, bundle manifests, run-ledger entries,
  and `bridges/local_r_pipeline_macbook/{R,R_bundle}` symlinks without opening
  the large H5ADs.

## Latest Run Verdict

- Date: `2026-04-24`
- Execution mode: `debug_massive`
- Scale mode: `massive`
- Checkpoint policy: `--scale-mode massive` (sets `mandatory_only` automatically)
- Prepare rerun: `no`
- Input: canonical prepared Zarr above.
- Result: clean rerun wrote `final_adata.h5ad`, `run_manifest.json`,
  `module_status.csv`, and `module_reconciliation.tsv`.
- Final object: `810218 x 30374`, `33G`.
- Module status: all requested modules `ok`.
- `pseudobulk_de`: `ok/completed` in the pipeline-native rerun.
- Reconciliation:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652/module_reconciliation.tsv`
- Remote R reporting:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652/r_plots/phase5_readable_20260424`

## Cleanup Truth

- User approved deleting old `results/` artifacts when they are reproduction/test
  outputs and each run is recorded.
- Pre-extended-run cleanup archive:
  `/home/zerlinshen/singlecell_factory/ops/cleanup_records/nc2024_pre_extended_real_run_20260424T051809Z`
- Deleted superseded result directories:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_LARGE_AUTO_20260423_035243`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222`
- These deleted directories must not be cited as currently existing canonical
  artifacts. Use the retained fresh stage-1 run and current clean rerun instead.
- Raw/reference data and the canonical prepared Zarr were not deleted.

## Most Important Carry-Forward Lessons

- Treat `large` as a capacity probe on this cohort, not the main completion lane.
- Use direct `massive` for actual full-cohort survivability, debugging, and
  recovery work.
- Use controller `large -> massive` only for orchestration validation.
- Do not rerun prepare unless the canonical prepared Zarr is invalid.
- Use `--scale-mode massive` for full-cohort runs; it automatically sets
  `checkpoint_policy=mandatory_only` to skip early-module AnnData checkpoints
  unless there is a specific need to preserve full per-module AnnData checkpoints.
- For reproduce-stage NC2024 work, interpret "all modules" as "all eligible
  modules", not "every imaginable module regardless of modality or contrast
  contract".
- Remote R plotting/reporting is remote-side by default on `ubuntu-tail` with:
  `/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript`
  and
  `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh`.
- The previous `/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript` runtime is
  retained as rollback, but it lacks R `arrow` and is not the default for v2
  parquet bundle plotting.
- The bridge folder name still says `local_r_pipeline_macbook`, but current
  operation is remote-first; do not describe the Mac-local R pipeline as active.
- Python still produces module-native QC/diagnostic plots. R is preferred only
  where it can render materially better publication-style figures from the same
  validated artifacts.
- Treat compact R bundle values as plotting/reporting handoff data, not DE-ready
  counts or new quantitative-expression evidence without full-object validation.
- `rna_velocity` without spliced/unspliced layers or loom/BAM+GTF should be
  treated as `skipped_missing_splicing_modality`, not as a surprising failure.
- `pseudobulk_de` should be treated as confirmatory only when an explicit
  contrast contract is provided; exploratory `group_vs_rest` output must be
  explicitly enabled.
- The winning implementation pattern for full-cohort pseudobulk is:
  keep lazy inputs through the general pipeline and materialize the counts layer
  only inside pseudobulk aggregation.
- Do not enable `cnv_inference`/`evolution` on the full cohort until a scale-safe
  CNV lane exists.
- Do not silently upgrade recovered modules: preserve original status files and
  cite recovery outputs separately.
- Use the local `reproduce-run-retention` skill when deciding whether a prior
  reproduce/test run can be deleted after evidence capture.
- Remote Codex / Claude Code agents should have the same Shenxin workflow skill
  coverage. Codex project skills live in `.codex/skills/` using standard Codex
  project skill management; Claude project skills live in `.claude/skills/`.
  `codex_skills/` remains only as a legacy compatibility mirror for historical
  project-local skills.
- For compact bundle handoff, current v2 parquet manifests validate against the
  new singlecell-to-multiomics contract. The `r_multiomics_arrow` environment can
  source `R_bundle/io_bundle.R`, validate NC2024 v2 manifests, read both current
  v2 bundles through `read_bundle()`, and plot them through
  `scripts/plot_remote_bundle_large.R`.
- `bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh` is now
  v2-aware: it reuses `bundle_manifest.json` parquet bundles when source
  `final_adata.h5ad` size/mtime and requested markers/obs/obsm are compatible,
  then delegates plotting to
  `/home/zerlinshen/r_multiomics_factory/scripts/plot_remote_bundle_large.R`.
- Remote report/figure outputs intended for Mac-side review must be pulled back
  as lightweight review copies under
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/human_review/figures/YYYY-MM-DD-<slug>/`
  with `README.md`, `TRANSFER_MANIFEST.tsv`, `checksums.sha256`, and provenance
  logs. Do not pull bulky primary compute objects by default.

## Current Scientific/Reporting State

- Phase-1 fresh massive reproduction package exists and remains the baseline
  package for the eight-module stage-1 reproduction.
- Phase-3 CSS fidelity audit remains the current controlled estimate of
  clean-vs-CSS precision:
  - clean-vs-CSS ARI `0.4190`, NMI `0.6540`
  - CSS-vs-full-massive Leiden ARI `0.5479`, NMI `0.7220`
  - exact cluster identity is weaker than broad biological interpretation.
- The clean rerun now reproduces the same broad module set but with
  `pseudobulk_de` pipeline-native success rather than a separate recovery-only
  artifact.
- Phase-5 package now exists for the clean rerun.
- The Phase-4 package now has a readable deck-style PDF plus a preserved
  `first_pass_unreadable` backup inside the same package.
- The Phase-4 package also now carries a second-pass readable R rerender with
  larger UMAP exports and split marker-dot panels.
- ELF3 should currently be framed primarily as tumor epithelial / tumor-enriched.
  In recovered pseudobulk DE:
  - `Tumor epithelial_vs_rest`: log2FC `6.24528665788722`,
    padj `3.1724278254870602e-24`
  - `Myeloid/Macro_vs_rest`: log2FC `-1.1154218654114758`,
    padj `8.227643967172654e-09`

## Read Before Next Remote Run

- `journal/2026-04-23-nc2024-stage1-recovery.md`
- `journal/2026-04-24-nc2024-fresh-massive-rerun.md`
- `journal/2026-04-24-nc2024-h5ad-schema-inspection.md`
- `journal/2026-04-24-nc2024-phase1-pdf-package.md`
- `journal/2026-04-24-nc2024-checkpoint-hardening.md`
- `journal/2026-04-24-nc2024-css-fidelity-audit.md`
- `journal/2026-04-24-nc2024-r-bridge-performance-hardening.md`
- `journal/2026-04-24-nc2024-remote-r-governance-correction.md`
- `journal/2026-04-24-nc2024-r-environment-hardening.md`
- `journal/2026-04-24-nc2024-bridge-code-review-hardening.md`
- `journal/2026-04-24-nc2024-artifact-inventory-cleanup.md`
- `journal/2026-04-24-nc2024-extended-full-cohort-real-run.md`
- `journal/2026-04-24-nc2024-phase4-report-pseudobulk-contract-hardening.md`
- `journal/2026-04-24-nc2024-clean-rerun-three-pass-resolution.md`
- `journal/2026-04-30-nc2024-architecture-contract-optimization.md`
- `journal/2026-04-30-nc2024-r-arrow-environment-validation.md`
- `journal/2026-04-30-nc2024-remote-to-mac-handoff-protocol.md`

## Current Open Operational Improvement

- Promote post-run pseudobulk recovery into a first-class pipeline recovery
  surface so future manifests can represent both original failure and recovered
  evidence cleanly.
- Add scale-safe CNV/evolution implementation before attempting full-cohort
  CNV/evolution.
- Automate the new remote-to-Mac handoff manifest/checksum generation if this
  pattern repeats across multiple figure/report packages.
- Add direct subprocess test coverage for the remote R wrapper.
- Add a first-class recovery manifest so pipeline failure and recovery can be
  linked more formally than a separate recovery directory.
- Add direct subprocess test coverage for the `r_multiomics_arrow` remote R
  wrapper path.

## Latest Hook-Confirmed Success

- Timestamp: `2026-04-24T20:15:00+08:00` final clean rerun artifacts observed.
- Success evidence:
  - current clean rerun final H5AD exists and is `33G`
  - current clean rerun `run_manifest.json` exists
  - current clean rerun `module_status.csv` exists
  - current clean rerun `module_reconciliation.tsv` exists
  - current clean rerun readable R report directory exists
  - current local Phase-5 package exists


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 06:05 CST

- Mainline remains honest reproduction of the Trevino Cell paper figures from the earliest public/computable inputs available locally; not the NC2024/NSCLC line.
- Latest completed panel: Figure 2E GPC enrichment, status `GENERATED_LOCAL_GPC_GENESET_ENRICHMENT_NOT_TOPGO_PARITY`.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2e_public_local_gpc_enrichment.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2E_gpc_enrichment`.
- Metrics: count-aligned GPCs 185 / paper 185; top TF enrichment `GO_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY` BH q=2.1e-10.
- Fresh verification: py_compile passed for Figure 2B-2E scripts; `git diff --check` passed; Trevino loader tests `8 passed in 0.78s`.
- Boundary: Figure 2E is local Fisher/BH enrichment from public resources, not exact topGO/Table-S2 parity. Continue with Figure 2F-2I multiome resource audit/panels.


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 06:20 CST

- Figure 2F/2G/2H/2I public multiome evidence generated and verified.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure2fghi_public_multiome_evidence.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_2F_2I_multiome`.
- Statuses: `{'2F': 'REPRODUCED_PUBLIC_AUTHOR_MULTIOME_QC_SUMMARY', '2G': 'GENERATED_PUBLIC_AUTHOR_RNA_AND_BOUNDED_ATAC_MULTIOME_PROJECTION', '2H': 'METHOD_RESOURCE_GAP_MULTIOME_LINKAGE_VENN_NOT_REPRODUCED', '2I': 'GENERATED_PUBLIC_MULTIOME_GPC_CORRELATION_CORRESPONDENCE'}`.
- Fresh verification: py_compile passed for Figure 2B-2I scripts; `git diff --check` passed; Trevino loader tests `8 passed in 0.78s`.
- Boundary: Figure 2H is a resource-gap audit, not a reproduced Venn; Figure 2G ATAC is bounded marker-NN projection, not exact author peak-LSI/uwot parity.


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 06:50 CST

- Figure 4/5 public author FCM panels generated/audited and verified.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure4_5_public_author_fcm_evidence.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_4_5_fcm_public_panels`.
- Statuses: `{'4A': 'GENERATED_AUTHOR_FCM_PUBLIC_PSEUDOBULK_OVERVIEW', '4B': 'REPRODUCED_AUTHOR_FCM_MODULE_HEATMAP_FROM_CENTERS', '4C': 'REPRODUCED_AUTHOR_FCM_SELECTED_GENE_HEATMAP', '4D': 'GENERATED_AUTHOR_FCM_MODULE_SCORE_UMAPS', '4E': 'REPRODUCED_AUTHOR_FCM_CENTROID_JACCARD_NETWORK', '4F': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4G': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4H': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4I': 'IMAGE_PANEL_NOT_COMPUTATIONAL_REPRODUCTION_LOCAL_PDF_ONLY', '4J': 'REPRODUCED_AUTHOR_FCM_GENE_MEMBERSHIP_EXPRESSION', '4K': 'IMAGE_PANEL_NOT_COMPUTATIONAL_REPRODUCTION_LOCAL_PDF_ONLY', '5A': 'REPRODUCED_AUTHOR_FCM_ASTRO_GENE_MEMBERSHIP_EXPRESSION', '5B': 'GENERATED_AUTHOR_FCM_LINKED_PEAK_MOTIF_ENRICHMENT_M13_VS_M14', '5C': 'GENERATED_AUTHOR_FCM_AQP4_POSITIVE_RECLUSTERING_LOCAL_RULE', '5D': 'GENERATED_PUBLIC_FCM_PSEUDOBULK_DE_NOT_DESEQ2_PARITY', '5E': 'METHOD_RESOURCE_GAP_BHADURI_PRIMARY10X_NOT_STAGED', '5F': 'METHOD_RESOURCE_GAP_BHADURI_PRIMARY10X_NOT_STAGED'}`.
- Fresh verification: FCM render command passed; py_compile passed; `git diff --check` passed; 14 PNGs verified; Trevino loader tests `8 passed in 0.77s`.
- Boundary: 4I/4K are image-only IHC panels, not recomputed matrix outputs; 5D is a local pseudo-bulk Welch-test proxy, not DESeq2/Table S5 parity; 5E/5F need Bhaduri/UCSC `primary10X` staging before reproduction can be claimed.
- Next: Figure 6 chromatin-state/GPC branch resource audit/panels.


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 07:05 CST

- Figure 6 public FCM/GPC branch panels generated/audited and verified.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure6_public_fcm_gpc_branch_evidence.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_6_fcm_gpc_branch`.
- Statuses: `{'6A': 'REPRODUCED_PUBLIC_FCM_CELL_CYCLE_MODULE_CORRELATION', '6B': 'GENERATED_PUBLIC_FCM_ATAC_PROJECTION_SCHEMATIC', '6C': 'REPRODUCED_AUTHOR_FCM_ATAC_BRANCH_PROJECTION_FROM_GA_ANCHORS', '6D': 'GENERATED_PUBLIC_FCM_BRANCH_GENE_ACTIVITY_HEATMAP_LOCAL_SPECIFICITY', '6E': 'GENERATED_PUBLIC_FCM_GPC_TF_MOTIF_GENE_BRANCH_HEATMAP', '6F': 'REPRODUCED_AUTHOR_FCM_GPC_ANCHOR_REPROJECTION_BRANCH_ARROWS', '6G': 'GENERATED_BOUNDED_MULTIOME_RNA_TO_FCM_NN_PROJECTION_COLORED_BY_ATAC_CLUSTER'}`.
- Fresh verification: Figure 6 render command passed; py_compile passed; `git diff --check` passed; 7 PNGs verified; Trevino loader tests `8 passed in 0.85s`.
- Boundary: 6D/6E use local public-branch specificity/means; 6G is bounded NN multiome RNA→FCM projection, not exact author uwot parity.
- Next: Figure 7 disease/BPNet/ASD resource audit.


## Current Wave5 Cell/Trevino all-figure reproduction update — 2026-05-18 07:06 CST

- Figure 7 public S6/BPNet-ASD panels generated/audited and verified. Main Trevino Cell Figures 1-7 now have public-resource evidence/outputs or explicit resource-gap audits.
- Evidence: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/figure7_public_asd_bpnet_evidence.json`; outputs: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figures/public_matrix/figure_7_asd_bpnet_public_panels`.
- Statuses: `{'7A': 'GENERATED_FROM_METHODS_AND_PUBLIC_TABLE_S6_COUNTS', '7B': 'GENERATED_S6D_CLUSTER_HIGH_EFFECT_SUMMARY_ON_PUBLIC_ATAC_UMAP_NOT_EXACT_AUTHOR_FISHER_DENOMINATORS', '7C': 'GENERATED_S6_COUNT_AUDIT_EXACT_OR_1_909_DENOMINATORS_NOT_PUBLIC_IN_S6', '7D': 'GENERATED_SFARI_AUDIT_CAPTION_BENCHMARK_PLUS_RECOMPUTED_COUNTS_NOT_EXACT_24_17_FROM_LOCAL_TABLES', '7E': 'REPRODUCED_FROM_PUBLIC_TABLE_S6E_MOTIF_OVERLAP_EXCESS', '7F': 'GENERATED_NFIA_S6_LOCUS_SUMMARY_NOT_EXACT_BPNET_LOGO_TRACK_PARITY', '7G': 'GENERATED_NPY_S6_LOCUS_SUMMARY_NOT_EXACT_BPNET_LOGO_TRACK_PARITY'}`.
- Fresh verification: Figure 7 render command passed; `python -m py_compile scripts/dev/wave5_render_public_figure7_asd_bpnet.py` passed; `git diff --check` passed; 7 PNGs verified; Trevino loader tests `8 passed in 0.78s`.
- Boundary: S6D duplicated `Cluster GluN6` yields 261/232 caption-like high-effect case/control vs paper 262/232; exact Figure 7B/7C Fisher denominators and 7F/7G BPNet logo/ref-alt tracks are not deposited in public S6/local Brain_ASD checkout.
- Next: final quality gate/code-review and consolidated completion handoff; no additional long compute phase is currently queued.


## Current Wave5 Cell/Trevino all-figure reproduction completion — 2026-05-18 07:13 CST

- Completed conditional public-resource all-main-figure reproduction/audit.
- Completion report: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/TREVINO_ALL_FIGURES_COMPLETION_REPORT.md`.
- Final quality gate: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88/python/figure_reproduction_evidence/all_figures_final_quality_gate_review.json`; code review `APPROVE/CLEAR`; scientific decision `conditional` due documented public-resource gaps.
- Evidence audit: 17 JSON files, 53 panel statuses, 0 missing referenced paths, 0 hash mismatches, 0 unexpected false checks.
- Pipeline viability: public RNA pipeline run has all modules `ok` and final AnnData under `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/wave5_trevino_public_rna_20260518_040458/final_adata.h5ad`.
- Boundary: not FASTQ/raw-fragment parity; Figure 7 exact BPNet denominators/model tracks and image-only panels remain honest resource gaps.
