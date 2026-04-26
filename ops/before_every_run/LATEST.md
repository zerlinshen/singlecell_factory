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

## Latest Run Verdict

- Date: `2026-04-24`
- Execution mode: `debug_massive`
- Scale mode: `massive`
- Checkpoint policy: `SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only`
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
- Keep `SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only` for full-cohort massive runs
  unless there is a specific need to preserve full per-module AnnData
  checkpoints.
- For reproduce-stage NC2024 work, interpret "all modules" as "all eligible
  modules", not "every imaginable module regardless of modality or contrast
  contract".
- Remote R plotting/reporting is remote-side by default on `ubuntu-tail` with:
  `/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript`
  and
  `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh`.
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

## Current Open Operational Improvement

- Promote post-run pseudobulk recovery into a first-class pipeline recovery
  surface so future manifests can represent both original failure and recovered
  evidence cleanly.
- Add scale-safe CNV/evolution implementation before attempting full-cohort
  CNV/evolution.
- Add package-level provenance checks for synced remote figures/tables.
- Add direct subprocess test coverage for the remote R wrapper.
- Add a first-class recovery manifest so pipeline failure and recovery can be
  linked more formally than a separate recovery directory.

## Latest Hook-Confirmed Success

- Timestamp: `2026-04-24T20:15:00+08:00` final clean rerun artifacts observed.
- Success evidence:
  - current clean rerun final H5AD exists and is `33G`
  - current clean rerun `run_manifest.json` exists
  - current clean rerun `module_status.csv` exists
  - current clean rerun `module_reconciliation.tsv` exists
  - current clean rerun readable R report directory exists
  - current local Phase-5 package exists
