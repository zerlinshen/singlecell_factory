# 2026-04-24 - NC2024 extended full-cohort real run

## Objective

- Run one approved full-cohort real lane with all eligible `singlecell_factory`
  and remote R-reporting modules after storage cleanup.
- Do not rerun prepare. Use the canonical prepared Zarr as the input.
- Record cleanup and every run artifact so deleting superseded `results/`
  directories does not compromise reproducibility.

## Starting Context

- Host: `ubuntu-tail`
- Remote working directory: `/home/zerlinshen/singlecell_factory`
- Canonical prepared input:
  `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- Retained fresh stage-1 evidence run:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
- Carry-forward lessons:
  - `large` is a capacity probe for this cohort, not the completion lane.
  - Use `execution_mode=debug_massive` for recovery/real-run work.
  - Avoid dense AnnData checkpoints on full-cohort massive runs.
  - Remote R plotting is the active reporting surface; the local Mac is for
    organization, review, comparison, and PDF packaging.

## What Was Run

- Cleanup first, then one direct `massive` extended full-cohort run.
- Execution mode: `debug_massive`
- Scale mode: `massive`
- Checkpoint policy: `SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only`
- Launch script:
  `/home/zerlinshen/singlecell_factory/scripts/run_nc2024_extended_full_cohort_20260424.sh`
- Launch log:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO.launch.log`
- Run directory:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003`
- Requested optional modules:
  `clustering,cell_cycle,batch_correction,differential_expression,annotation,trajectory,pseudo_velocity,pathway_analysis,cell_communication,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,pseudobulk_de,cell_fate,composition,metacell,paper_repro,validate_cbioportal`
- Deferred modules:
  - `rna_velocity`: no spliced/unspliced layers and no loom/BAM/GTF source.
  - `cnv_inference`: current full-cohort implementation can densify CNV
    matrices; defer to a scale-safe CNV lane.
  - `evolution`: depends on CNV and was deferred with CNV.

## Outcome

- `partial-success with recovered pseudobulk`
- The run wrote a full final object and all expected core manifests.
- `module_status.csv` has `20 ok` rows and one original pipeline-level failed row:
  `pseudobulk_de`.
- `pseudobulk_de` was recovered post-run from `final_adata.h5ad`.
- Remote R plotting bundle/report was generated.
- Local Phase-4 PDF reproduction package was generated.

## Evidence

- Final object:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/final_adata.h5ad`
  - size: `33G`
  - shape: `810218 x 30374`
- Manifest:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/run_manifest.json`
- Module status:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/module_status.csv`
- Pseudobulk recovery:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/pseudobulk_de_recovery`
- Remote R report:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/r_plots/extended_real_run_main_20260424`
- Local Phase-4 package:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase4_extended_full_cohort_real_run_20260424`
- Local agent run ledger:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/agent_runs/2026-04-24-full-cohort-extended-real-run/RUN.md`

## Problems Encountered

- Storage pressure existed before the run because old full-cohort result
  directories and bulky per-module checkpoints occupied too much space.
- The original pipeline `pseudobulk_de` row failed with:
  `_cs_matrix.sum() got an unexpected keyword argument 'keepdims'`.
- `pathway_analysis` fell back because the installed `decoupler` did not expose
  the expected PROGENy helper.
- `composition` fell back to chi-squared because pertpy/scCODA was unavailable.
- `metacell` fell back to MiniBatchKMeans because SEACells was unavailable.

## What Was Changed

- Storage cleanup archived metadata and deleted superseded old result
  directories:
  - cleanup archive:
    `/home/zerlinshen/singlecell_factory/ops/cleanup_records/nc2024_pre_extended_real_run_20260424T051809Z`
  - deleted:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_LARGE_AUTO_20260423_035243`
  - deleted:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552`
  - deleted:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222`
- `workflow/modular/context.py` now supports
  `SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only`, which keeps JSON checkpoint
  sidecars without writing massive per-module AnnData checkpoints.
- `workflow/modular/modules/pseudobulk_de.py` was patched for current
  `pydeseq2` contrast API drift while preserving fallback behavior.
- Remote README was updated to describe the current remote R-first reporting
  contract, cleanup policy, extended run, and recovered pseudobulk boundary.
- Local reproduction README and Phase-4 package were updated with the new
  evidence package.

## Resolution

- Disk was freed without deleting raw/reference data or the retained fresh
  stage-1 run.
- The extended run completed far enough to write a full final AnnData object.
- A synthetic `pydeseq2` smoke test passed after patching the contrast API.
- Post-run pseudobulk recovery completed from `final_adata.h5ad` in about
  `149` seconds and produced recovered DE outputs.
- Remote R reporting produced PNG/PDF outputs from a compact bundle rather than
  forcing full-object Seurat conversion.

## Cautions for the Next Run

- Do not treat the old 2026-04-23 direct/controller result directories as
  available; they were deliberately deleted after archival.
- Do not rewrite the original `module_status.csv` to hide the pipeline-level
  `pseudobulk_de` failure. Use the recovery directory for recovered evidence.
- Do not rerun prepare unless the canonical prepared Zarr is invalid.
- Do not use `large` as the completion lane for this cohort.
- Do not enable CNV/evolution on the full cohort until the CNV implementation is
  made scale-safe.
- If new modules or more full-cohort reruns are needed, keep metadata-only
  checkpoint policy enabled unless there is a concrete reason to preserve full
  intermediate AnnData objects.

## Improvement Ideas

- Promote post-run recovery modules into the official pipeline wrapper so the
  manifest can record both original failure and recovered evidence explicitly.
- Add a scale-safe CNV/evolution lane that avoids cell-by-gene dense matrices.
- Add direct pytest/subprocess coverage for the remote R wrapper instead of only
  content assertions plus shell smoke execution.
- Add package-level provenance checks for synced remote figures and tables.

## Classification

- `canonical`
  - canonical prepared Zarr:
    `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
  - current extended full-cohort run:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003`
  - retained fresh stage-1 evidence run:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
  - Phase-4 local package:
    `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase4_extended_full_cohort_real_run_20260424`
- `evidence-only`
  - cleanup archive:
    `/home/zerlinshen/singlecell_factory/ops/cleanup_records/nc2024_pre_extended_real_run_20260424T051809Z`
  - launch log:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO.launch.log`
- `superseded`
  - deleted 2026-04-23 stage-1 large/direct/controller result directories listed
    above.
- `failed exploratory`
  - none newly preserved as a failed run directory; the original pseudobulk
    failure is preserved in the canonical extended run status file.
