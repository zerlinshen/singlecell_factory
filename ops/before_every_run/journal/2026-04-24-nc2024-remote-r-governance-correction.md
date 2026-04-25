# 2026-04-24 NC2024 Remote R Governance Correction

## Task

Correct the NC2024 single-cell R plotting contract and workspace governance
after confirming that local R plotting has been removed and active plotting is
remote-side.

## Objective

- Stop describing local/Mac R plotting as the current bridge.
- Verify the remote R bundle plotting wrapper.
- Add a result-first human entrance and a separate agent run ledger.
- Preserve all existing canonical artifacts without moving or deleting them.

## What Was Attempted

- Updated remote `singlecell_factory` README and bridge docs to state that
  NC2024-scale R plotting/reporting is remote-side by default.
- Added remote wrapper:
  `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh`
- Reworked legacy wrapper:
  `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/pull_and_plot_remote_result.sh`
  to delegate to remote plotting rather than pulling data for local R.
- Ran remote smoke plotting against the current full-cohort `r_bundle`.
- Added local governance entries:
  `00_HUMAN_START_HERE.md`, `human_review/`, and `agent_runs/`.
- Added a reusable skill:
  `/Users/zerlinshen/.codex/skills/reproduction-workspace-governance`

## What Succeeded

- Remote R runtime exists and is usable:
  `/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript`
- Remote R parse succeeded for:
  `R/remote_bundle_manifest.R`
  and `scripts/plot_remote_bundle_large.R`
- Remote smoke output succeeded:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_after_patch`
- Standalone remote bundle integrity validation succeeded with byte-size and
  SHA256 checks against `obs`, `X_umap`, `X_pca`, and `marker_expr`.
- The large plotting smoke then succeeded on the files consumed by
  `plot_remote_bundle_large.R` (`obs`, `X_umap`, and `marker_expr`):
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_integrity`
- Negative integrity test succeeded: a temporary fake bundle with an incorrect
  `obs` SHA256 failed with `Bundle file obs SHA256 mismatch`.
- Safe bundle reuse succeeded after stripping CR characters from manifest TSV
  values during shell comparisons:
  `Reusing validated R bundle: .../r_bundle`
  followed by:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_reuse_fixed`
- Corrupted-bundle fallback succeeded after fixing cwd coupling in the wrapper:
  a deliberately damaged temp bundle failed with a byte-size mismatch, triggered
  re-export, and then plotted successfully:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_fallback_reexport`
- No new full-cohort stage-1 result directory was created.

## What Failed Or Was Corrected

- Initial smoke revealed that `validate_bundle_cell_alignment()` returned all
  cell IDs at top level, causing huge console output.
- Fix:
  `validate_bundle_cell_alignment()` now returns `invisible(ref)`.
- Marker dot plot subtitles now state whether the plot is sampled and how many
  cells were shown, preventing visual previews from being over-read as full
  statistics.
- `R/remote_bundle_manifest.R` now validates manifest-backed byte sizes and
  SHA256 hashes before plotting.
- `scripts/run_remote_bundle_plot.sh` now reuses an existing bundle only when
  source paths match, the manifest is newer than `final_adata.h5ad`, and the R
  integrity validator passes. Set `FORCE_R_BUNDLE_EXPORT=1` to force a fresh export.
- The wrapper's integrity-fail fallback now runs validation in a subshell and
  re-enters `$FACTORY_ROOT` before exporting, so a failed reuse check refreshes
  the bundle instead of looking for the exporter under the bridge directory.

## Artifact Classification

- Canonical:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_bundle`
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh`
- Evidence-only:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_after_patch`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_integrity`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_reuse_fixed`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_fallback_reexport`
- Superseded:
  - older language that describes local/Mac R plotting as the active default
- Failed exploratory:
  - none requiring cleanup

## Remaining Risks

- The bridge folder name `local_r_pipeline_macbook` is historical and can still
  confuse new agents; current docs explicitly call this out.
- The governance pass added safe indexes and run records, but did not physically
  migrate old artifacts. A destructive cleanup or migration should be a separate
  planned lane.
- Remote README/script updates and remote R smoke completed before a later SSH
  connectivity timeout.
- After SSH recovered, remote `ops/before_every_run/LATEST.md` and this journal
  entry were synced back to:
  `/home/zerlinshen/singlecell_factory/ops/before_every_run/`

## Next Operator Notes

- Use remote R plotting by default for NC2024 single-cell factory outputs.
- Do not rerun full cohort just to refresh plot/report packaging.
- Start humans at `00_HUMAN_START_HERE.md` and agents at `agent_runs/README.md`
  plus `before-every-run/LATEST.md`.

## Ralph Re-verification: 2026-04-24T03:51:20Z

- Mode: `execution_mode=verification_only`; no pipeline, prepare, or full-cohort
  rerun was launched.
- Remote artifact existence was rechecked for the canonical prepared input,
  fresh massive run, `final_adata.h5ad`, `r_bundle`, `bundle_manifest.json`,
  remote plotting wrapper, R manifest helper, and remote `LATEST.md`.
- `bash -n` passed for:
  `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh`
- R-side integrity validation passed for:
  `obs`, `X_umap`, `X_pca`, and `marker_expr` via
  `validate_remote_bundle(..., allow_missing_manifest=FALSE)`.
- Wrapper inspection confirmed:
  `FORCE_R_BUNDLE_EXPORT`, CR-stripped manifest TSV reads, subshell R integrity
  validation, and `cd "$FACTORY_ROOT"` before fallback export.
- Existing evidence-only plot directories were rechecked and still contain the
  expected PNG/PDF outputs:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_integrity`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_reuse_fixed`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/governance_smoke_fallback_reexport`
- Full-cohort result directory listing remained unchanged; no new
  `NC2024_NSCLC_FULL_COHORT_*` directory appeared.
