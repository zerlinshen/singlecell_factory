# NC2024 Pre-Clean Rerun Cleanup Archive

Timestamp:
- `2026-04-24T09:53:18Z`

Purpose:
- record the pre-rerun cleanup boundary before launching the clean
  full-cohort rerun
- preserve current canonical extended-run provenance
- permit deletion only of superseded scratch/transient artifacts

Current canonical run retained:
- `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003`

Protected keep-list:
- current canonical run dir above
- current launch log:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO.launch.log`
- current `pseudobulk_de_recovery/`
- current readable rerender:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/r_plots/extended_real_run_readable_20260424`
- current local Phase-4 package
- canonical prepared input:
  `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`

Delete-eligible before rerun:
- stale PID file for the superseded prior project name:
  `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO.pid`

Not delete-eligible before rerun:
- the current canonical extended run itself
- raw data
- prepared canonical input
- current local package

Archive contents:
- this README
- `current_canonical_run.file_inventory.tsv`
- `deletions.tsv`
