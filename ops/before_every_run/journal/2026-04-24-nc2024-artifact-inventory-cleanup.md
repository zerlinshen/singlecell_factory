# 2026-04-24 NC2024 Artifact Inventory Cleanup

## Task

Answer whether the workspace should delete intermediates and rerun the full
cohort by first classifying current artifacts.

## Objective

- Preserve canonical evidence.
- Remove only disposable temp artifacts.
- Avoid another full-cohort run unless a concrete reproducibility gap appears.

## Execution Mode

- `execution_mode=verification_cleanup`
- No prepare, controller, benchmark, or full-cohort stage-1 run was launched.

## What Was Attempted

- Read local and remote `before-every-run/LATEST.md`.
- Inventoried remote NC2024 result directories and fresh-run `r_plots`.
- Checked exact `/tmp` evidence bundle/log paths from the bridge hardening lane.
- Verified canonical prepared input, prior successful runs, fresh
  `final_adata.h5ad`, fresh `r_bundle`, and fallback proof log still exist.

## Cleanup Performed

Deleted only disposable `/tmp` artifacts:

- `/tmp/nc2024_bridge_review_bundle_20260424` (`126M`)
- `/tmp/nc2024_bridge_custom_contract_bundle_20260424` (`100M`)
- `/tmp/nc2024_bridge_custom_contract_20260424.log` (`4.0K`)
- `/tmp/nc2024_bridge_custom_contract_reuse_20260424.log` (`4.0K`)
- `/tmp/nc2024_bridge_canonical_reuse_after_contract_fix_20260424.log` (`4.0K`)
- `/tmp/nc2024_remote_parse_after_deslop.log` (`8.0K`)

Approximate reclaimed space: `226M`.

## Artifact Classification

- `canonical`:
  - `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_WITH_FALLBACK.launch.log`
- `reusable cache`:
  - fresh-run `r_bundle`
  - fresh-run schema inspection / validation outputs
- `evidence-only preserved`:
  - cited bridge/R smoke plot directories under fresh-run `r_plots`
  - CSS fidelity and module-specific small analysis outputs
- `superseded temp deleted`:
  - exact `/tmp` paths listed above

## Validation

- All deleted temp paths were absent after cleanup.
- Canonical existence gate passed.
- No new `NC2024_NSCLC_FULL_COHORT*` result directory appeared.

## Remaining Risks

- The three 84-85G successful full-cohort directories remain intentionally
  preserved because current run memory marks them as canonical/prior-canonical
  evidence. Do not delete them without a separate archival policy.
- Some small evidence-only plot directories remain because they are cited in
  README/run memory and are cheap to keep.

## Next Operator Notes

- Do not do a full-cohort real rerun for cleanliness alone.
- If space pressure becomes severe, propose an explicit archival policy for
  older prior-canonical full-cohort successes rather than deleting them ad hoc.
- The next productive lane remains targeted ELF3 Myeloid/Macro or missed-result
  analysis against the existing fresh `final_adata.h5ad` and `r_bundle`.
