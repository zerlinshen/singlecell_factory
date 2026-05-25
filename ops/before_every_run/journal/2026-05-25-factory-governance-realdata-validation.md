# Factory Governance Real-Data Validation - 2026-05-25

## Scope

- Task: G004 of the factory-governance ultragoal, after Claude critique that
  structure-only gates were being overrepresented as scientific validation.
- Execution mode: controller_validation with bounded real-data proof.
- Host: `/home/zerlinshen` on `ubuntu-tail`.

## Failed Exploratory Run

- Attempted HGMM modular CLI rerun:
  `/home/zerlinshen/projects/hgmm-smoke/runs/2026-05-24T1930Z-0321773`.
- Command used `SC_REQUIRE_PROJECT_ROOT=1`, `--checkpoint`,
  `--parallel-workers 1`, `--gpu-mode off`, and `--allow-dirty`.
- Outcome: the process created only the project scaffold and hung for more than
  three hours with no produced files under `python/`, `python/bundle/`, `r/`, or
  `logs/`.
- Action: process terminated with SIGTERM; added a run `README.md` so the suite
  run-bundle gate treats it as an explicit failed exploratory scaffold.
- Do not cite this directory as validation evidence.

## Successful Bounded Proof

- New governed project run:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-24T2347Z-0321773`.
- Source: retained Round9 LUSC consensus OR lane:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1756Z-9907a7d/python/lusc_ps01_round9_consensus_scdbl_or_20260522_015643/final_adata.h5ad`.
- Source shape: `65669` cells by `17267` genes with `X_pca` and `X_umap`.
- Export: `singlecell_factory/scripts/export_singlecell_r_bundle.py` wrote a
  v2.1 compact bundle with `5000` deterministic cells, `17` marker genes,
  metadata, `X_pca`, and `X_umap`.
- R render: default `r_multiomics` env failed because package `arrow` was
  missing; retry with `r_multiomics_arrow` completed
  `scripts/plot_remote_bundle.R --project-root ... --run-id ...`.
- Output: `7` PNG figures and `1` PDF under the governed project run `r/`
  directory.
- Mechanical QA: `7/7` PNGs nonblank, all dimensions at least `3200x1280`.
- Evidence:
  - `runs/2026-05-24T2347Z-0321773/manifest.json`
  - `runs/2026-05-24T2347Z-0321773/evidence/round9_cross_factory_bundle_plot_validation.md`
  - `ledger/human_review/2026-05-25-cross-factory-bundle-plot-validation.md`

## Interpretation

This proof validates the current real-data cross-factory handoff from a retained
Python AnnData through the R bundle reader and plotting helpers. It is stronger
than the suite structural gate because it uses real Round9 data and current
code paths.

It remains bounded. It is not full raw-input pipeline reproduction, not NG2025
end-to-end validation, and not evidence for new quantitative biological claims.
The R environment split remains a release-readiness risk until one canonical env
can run the parquet-bundle plotting path without retry.

## Verification

- `python scripts/export_singlecell_r_bundle.py ...` completed under timeout.
- `Rscript scripts/plot_remote_bundle.R ...` completed under timeout in
  `r_multiomics_arrow`.
- PIL PNG QA passed for all `7` rendered PNGs.
- `python scripts/validate_run_bundle.py` passed with warnings only; the new
  Round9 validation run is governed by `manifest.json`.
- `bash scripts/run_all_gates.sh` passed `8/8`.
