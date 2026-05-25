# Yost Exact-Panel Heatmap/Contact Package - 2026-05-23

## Scope

- Task: continue Yost 2025 Nature Genetics exact-panel figure reproduction from
  `/home/zerlinshen/projects/ng2025-3d-genome/NEXT_AGENT_YOST_EXACT_PANEL_REPRODUCTION.md`.
- Execution mode: `controller_validation` for plotting/reporting only.
- Project root: `/home/zerlinshen/projects/ng2025-3d-genome`.
- Output package:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T001058`.

## Objective

Replace the prior author-code-informed summary/proxy figure lane with a
panel-status-governed package for supported Yost Fig. 1 heatmaps/contact/track
panels, while keeping STK11/LKB1 separate from original paper-panel claims.

## What Was Attempted

- Added renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_exact_panel_heatmaps.R`.
- Used the author R heatmap route from `Plot_heatmaps.R`:
  sample-sample Pearson correlations, pheatmap dendrogram ordering,
  TCGA cancer-type annotations, and Yost TCGA palettes.
- Used the author track/contact route from `plotting_functions.R` where
  available: `strawr` `.hic` extraction, stacked bigWig tracks, and loop arcs.
- Copied original paper figure crops into `00_original_references/` and created
  a Fig. 1 comparison sheet.

## What Succeeded

- Produced exact/status-governed outputs:
  - Fig. 1f H3K27ac 1D-signal sample-correlation heatmap:
    `EXACT_AUTHOR_CODE_DATA_AVAILABLE`.
  - Fig. 1f FitHiChIP loop-interaction sample-correlation heatmap:
    `EXACT_AUTHOR_CODE_DATA_AVAILABLE`.
  - Fig. 1f LUAD/LUSC CALDER comp-rank heatmap: `LUNG_SUBSET_EXACT`.
  - Fig. 1c LUAD/LUSC chr8 HiChIP contact maps: `LUNG_SUBSET_EXACT`.
  - Fig. 1d/e-style MYC/PVT1 tracks/loop arcs: `AUTHOR_CODE_ADAPTED`.
  - STK11/LKB1 track/loop extension: `STK11_EXTENSION_NOT_ORIGINAL_PANEL`,
    under `04_STK11_LKB1_extension/`.
- Marked full-cohort CALDER Fig. 1f as `UNSUPPORTED_INPUT_GAP` because only
  LUAD/LUSC CALDER outputs are currently confirmed.
- Wrote `PANEL_STATUS.tsv`, `PANEL_SOURCE_MAP.tsv`, `COMPARISON_MAP.tsv`,
  `QC_REPORT.md`, and `logs/QC_METRICS.tsv`.

## What Failed Or Needed Recovery

- First render produced contact-map `geom_raster` uneven-interval warnings.
  The renderer was revised to use `geom_tile` for contact matrices and rerun.
- An intermediate status table included the comparison sheet as a panel row.
  The renderer was revised so comparison sheets are recorded in
  `COMPARISON_MAP.tsv`, not `PANEL_STATUS.tsv`.

## Verification

- Final renderer completed:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T001058`.
- Output counts: 18 PNG files including original references, 7 SVG, 7 PDF, and
  7 TIFF outputs.
- R PNG QA passed for all 18 PNG files; no blank PNG detected.
- Filename contamination check found 0 Miao/Qiaowei outputs.
- STK11 isolation check: 0 STK11 files under `01_exact_reproduced_panels` and
  `02_lung_subset_exact_panels`; 4 files under `04_STK11_LKB1_extension`.
- `git diff --check` passed for
  `plotting_factory/r/report/render_yost_exact_panel_heatmaps.R`.

## Artifact Classification

- Canonical:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T001058`.
- Evidence-only / provenance:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T000054`
  and
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T000529`.
- Superseded proxy:
  `/home/zerlinshen/projects/reproductions/figure_packages/3d_lung_yost_author_code_figures_20260522`.
- Canonical renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_exact_panel_heatmaps.R`.

## Remaining Risks

- Full-paper-scale CALDER subcompartment Fig. 1f is still not exact because
  all-cancer CALDER `all_sub_compartments.tsv` files are not staged.
- The MYC/PVT1 track/loop panel is author-code adapted, not exact EIS parity,
  because EIS-specific input/state is not claimed.
- Fig. 2-5 exact panel work remains future scope and should continue from the
  same status vocabulary rather than reviving summary/proxy figures.

## Next Operator Notes

- Do not use the old author-code-informed package as proof of panel-exact
  reproduction.
- Use `PANEL_STATUS.tsv` and `PANEL_SOURCE_MAP.tsv` as the source of truth for
  what is exact, lung-subset exact, adapted, input-gapped, or STK11 extension.
- Keep STK11/LKB1 extension outputs separate unless a later retention decision
  explicitly reclassifies them.
- When full-cohort CALDER outputs are staged, rerun only the CALDER heatmap row
  or rerun the renderer with the expanded CALDER source list.
