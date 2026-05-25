# Yost STK11/LKB1 Individual Single-Script Figure Package

Date/time: 2026-05-23

## Objective

Correct the STK11/LKB1 storyline package so every final figure is generated as
its own data-driven figure, not copied from an older rendered package and not
cropped from a larger generated figure.

## What Was Attempted

- Replaced the old storyline handoff route with
  `render_yost_stk11_storyline_individual_figures.R`.
- Added one R entrypoint per final output figure under
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/yost_stk11_individual/`.
- Added a real-data NSCLC support-table extractor at
  `/home/zerlinshen/projects/ng2025-3d-genome/scripts/extract_yost_stk11_nsclc_singlecell_support.py`.
- Re-rendered 10 existing storyline figures plus 4 NSCLC support/audit figures.

## What Succeeded

- Canonical remote package:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.
- Output count: 14 PNG, 14 SVG, 14 PDF, and 14 TIFF files.
- `SCRIPT_MANIFEST.tsv`: all 14 scripts exited 0.
- `PNG_QA.tsv`: all 14 PNGs are nonblank and high resolution.
- `FIGURE_STATUS.tsv`: all rows record
  `independent_data_render_no_copy_no_crop`.
- Visual spot checks confirmed that Fig. 1f titles are complete on separate
  lines, STK11/LKB1 tracks/arcs has no in-image extension title, and the
  single-cell/metadata-audit figures are legible.

## What Failed Or Was Corrected

- The old package
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_single_figures_20260523T010004`
  was usable visually but did not satisfy the stricter one-script-per-figure
  rule because it copied exact-panel rendered outputs into the storyline
  package.
- Initial remote driver execution needed a path-with-spaces fix; the driver now
  runs individual scripts from the script directory.
- Initial heatmap annotation colors were narrowed too much; the common palette
  now retains the full TCGA cancer-type palette needed by the all-cancer Fig. 1f
  heatmaps.

## What Changed

- New current renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_stk11_storyline_individual_figures.R`.
- New individual script directory:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/yost_stk11_individual/`.
- Deprecated compatibility wrapper:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_stk11_storyline_single_figures.R`.
- Updated project README and retention policy to name the individual package as
  the current handoff and to mark the old 10-figure package as superseded.

## Artifact Classification

- Canonical: `yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.
- Superseded: `yost_stk11_storyline_single_figures_20260523T010004` for final
  handoff.
- Evidence-only: earlier exact-panel package remains useful for comparison and
  context, but not as a rendered-output source for this storyline handoff.
- Failed exploratory: none retained in this package. The scratch
  `debug_individual_01` test directory can be deleted if space cleanup is
  needed.

## Residual Risks

- CN and methylation tracks are still proxy summaries, not allele-specific LOH
  or validated promoter-silencing calls.
- STK11-mut versus WT LUAD comparisons are exploratory because the staged split
  is `n=2` versus `n=2`.
- WCH and Round9 single-cell resources support cell-state/checkpoint context,
  but current remote NSCLC resources do not contain audited STK11 genotype,
  therapy-response, or survival labels. Clinical efficacy remains unsupported.

## Next Operator Notes

- Use the individual renderer and `render_XX_*.R` scripts for final figure
  handoff.
- Do not resurrect rendered-output copying, contact-sheet cropping, or
  comparison-sheet delivery for the STK11/LKB1 storyline.
- Mac copies are archives only. Remote project and factory records remain the
  source of truth.
