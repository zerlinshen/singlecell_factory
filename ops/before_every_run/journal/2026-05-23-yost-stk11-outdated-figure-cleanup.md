# Yost STK11/LKB1 Outdated Figure Cleanup

Date/time: 2026-05-23

## Objective

Preserve only the final 14 Yost/STK11 individual figures and delete outdated
generated figure packages locally and remotely.

## Retained Source Of Truth

- Remote:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.
- Local archive:
  `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.

## Deleted Remote Generated Figure Packages

- `/home/zerlinshen/projects/reproductions/figure_packages/3d_lung_author_style_single_figures_20260522`
- `/home/zerlinshen/projects/reproductions/figure_packages/3d_lung_stk11_ppt_20260521`
- `/home/zerlinshen/projects/reproductions/figure_packages/3d_lung_stk11_single_figures_20260522`
- `/home/zerlinshen/projects/reproductions/figure_packages/3d_lung_yost_author_code_figures_20260522`
- `/home/zerlinshen/projects/reproductions/figure_packages/debug_individual_01`
- `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T000054`
- `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T000529`
- `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T001058`
- `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_single_figures_20260523T010004`

## Deleted Local Generated Figure Packages

- `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/3d_lung_author_style_single_figures_20260522`
- `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/3d_lung_stk11_ppt_20260521`
- `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/3d_lung_stk11_single_figures_20260522`
- `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/3d_lung_yost_author_code_figures_20260522`
- `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/yost_2025_exact_panel_heatmaps_20260523T001058`
- `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/yost_stk11_storyline_single_figures_20260523T010004`

## Verification

- Local `figures/` now contains only
  `yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.
- Remote `figure_packages/` Yost/3D-lung listing now contains only
  `yost_stk11_storyline_individual_figures_20260523T_single_script_per_figure`.
- Final package still contains 14 PNG, 14 SVG, 14 PDF, and 14 TIFF outputs.
- Final package size is 16M locally and 16M remotely.

## Preserved Boundary

Raw data, prepared inputs, source code, original references, project ledgers,
README files, and final package metadata were not deleted.

## Artifact Classification

- Canonical: final 14-figure individual package.
- Superseded and deleted: old generated figure packages listed above.
- Evidence-only preserved outside deletion scope: scripts, project records, and
  run/source tables.
- Failed exploratory deleted: `debug_individual_01`.

## Next Operator Notes

Do not expect old generated figure packages to exist. Rerun the corresponding
plotting scripts if an old exact-panel or comparison package is needed for a new
review branch.
