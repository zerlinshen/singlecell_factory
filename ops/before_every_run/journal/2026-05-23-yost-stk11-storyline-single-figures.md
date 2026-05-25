# Yost STK11/LKB1 Storyline Single-Figure Package - 2026-05-23

## Scope

- Task: execute the approved STK11/LKB1 storyline plan on the remote host, using only real NG2025/Yost project data and real run outputs.
- Execution mode: `controller_validation` for plotting/reporting only.
- Project root: `/home/zerlinshen/projects/ng2025-3d-genome`.
- Renderer: `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_stk11_storyline_single_figures.R`.
- Output package: `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_single_figures_20260523T010004`.

## Objective

Create the closest complete STK11/LKB1 storyline currently supported by existing data, while keeping every output as a standalone Yost-style figure and preserving the boundary between original Yost panel correspondences and STK11/LKB1 extension figures.

## What Was Attempted

- Added a remote R renderer in `plotting_factory/r/report/`.
- Used the reviewed exact-panel package as the source for direct Fig. 1f context panels and the STK11/LKB1 locus track/arc extension.
- Used real T10/T11 run tables for STK11 alteration, local activity, STK11/KEAP1/NFE2L2 module, downstream pathway, TF motif, and CD274 immune-context figures.
- Updated the exact-panel renderer so the CALDER panel title is not clipped and the STK11/LKB1 extension panel has no non-original title.

## What Succeeded

- Produced 10 standalone figures, each exported as PNG, SVG, PDF, and TIFF.
- Wrote `FIGURE_STATUS.tsv`, `STORY_STATUS.tsv`, `DATA_GAPS.tsv`, `PNG_QA.tsv`, and `README.md`.
- Kept original Yost titles only for direct Fig. 1f correspondences.
- Kept STK11/LKB1 extension figures titleless in the image surface and separated under `01_stk11_lkb1_storyline/`.

## Verification

- `Rscript r/report/render_yost_stk11_storyline_single_figures.R` completed on the remote host.
- Output counts: 10 PNG, 10 SVG, 10 PDF, and 10 TIFF.
- PNG QA passed for all 10 PNGs: each image had width/height greater than 500 pixels, non-zero pixel variance, and file size greater than 10 KB.
- `git diff --check` passed for:
  - `r/report/render_yost_exact_panel_heatmaps.R`
  - `r/report/render_yost_stk11_storyline_single_figures.R`
  - `r/report/README.md`
- Visual spot-check fixed the prior cropped CALDER title and removed the STK11 extension title.

## Artifact Classification

- Canonical storyline package:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_stk11_storyline_single_figures_20260523T010004`.
- Supporting exact-panel package:
  `/home/zerlinshen/projects/reproductions/figure_packages/yost_2025_exact_panel_heatmaps_20260523T001058`.
- Renderer source:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_stk11_storyline_single_figures.R`.

## Remaining Risks

- STK11 relative copy-number values remain proxy evidence, not allele-specific LOH.
- STK11 methylation values remain probe/beta proxy evidence, not validated promoter silencing.
- STK11-mut versus WT LUAD pathway expression remains exploratory because the staged LUAD comparison is `n=2` versus `n=2`.
- Bulk ATAC motif context is LUAD/LUSC context evidence, not STK11-genotype-specific proof.
- External single-cell support remains gap-labeled until lung/NSCLC metadata and STK11 genotype or suitable cell-state metadata are audited.

## Next Operator Notes

- Do not treat the storyline package as paper-panel reproduction. It is an extension story package anchored by exact-panel context and real STK11/LKB1 run outputs.
- Do not merge STK11/LKB1 extension figures into original-panel reproduction folders.
- If suitable external lung single-cell data become available, add it as a separate audited figure and keep genotype claims out unless genotype metadata are present.
