# Yost author-code-informed figure package - 2026-05-22

## Scope

- Replaced the mixed lung 3D-genome figure-generation lane with a Yost-only
  author-code-informed renderer.
- Current canonical remote package:
  `/home/zerlinshen/projects/reproductions/figure_packages/3d_lung_yost_author_code_figures_20260522`.
- Current Mac handoff copy:
  `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/3d_lung_yost_author_code_figures_20260522`.

## Actions

- Added R renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/r/report/render_yost_author_style_figures.R`.
- Added/updated plotting governance docs:
  - `plotting_factory/r/report/README.md`
  - `plotting_factory/python/report/README.md`
  - `plotting_factory/README.md`
  - `plotting_factory/AGENTS.md`
  - `plotting_factory/figure_memory/lung_3d_genome_2025/README.md`
  - `plotting_factory/contracts/lung_3d_author_style_profiles.json`
- Updated the author-style contract so Yost is the only current generated
  author-code route; Miao Liu and Qiaowei Liu are record-only/source-memory
  routes unless explicitly reopened.

## Output

- 6 Yost-only figure records in `FIGURE_MAP.tsv`.
- Formats per figure: PNG, SVG, PDF, TIFF.
- Groups:
  - `01_Yost2025_NatGenet_author_code_informed_reproduced_figures`
  - `02_Yost2025_NatGenet_STK11_LKB1_deep_dive`
- Figure mapping:
  - Yost Fig. 1-style HiChIP/CALDER2 landscape
  - Yost Fig. 2-style differential-loop boundary
  - Yost Fig. 3-style lineage/accessibility support
  - Yost Fig. 4-style gene-loop/RNA coupling
  - Yost Fig. 5-style CNV-aware loop interpretation
  - STK11/LKB1 exploratory extension, explicitly not an original paper figure

## Verification

- Renderer completed with R 4.5.3 from
  `/home/zerlinshen/conda/envs/r_multiomics/bin`.
- Remote package contains `6` PNG, `6` SVG, `6` PDF, `6` TIFF, `FIGURE_MAP.tsv`,
  and `README.md`.
- Filename check confirmed no Qiaowei or Miao outputs in the Yost-only package.
- Author-style profile validator:
  `OK profiles=3 groups=2 record_only=2`; runtime check found Yost `Rscript`
  and `python3`, while Miao/Qiaowei are skipped as record-only.
- `git diff --check` passed for the renderer, contract, and README/AGENTS
  governance files.
- Manual visual QA was performed on the key PNG outputs; the second R render
  removed the stray ggplot text-legend glyph and edge-label clipping.

## Cautions

- The STK11/LKB1 figure is a project-specific extension using Yost plotting
  grammar; do not describe it as an original Yost paper figure.
- The legacy Python/matplotlib renderer remains a processed-route summary and
  must not be described as author-code parity.
- Local Mac copies remain handoff artifacts; the remote project package is the
  source of truth.
