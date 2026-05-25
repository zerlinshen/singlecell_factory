# Lung 3D Genome Single-Figure Redesign

Date: 2026-05-22

## Scope

Redesigned the lung 3D-genome reproduction figure handoff from a small set of
combined PPT-style panels into individual publication-style figures mapped to
original paper figures or explicit extension/boundary figures.

This was a plotting/reporting run, not a new biological re-analysis.

## Preflight

- Read current before-every-run memory and project governance records before
  remote mutation.
- Preserved the suite contract:
  - factories are tools
  - `/home/zerlinshen/projects` owns evidence, figures, and conclusions
  - `plotting_factory` owns presentation-only rendering
  - the Mac workspace owns coordination/handoff copies only

## Actions

- Added plotting renderer:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/python/report/render_lung_3d_single_figures.py`
- Added plotting renderer README:
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/python/report/README.md`
- Used isolated render environment:
  `/home/zerlinshen/.venvs/plotting_factory_figures`
- Generated project-owned figure package:
  `/home/zerlinshen/projects/reproductions/figure_packages/3d_lung_stk11_single_figures_20260522`
- Synced handoff copy to Mac:
  `/Users/zerlinshen/Downloads/1. Codex/remote_handoffs/3d_genome_reproductions_2026-05-22/figures/3d_lung_stk11_single_figures_20260522`

## Outputs

- 12 standalone figure records in `FIGURE_MAP.tsv`.
- 12 PNG, 12 SVG, and 12 PDF files.
- The figures cover:
  - Yost 2025 Nat Genet Fig. 1-5 style reproductions
  - Yost strict scATAC boundary and STK11/LKB1 exploratory extension
  - Miao Liu 2025 Nat Genet Fig. 2-style processed single-cell 3D genome
    progression and STK11 boundary check
  - Qiaowei Liu 2025 Nat Commun Fig. 2, Fig. 3, and Fig. 6 style processed or
    proxy reproductions

## Scientific Conclusions

- Yost 2025 remains partially reproduced at LUAD/LUSC MVS processed scope:
  HiChIP loop landscape, CALDER2 A/B rendering, CNV-aware loop interpretation,
  and lineage motif/accessibility support are computable. Strict cell-type
  scATAC/chromVAR reproduction remains not reproduced because the required raw
  ArchR/count route and LUSC matched public scATAC are absent or resource-gated.
- Miao Liu 2025 is the true single-cell 3D genome paper in this bundle.
  Processed 4DN/source-data evidence supports the central progression claim, but
  raw WES/SRA/4DN parity was not attempted.
- Qiaowei Liu 2025 supports a TP63-MYC chromatin-loop mechanism at
  processed/proxy level. It is not a single-cell paper and is not STK11-centric.
- STK11/LKB1 should be framed as a Yost LUAD exploratory extension only. Miao is
  a boundary/negative STK11 result, and Qiaowei should not be cited as STK11
  evidence.

## Verification

- `python -m py_compile` passed for the renderer.
- Renderer completed and wrote 12 figures.
- PNG QA passed for 12/12 figures: each PNG was nonblank and larger than
  1600 x 1100 pixels.
- `git diff --check -- python/report/render_lung_3d_single_figures.py` passed.
- Manual visual QA inspected key figures after iterative layout fixes:
  Yost Fig. 1, Yost Fig. 4, Miao Liu Fig. 2, and Qiaowei Liu Fig. 2.

## Cautions

- Do not describe this package as full raw FASTQ/WES/Micro-C/HiChIP parity.
- Do not treat the Mac handoff directory as the source of truth.
- Do not move scientific outputs into `plotting_factory`; new figure outputs
  belong under `/home/zerlinshen/projects/...`.
