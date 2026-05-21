# Round9 LUSC closure - 2026-05-22

## Objective

Close the Round9 single-cell comparison slice with a second-opinion doublet
lane, downstream consequence comparison, and vector figure QA rather than
judging the pipeline from the baseline Scrublet warning alone.

## Execution Mode

`benchmark`

## Runs

- Baseline refresh:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1744Z-9907a7d/python/lusc_ps01_round9_baseline_refresh_20260522_014447`
- scDblFinder:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1749Z-9907a7d/python/lusc_ps01_round9_scdblfinder_20260522_014947`
- Consensus OR:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1756Z-9907a7d/python/lusc_ps01_round9_consensus_scdbl_or_20260522_015643`

All runs used project-root output, `SC_REQUIRE_PROJECT_ROOT=1`,
`--allow-dirty`, `--gpu-mode off`, `--parallel-workers 1`, and the same broad
downstream optional module set.

## Succeeded

- scDblFinder backend was wired through the modular CLI and R mtx bridge.
- Consensus supports `scrublet_scdblfinder` and records pair, logic, agreement,
  pairwise agreement, and per-backend rates.
- All three LUSC lanes completed with `14` modules `ok`.
- `pseudobulk_de` truthfully skipped because no explicit contrast contract was
  provided.
- Baseline Scrublet call rate remained implausibly low: `21` doublets
  (`0.03%`, undercall ratio `206.9x`).
- scDblFinder rescued the call rate: `6737` doublets (`9.30%`).
- Consensus OR rescued the call rate: `6751` doublets (`9.32%`), agreement
  `0.9069`, final retained cells `65669`.
- Downstream metrics remained stable enough for a conditional candidate:
  consensus `35` clusters versus baseline `36`, `0.0%` annotation unknown rate,
  and `10488` significant DE genes versus baseline `10771`.

## Reports And Figures

- Comparison JSON:
  `/home/zerlinshen/projects/round9-singlecell-comparison/reports/round9_lusc_doublet_downstream_comparison.json`
- Comparison Markdown:
  `/home/zerlinshen/projects/round9-singlecell-comparison/reports/round9_lusc_doublet_downstream_comparison.md`
- Curated figure bundle:
  `/home/zerlinshen/projects/round9-singlecell-comparison/reports/figures/round9_lusc_closure`
- Figure QA:
  `/home/zerlinshen/projects/round9-singlecell-comparison/reports/figures/round9_lusc_closure/figure_validation.json`
- Human conclusion:
  `/home/zerlinshen/projects/round9-singlecell-comparison/ledger/human_review/2026-05-22-round9-lusc-closure-human-conclusion.md`

Figure mechanical QA passed with `9` valid figures, `6` valid vector outputs,
`3` rasters, and `0` warnings.

## Decision

Keep Scrublet as the global default because hgmm ground truth still supports it.
For heterogeneous tissue/tumor atlas data where the Scrublet under-call
diagnostic fires, use scDblFinder or consensus OR with `scrublet_scdblfinder`
as a conditional second-opinion lane. Do not promote this to a universal
default without another different-shape dataset or a stricter data-shape
router.

## Remaining Risks

- Manual visual review is still required before treating the curated vector
  bundle as publication-ready.
- `pseudobulk_de` remains unavailable for confirmatory claims until a real
  contrast spec is supplied.
- Broader default changes remain out of scope for this round.
