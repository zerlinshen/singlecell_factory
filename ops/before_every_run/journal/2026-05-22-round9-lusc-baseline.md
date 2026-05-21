# Round9 LUSC baseline run - 2026-05-22

## Objective

Start executing the Round9 full single-cell comparative optimization plan on a
second real data shape, with downstream modules and figure QA included.

## Attempted

- Read `ops/before_every_run/LATEST.md` and the recent Cell/Trevino module-smoke
  journal before launching.
- Execution mode: `benchmark`.
- Launched a project-root run under
  `/home/zerlinshen/projects/round9-singlecell-comparison`.
- Input boundary:
  `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/data/raw/wch_lung_cancer_atlas/prepared_input.h5ad`.
- Command used `SC_REQUIRE_PROJECT_ROOT=1`, `--allow-dirty`, `--gpu-mode off`,
  and requested upstream plus downstream modules:
  `clustering,annotation,differential_expression,pathway_analysis,composition,gene_regulatory_network,cell_communication,immune_phenotyping,tumor_microenvironment,pseudobulk_de,metacell`.

## Succeeded

- Run exit status: `0`.
- Run path:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1718Z-9907a7d/python/lusc_ps01_round9_baseline_20260522_011822`.
- `module_status.csv` records `ok` for:
  `cellranger`, `qc`, `ambient_correction`, `doublet_detection`, `clustering`,
  `annotation`, `differential_expression`, `gene_regulatory_network`,
  `metacell`, `cell_communication`, `composition`, `immune_phenotyping`,
  `tumor_microenvironment`, and `pathway_analysis`.
- `pseudobulk_de` truthfully skipped:
  `No explicit confirmatory pseudobulk contrast was provided.`
- Final AnnData shape: `72399 x 17267`.
- Figure validator found `30` valid PNGs and no corrupt files.

## Failed Or Risky

- Scrublet under-called doublets: `0.03%` call rate, `21` calls, undercall
  ratio `206.9x` versus expected `6%`.
- Candidate manuscript vector gate fails as expected because all `30` figures
  are PNG and there are `0` PDF/SVG vector outputs.
- Current run manifest did not include ambient skipped/no-trigger metadata even
  though AnnData `uns["ambient_correction"]` did. The module was patched after
  the run so future skip/dry-run/disabled paths expose ambient decision fields
  in manifest metadata.
- Scanpy DE emitted many DataFrame fragmentation performance warnings. The run
  completed, but DE is the dominant runtime (`147.877` seconds of `255.114`).

## Artifact Classification

- LUSC baseline run: `evidence-only`.
- Project root:
  `/home/zerlinshen/projects/round9-singlecell-comparison`.
- Human conclusion:
  `/home/zerlinshen/projects/round9-singlecell-comparison/ledger/human_review/2026-05-22-round9-lusc-baseline-human-conclusion.md`.
- Run evidence note:
  `/home/zerlinshen/projects/round9-singlecell-comparison/runs/2026-05-21T1718Z-9907a7d/evidence/round9_lusc_baseline_conclusion_20260522.md`.
- Figure QA artifacts:
  `/home/zerlinshen/.omx/state/round9-lusc-baseline-figure-validation-20260521T172300Z.json`
  and `.md`.

## Next Operator Reminders

- Do not promote Scrublet as sufficient for this tumor data shape.
- Next lane should rerun the same LUSC input with a doublet second-opinion or
  consensus strategy, then compare downstream biology and figure quality.
- Add an explicit pseudobulk contrast spec before expecting `pseudobulk_de` to
  provide confirmatory DE.

