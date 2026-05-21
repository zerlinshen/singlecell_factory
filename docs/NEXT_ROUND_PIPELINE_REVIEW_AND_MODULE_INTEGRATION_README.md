# Next Round Pipeline Review And Module Integration README

Date: 2026-05-22

This file is the handoff surface for the next session that will continue
method discovery, GitHub module intake, downstream optimization, and
publication-quality figure review across the three-factory pipeline.

## Current Verdict

Round9 validated that a single benchmark is not enough to replace a global
default. On LUSC PS01, Scrublet under-called doublets and scDblFinder /
consensus rescued a plausible tissue-atlas failure mode. On the earlier hgmm
ground-truth benchmark, Scrublet remains the stronger default. Therefore:

- Keep `scrublet` as the global default for now.
- Use `scdblfinder` or `consensus` as a second-opinion lane when diagnostics
  show under-calling or when tumor/tissue data shape matches the LUSC pattern.
- Promote a new default only after either a second different-shape dataset
  confirms it or a documented data-shape conditional activates it.
- Compare downstream consequences, not only upstream QC. A backend that changes
  doublet calls must be evaluated through annotation, DE, pathway, cell-cell
  communication, GRN, metacell, and figure outputs.

## Canonical Paths

Suite root:

`/home/zerlinshen/Bioinformatics Research Pipeline`

Factory repos:

- `singlecell_factory`: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory`
- `r_multiomics_factory`: `/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory`
- `plotting_factory`: `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory`

Current Round9 project evidence:

- Project root: `/home/zerlinshen/projects/round9-singlecell-comparison`
- Comparison report: `/home/zerlinshen/projects/round9-singlecell-comparison/reports/round9_lusc_doublet_downstream_comparison.md`
- Comparison JSON: `/home/zerlinshen/projects/round9-singlecell-comparison/reports/round9_lusc_doublet_downstream_comparison.json`
- Figure bundle: `/home/zerlinshen/projects/round9-singlecell-comparison/reports/figures/round9_lusc_closure`
- Human conclusion: `/home/zerlinshen/projects/round9-singlecell-comparison/ledger/human_review/2026-05-22-round9-lusc-closure-human-conclusion.md`

Temporary paper-grounded references for 3D genome work:

- Liu 2025 TP63/MYC loops raw mirror: `/home/zerlinshen/data/raw/liu_2025_natcomm_tp63_myc_loops`
- Liu 2025 TP63/MYC loops published figures: `/home/zerlinshen/data/external/published_figures/liu_2025_natcomm_tp63_myc_loops`
- Liu 2025 KRAS single-cell 3D genome raw mirror: `/home/zerlinshen/data/raw/liu_2025_natgenet_kras_single_cell_3d_genome`
- Liu 2025 KRAS single-cell 3D genome published figures: `/home/zerlinshen/data/external/published_figures/liu_2025_natgenet_kras_single_cell_3d_genome`

## Owner By Primary Output

Use this rule before moving code:

- `singlecell_factory` owns Python upstream execution, AnnData truth, modular
  CLI/config wiring, manifests, and cross-factory validation.
- `r_multiomics_factory` owns R-heavy primary analysis when the R method creates
  primary objects, statistics, or biological interpretation.
- `plotting_factory` owns presentation-only visual rendering, layout, palettes,
  figure schemas, and multi-panel figure composition. It must not own biological
  conclusions.

If a method produces a primary biological statistic, integrate it first in the
scientific owner repo. If it only renders already-computed tables, integrate it
in `plotting_factory`.

## GitHub Method Intake Ladder

For each external GitHub method or module:

1. Record the upstream repository, commit SHA, license, method paper, and
   minimum runtime dependencies.
2. Reproduce a small upstream example or a minimal smoke input without changing
   our pipeline.
3. Map the method to one pipeline owner:
   `singlecell_factory`, `r_multiomics_factory`, or `plotting_factory`.
4. Write a wrapper with an explicit input/output contract. Avoid direct
   cross-repo imports except through approved bridges or subprocess drivers.
5. Add CLI/config wiring only after the wrapper contract is stable.
6. Add tests that cover failure behavior, output shape, and provenance metadata.
7. Run at least one production smoke lane under `/home/zerlinshen/projects`.
8. Compare downstream consequences against the current baseline.
9. Update docs, run memory, and handoff records before committing.

Do not use `install_github` inside production pipeline execution. Pin and audit
external source locations first, then call local clones or installed env
packages through explicit paths and fail-loud checks.

## 3D Genome Paper-Repo Intake Strategy

The single-cell pipeline has already had many rounds of repair. The next
highest-value gap is the 3D genome side. Treat it differently from ordinary
single-cell module tuning: many 3D genome methods are tightly coupled to a
specific paper repository, genome build, contact-matrix format, normalization,
resolution, and figure recipe. The right first move is often to find the
paper's own repository or a widely adopted implementation, reproduce the
paper-level output, and only then decide what belongs in our reusable pipeline.

Use a paper-repo-first ladder:

1. Pick one 3D genome paper target and define the exact temporary ground truth:
   published loop/TAD/compartment calls, cell-state pattern, enhancer-promoter
   linkage, 3D genome-associated gene program, or figure panel.
2. Locate the paper repository, release tag or commit SHA, license, data
   accession, genome build, restriction enzyme / assay type, matrix format
   (`.hic`, `.cool`, `.mcool`, sparse TSV), bin resolution, and normalization
   method.
3. Reproduce the upstream repository outside our factories first. Keep logs,
   dependency versions, command lines, and a minimal successful output bundle.
4. Compare reproduced outputs against the published figure or supplementary
   table before writing pipeline wrappers. The paper output is a temporary
   ground truth for method evaluation, not an absolute biological truth.
5. Only after reproduction succeeds, split ownership:
   - primary contact-map processing, loop/TAD/compartment calling, statistical
     calls, or 3D genome interpretation belong in the scientific factory;
   - single-cell links and AnnData handoff belong in `singlecell_factory`;
   - figure rendering belongs in `plotting_factory`.
6. Wrap the method with explicit contracts and fail-loud dependency checks.
   Do not paste a whole paper script into production or silently change
   normalization, resolution, genome build, or filtering defaults.
7. Run an end-to-end smoke lane after the 3D genome method is wrapped, because
   the whole workflow matters: 3D genome evidence must still connect cleanly to
   single-cell state, downstream biological interpretation, and final figures.

Likely pitfalls to guard against:

- Paper repositories may be stale, untagged, partially manual, or missing
  private intermediate files.
- Published figures may use preprocessed matrices rather than raw FASTQ-level
  processing; record which level is being reproduced.
- Loop/TAD/compartment results are sensitive to genome build, bin size,
  normalization, blacklists, mappability, and sparse-matrix filtering.
- `.hic` and `.cool/.mcool` workflows are not interchangeable unless conversion
  and balancing are validated.
- Different tools optimize different endpoints: loop callers, TAD callers,
  compartment callers, single-cell 3D embedding, and enhancer-gene linkage
  should not be ranked by one metric.
- A pretty reproduced panel is not enough. Require a numeric comparison table
  plus manuscript-grade figure QA.

## Evaluation Matrix

Every candidate method should be judged on four layers:

1. Upstream module behavior:
   - runtime, memory, failure mode, fallback behavior
   - call rates or primary statistics
   - metadata and manifest provenance

2. Downstream biological stability:
   - retained cells and clusters
   - annotation unknown rate
   - marker consistency
   - DE gene count and pseudobulk limitation
   - pathway and ligand-receptor overlap
   - GRN / TF / metacell consequences when enabled

3. Cross-dataset generality:
   - one tumor/tissue dataset is not enough for a default change
   - require a second different-shape dataset or a data-shape conditional
   - document both winners and losers; methods have different valid domains

4. 3D genome paper-grounded validity:
   - reproduction distance from published loops, TADs, compartments, or
     enhancer-gene links
   - preservation of the paper's main cell-state or regulatory conclusion
   - sensitivity to genome build, bin resolution, normalization, and input
     matrix format
   - handoff quality into single-cell interpretation and final figure panels

5. Figure and publication quality:
   - vector SVG/PDF outputs for manuscript panels
   - no debug PNGs as final panels
   - labels readable at manuscript scale
   - no label/legend/panel overlap
   - color palettes distinguish methods and remain print-safe
   - figure conclusion matches the comparison report

## Figure Policy

The final principle is simple: choose the path that produces the most beautiful,
truthful, publication-grade figure. Do not choose a plotting path mechanically.

- Use `$nature-figure` whenever it can improve the final manuscript-style
  figure, including multi-panel design, published-figure comparison, visual QA,
  palette/layout polishing, SVG/PDF export review, and manuscript-scale label
  checks.
- Do not rewrite a validated R analysis path only to force it through a skill.
  Some R-based analyses naturally emit their own ggplot2, patchwork,
  ComplexHeatmap, or method-native plots. Keep those native outputs when they
  are part of the analysis contract, then polish or QA them when feasible.
- If `$nature-figure` is invoked, choose the backend by context: use R for
  plotting_factory-style final figures and R-native outputs; use Python only
  when continuing existing Python comparison scripts or matplotlib/seaborn
  figure code.
- Debug PNGs are never final manuscript panels. A final panel must be reviewed
  for readability, overlap, color, vector/raster export quality, and whether
  it supports the biological conclusion.
- The stop criterion for figures is not "a plot exists." Stop only when the
  figure is the best available presentation of the evidence under the current
  data, method, and journal-style constraints.

## Current Modules To Continue Reviewing

Already integrated or partially integrated:

- Doublet second-opinion lanes in `singlecell_factory`:
  `scrublet`, `doubletfinder`, `scdblfinder`, and `consensus`.
- Ambient correction mtx bridge for production DecontX.
- Post-run comparison and figure validators in `singlecell_factory/scripts`.
- R-side method drivers in `r_multiomics_factory/scripts`:
  `decontx_run.R`, `doubletfinder_run.R`, `scdblfinder_run.R`, `soupx_run.R`,
  and `nichenet_run.R`.
- Presentation helpers in `plotting_factory/r` for palettes, layout, DE plots,
  heatmaps, 3D dim plots, sankey plots, and report-level figure packages.

Still needs a deliberate next-session gap matrix:

- Milty Omix / broader multi-omics modules not yet wired into the modular CLI.
- 3D genome paper-repo intake: choose the strongest paper repositories and
  reproduce them before converting anything into reusable modules.
- NicheNet regulatory layer versus existing LIANA ligand-receptor layer.
- SoupX versus DecontX under raw-barcode and filtered-only input shapes.
- Batch integration backends: Harmony/RPCA/CCA versus fastMNN/MNN/ComBat.
- Reference projection and label transfer methods: Seurat anchors, Symphony,
  SingleR, and scmap.
- Figure factory migration: move reusable pure rendering from one-off scripts
  into stable `plotting_factory` APIs and schemas.

## Recommended Next Session Plan

1. Build a gap matrix with rows = candidate methods and columns = owner repo,
   upstream SHA/license, input shape, output contract, dependencies, tests,
   production evidence, and figure needs.

2. Pick two data shapes for confirmation:
   - LUSC/tumor tissue atlas shape from Round9.
   - A different-shape benchmark, preferably hgmm or a public brain/development
     dataset where the expected doublet and ambient profiles differ.
   - For 3D genome work, add one paper-grounded dataset from the Liu 2025
     references above and score reproduction against its published figures or
     supplementary tables.

3. Run paired lanes through the same downstream module set:
   - baseline
   - candidate backend
   - consensus or conditional router

4. For 3D genome candidates, first reproduce the paper repository in isolation,
   then wrap only the stable method boundary. Keep the upstream clone, command
   log, dependency record, and reproduced figure/table comparison as evidence.

5. Generate one comparison report per method family, not only one upstream
   metric table.

6. Generate a figure bundle for every accepted comparison using vector outputs
   plus a manual visual-review checklist. Use `$nature-figure` whenever it can
   improve the final appearance or QA, but keep R-native analysis plots in R
   when that is the correct method contract.

7. Only then decide whether a method becomes:
   - default
   - data-shape conditional
   - optional second-opinion lane
   - rejected / not worth integrating

## Verification Commands

Singlecell targeted gate:

```bash
cd "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory"
pytest -q tests/test_compare_singlecell_runs.py tests/test_doublet_detection_backends.py tests/test_ambient_correction.py tests/test_validate_figure_outputs.py --no-cov
python scripts/validate_figure_outputs.py --run-dir /home/zerlinshen/projects/round9-singlecell-comparison/reports/figures/round9_lusc_closure --require-vector --json
python scripts/validate_project_governance.py /home/zerlinshen/projects/round9-singlecell-comparison --json
git diff --check
```

R method syntax gate:

```bash
cd "/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory"
/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript -e 'files <- c("R/batch_integration_module.R", "R/projection_module.R", "scripts/decontx_run.R", "scripts/doubletfinder_run.R", "scripts/scdblfinder_run.R", "scripts/soupx_run.R", "scripts/nichenet_run.R"); for (f in files) parse(file=f); cat("parsed", length(files), "r_multiomics files\n")'
git diff --check
```

Plotting syntax gate:

```bash
cd "/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory"
/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript -e 'files <- c("r/core/layout_helpers.R", "r/core/palette_engine.R", "r/general/de_plots.R", "r/general/dim_plots_3d.R", "r/general/heatmaps.R", "r/general/sankey_plots.R", "r/report/lung_3d_reproduction_ppt_figures.R"); for (f in files) parse(file=f); cat("parsed", length(files), "plotting files\n")'
git diff --check
```

## Stop Condition For Next Round

Do not stop at "the module runs." Stop only when:

- the method has a documented owner and input/output contract;
- at least one production lane has artifact evidence under `/home/zerlinshen/projects`;
- downstream consequences are compared against baseline;
- figure outputs pass mechanical QA and have a manual-review checklist;
- docs and run memory are updated;
- any default change satisfies the two-dataset or data-shape-conditional rule.
