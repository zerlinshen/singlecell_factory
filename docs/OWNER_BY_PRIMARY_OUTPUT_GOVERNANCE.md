# Owner-by-primary-output governance

This document defines which factory owns scientific truth when a workflow crosses
`singlecell_factory`, `r_multiomics_factory`, and `plotting_factory`.


The public umbrella name for this governed three-repository suite is
**Bioinformatics Research Pipeline**. The name combines the repositories at the
architecture/product level only; it does not collapse ownership boundaries or
make `plotting_factory` responsible for biological conclusions.

The canonical local suite root is
`/home/zerlinshen/Bioinformatics Research Pipeline`. The historical
`/home/zerlinshen/<repo>` factory paths are retired; do not recreate symlinks or
write new validators that depend on those aliases.

## Rule

Use **owner-by-primary-output**:

1. **`singlecell_factory` is the default global control plane** for Python-heavy
   single-cell upstream work, AnnData truth, bundle export, cross-factory final
   validation, and governance records.
2. **`r_multiomics_factory` owns R-heavy/spatial primary scientific truth** when
   preprocessing, statistical analysis, biological interpretation, and primary
   data objects are produced by R workflows. Examples include Visium / Space
   Ranger outputs consumed directly by Seurat, SpatialExperiment, Giotto,
   BayesSpace, or other R-first spatial transcriptomics workflows.
3. **`plotting_factory` is a presentation-only plotting surface**. It owns plot
   helpers, color palettes, theme/layout details, rendering contracts, and figure
   schemas. `plotting_factory` must not own biological conclusions.

## Practical routing

| Workflow shape | Local scientific owner | Global/reporting relationship |
| --- | --- | --- |
| Cell Ranger / Scanpy / AnnData preprocessing, then R plotting | `singlecell_factory` | Reports directly from `singlecell_factory` governance records. |
| R-heavy spatial or multi-omics workflow where R creates primary objects and conclusions | `r_multiomics_factory` | `r_multiomics_factory` keeps local run governance and reports summary/evidence to global cross-factory records. |
| Pure color, layout, theme, or rendering change | `plotting_factory` | No biological claim ownership; validate figure rendering and keep upstream gates aligned. |

## Required evidence by owner

### `singlecell_factory`

- run manifest / module status / project-root ledger;
- bundle manifest and provenance;
- cross-factory validator output;
- final governance report under `ops/governance_records/`.

### `r_multiomics_factory`

For R-heavy/spatial ownership, record at minimum:

- project/run manifest;
- R `sessionInfo()` or equivalent package lock evidence;
- input provenance for Space Ranger / 10x / public matrices;
- primary R object paths, such as Seurat or SpatialExperiment objects;
- module/status table for R preprocessing, statistics, and interpretation steps;
- expected figures and claim boundaries;
- bridge summary back to `singlecell_factory` global governance when the project
  participates in cross-factory reporting.

### `plotting_factory`

- figure schema / style contract;
- rendering smoke or parity evidence;
- explicit upstream owner for any biological interpretation shown by the figure.

## Validation

`singlecell_factory/scripts/validate_ci_governance.py` checks that all three
repos keep this ownership language synchronized. Current real-data effectiveness
is validated through the NC2024 + Cell/Trevino gates in
`ops/governance_records/2026-05-18-two-real-dataset-final-validation/`.
