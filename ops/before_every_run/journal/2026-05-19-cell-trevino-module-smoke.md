# Cell/Trevino module-smoke run - 2026-05-19

## Objective

Test additional `singlecell_factory` RNA modules on real Cell/Trevino-derived data while preserving the public-resource truth boundary and current Cell/Trevino retention policy.

## Attempted

- Built a derived smoke input from `/home/zerlinshen/projects/wave5-trevino/inputs/gse162170_public_rna_prepared/prepared_input.h5ad`.
- Output input: `/home/zerlinshen/projects/wave5-trevino/inputs/module_smoke_20260519_cell_rna_3000/prepared_input.h5ad`.
- Input shape: `3000 x 33355`, stratified across `pcw16`, `pcw20`, `pcw21`, and `pcw24`; `8` sample groups; `23` original Seurat clusters.
- Custom marker JSON: `/home/zerlinshen/projects/wave5-trevino/configs/module_smoke_20260519/cell_brain_development_markers.json`.
- Custom signature JSON: `/home/zerlinshen/projects/wave5-trevino/configs/module_smoke_20260519/cell_brain_development_signatures.json`.
- Launched project-root run `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-19T0956Z-82f1964` with `--allow-dirty`, `--checkpoint`, CPU clustering/DE (`--gpu-mode off`), relaxed QC thresholds for prepared matrices, and exploratory pseudobulk mode.
- Full launch command/logs: `/home/zerlinshen/projects/wave5-trevino/ledger/run_records/2026-05-19T0956Z-82f1964/` and copied into the run `python/logs/` directory.

## Succeeded

- Run exit status: `0`.
- All module statuses were `ok`:
  - `cellranger`, `qc`, `doublet_detection`, `clustering`
  - `annotation`, `cell_cycle`, `differential_expression`
  - `gene_signature_scoring`, `metacell`, `trajectory`
  - `composition`, `pathway_analysis`, `pseudobulk_de`, `cell_fate`
- Important outputs include:
  - `annotation/cell_type_annotation.csv`
  - `cell_cycle/cell_cycle_scores.csv`
  - `differential_expression/marker_genes.csv`
  - `gene_signature_scoring/gene_signature_scores.csv`
  - `metacell/metacells.h5ad`
  - `trajectory/dpt_pseudotime.csv`
  - `composition/composition_counts.csv`
  - `pathway_analysis/pathway_enrichment.csv`
  - `pseudobulk_de/pseudobulk_de_results.csv`
  - `cell_fate/fate_probabilities.csv`

## Metrics

- Raw input: `3000` cells, `33355` genes.
- After QC: `3000` cells, `22556` genes.
- Doublets detected/removed: `3`; final cells `2997`.
- Clusters: `19`.
- Cell-cycle phase counts: G1 `2400`, S `414`, G2M `183`.
- DE significant genes: `2279`.
- Gene signatures scored: `15` total, including Cell/Trevino smoke signatures and built-ins.
- Metacells: `50`, median size `47`.
- Composition groups: `8` sample groups, `5` annotated cell types.
- Exploratory pseudobulk status: completed.
- Cell-fate terminal states: `5`.
- Pipeline wall time: about `58` seconds.

## Expected fallbacks and cautions

- RAPIDS scrublet rejected dtype and doublet detection used CPU scrublet fallback.
- `SEACells` was unavailable, so metacell used MiniBatchKMeans fallback.
- `pertpy` was unavailable, so composition used fallback tests.
- `gseapy` and `decoupler` were unavailable, so pathway analysis used the built-in fallback.
- `pydeseq2` was unavailable, so pseudobulk used Mann-Whitney fallback.
- `CellRank` was unavailable, so cell-fate used manual fallback.
- These fallbacks are acceptable for module-interface smoke evidence, but not for new Cell/Trevino biological claims.

## Validation

- `scripts/validate_project_governance.py --json /home/zerlinshen/projects/wave5-trevino` returned overall `pass` after the run.
- Warnings are existing advisory project.yaml reproduction fields and legacy run layout warnings; the new run had no findings.

## Artifact classification

- Canonical human-facing Cell/Trevino source of truth: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88`.
- Canonical linked pipeline run: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88`.
- New module-smoke run: `evidence-only`.
- New derived input/configs: `evidence-only` module-smoke artifacts.
- Existing Cell/Trevino runs: untouched; do not delete without a separate Cell cleanup inventory.

## Next operator reminders

- Do not interpret the module-smoke annotation, pathway, pseudobulk, or cell-fate outputs as manuscript biology.
- If exact Cell/Trevino reproduction is requested, stage missing FASTQ/fragments/BPNet resources first and update the project policy.
- If module dependency completeness matters, consider installing/validating optional backends (`SEACells`, `pertpy`, `gseapy`, `decoupler`, `pydeseq2`, `CellRank`) in a separate environment task.
