---
name: validate-analysis-outputs
description: Validate scientific and structural quality of pipeline outputs. Use after a run finishes to check module completeness, key plots/tables, and metadata sanity before reporting results.
argument-hint: [run-dir]
---

Validate run outputs using manifest-driven checks and module artifacts.

1. Resolve target run.
- Use provided run dir or newest folder under `results/`.

2. Perform structural checks.
- Ensure `run_manifest.json`, `module_status.csv`, and `final_adata.h5ad` exist.
- Ensure each `ok` module has a corresponding subdirectory.
- Flag any failed modules immediately.

3. Perform high-signal content checks.
- QC: pre/post filter plots exist; `cells_after_qc` less than `raw_cells`.
- Clustering: `umap_leiden.png` exists; `n_clusters >= 2`.
- Annotation: `cell_type_annotation.csv` and `umap_cell_type.png` exist.
- Differential expression: `marker_genes.csv` exists and is non-empty.
- Trajectory (if enabled): `dpt_pseudotime.csv` exists and has non-null values.
- CNV/evolution (if enabled): required CSV/PNG outputs exist.

4. Cross-check metadata consistency.
- Compare `module_status.csv` order against dependency order expectations.
- Confirm key metadata fields are numeric and in plausible ranges.
- Highlight missing metadata for modules that reported `ok`.

5. Output a short validation verdict.
- PASS with caveats, or FAIL with blocking issues.
- Include exact paths for missing or suspect artifacts.
- Recommend rerun scope (full, partial, or `--resume-from`).
