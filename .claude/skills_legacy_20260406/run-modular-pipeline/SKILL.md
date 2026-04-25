---
name: run-modular-pipeline
description: Run this repository's modular scRNA-seq pipeline end-to-end or as a focused subset. Use when asked to execute analysis runs, choose optional modules, tune CLI parameters, or resume from checkpoints.
argument-hint: [project-name] [sample-root] [optional-modules]
---

Run the modular workflow with reproducible commands and safe defaults.

1. Validate inputs before execution.
- Confirm `sample-root` exists and contains `outs/filtered_feature_bc_matrix`.
- If `--outs-dir` is passed, verify it exists.
- Prefer `--checkpoint` for runs longer than a few modules.

2. Build the command from this template.
```bash
python -m workflow.modular.cli \
  --project <project-name> \
  --sample-root <sample-root> \
  --optional-modules <comma-separated-modules> \
  --checkpoint
```

3. Use module subsets intentionally.
- Fast structural run: `clustering,differential_expression,annotation`.
- Trajectory-focused run: add `trajectory,pseudo_velocity`.
- Tumor/TME run: add `cnv_inference,immune_phenotyping,tumor_microenvironment,evolution`.

4. Choose parallelism carefully.
- Start with `--parallel-workers 1` for debugging.
- Increase to `2-4` only after a clean sequential run.
- Remember mutating stages (`batch_correction`, `rna_velocity`) are forced sequentially by pipeline logic.

5. Verify completion immediately.
- Read `<run_dir>/module_status.csv` and fail the run if any `status != ok`.
- Read `<run_dir>/run_manifest.json` and report key metadata (`n_clusters`, cells after QC, DE counts, CNV summary, etc.).

6. On failure, resume instead of rerunning everything.
- Use `--resume-from <failed-module>` with the same `--project` and module list.
- If checkpoint files are missing, rerun from start with `--checkpoint`.
