# 2026-05-19 two-dataset full simulation and DE/checkpoint patch

## Runs

- Cell/Trevino full public RNA controller-validation run completed: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-19T1005Z-82f1964`.
  - Full RNA input: `57868 x 33355`; after QC `56879 x 25519`; final after doublet removal `55568` cells.
  - Modules all ok through clustering, annotation, cell_cycle, DE, signatures, metacell, trajectory, composition, pathway, pseudobulk, cell_fate.
  - Key metrics: 25 clusters; annotation unknown pct 0.0; 4955 significant DE genes; 15 signatures; 740 metacells; 7 terminal states; wall about 185.836 s.
  - Full simulation figure/result comparison: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-19T1005Z-82f1964/evidence/full_sim_annotation_results/`.
  - Label comparison vs Trevino seurat clusters: ARI 0.567526, NMI 0.750967; figure marker panel present with no missing markers.
- NC2024 full tumor no-checkpoint run reached main DE evidence, then was terminated intentionally at status 143 to stop a redundant old-code DE rerun: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-19T1010Z-82f1964`.
  - Full tumor subset: 877023 cells from full cohort input under `singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort`.
  - Completed evidence before termination: doublet plot, CSS clustering UMAP, annotation CSV/plots, batch correction before/after UMAPs, context validation JSON, and main DE marker CSVs.
  - Main DE outputs: `marker_genes_all.csv` 1901 lines, `marker_genes.csv` 1893 lines, `marker_top5_by_cluster.csv` 96 lines.
  - Observed failure family: `rapids-singlecell 0.13.4` lacks `rapids_singlecell.tl.rank_genes_groups`, so Wilcoxon DE falls back to CPU despite GPU clustering being available.
  - Observed full-scale pressure: CPU DE took about 47 minutes for main DE, triggered MemoryWatchdog warnings, then substate DE raised `sparse subset too large to densify safely` and old `run()` logic retried the whole DE implementation.

## Code Changes Applied

- `workflow/modular/modules/differential_expression.py`
  - API-gates RAPIDS DE with `hasattr(rsc.tl, "rank_genes_groups")` and records `de_rapids_singlecell_version` plus a clear fallback reason.
  - Uses existing sparse Welch fallback for `scale_mode=massive` sparse CPU DE; records `de_backend=cpu_sparse` and `de_test_actually_used=sparse_welch_fallback`.
  - Makes sparse/dense Welch fallback correction-aware (`benjamini-hochberg` or `bonferroni`) and records `de_correction_actually_used`.
  - Re-raises `MemoryGuardError` after recording `mem_warnings`, so the orchestrator does not retry and also cannot falsely finalize DE as `ok`.
  - Skips oversized substate DE with metadata `de_substate_skip_reasons` instead of raising and triggering a full retry.
- `workflow/modular/context.py`
  - Handles AnnData `Dataset2D` checkpoint write failures by skipping full AnnData checkpoint, removing partial checkpoints, and keeping JSON sidecars.
  - Supports metadata/sidecar/no-adata checkpoint policy values via `SC_CHECKPOINT_POLICY` and legacy `SCF_MASSIVE_CHECKPOINT_POLICY`, including `metadata_only`.
- Tests added for RAPIDS API fallback, massive sparse DE routing, Bonferroni fallback correction, MemoryGuard no-retry/non-ok behavior, oversized substate skip, Dataset2D checkpoint skip, and checkpoint policy aliases/precedence.
- `README.md`
  - Documents module-specific GPU behavior, RAPIDS DE fallback, sparse Welch correction metadata, checkpoint policy aliases, and the non-resumable nature of sidecar-only checkpoints.

## Verification

- `py_compile` passed for changed source and tests.
- `git diff --check` passed for changed files.
- Targeted pytest passed with `CONDA_PREFIX=/home/zerlinshen/conda/envs/sc_gpu_stable PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 -o addopts=''`: 9 passed, 2 expected sparse numeric warnings.
- Real environment smoke passed: actual `rapids-singlecell 0.13.4` recorded missing `tl.rank_genes_groups`, fell back to `cpu_sparse`, used `sparse_welch`, honored `de_correction=bonferroni`, and wrote `marker_genes_all.csv`.

## Retention / Next Run Guidance

- Preserve Cell full run `2026-05-19T1005Z-82f1964` as evidence-only controller-validation simulation.
- Preserve NC partial run `2026-05-19T1010Z-82f1964` logs and lightweight outputs as failure-family evidence; it is not a final source-of-truth run and has no final manifest/final AnnData.
- Future NC full rerun should use patched code, preferably with checkpoint enabled; Dataset2D checkpoint will keep sidecar evidence instead of crashing. Sidecar-only checkpoint records provenance and is not a resumable AnnData checkpoint.
- Do not interpret RAPIDS DE fallback as GPU clustering failure. GPU clustering and GPU batch postprocessing are separate; only DE method availability was missing in the installed RAPIDS version.
