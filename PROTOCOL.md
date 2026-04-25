# singlecell_factory Beginner Protocol (Complete, Practical)

This protocol is a beginner-friendly, end-to-end guide for running the current modular scRNA-seq pipeline.

It is aligned with the current codebase (`workflow/modular/*`) and CLI (`python -m workflow.modular.cli`).


## Bridge Architecture

The canonical R source of truth is `multiomics_r_factory/`:

- `multiomics_r_factory/R/` — Seurat-based high-level analysis modules (12 files)
- `multiomics_r_factory/R_bundle/` — bundle-path-specific helpers (`remote_bundle_manifest.R`)

The bridge at `bridges/local_r_pipeline_macbook/` references these via symlinks:

- `bridges/local_r_pipeline_macbook/R` → `../../../multiomics_r_factory/R`
- `bridges/local_r_pipeline_macbook/R_bundle` → `../../../multiomics_r_factory/R_bundle`

**Rule**: `bridges/.../R` and `bridges/.../R_bundle` must always be symlinks, never real directories.
Verify with: `bash scripts/ci/check_bridge_symlink.sh`

Do NOT place real R files under `bridges/local_r_pipeline_macbook/R/` or `bridges/local_r_pipeline_macbook/R_bundle/`.
Edit R sources in `multiomics_r_factory/` only.

---

## 0. Before Every Meaningful Remote Run

Read first:
- `/home/zerlinshen/singlecell_factory/ops/before_every_run/LATEST.md`
- the newest relevant entry under `/home/zerlinshen/singlecell_factory/ops/before_every_run/journal/`

Operational rule for the NC2024 full cohort:
- `execution_mode = debug_massive`
  - use direct `massive`
  - stop at first failing module
  - patch only that module
- `execution_mode = controller_validation`
  - use the official orchestration path `large -> massive`
  - treat `large` as a capacity probe, not the main completion lane

After a meaningful remote run finishes or fails:
- update `/home/zerlinshen/singlecell_factory/ops/before_every_run/journal/`
- update `/home/zerlinshen/singlecell_factory/ops/before_every_run/LATEST.md`
- keep one evidence-rich failed run if it explains the winning fix

Targeted evidence lane for NC2024 subtype checkpoint closure:
- If stage-1 baseline is already proven, do not rerun the full cohort just to inspect subtype checkpoint hierarchy.
- Prefer the rerunnable audit script over existing subtype LIANA outputs:
  - `python scripts/audit_nc2024_subtype_checkpoint_pairs.py --luad-raw-csv ... --lusc-raw-csv ... --luad-top20-csv ... --lusc-top20-csv ... --output-dir ...`
- After the raw pair audit, generate the higher-level hierarchy summary:
  - `python scripts/summarize_nc2024_checkpoint_hierarchy.py --audit-dir <checkpoint-audit-run-dir>`
- For the standard verification lane, prefer:
  - `bash scripts/verify_nc2024_subtype_checkpoint_audit.sh`
- Canonical example inputs:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/luad/cell_communication/cell_communication_liana.csv`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/lusc/cell_communication/cell_communication_liana.csv`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_adenocarcinoma_checkpoint_top20.csv`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_squamous_cell_carcinoma_checkpoint_top20.csv`
- Use this lane to separate:
  - raw subtype LIANA pair presence
  - current paper-facing top20 visibility
- Keep the checkpoint verdict `partial` unless the paper-facing hierarchy itself is reproduced, not merely because the raw pairs are present.


---

## 1. What You Will Do

You will:
1. Set up the environment.
2. Validate your 10X input directory.
3. Run a minimal analysis first.
4. Run a full analysis.
5. Read outputs (`run_manifest.json`, `module_status.csv`, figures/tables).
6. Recover from failures with checkpoint/resume.
7. Tune parameters safely.

Pipeline structure:
- Mandatory modules (always): `cellranger -> qc -> doublet_detection`
- Optional modules (22): selected by `--optional-modules` with auto dependency resolution.

---

## 2. Prerequisites

### 2.1 Software

```bash
cd /home/zerlinshen/singlecell_factory
conda env create -f environment.yml
conda activate sc10x
export MPLCONFIGDIR=$PWD/.mplconfig
export NUMBA_CACHE_DIR=/tmp/numba_cache
```

Recommended for stable plotting/JIT behavior:
- `MPLCONFIGDIR`
- `NUMBA_CACHE_DIR`

### 2.2 Input Data Format

Your sample root must contain:

```text
data/raw/<your_sample>/outs/filtered_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

If this folder is missing, either:
- point `--outs-dir` to the real matrix folder, or
- provide Cell Ranger inputs and allow rerun (`--fastq-dir`, `--transcriptome-dir`).

---

## 3. Quick Validation Before First Run

```bash
test -d data/raw/<your_sample>/outs/filtered_feature_bc_matrix && echo "matrix found"
```

Optional (recommended):
```bash
ls data/raw/<your_sample>/outs/filtered_feature_bc_matrix
```

---

## 4. First Run (Minimal, Safe)

Run the smallest useful workflow first:

```bash
python -m workflow.modular.cli \
  --project demo_minimal \
  --sample-root data/raw/<your_sample> \
  --optional-modules clustering
```

What this gives you:
- QC-filtered and de-doubleted data
- PCA/UMAP/Leiden clusters
- baseline figures to confirm data quality

---

## 5. Full Run (Local, No Network Dependency)

```bash
python -m workflow.modular.cli \
  --project demo_full_local \
  --sample-root data/raw/<your_sample> \
  --optional-modules clustering,cell_cycle,batch_correction,differential_expression,annotation,trajectory,pseudo_velocity,rna_velocity,cnv_inference,pathway_analysis,cell_communication,gene_regulatory_network,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,evolution,pseudobulk_de,cell_fate,composition,metacell \
  --velocity-bam data/raw/<your_sample>/outs/possorted_genome_bam.bam \
  --transcriptome-dir /path/to/refdata-gex-GRCh38-2024-A
```

Notes:
- `rna_velocity` needs either `--velocity-loom` OR (`--velocity-bam` + resolvable GTF via `--velocity-gtf` or `--transcriptome-dir`).
- `validate_cbioportal` is not included above (avoids network dependency).

---

## 6. Module Catalog (Current Pipeline)

Mandatory:
- `cellranger`
- `qc`
- `doublet_detection`

Optional (22):
- `clustering`
- `cell_cycle`
- `batch_correction`
- `differential_expression`
- `annotation`
- `trajectory`
- `pseudo_velocity`
- `rna_velocity`
- `cnv_inference`
- `pathway_analysis`
- `cell_communication`
- `gene_regulatory_network`
- `validate_cbioportal`
- `immune_phenotyping`
- `tumor_microenvironment`
- `gene_signature_scoring`
- `evolution`
- `pseudobulk_de`
- `cell_fate`
- `composition`
- `metacell`
- `paper_repro`

Dependency handling is automatic: if you request a downstream module, upstream modules are auto-included.

---

## 7. Where Results Go

Default output root:
- `/home/zerlinshen/singlecell_factory/results`

Run directory format:

```text
<output-dir>/<project>_<timestamp>/
├── final_adata.h5ad
├── run_manifest.json
├── module_status.csv
├── .checkpoints/              # only if --checkpoint
└── <module_name>/             # figures/tables per module
```

---

## 8. How To Read Success/Failure Correctly

Use `module_status.csv` and `run_manifest.json`.

Normalized status values:
- `ok`: module completed successfully.
- `skipped`: module intentionally skipped (for example: missing prerequisite keys, single-batch/single-sample conditions, or unavailable input-specific requirements).
- `failed`: module error.

Important:
- Mandatory modules must be `ok`.
- Optional modules can be `ok` or `skipped` with a meaningful reason.

---

## 9. Crash Recovery (Checkpoint + Resume)

Enable checkpoints:

```bash
python -m workflow.modular.cli \
  --project demo_resume \
  --sample-root data/raw/<your_sample> \
  --optional-modules clustering,differential_expression,annotation \
  --checkpoint
```

Resume:

```bash
python -m workflow.modular.cli \
  --project demo_resume \
  --sample-root data/raw/<your_sample> \
  --optional-modules clustering,differential_expression,annotation \
  --checkpoint \
  --resume-from annotation
```

Current resume behavior:
- Reuses the latest run directory for that `--project` containing `.checkpoints`.
- Finds checkpoint nearest before `--resume-from` by searching backward in execution order.
- Raises explicit `FileNotFoundError` if no suitable checkpoint exists.

---

## 9A. NC2024 Full-Cohort Remote Reproduction

For the current `E-MTAB-13526` reproduction lane, do not use the older tumor-only scripts.

Use this owned full-cohort path instead:

```bash
cd /home/zerlinshen/singlecell_factory
bash scripts/run_emtab13526_full_cohort_with_fallback.sh
```

This controller will:
1. write a preflight inventory under `results/`
2. delete stale derived artifacts from earlier failed attempts
3. build or validate `data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
4. run one `large` stage-1 probe
5. auto-fallback once to `massive` only if the `large` failure is clearly capacity-related

Owned scripts:
- `scripts/prepare_emtab13526_full_cohort_zarr.py`
- `scripts/run_emtab13526_full_cohort_stage1.sh`
- `scripts/run_emtab13526_full_cohort_with_fallback.sh`

Expected prepared-input artifacts:
- `data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- `data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.summary.json`
- `data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.ready`

The prepare contract is:
- one unsplit object for all `81` samples
- zarr parts only
- write `prepared_input.summary.json` with sample-level retained barcode counts
- CSR sparse storage
- require `prepared_input.summary.json` and `X.shape[0] == retained_barcodes_total` before trusting reuse or writing `prepared_input.ready`
- required merged `obs` columns:
  - `sample`
  - `patient`
  - `batch`
  - `disease`
  - `condition`
  - `sorting`
  - `sampling_site`
  - `sex`
  - `original_source_name`
  - `tumor_type`

Do not judge success only by directory existence or a stale ready sentinel; reuse must pass the summary-backed shape check.

For a successful run, verify:
- `run_manifest.json`
- `module_status.csv`
- `module_status.csv` has `ok` for:
  - `cellranger`
  - `qc`
  - `doublet_detection`
  - `clustering`
  - `annotation`
  - `composition`
  - `immune_phenotyping`
  - `tumor_microenvironment`
- required result tables exist:
  - `annotation/cell_type_annotation.csv`
  - `composition/composition_proportions.csv`
  - `immune_phenotyping/immune_phenotyping.csv`
  - `tumor_microenvironment/tme_scores_per_cell.csv`

Final controller states:
- `FINAL_STATUS=SUCCESS_LARGE`
- `FINAL_STATUS=SUCCESS_MASSIVE`
- `FINAL_STATUS=STOP_NO_FALLBACK`
- `FINAL_STATUS=STOP_AFTER_MASSIVE_FAILURE`

Current verified state for this lane:
- The frozen full-cohort prepared input validates at `884050 x 33538` with `81` samples and `retention_fraction = 0.001642391120829157`.
- Direct solo debug baseline: `results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222` has all stage-1 modules `ok`.
- Controller validation baseline: `results/NC2024_NSCLC_FULL_COHORT_STAGE1_WITH_FALLBACK.launch.log` now records truthful promotion from `large` to `massive`, and the promoted run `results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552` completes with all stage-1 modules `ok`.
- Operational rule: use direct `massive` for blocker-by-blocker debugging; use the controller when validating end-to-end orchestration.

Interpretation note:
- `large` is only a probe here.
- `massive` is the scale-protective fallback and should be expected when eager full-object loading or standard clustering paths exceed memory.

---

## 10. Beginner Parameter Cheat Sheet

### 10.1 High-impact, safe-to-change

- `--leiden-resolution` (default `0.8`)
  - higher: more/smaller clusters
  - lower: fewer/larger clusters

- `--n-pcs` (default `40`)
  - lower for speed, higher for complex datasets

- `--de-method` (default `wilcoxon`)
  - keep `wilcoxon` for most biological use cases

- `--batch-method` (default `harmony`)
  - only use when true multi-batch effect exists

### 10.2 QC thresholds

Defaults:
- `--min-genes 200`
- `--max-genes 7000`
- `--min-counts 500`
- `--max-counts 50000`
- `--max-mito-pct 20`
- `--max-ribo-pct 50`
- `--min-cells 3`

Change carefully and rerun QC/clustering for sanity.

---

## 11. Parallelism and Reliability Notes

- `--parallel-workers > 1` enables tiered parallel execution for safe appending modules.
- Structurally mutating modules run sequentially.
- Memory safety guard can reduce effective parallel worker count.
- Runtime telemetry is recorded in `run_manifest.json -> metadata.module_runtime_sec` and `metadata.pipeline_wall_seconds`.

---

## 12. RNA Velocity (Common Confusion)

To enable `rna_velocity`, provide at least one of:
1. `--velocity-loom`
2. `--velocity-bam` plus GTF resolution (`--velocity-gtf` or `--transcriptome-dir`)

Helpful details:
- Velocity extraction cache default: `/tmp/singlecell_factory_velocity_cache`
- Override cache directory:

```bash
export SCF_VELOCITY_CACHE_DIR=/path/to/fast_ssd_cache
```

---

## 13. Pseudobulk DE (Accuracy-Critical)

`pseudobulk_de` uses raw UMI counts from `adata.layers["counts"]`.

Current behavior:
- If counts layer is missing or invalid, module is skipped with explicit status.
- Single-sample input will skip pseudobulk comparison (expected behavior).

---

## 14. Basic Interpretation Checklist

After a run:
1. Check `module_status.csv` first.
2. Inspect `qc/qc_violin_post_filter.png`.
3. Inspect `clustering/umap_leiden.png`.
4. If annotation enabled, inspect `annotation/umap_cell_type.png`.
5. If DE enabled, inspect `differential_expression/marker_genes.csv` and volcano/heatmap outputs.
6. Open `run_manifest.json` for backends, runtimes, skip reasons, and metadata.

---

## 15. Troubleshooting (Beginner FAQ)

### 15.1 "Cell Ranger output not found"

- Verify `--sample-root` and `--outs-dir`.
- Confirm `filtered_feature_bc_matrix` exists.

### 15.2 "No checkpoint directory found"

- You must run once with `--checkpoint` before using `--resume-from`.

### 15.3 "Optional module skipped"

- This is often expected.
- Read skip message in `module_status.csv` / `run_manifest.json`.

### 15.4 "GPU not used"

- Pipeline auto-detects GPU backends.
- If unavailable or failure occurs, it falls back to CPU and records backend metadata.

### 15.5 "RNA velocity failed"

- Usually missing loom/BAM/GTF requirements.
- Re-run with valid `--velocity-bam` and `--transcriptome-dir` (or `--velocity-gtf`).

---

## 16. Reproducibility Protocol (Recommended)

For every production run:
1. Keep command line in a `run.sh` file.
2. Keep `run_manifest.json`, `module_status.csv`, and `final_adata.h5ad` together.
3. Record environment (`conda list > conda_env_export.txt`).
4. If using network-dependent modules (`validate_cbioportal`), note run date and connectivity.

---

## 17. Copy-Paste Templates

### 17.1 Local full run with checkpointing

```bash
python -m workflow.modular.cli \
  --project my_project_full \
  --sample-root data/raw/<your_sample> \
  --optional-modules clustering,cell_cycle,batch_correction,differential_expression,annotation,trajectory,pseudo_velocity,rna_velocity,cnv_inference,pathway_analysis,cell_communication,gene_regulatory_network,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,evolution,pseudobulk_de,cell_fate,composition,metacell \
  --velocity-bam data/raw/<your_sample>/outs/possorted_genome_bam.bam \
  --transcriptome-dir /path/to/ref \
  --checkpoint
```

### 17.2 Resume template

```bash
python -m workflow.modular.cli \
  --project my_project_full \
  --sample-root data/raw/<your_sample> \
  --optional-modules clustering,cell_cycle,batch_correction,differential_expression,annotation,trajectory,pseudo_velocity,rna_velocity,cnv_inference,pathway_analysis,cell_communication,gene_regulatory_network,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,evolution,pseudobulk_de,cell_fate,composition,metacell \
  --checkpoint \
  --resume-from <module_name>
```

---

## 18. Final Notes

- Start small (`clustering` only), verify quality, then scale to full modules.
- Treat `skipped` as informative, not automatically bad.
- Use checkpoint/resume for long runs.
- Prefer stable defaults unless you have a concrete biological reason to change parameters.

---

## Run Audit Ledger

Every modular pipeline run automatically writes a JSON audit record to `ops/run_ledger/<project>_<timestamp>.json` (schema version `run_ledger_v1`). The ledger is **observational only** — a failure to write it never aborts the pipeline run.

### What is captured

| Field | Description |
|---|---|
| `schema_version` | Always `"run_ledger_v1"` |
| `timestamp_utc` / `end_timestamp_utc` | ISO-8601 UTC start and end times |
| `project` | Project name passed via `--project` |
| `cli_args` | Full `sys.argv` at launch |
| `env` | Filtered env vars: `SC_*`, `SCF_*`, `CUDA_*`, `CONDA_*` prefixes + small standard allowlist |
| `git_sha` / `git_branch` / `git_dirty` | Repository state at launch |
| `conda_env` | `$CONDA_DEFAULT_ENV` |
| `python_version` | Full Python version string |
| `hostname` | Machine hostname |
| `sample_root` | Dataset root from config |
| `optional_modules` | Modules requested for the run |
| `module_results` | Per-module: name, status, message, wall_seconds, rss_peak_bytes |
| `total_wall_seconds` | End-to-end pipeline wall time |
| `peak_rss_bytes` | Process RSS peak sampled at 1-second intervals via psutil |
| `final_adata_sha256` | SHA-256 of `final_adata.h5ad` (empty string if not produced) |

### Implementation

- `workflow/modular/_run_ledger.py` — `RunLedger` class with `record_start()`, `record_module()`, `record_end()`, `write()`.
- Instantiated in `workflow/modular/cli.py` `main()` before `run_pipeline()`.
- Per-module hooks fire in `workflow/modular/pipeline.py` via `_ledger_record_module()` after every module in both sequential and parallel execution paths.
- All ledger calls are wrapped in `try/except`; warnings are logged but the run continues.

### Ledger output location

```
ops/run_ledger/<project>_<YYYYMMDDTHHMMSS>.json
```

Ledger files are tracked by git (see `.gitignore` exception `!ops/run_ledger/*.json`).

### Retention policy

Keep all ledger files indefinitely. Records are append-only and are never auto-pruned. They serve as the primary audit trail for production runs including the NC2024 sparse-exact 900k cohort (Phase 6).
