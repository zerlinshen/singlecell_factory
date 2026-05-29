# singlecell_factory Beginner Protocol (Complete, Practical)

## Authority / Read This First

AI agents should start with `AI_AGENT_PROTOCOL.md` before using this file. This
file remains the deep operational guide for running and recovering the modular
pipeline; it is not deprecated.

This protocol is a beginner-friendly, end-to-end guide for running the current modular scRNA-seq pipeline.

It is aligned with the current codebase (`workflow/modular/*`) and CLI (`python -m workflow.modular.cli`).

**Architecture note (2026-05+):** As of 2026-05, all pipeline runs must supply
`--project-root <path>` to write outputs outside the factory tree. See the
"Architecture (2026-05+)" section in `README.md` and the full plan at
`/home/zerlinshen/.omc/plans/factory-project-separation.md`. Bootstrap new
projects with `/home/zerlinshen/projects-bootstrap/omc-new-project <project-id>`.

**Round-1a governance (2026-05):** `README.md` now documents two new subsections
under "Architecture (2026-05+)": `### Environment Switches` (covering the
`SC_REQUIRE_PROJECT_ROOT` fail-fast gate) and `### R-factory SHA fields`
(covering dual-SHA provenance in `run_manifest.json` and `bundle/provenance.json`).
See `/home/zerlinshen/.omc/plans/factories-optimization-round1.md` for the
full Round-1a governance plan. The sibling repos `r_multiomics_factory` (R-native
analysis; renamed from `multiomics_r_factory` 2026-05-18) and `plotting_factory`
(dual-language plotting, introduced 2026-05-18; `python/` + `r/` subtrees +
`schema/` + vendored `contracts/figure_bundle_schema.yaml`) carry matching
documentation for the R-side contracts (renv bootstrap, schema hard-error,
SHA mismatch warning). See `ops/governance_records/2026-05-18-three-factory-trifurcation/ADR.md`
for the trifurcation decision record.

**Current validation ledger (2026-05-27):** The latest integration-biology and
multiomics validation artifacts are under
`/home/zerlinshen/projects/pipeline-validation-20260527/`, with repo-native
progress tracked in `.omx/ultragoal/ledger-integration-biology-multiomics-validation-20260527.jsonl`.
Use those reports and `ops/before_every_run/LATEST.md` before claiming marker
retention, mixing, rare-population preservation, annotation/DE, or 3D-genome
readiness. The current post-review state includes explicit marker/embedding
negative controls, G005 raw-count-compatible sample exclusions, and the G006
H3K27ac HiChIP low-information boundary. The final registered verdict is
`PASS_SUPPORTED_NOT_FINAL`.

### Remote Factory/Project Governance

Remote governance is implemented in `singlecell_factory` before any local Mac
sync. Use `scripts/validate_project_governance.py` to inspect project roots
read-only, and write governance reports only to factory control-plane locations
such as `ops/governance_records/`. Do not write validation reports into
canonical project run outputs.

The project run contract is documented in
`docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md` and
`contracts/project_run_contract.yaml`. It keeps root `manifest.json` separate
from producer-native Python/R/bundle provenance.


## Bridge Architecture

The canonical R source of truth is `r_multiomics_factory/`:

- `r_multiomics_factory/R/` — Seurat-based high-level analysis modules (12 files)
- `r_multiomics_factory/R_bundle/` — bundle-path-specific helpers (`remote_bundle_manifest.R`)

The bridge at `bridges/local_r_pipeline_macbook/` references these via symlinks:

- `bridges/local_r_pipeline_macbook/R` → `../../../r_multiomics_factory/R`
- `bridges/local_r_pipeline_macbook/R_bundle` → `../../../r_multiomics_factory/R_bundle`

**Rule**: `bridges/.../R` and `bridges/.../R_bundle` must always be symlinks, never real directories.
Verify with: `bash scripts/ci/check_bridge_symlink.sh`

Do NOT place real R files under `bridges/local_r_pipeline_macbook/R/` or `bridges/local_r_pipeline_macbook/R_bundle/`.
Edit R sources in `r_multiomics_factory/` only.

### Bundle export CLI (v2.1 + multi-modal flags)

Default schema is `singlecell_r_bundle_v2.1`. Force legacy with
`--schema-version v2`. Modality emission is opt-in:

```
python scripts/export_singlecell_r_bundle.py \
  --adata <run_dir>/final_adata.h5ad \
  --output <run_dir>/r_bundle \
  --schema-version v2.1 \
  --include-protein --protein-obsm-key protein_clr \
  --include-spatial --spatial-obsm-key spatial \
  --include-multimodal-obsm --multimodal-obsm-keys X_wnn X_mofa
```

Spatial library image paths are STRING-only (`--no-spatial-image-paths`
suppresses them); image bytes are never serialized into the bundle. The R
reader skips unknown extension keys with a `[singlecell_r_bundle_v2.1] skipping
unknown extension '<name>'` message for forward-compat.

### RSCRIPT_BIN env override

`tests/conftest.py` honors the `RSCRIPT_BIN` environment variable for
Python-to-R subprocess parity tests (`tests/test_python_r_parity.py`,
`tests/test_r_bundle_contract.py`). Default resolution remains the renv-pinned
`r_multiomics_arrow` env for production parity, while `r_multiomics` is also
parquet-capable as of 2026-05-29 (`r-arrow=24.0.0`) and can be selected
explicitly. Override when running tests against a non-default Rscript:

```
RSCRIPT_BIN=/path/to/Rscript pytest -q tests/test_python_r_parity.py
```

### Focused Python verification without coverage false failures

The root `pyproject.toml` intentionally enables coverage and a global
`fail_under` for full-suite pytest runs. For narrow contract or smoke checks,
append `--no-cov` so the focused lane fails only on the selected tests, not on
unrelated whole-project coverage scope:

```bash
pytest -q tests/test_bundle_sha_python_r_parity.py tests/test_r_bundle_contract.py tests/test_cli_project_root_smoke.py --no-cov
```

Do not remove the global coverage addopts or lower `fail_under`; full coverage
remains the separate broad regression gate via `pytest -q`.

---

## 0. Before Every Meaningful Remote Run

Read first:
- `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/LATEST.md`
- the newest relevant entry under `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/journal/`

Operational rule for the NC2024 full cohort:
- `execution_mode = debug_massive`
  - use direct `massive`
  - stop at first failing module
  - patch only that module
- `execution_mode = controller_validation`
  - use the official orchestration path `large -> massive`
  - treat `large` as a capacity probe, not the main completion lane

After a meaningful remote run finishes or fails:
- update `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/journal/`
- update `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/LATEST.md`
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
  - `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/luad/cell_communication/cell_communication_liana.csv`
  - `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/lusc/cell_communication/cell_communication_liana.csv`
  - `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_adenocarcinoma_checkpoint_top20.csv`
  - `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_squamous_cell_carcinoma_checkpoint_top20.csv`
- Use this lane to separate:
  - raw subtype LIANA pair presence
  - current paper-facing top20 visibility
- Keep the checkpoint verdict `partial` unless the paper-facing hierarchy itself is reproduced, not merely because the raw pairs are present.


---

## 0A. Opt-in environment variables (Principle 9 / F-3 pattern)

The audit Round-2 closing added explicit opt-in gates for previously silent
algorithmic compromises. None of these should ever be set in a
publication-claim run unless the operator has explicitly verified the
downstream consequences and recorded them in the run ledger.

| Variable | What it permits | Default | Recorded as |
|---|---|---|---|
| `SC_ALLOW_WELCH_FALLBACK=1` | `differential_expression` module: Welch t-test fallback when sparse-CPU DE path is missing | unset → strict mode raises | `ctx.metadata["de_welch_opt_in_acknowledged"] = True` |
| `SC_ALLOW_DOUBLET_NULL_FALLBACK=1` | `doublet_detection` module: all-singlets fallback when **both** rsc.scrublet AND CPU scrublet raise (failure-based path only; the tiny-dataset degenerate path remains exempt) | unset → strict mode raises | `ctx.metadata["doublet_null_fallback_opt_in_acknowledged"] = True`; `doublet_method_actually_used = "fallback_all_singlets_opt_in"` |
| `SC_ALLOW_CELLRANK_FALLBACK=1` | `cell_fate` module: connectivity-diffusion fallback when CellRank raises a non-ImportError exception (ImportError fallback always allowed and stamped) | unset → strict mode raises | `ctx.metadata["cell_fate_fallback_opt_in_acknowledged"] = True`; `adata.uns["cell_fate"]["engine"] = "fallback_connectivity_diffusion"` |
| `SC_ALLOW_BATCH_BACKEND_SKIP=1` | `batch_correction` module: silent skip when an explicitly-selected scvi/mnn/fastmnn backend fails (Harmony unchanged — always fail-loud) | unset → strict mode raises | `ctx.metadata["batch_correction_skip_opt_in_acknowledged"] = True`; `batch_correction_method_actually_used = "none_opt_in_skip"` |
| `SC_AMBIENT_TRIGGERS_DISABLE=1` | `ambient_correction` module: force-skip DecontX irrespective of QC triggers (paper-faithful reproduction mode) | unset → triggered-on policy | `adata.uns["ambient_correction"]["decision"] = "force_skip_env_disabled"` |
| `SC_REQUIRE_PROJECT_ROOT=1` | Pipeline entry: fail-fast on missing `--project-root` instead of falling back to legacy `output/` | unset → DeprecationWarning + legacy fallback | controller error exit |

A run that sets any of `SC_ALLOW_*_FALLBACK` should record:
- the env var actually set
- the resulting `*_actually_used` field from `run_manifest.json`
- the reviewer's explicit acknowledgement in `ops/before_every_run/`

---

## 0B. Tier-3 methodology dependencies (2026-05-22 ralplan close)

These conda envs and reference assets are required when running the
upgraded canonical modules. All are local, no further downloads needed.

| Module | Env / asset | Path |
|---|---|---|
| `chromvar` full mode (pychromvar) | `sc10x_methods` conda env (pychromvar 0.0.4, biopython 1.87) | `/home/zerlinshen/conda/envs/sc10x_methods/bin/python` |
| `chromvar` PWMs | JASPAR2024 CORE non-redundant | `singlecell_factory/data/references/JASPAR2024/JASPAR2024_CORE_non-redundant_pfms_meme.txt` (sha256 stored alongside) |
| `chromvar` genome | GRCh38 (Cell Ranger refdata) | `singlecell_factory/ref/reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa` |
| `multimodal_integration` MOFA | `sc10x_methods` (muon 0.1.7, mofapy2 0.7.4) | same env |
| `multimodal_integration` WNN | `r_multiomics` (Seurat 5.4.0) | `/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript` |
| `spatial_neighborhoods` | `sc10x_methods` (squidpy 1.8.1) | same env |
| `projection_module.R` | `r_multiomics` (Seurat 5.4.0, SingleR 2.12.0, celldex 1.20.0, harmony 1.2.4) | same Rscript |
| `projection_module.R` symphony backend | optional install (`remotes::install_github("immunogenomics/symphony")`) | pending GitHub rate-limit window |

Smoke tests (all exit 0 on 2026-05-22):

```bash
# chromVAR pychromvar
/home/zerlinshen/conda/envs/sc10x_methods/bin/python \
    "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/scripts/test_chromvar_pychromvar_smoke.py"

# multimodal MOFA
/home/zerlinshen/conda/envs/sc10x_methods/bin/python \
    "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/scripts/test_multimodal_wnn_mofa_smoke.py"

# spatial neighborhoods (Visium V1 Mouse Brain)
/home/zerlinshen/conda/envs/sc10x_methods/bin/python \
    "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/scripts/test_spatial_neighborhoods_visium_smoke.py"

# projection (R)
/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript \
    "/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory/scripts/test_projection_module_smoke.R"
```

Evidence record:
`governance/external_references/audits/tier3_methodology_uplift_2026-05-22.md`.

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
- Mandatory modules (always): `cellranger -> qc -> ambient_correction -> doublet_detection`
- Optional modules (22): selected by `--optional-modules` with auto dependency resolution.

Doublet backend discipline:
- `--doublet-backend scrublet` remains the global default because hgmm species-mix ground truth still favors Scrublet.
- Do not replace that default from one tissue benchmark. If Scrublet under-calls on heterogeneous tumor/tissue data, rerun a conditional second-opinion lane with `--doublet-backend scdblfinder` or `--doublet-backend consensus --doublet-consensus-pair scrublet_scdblfinder --doublet-consensus-logic or`.
- Treat consensus/scDblFinder recommendations as data-shape conditional until a second different-shape dataset supports a broader default change.

---

## 2. Prerequisites

### 2.1 Software

```bash
cd /home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory
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

The machine-readable hierarchy lives in
`workflow/modular/module_catalog.py`. It defines:

- module dependencies
- architectural layer
- modality/owner tags
- whether outputs are intended for R/report bundle consumption

`workflow/modular/pipeline.py` derives `MODULE_DEPENDENCIES` from that catalog
for backward compatibility. CLI help derives its optional-module list from the
same source. When adding or moving a module, update the catalog first, then the
implementation/registry/tests.

Mandatory:
- `cellranger`
- `qc`
- `ambient_correction`
- `doublet_detection`

Optional (23):
- `clustering`
- `cell_cycle`
- `integration_select` (opt-in via `--select-integration`; per-run discovery gate that sets `--batch-method` automatically, routes a Harmony pick to `harmony_backend=direct`, and fails loud on degraded candidates or a non-firing shuffle control)
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
- `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/results`

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

## 9. Paper Reproduction Ladder

Paper reproduction is faithful-first, then context-optimized.

1. Clone or stage the upstream paper repository, scripts, and supplementary methods when available.
2. Pin repository URL, commit/tag, DOI, data accessions, license, environment files, and raw/processed input boundaries.
3. Run raw-data reproduction first when public raw data exists and host capacity allows it.
4. If raw data is missing or impractical, use the earliest public computable input and label the boundary explicitly.
5. Reproduce data objects and figures separately. Data-object parity alone is not figure parity; figure parity alone is not object-level reproducibility.
6. Record figure/claim outcomes as exact, approximate, proxy, unsupported, or resource gap.
7. Map each paper method to factory capability:
   - existing module
   - parameter/config change
   - new reusable module
   - paper-specific script that should not enter the factory
8. Add or update modules only when the method should be reusable across projects.
9. Keep paper-faithful parameters separate from context-optimized defaults.
10. After faithful reproduction, optimize for our context: biological question, wet-lab decision, cohort scale, modality mix, memory limits, and downstream hypotheses.

Required project-policy fields for reproduction work:

- `upstream_repository`
- `raw_data_reproduction`
- `data_object_reproduction`
- `figure_reproduction`
- `module_gap_decisions`
- `context_optimization_decisions`

Use `workflow/modular/modules/paper_repro` or the project ledger for claim evidence; use `develop-and-integrate-module` / `singlecell-factory-module-delivery` when a paper method becomes a reusable factory module.

---

## 9A. NC2024 Full-Cohort Remote Reproduction

For the current `E-MTAB-13526` reproduction lane, do not use the older tumor-only scripts.

Use this owned full-cohort path instead:

```bash
cd /home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory
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
  - use `--leiden-resolution-sweep 0.5,0.8,1.0` on first-pass tuning runs to write a diagnostic table without changing final labels
  - if batch correction reruns Leiden, inspect `batch_correction/leiden_resolution_sweep.csv` because it reflects the corrected graph

- `--n-pcs` (default `40`)
  - lower for speed, higher for complex datasets

- `--de-method` (default `wilcoxon`)
  - keep `wilcoxon` for most biological use cases

- `--batch-method` (default `harmony`)
  - only use when true multi-batch effect exists
  - inspect `batch_correction/batch_mixing_metrics.json` with the before/after UMAPs; better correction should reduce same-batch neighbor fraction and increase normalized batch entropy without erasing real cell-type separation

- `--harmony-backend` (default `auto`; choices `auto`/`cpu`/`gpu`/`direct`)
  - `direct` calls the canonical `harmonypy.run_harmony` and stores `Z_corr.T` as `(n_cells, n_pcs)` — the proven-working path when the rapids (GPU CUBLAS) and scanpy-external (wrong-shape) wrappers fail on the host env
  - Principle 9 preserved: non-convergence raises unless `SC_ALLOW_HARMONY_NON_CONVERGENCE=1`

- `--select-integration` (opt-in; default off)
  - runs the `integration_select` discovery gate after clustering and BEFORE `batch_correction`; it scores required baseline/Harmony/scVI/shuffle candidates and SETS `cfg.batch.method` for THIS dataset (a Harmony pick also sets `harmony_backend=direct`)
  - fails loud if the required candidate set is degraded or the shuffle falsifiability control is absent/non-firing; extreme-theta is reported but not gating
  - expensive scVI seed sweep, so it is NOT default-on; the per-run cost is bounded by a default-ON cache keyed on baseline embedding, batch-label state, scVI config/training params, gate params, and gate code version
  - outputs (`integration_recommendation.json`, `integration_scoreboard.csv`, `integration_audit.md`) land under the run's `integration_select/` module dir, never inside the factory tree
  - tune the mixing gate with `--integration-margin-mix` (default 0.05)

### 10.2 QC thresholds

Defaults:
- `--min-genes 200`
- `--max-genes 7000`
- `--min-counts 500`
- `--max-counts 50000`
- `--max-mito-pct 20`
- `--max-ribo-pct 50`
- `--min-cells 3`

The QC defaults are fixed and configurable, not adaptive. Each run writes
`qc/qc_threshold_audit.json` with threshold values, metric quantiles, and
filter-failure counts. Do not promote new QC defaults from one benchmark alone;
use either a second dataset with a different shape or an explicit data-shape
conditional. Change carefully and rerun QC/clustering for sanity.

### 10.3 Annotation before CNV

`cnv_inference` depends on `annotation` in the module DAG. The CNV module writes
`cnv_inference/cnv_annotation_qc.json`, and annotation writes
`annotation/epithelial_marker_qc.json` plus EPICAM/KRT8/KRT18 summaries when
those genes are present. Treat epithelial subcluster and CNV conclusions as
blocked until these marker and annotation checks are plausible.

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
