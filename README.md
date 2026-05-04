# singlecell_factory v5.0 — Modular scRNA-seq Pipeline

## Authority / Read This First

AI agents must start with [AI_AGENT_PROTOCOL.md](AI_AGENT_PROTOCOL.md). That
file is the onboarding index; `AGENTS.md` / `CLAUDE.md` remain runtime-specific
authorities, and `PROTOCOL.md` remains the deep operational guide.

A comprehensive, production-ready single-cell RNA-seq analysis framework with **mandatory QC + 22 optional analysis modules + automatic dependency resolution + GPU acceleration + categorized output**.

Designed for 10X Genomics datasets. Tested on lung squamous cell carcinoma (LUSC) 3K cells.

Beginner entrypoint: see [PROTOCOL.md](PROTOCOL.md) for a complete step-by-step guide.

## Generated overview

Self-contained technical PDF + Markdown source covering both `singlecell_factory`
and `multiomics_r_factory`: see [docs/FACTORIES_OVERVIEW.pdf](docs/FACTORIES_OVERVIEW.pdf)
(rendered) and [docs/FACTORIES_OVERVIEW.md](docs/FACTORIES_OVERVIEW.md) (diff-friendly source).
Regenerate with `python scripts/generate_factories_report.py`.

## Publication & Reproducibility Documentation (NC2024)

| Document | Purpose |
|---|---|
| [ops/nc2024_methodology_audit/AUDIT_2026-04-26_v2.md](ops/nc2024_methodology_audit/AUDIT_2026-04-26_v2.md) | Parameter table, paper alignment, deliberate differences, reproducibility manifest |
| [ops/nc2024_methodology_audit/ALIGNMENT_REPORT_v2_2026-04-26.md](ops/nc2024_methodology_audit/ALIGNMENT_REPORT_v2_2026-04-26.md) | Cell type proportions vs Sanchez-Mejias 2024; v1 (broken) → v2 (fixed) cluster annotations |
| [ops/nc2024_methodology_audit/SMALL_REAL_VALIDATION_2026-04-26.md](ops/nc2024_methodology_audit/SMALL_REAL_VALIDATION_2026-04-26.md) | 100k staircase validation gate for the cluster_voting annotation fix |
| [ops/nc2024_methodology_audit/SEGFAULT_TRACE_2026-04-26.md](ops/nc2024_methodology_audit/SEGFAULT_TRACE_2026-04-26.md) | Root-cause + fix for the post-completion C-extension teardown segfault |
| [ops/nc2024_methodology_audit/PAPER_REPRO_REPORT_2026-04-27.md](ops/nc2024_methodology_audit/PAPER_REPRO_REPORT_2026-04-27.md) | Biology-level reproduction of the four core Sanchez-Mejias 2024 findings on the v2 cohort (verdict: 3 PASS / 1 PARTIAL) |
| [ops/MAC_PULL_RECIPE_2026-04-26.md](ops/MAC_PULL_RECIPE_2026-04-26.md) | Tailscale scp command + R load command for downstream plotting on Mac |
| [docs/PUBLICATION_READY.md](docs/PUBLICATION_READY.md) | Methods section template and citation patterns for manuscript drafting |

Canonical NC2024 outputs are the **v2** runs at `results/nc2024_tumor_20260426_v2/` and `results/nc2024_bh_20260426_v2/`. The matching v1 directories (without `_v2`) shipped with a known annotation labeling bug and must not be cited.

### Engineering Disciplines (Phase 7+)

**Capability flags** replace the multivalent `--scale-mode` (the old preset is still accepted as a backward-compat bundle):

- `--lazy-read {auto,true,false}` (auto = file > 5GB)
- `--doublet-strategy {auto,grouped,whole,skip}` (auto = grouped when n_obs ≥ 100k and a `sample` column is present)
- `--clustering-engine {auto,sparse_exact,css,gpu}`
- `--checkpoint-policy {full,mandatory_only,metadata_only}`
- `--annotation-strategy {cluster_voting,cell_argmax}` (default `cluster_voting`; `cell_argmax` retained as a fallback only — it drifts on >100k cohorts)

**Staircase test discipline** — all pipeline changes must be validated at increasing scale before being considered production-ready:

| Tier | Dataset | Cells (approx) | Marker |
|---|---|---|---|
| nano | synthetic CSR (conftest fixture) | ~5k | `pytest -m nano` |
| small_real | NC2024 subset (10 samples) | ~100k | `pytest -m small_real` |
| medium_real | NC2024 subset (50 samples) | ~500k | `pytest -m medium_real` |
| full_real | NC2024 full cohort | ~884k | `pytest -m full_real` |

A change that passes only nano/small_real is not cleared for full_real runs. Build fixtures with `scripts/build_staircase_fixtures.py`.

**Densify policy** — the pipeline enforces a grep-ban on unmarked `toarray()` / `todense()` calls in CI (`tests/test_densify_audit.py`). Any deliberate densification must carry a `# densify-allowed: <reason>` annotation on the same line, or route through `workflow/modular/_densify_policy.py:plan_densify()` which returns a `{GO, CHUNK, ABORT}` decision based on free memory and configured caps. This prevents silent memory explosions at 884k-cell scale. The Phase C modality modules (`protein_adt`, `spatial_neighborhoods`, `multimodal_integration`) carry `# densify-allowed: <reason>` annotations at every dense intermediate (protein panels are O(100) features, WNN UMAP is O(n_cells x 2), spatial neighborhood means are O(n_genes)).

**MemoryEnforcer cooperative abort** — `workflow/modular/_mem_guard.py` replaces the earlier observational MemoryGuard. A pre-flight RSS budget check + watchdog Event signals modules to abort at the next chunk boundary (raising `SkipModule`), rather than letting the kernel SIGKILL Python at the OOM threshold. Enable with `SC_MEM_GUARD=on SC_MEM_WATCHDOG=on`. Modules must not catch this exception.

**Shutdown cleanup** — `workflow/modular/_shutdown.py` runs explicit cupy / torch / zarr cleanup at interpreter exit, eliminating the post-completion C-extension teardown segfault that previously affected exit codes (run outputs were intact but `set -e` propagated exit 139 and skipped downstream stages). See `ops/nc2024_methodology_audit/SEGFAULT_TRACE_2026-04-26.md`.

**Pre-commit doc-sync gate** — `.githooks/pre-commit` (activated via `git config core.hooksPath .githooks`) runs two stages before every commit: (1) the existing reference manager that keeps the README citation list in sync with module-level `__references__` blocks, then (2) the global `repo-doc-sync` drift detector that validates README + `AGENTS.md` + `AI_AGENT_PROTOCOL.md` against the current canonical state (latest run dir, latest ledger, latest audit doc, missing DOIs, stale `Current State (YYYY-MM-DD)` blocks, missing v1 do-not-cite when v2 exists). The hook blocks commits when drift is found. The detector is installed globally at `~/.claude/skills/repo-doc-sync/` (also symlinked into `~/.codex/skills/` and `~/.kimi/skills/`). Bypass with `git commit --no-verify` only when drift is intentional and documented.

**Module catalog contract** — `workflow/modular/module_catalog.py` is the
single source for module dependency, architectural layer, modality, owner, and
R-bridge readiness metadata. `workflow/modular/pipeline.py` still exposes the
legacy `MODULE_DEPENDENCIES` shape for compatibility, but it is derived from the
catalog. CLI help also reads the catalog, so adding a module now requires one
catalog edit plus the normal implementation/registry/tests instead of scattered
README/CLI/DAG string updates.

**Singlecell → multiomics bridge contract** — compact bundle exports keep
`singlecell_factory` as the upstream truth and `multiomics_r_factory/R_bundle/`
as the downstream reader. v2 bundles now require manifest-backed file records
for parquet payloads and, when marker expression is exported as `mtx.gz`, the
barcode/gene sidecars are recorded with byte sizes and SHA256 values. This keeps
large R plotting/reporting handoff separate from full-object computation and
prevents sidecar drift.

**NC2024 architecture validation** — run
`python3 scripts/validate_nc2024_architecture_contract.py` for a no-rerun
controller-validation smoke. It checks current v2 run directories, r_bundle
manifests, run-ledger entries, and bridge symlinks without opening the 33G H5ADs.
Use `--verify-sha` when you specifically want bundle file hashing.

## Current Operational Defaults For NC2024-Style Full-Cohort Runs

- Read run memory first:
  - `/home/zerlinshen/singlecell_factory/ops/before_every_run/LATEST.md`
- Both operating modes are valid:
  - operate directly in this remote repo for compute, pipeline execution,
    run triage, and remote R/reporting work
  - orchestrate from the Mac over SSH when the task is coordination, review,
    handoff, or report-packaging oriented
- In both modes, this remote repo remains the run-truth surface. Prefer
  `run_manifest.json`, `module_status.csv`, `ops/before_every_run/LATEST.md`,
  and `ops/run_ledger/` over local summaries when deciding canonical state.
- For the NC2024 full cohort, treat `large` as a capacity probe rather than the main completion lane.
- Use direct `massive` for debug and recovery work.
- Use controller `large -> massive` only for orchestration validation.
- The canonical prepared input is:
  - `/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr`
- Historical stage-1 direct/controller success runs from `2026-04-23` were
  superseded for storage governance and then deleted after metadata archival:
  - archive:
    `/home/zerlinshen/singlecell_factory/ops/cleanup_records/nc2024_pre_extended_real_run_20260424T051809Z`
  - deleted direct success:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222`
  - deleted controller-fallback success:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552`
- The retained fresh stage-1 evidence run is:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
- The current clean full-cohort rerun source of truth is:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652`
  - launch log:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO.launch.log`
  - result summary: `810218 x 30374`, `33G` final H5AD, all requested
    modules `ok`, and pipeline-native `pseudobulk_de = ok/completed`
  - row-level reconciliation:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652/module_reconciliation.tsv`
  - remote R report bundle:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO_20260424_193652/r_plots/phase5_readable_20260424`
- The earlier extended full-cohort real-run is retained as predecessor
  evidence, not the current source of truth:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003`
  - launch log:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO.launch.log`
  - final object:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/final_adata.h5ad`
  - result summary: `810218 x 30374`, `33G` final H5AD, `20 ok` modules plus
    one original `pseudobulk_de` failure row recovered post-run from
    `final_adata.h5ad`
  - recovery outputs:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/pseudobulk_de_recovery`
  - remote R report bundle:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO_20260424_132003/r_plots/extended_real_run_main_20260424`
- `2026-04-25` methodology / optimization status:
  - authoritative audit:
    `/home/zerlinshen/singlecell_factory/ops/nc2024_methodology_audit/AUDIT_2026-04-25.md`
  - paper-aligned launcher:
    `/home/zerlinshen/singlecell_factory/scripts/run_nc2024_paper_aligned_20260425.sh`
  - sparse-exact exploratory launcher:
    `/home/zerlinshen/singlecell_factory/scripts/run_nc2024_full_cohort_sparse_exact_20260425.sh`
  - important changes: `sparse_exact` clustering path, densify guards,
    paper-aligned Scrublet/Harmony/Leiden/DE defaults, tumor-vs-background /
    healthy `--cohort-subset`, R bundle v2 subprocess contract tests, and
    staircase real-data fixtures
  - observed `NC2024_NSCLC_FULL_COHORT_SPARSE_EXACT_REAL_AUTO_*` result
    directories from `2026-04-25` are not canonical successful runs unless a
    later audit finds `final_adata.h5ad`, `run_manifest.json`, and
    `module_status.csv`; the current inspection saw only early mandatory
    outputs/checkpoints and a zero-byte sparse-exact launch log
- The local machine is now treated as an organization/review surface for processed outputs, not as a maintained local R pipeline.
- For reproduce-stage NC2024 full-cohort work, interpret "all modules" as
  "all eligible modules", not "every imaginable module regardless of modality
  or contrast contract".
- Current eligibility guardrails:
  - `rna_velocity` is eligible only when true splicing modality is available
    (`spliced/unspliced` layers, loom, or BAM+GTF).
  - `pseudobulk_de` is eligible only when raw counts exist and either:
    - an explicit confirmatory contrast contract is provided, or
    - exploratory `group_vs_rest` is explicitly enabled.
- Plotting boundary:
  - Python still emits module-native QC, diagnostic, and exact pipeline figures.
  - R is preferred for publication-style large-cohort figures where it is
    better: rasterized UMAPs, dot plots, heatmaps, composition panels, and
    presentation-ready figure boards.
- R plotting/reporting for large NC2024 outputs is remote-side:
  - R runtime: `/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript`
  - legacy R runtime retained for rollback: `/home/zerlinshen/conda/envs/r_multiomics/bin/Rscript`
  - bridge scripts: `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/scripts/`
  - historical note: the bridge folder name still says `local_r_pipeline_macbook`, but current operation is remote-first.
  - preferred command:
    ```bash
    bash bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh \
      /home/zerlinshen/singlecell_factory/results/<run> \
      /home/zerlinshen/singlecell_factory/results/<run>/r_plots/main \
      cell_type leiden
    ```
  - bridge controls for prettier/high-throughput remote plots:
    ```bash
    R_PLOT_THREADS=8 \
    R_BUNDLE_MARKERS=ELF3,EPCAM,KRT8,KRT18,PTPRC,CD3E,LYZ,MS4A1,NKG7 \
    bash bridges/local_r_pipeline_macbook/scripts/run_remote_bundle_plot.sh \
      /home/zerlinshen/singlecell_factory/results/<run> \
      /home/zerlinshen/singlecell_factory/results/<run>/r_plots/main \
      cell_type leiden 200000
    ```
  - reuse guard: the wrapper reuses a bundle only when source paths match,
    `final_adata.h5ad` size/mtime metadata match the manifest, optional
    `R_BUNDLE_MARKERS` / effective `R_BUNDLE_OBS_COLS` / effective
    `R_BUNDLE_OBSM` requests match, and the R validator passes file SHA256
    checks.
  - plotting entry point: `PLOT_SCRIPT` defaults to
    `/home/zerlinshen/multiomics_r_factory/scripts/plot_remote_bundle_large.R`,
    the upstream v1/v2-aware plotting script.
  - parameter contract: `R_BUNDLE_OBS_COLS` and `R_BUNDLE_OBSM` are additive,
    not replacement controls. The wrapper always keeps selected `group_by`,
    selected `cluster_by`, `X_umap`, and `X_pca`; ROI marker-dot plots use every
    numeric marker column exported in the marker-expression payload.
  - plotting guard: R UMAPs use rasterized points through `ggrastr` when
    available, so R can draw publication-style large-cohort figures without
    repeatedly loading the full expression matrix.
- Remote R environment hardening:
  - default v2 parquet plotting environment:
    `/home/zerlinshen/conda/envs/r_multiomics_arrow`
  - legacy rollback environment:
    `/home/zerlinshen/conda/envs/r_multiomics`
  - `r_multiomics_arrow` is a clone of `r_multiomics` with conda-forge
    `r-arrow 24.0.0` / `libarrow 24.0.0`; it was validated on `2026-04-30`
    against both NC2024 v2 bundles through `read_bundle()` and
    `scripts/plot_remote_bundle_large.R`.
  - installed and validated for NC2024 reporting:
    `Seurat 5.4.0`, `SeuratObject 5.4.0`, `readr`, `ggrastr`,
    `scattermore`, `pheatmap`, `ComplexHeatmap`, `circlize`, `hdf5r`,
    `harmony`, `BiocManager`, `R.utils`, `zellkonverter`, `remotes`, and
    `arrow`
  - validation artifact:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/r_env_dependency_smoke_20260424`
  - validation used the existing manifest-backed `r_bundle`, sampled `100000`
    cells for raster UMAP, generated `pheatmap` and `ComplexHeatmap` outputs,
    and opened `final_adata.h5ad` via `hdf5r` without converting the full object
    into Seurat.
  - `SeuratDisk` is intentionally not installed in this environment: the
    conda-forge package currently conflicts with R 4.5 / `zellkonverter`
    through old `spatstat` requirements. Use `zellkonverter`/`hdf5r` for small
    H5AD bridge checks and the compact bundle for NC2024-scale plotting.
  - bridge code-review hardening from `2026-04-24T04:10:07Z` produced:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_pretty_export_20260424`
    and verified reuse at:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_pretty_reuse_20260424`
  - Ralph follow-up at `2026-04-24T04:15:10Z` refreshed the canonical
    `r_bundle` itself with the stricter manifest keys and then verified default
    wrapper reuse:
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_refresh_20260424`
    and
    `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/r_plots/bridge_review_canonical_reuse_20260424`
- Artifact inventory cleanup from `2026-04-24T04:50:00Z` used
  `execution_mode=verification_cleanup`, not a pipeline rerun:
  - preserved canonical prepared input, prior successful direct/controller
    runs, the fresh `massive` run, fresh `final_adata.h5ad`, fresh `r_bundle`,
    and the fallback proof log
  - deleted only disposable `/tmp` bridge bundles/logs:
    `/tmp/nc2024_bridge_review_bundle_20260424`,
    `/tmp/nc2024_bridge_custom_contract_bundle_20260424`,
    `/tmp/nc2024_bridge_custom_contract_20260424.log`,
    `/tmp/nc2024_bridge_custom_contract_reuse_20260424.log`,
    `/tmp/nc2024_bridge_canonical_reuse_after_contract_fix_20260424.log`, and
    `/tmp/nc2024_remote_parse_after_deslop.log`
  - reclaimed about `226M`
  - no prepare, controller, benchmark, or full-cohort stage-1 run was launched
  - cleanliness alone is not a reason to rerun NC2024 full cohort; require a
    concrete reproducibility gap or acceptance criterion first
- Pre-extended-run cleanup from `2026-04-24T05:18:09Z` used
  `execution_mode=storage_governance_cleanup` before the approved full-cohort
  extended real run:
  - preserved raw/reference data, canonical prepared Zarr, fresh retained
    stage-1 run, cleanup metadata, launch logs, final H5ADs, manifests, and
    local-synced reproduction packages
  - deleted old superseded 2026-04-23 result directories and bulky checkpoint
    AnnData objects whose metadata sidecars had already been archived
  - reclaimed enough disk for the extended run without compromising the run
    ledger
  - deletion of test/reproduction result directories is acceptable only when
    each run is recorded and the current canonical artifacts are explicitly
    named in `ops/before_every_run/LATEST.md`
- Extended real-run execution from `2026-04-24T13:20:03+08:00` used
  `execution_mode=debug_massive`, direct `massive`, and
  `SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only`:
  - no prepare rerun was performed; input was the canonical prepared Zarr
  - `pydeseq2` API drift in `pseudobulk_de` was patched after the original
    pipeline row failed; a synthetic contrast smoke test and post-run recovery
    from `final_adata.h5ad` passed
  - original `module_status.csv` remains unedited and still records
    `pseudobulk_de` as failed; use the recovery directory above for recovered
    pseudobulk evidence
  - known fallbacks are evidence, not hidden failures: `pathway_analysis` used
    fallback when decoupler PROGENy was unavailable, `composition` used
    chi-squared fallback when pertpy/scCODA was unavailable, and `metacell`
    used MiniBatchKMeans fallback when SEACells was unavailable

## NC2024 Targeted Post-Baseline Evidence Lanes

Once the NC2024 full-cohort stage-1 baseline is already proven, prefer targeted evidence lanes over repeating stage-1.

- Do not rerun full-cohort stage-1 just to inspect subtype checkpoint hierarchy.
- Use the existing subtype communication outputs as inputs when the question is:
  - which paper-relevant checkpoint pairs exist in LUAD/LUSC raw subtype LIANA outputs
  - which of those pairs are hidden by the current checkpoint top20 summary surface
- Rerunnable checkpoint audit script:
  - `python scripts/audit_nc2024_subtype_checkpoint_pairs.py --luad-raw-csv ... --lusc-raw-csv ... --luad-top20-csv ... --lusc-top20-csv ... --output-dir ...`
- Higher-level checkpoint hierarchy summary:
  - `python scripts/summarize_nc2024_checkpoint_hierarchy.py --audit-dir <checkpoint-audit-run-dir>`
- One-command verification + regeneration wrapper:
  - `bash scripts/verify_nc2024_subtype_checkpoint_audit.sh`
- Current canonical example inputs:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/luad/cell_communication/cell_communication_liana.csv`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_CELLCOMM_AUTO_20260423_073528/lusc/cell_communication/cell_communication_liana.csv`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_adenocarcinoma_checkpoint_top20.csv`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_SUBTYPE_LR_FOCUS_AUTO_20260423_073900/lung_squamous_cell_carcinoma_checkpoint_top20.csv`
- This lane is intended to answer a paper-facing summary question, not to replace the full baseline or to silently upgrade the checkpoint verdict to `matched`.


## Features

- **3 mandatory modules** (cellranger, QC, doublet detection) ensure data quality baseline
- **22 optional analysis modules** covering the full scRNA-seq workflow
- Automatic topological dependency resolution — just list what you want, dependencies are auto-included
- **GPU acceleration** — auto-detected rapids-singlecell backend for clustering, batch post-processing, DE ranking, and evolution clone markers
- **Categorized output** — each module's figures and tables in its own subfolder
- **Checkpoint & resume** — zarr-accelerated checkpoints with h5ad fallback; `--resume-from` reuses the latest checkpointed run for the same project and backtracks to the nearest available checkpoint
- **Massive checkpoint policy** — set `SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only` with `--checkpoint` to keep per-module JSON sidecars without writing 10GB+ AnnData checkpoints for every full-cohort module.
- **Large-dataset modes** — `--scale-mode large|massive` selects safer defaults for 100k+ and several-hundred-thousand+ cell datasets without changing the underlying biological model family
- **Parallel execution** — thread-safe parallel tiers with cost-aware scheduling, merge-back safety warnings for structural mutations
- **Module contracts** — `requires_keys` / `provides_keys` declarations enable pre-flight validation; missing upstream data skips optional modules gracefully instead of crashing
- **Normalized status tracking** — module results are recorded as `ok` / `skipped` / `failed` in both `run_manifest.json` and `module_status.csv`
- **Module runtime telemetry** — per-module wall-time automatically stored in `run_manifest.json`
- **Raw-count integrity for pseudobulk** — `cellranger` stores raw UMI matrix in `adata.layers["counts"]`; `pseudobulk_de` consumes this layer only
- **Reference-aware annotation (optional)** — KNN label transfer from reference `h5ad` can override low-certainty marker labels
- Multi-backend support: each module auto-detects the best available tool
- Validated against cBioPortal mutation data
- Engineering principle: **accuracy and reproducibility first, performance second**

## Architecture

### Mandatory Modules (always run)

| Module | Function |
|---|---|
| `cellranger` | Load `filtered_feature_bc_matrix`, unified data entry |
| `qc` | Cell/gene filtering (mito/ribo/hemoglobin/n_genes/UMI), QC visualizations |
| `doublet_detection` | Scrublet doublet detection and removal |

### Optional Modules (auto-dependency resolution)

| Module | Function | Depends on |
|---|---|---|
| `clustering` | Normalization, HVG, PCA, UMAP, Leiden clustering | doublet_detection |
| `cell_cycle` | Cell cycle scoring (S/G2M), optional regression | clustering |
| `batch_correction` | Multi-sample batch correction (Harmony/BBKNN/Combat/Scanorama/scVI/MNN/fastMNN-style) | clustering |
| `differential_expression` | Cluster marker genes (`wilcoxon` default, configurable), significance filtering | clustering |
| `annotation` | Marker-based cell type annotation with confidence scores | clustering |
| `trajectory` | PAGA trajectory graph + DPT pseudotime + gene expression dynamics | clustering |
| `pseudo_velocity` | Pseudo-RNA velocity with arrow/stream plots (vectorized) | trajectory |
| `rna_velocity` | Real RNA velocity (scVelo: stochastic/dynamical) | clustering |
| `cnv_inference` | Expression-based CNV inference (infercnvpy / vectorized sliding window) | clustering |
| `pathway_analysis` | Gene set enrichment (gseapy/decoupler/built-in Hallmark, BH FDR) | differential_expression |
| `cell_communication` | Ligand-receptor interactions (LIANA / manual L-R scoring) | annotation |
| `gene_regulatory_network` | TF activity inference (decoupler+DoRothEA / manual TF-target) | clustering |
| `validate_cbioportal` | Cross-reference DE genes with cBioPortal mutation data | differential_expression |
| `immune_phenotyping` | 15 immune subtypes + exhaustion/cytotoxicity/activation scores | annotation |
| `tumor_microenvironment` | TME scoring (CYT/TIS/IFN-gamma/ESTIMATE) + checkpoint profiling | annotation |
| `gene_signature_scoring` | 10 built-in cancer signatures + custom JSON signatures | clustering |
| `evolution` | CNV-based clonal clustering, phylogenetic dendrogram, pseudotime-ordered evolution | cnv_inference + trajectory |
| `pseudobulk_de` | Pseudobulk DE (pydeseq2 / Mann-Whitney fallback) — statistically proper multi-sample comparison | differential_expression |
| `cell_fate` | Probabilistic cell fate mapping (CellRank / diffusion-based fallback) | trajectory |
| `composition` | Differential cell type composition analysis (pertpy/scCODA / chi-squared fallback) | annotation |
| `metacell` | Metacell aggregation (SEACells / MiniBatchKMeans fallback) — noise reduction for large datasets | clustering |
| `paper_repro` | Paper-driven reproduction ledger: track paper/repo/commit/license and validate figure parity against pipeline outputs | clustering |

### Module Dependency DAG

```
cellranger -> qc -> doublet_detection -> clustering -+-> differential_expression -+-> pathway_analysis
                                                      |                            +-> validate_cbioportal
                                                      |                            +-> pseudobulk_de
                                                      +-> annotation -+-> cell_communication
                                                      |               +-> immune_phenotyping
                                                      |               +-> tumor_microenvironment
                                                      |               +-> composition
                                                      +-> trajectory -+-> pseudo_velocity
                                                      |               +-> cell_fate
                                                      +-> cnv_inference --+--> evolution
                                                      |                        (requires both trajectory + cnv_inference)
                                                      +-> cell_cycle
                                                      +-> batch_correction
                                                      +-> rna_velocity
                                                      +-> gene_regulatory_network
                                                      +-> gene_signature_scoring
                                                      +-> metacell
                                                      +-> paper_repro
                                                      +-> protein_adt
                                                      +-> spatial_ingest --+--> spatial_neighborhoods
                                                      +-> multimodal_integration   [EXPERIMENTAL]
```

Standalone flowchart artifact (generated from `workflow/modular/pipeline.py`):

- SVG: `docs/module_dependency_flow.svg`
- PNG: `docs/module_dependency_flow.png`
- Regenerate:
  ```bash
  MPLCONFIGDIR=$PWD/.mplconfig NUMBA_CACHE_DIR=/tmp/numba_cache \
  PYTHONPATH=. python scripts/generate_module_flowchart.py
  ```

![Module Dependency Flow](docs/module_dependency_flow.png)

## Project Structure

```
singlecell_factory/
├── workflow/modular/
│   ├── cli.py              # CLI entry point
│   ├── config.py           # Dataclass configuration
│   ├── context.py          # Runtime context (per-module output dirs, checkpointing)
│   ├── pipeline.py         # Pipeline orchestration + dependency resolution
│   ├── perf_baseline.py    # Performance baselines
│   └── modules/            # 25 analysis modules
├── data/raw/               # Input datasets
├── results/                # Pipeline output (each run = timestamped folder with analysis subfolders)
├── tests/                  # Test suite
├── reports/                # Analysis reports
├── environment.yml         # Conda environment
└── pyproject.toml          # Python project config
```

## Installation

```bash
cd /home/zerlinshen/singlecell_factory
conda env create -f environment.yml
conda activate sc10x
export MPLCONFIGDIR=$PWD/.mplconfig
export NUMBA_CACHE_DIR=/tmp/numba_cache
```

`MPLCONFIGDIR` and `NUMBA_CACHE_DIR` are strongly recommended for stable `scanpy/scvelo` startup in some Conda environments.

Optional GPU environment (NVIDIA):

```bash
conda env create -f environment_gpu.yml
conda activate sc_gpu
export MPLCONFIGDIR=$PWD/.mplconfig
export NUMBA_CACHE_DIR=/tmp/numba_cache
```

## Dependency Requirements

- Core: `python>=3.10`, `scanpy`, `anndata`, `numpy`, `scipy`, `pandas`, `matplotlib`, `scikit-learn`
- Mandatory-module runtime: `scrublet`
- Optional backends (auto-detected at runtime):
  - `rapids-singlecell`, `cupy` (GPU clustering + batch post-processing + DE + evolution clone-marker ranking)
  - `harmonypy` / `bbknn` / `scanorama` / `scvi-tools` / `mnnpy` (batch correction)
  - `infercnvpy`, `pybiomart` (CNV)
  - `gseapy`, `decoupler` (pathway / TF activity)
  - `liana` (cell-cell communication)
  - `cellrank` (cell fate)
  - `pydeseq2` (pseudobulk DE)
  - `SEACells` (metacell)
  - `scvelo`, `pysam` (RNA velocity)

If a backend is missing, the corresponding module falls back to an implemented alternative whenever possible.

## Data Preparation (LUSC)

The LUSC 3K dataset should be placed at:

```
data/raw/lung_carcinoma_3k_count/outs/filtered_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

## Quick Start

### User-friendly entrypoint: `scfactory`

`scripts/scfactory.py` is a thin UX wrapper around the canonical
`python -m workflow.modular.cli` invocation. It auto-detects modality
(RNA-only / CITE-seq / spatial / multimodal) from an `.h5ad` file or a
sample-root directory, picks sensible optional modules, and can also
dispatch the v2.1 R bundle export. It does **not** change pipeline
behavior — pass `--optional-modules` to override the auto plan, or
`--dry-run` to preview the underlying CLI invocation.

```bash
# preview the planned modular invocation for an .h5ad input
python scripts/scfactory.py run path/to/sample.h5ad --dry-run

# run + export an R bundle with modality-aware defaults
python scripts/scfactory.py run data/raw/lung_carcinoma_3k_count \
    --project demo_run --bundle

# read-only environment health check (Rscript, bridges, deps, last run)
python scripts/scfactory.py doctor          # human-readable
python scripts/scfactory.py doctor --json   # machine-readable, exits non-zero on FAIL
```

#### Recipes (presets)

Recipes pre-package optional modules + capability-flag env vars + bundle
config for common workflows. Pass `--recipe NAME` to `scfactory run`; precedence is
`--optional-modules` > `--recipe` > auto-detect. Recipe `env` is applied
to the subprocess only (parent env untouched). Requires `pyyaml`
(only when `--recipe` / `--list-recipes` is used).

```bash
python scripts/scfactory.py run --list-recipes        # list all
python scripts/scfactory.py run sample.h5ad --recipe quick_explore --dry-run
```

Starter recipes in `recipes/`:

- `quick_explore` — RNA-only first look: clustering + DE, bundle on (most users start here).
- `nc2024_paper` — NSCLC paper-faithful repro: clustering + DE + annotation + trajectory + paper_repro, `scale-mode=massive`, `SC_CLUSTERING_ENGINE=sparse_exact`, bundle off by default.
- `cite_seq_full` — CITE-seq RNA + ADT (CLR) with bundle protein extension.
- `visium_neighborhoods` — Visium spatial: ingest + neighborhoods (squidpy) + bundle spatial extension.

### Execution Profiles (Human + AI)

Use one of the following copy-paste profiles directly.

Eligibility note:
- Only include `rna_velocity` when true splicing inputs are available.
- Only include `pseudobulk_de` as a confirmatory module when an explicit
  contrast contract is provided. Otherwise treat it as exploratory and enable it
  intentionally.

1. **Full local analysis (recommended, no external network dependency)**

```bash
python -m workflow.modular.cli \
  --project LUSC_full_local \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,cell_cycle,batch_correction,differential_expression,annotation,trajectory,pseudo_velocity,rna_velocity,cnv_inference,pathway_analysis,cell_communication,gene_regulatory_network,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,evolution,pseudobulk_de,cell_fate,composition,metacell \
  --velocity-bam data/raw/lung_carcinoma_3k_count/outs/possorted_genome_bam.bam \
  --transcriptome-dir ref/reference/refdata-gex-GRCh38-2024-A \
  --pseudobulk-exploratory-group-vs-rest
```

2. **Full analysis with online cancer-cohort validation**

```bash
python -m workflow.modular.cli \
  --project LUSC_full_online \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,cell_cycle,batch_correction,differential_expression,annotation,trajectory,pseudo_velocity,rna_velocity,cnv_inference,pathway_analysis,cell_communication,gene_regulatory_network,validate_cbioportal,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,evolution,pseudobulk_de,cell_fate,composition,metacell \
  --velocity-bam data/raw/lung_carcinoma_3k_count/outs/possorted_genome_bam.bam \
  --transcriptome-dir ref/reference/refdata-gex-GRCh38-2024-A \
  --pseudobulk-exploratory-group-vs-rest
```

3. **Fast baseline (no RNA velocity)**

```bash
python -m workflow.modular.cli \
  --project LUSC_fast \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression,annotation,trajectory,pseudo_velocity,cnv_inference,pathway_analysis
```

### Large Dataset Modes

Use `--scale-mode` to switch from convenience defaults to memory-safer presets:

| Mode | Intended scale | Default optional modules when not overridden | Key parameter shifts |
|---|---|---|---|
| `standard` | up to ~100k cells | `clustering,differential_expression,annotation,trajectory,pseudo_velocity` | full default behavior |
| `large` | ~100k-300k cells | `clustering,annotation,differential_expression` | fewer HVGs/PCs, slightly lower graph density |
| `massive` | ~300k to ~1M cells | `clustering` | clustering-first pass, reduced HVGs/PCs/neighbors, lighter memory footprint; auto-prefers CSS-style representation when sample labels are available |

Examples:

```bash
python -m workflow.modular.cli \
  --project NSCLC_100K_large \
  --sample-root data/raw/my_large_dataset \
  --scale-mode large \
  --checkpoint
```

```bash
python -m workflow.modular.cli \
  --project NSCLC_900K_massive \
  --sample-root data/raw/my_massive_dataset \
  --scale-mode massive \
  --checkpoint
```

Recommended interpretation:
- `~100k cells` is already a large dataset, but it is not an extreme scale for a 96 GB workstation if runs are staged sensibly.
- `300k+ cells` is a very large dataset and should usually start with `--scale-mode massive` or a staged subset-first workflow.
- `~1M cells` is an extreme scale for classic full-object scRNA workflows and should be treated as clustering-first, then subset/refine later.
- In `--scale-mode massive`, the pipeline now prefers a CSS-style clustering representation when `adata.obs["sample"]` is available; otherwise it automatically falls back to the lighter clustering-first route.

### Memory Pressure and Swap Guidance

When large runs approach RAM limits:

- Prefer reducing module scope first (`--scale-mode large` or `--scale-mode massive`) before changing hardware assumptions.
- A moderate SSD-backed swapfile (for example `30-40 GB`) can help prevent abrupt OOM kills and make checkpoint-heavy runs more forgiving.
- On this workstation, Harmony batch correction is intentionally routed to CPU in normal `auto/off` usage because the current `harmonypy` wrapper path was unstable under the previous GPU route; the core Harmony algorithm itself remains valid and was verified separately on both CPU and CUDA.
- Current GPU reality on this workstation is module-specific rather than globally broken: basic CUDA, PyTorch CUDA, CuPy, and Harmony core all work; the main unstable paths are RAPIDS PCA (`CUSOLVER_STATUS_INTERNAL_ERROR`) and RAPIDS DE (`CUBLAS_STATUS_NOT_INITIALIZED`) on real workloads.
- For clustering, the pipeline now supports a hybrid fallback path: CPU PCA followed by GPU neighbors/UMAP/Leiden when GPU PCA is the failing substep.
- Swap is a safety buffer, not real RAM: it may keep a run alive, but it can become much slower once the workflow starts paging heavily.
- For large but not extreme runs (for example around `100k` cells), swap can be a useful bridge before a RAM upgrade.
- For extreme runs (several hundred thousand to ~1M cells), swap alone is usually not enough; use staged execution, lighter module sets, subset-first refinement, and consider adding RAM.
- Put swap on fast NVMe storage when possible.
- If a run is repeatedly `Killed`, treat that as a sign to lower memory pressure first, then add swap, then consider hardware upgrades.

Practical order of operations for large datasets:
1. Switch to `--scale-mode large` or `--scale-mode massive`.
2. Reduce modules to the minimum biologically necessary first pass.
3. Add `30-40 GB` swap as an OOM safety net if RAM is tight.
4. Split downstream analysis into subset-specific second-stage runs.
5. Upgrade RAM when very large runs become routine rather than exceptional.

## NC2024 Full-Cohort Remote Lane

The NC2024 `E-MTAB-13526` reproduction lane now has a dedicated full-cohort remote path.
Use these owned scripts instead of the older tumor-only launchers:

- `scripts/prepare_emtab13526_full_cohort_zarr.py`
- `scripts/run_emtab13526_full_cohort_stage1.sh`
- `scripts/run_emtab13526_full_cohort_with_fallback.sh`

Behavior contract:
- build one unsplit `full_cohort/prepared_input.zarr` covering all `81` samples
- emit `full_cohort/prepared_input.summary.json` with per-sample retained-barcode accounting
- use zarr parts only, never `h5ad` parts
- validate CSR sparse storage and require `X.shape[0] == retained_barcodes_total` before writing `prepared_input.ready`
- treat stale `prepared_input.zarr` directories without the matching summary JSON as invalid for reuse
- run one `large` stage-1 probe first
- auto-promote once to `massive` only on explicit capacity failures

Primary human entrypoint:

```bash
bash scripts/run_emtab13526_full_cohort_with_fallback.sh
```

The controller owns:
- preflight inventory capture
- cleanup of stale derived artifacts
- prepare-if-missing or validate-only prepare checks
- summary-backed prepared-input reuse checks (directory + ready sentinel + `prepared_input.summary.json`)
- one `large` run
- one optional `massive` retry
- final terminal status logging

Final controller states:
- `FINAL_STATUS=SUCCESS_LARGE`
- `FINAL_STATUS=SUCCESS_MASSIVE`
- `FINAL_STATUS=STOP_NO_FALLBACK`
- `FINAL_STATUS=STOP_AFTER_MASSIVE_FAILURE`

Current verified NC2024 status:
- Frozen prepared input is validated at `884050 x 33538` with `81` samples and `retention_fraction = 0.001642391120829157`.
- A direct solo `massive` stage-1 run now succeeds end-to-end at `results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_CLUSTER_FIX_AUTO_20260423_031222`.
- A full controller validation also succeeds truthfully: `large` is capacity-killed after `qc`, controller logs `LARGE_FAILED_CAPACITY_PROMOTING_TO_MASSIVE`, and the promoted `massive` run succeeds at `results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_AUTO_20260423_035552` with `FINAL_STATUS=SUCCESS_MASSIVE`.
- Recommended usage now is: direct `massive` for module debugging, controller for production-like orchestration validation.

Success for either mode is defined by both:
- `run_manifest.json`
- `module_status.csv`

Required `module_status.csv == ok` modules:
- `cellranger`
- `qc`
- `doublet_detection`
- `clustering`
- `annotation`
- `composition`
- `immune_phenotyping`
- `tumor_microenvironment`

Required downstream artifacts:
- `annotation/cell_type_annotation.csv`
- `composition/composition_proportions.csv`
- `immune_phenotyping/immune_phenotyping.csv`
- `tumor_microenvironment/tme_scores_per_cell.csv`

Important caveat:
- `large` is only a bounded probe in this lane.
- The protected scale-aware fallback is `massive`, because lazy zarr loading, grouped scrublet, and CSS sparse/SVD protections are `massive`-only in the current codebase.

### Sample-wise / Batch-wise Staging for Massive Cohorts

This is a common and recommended strategy for very large datasets, especially multiplex or cohort-scale studies.

Important distinction:
- The reason to stage by sample/batch is not only to save RAM.
- It also improves data hygiene by applying QC and doublet handling before all cells are merged into one giant object.

Use this strategy when the dataset contains multiple samples, donors, batches, lanes, or multiplexed barcodes. For example, the 900k NSCLC example includes many per-sample entries in its aggregation sheet rather than one biologically homogeneous sample.

Practical note for this pipeline:
- CSS-style integration needs per-cell sample labels in `adata.obs["sample"]`.
- If an aggregated matrix lacks those labels, massive mode will not force CSS blindly; it will fall back to the simpler clustering-first path until sample mapping is restored.
- For public multiplex datasets that expose stable probe-barcode groups but not full per-cell sample assignment, a probe-group proxy can still be used as an engineering-stage CSS approximation for large-data staging. This is useful for RAM-constrained first-pass processing, but it should be described as a proxy strategy rather than a final publication-grade sample assignment.

Recommended workflow:
1. Run the same QC logic per sample or per multiplexed sample group.
2. Remove low-quality cells and doublets before full aggregation when possible.
3. Build a first-pass clustering or reference on a reduced representative object.
4. Integrate or map the remaining cells back onto that reference.
5. Run heavier downstream modules on biologically relevant subsets rather than the full giant object.

How this relates to batch effects:
- Batch effects are not created by QC itself.
- Batch effects usually come from differences in donor, library prep, lane, chemistry, run, or processing time.
- Applying the same QC standard across samples helps reduce technical noise and makes later integration more consistent, but it does not by itself remove batch effects.
- Proper integration or batch-aware modeling is still needed when multiple samples are combined.

Practical guidance for a 96 GB workstation:
- Around `100k` cells: often manageable with staged execution.
- Several hundred thousand cells: strongly prefer sample-wise staging plus `--scale-mode massive`.
- Near `1M` cells: do not treat as an ordinary single-object run; use clustering-first, reference mapping, and subset refinement.

### Relationship to multiomics_r_factory

`singlecell_factory` is the upstream analysis engine.

It is responsible for:
- raw and processed single-cell data handling
- checkpointed analysis runs
- module outputs and analysis artifacts
- producing the result directories that downstream reporting consumes

`multiomics_r_factory` is the downstream R/report workspace that depends on outputs generated here.

Repository links:
- `singlecell_factory`: <https://github.com/zerlinshen/singlecell_factory>
- `multiomics_r_factory`: <https://github.com/zerlinshen/multiomics_r_factory>

Practical dependency direction:
- `singlecell_factory` -> `multiomics_r_factory`

That means:
- large objects and primary analysis should originate here
- downstream R plotting/report work is remote-side by default, either through the bridge mirror under `singlecell_factory` or the broader remote `multiomics_r_factory`
- the local Mac is a review/organization surface, not the maintained R plotting runtime for NC2024-scale work
- the two remote workspaces should be treated as linked analysis/report layers rather than unrelated repositories

### Remote R Pipeline Bridge

A remote R plotting/report workflow is kept in two forms:
- Bridge mirror inside singlecell_factory:
  - `/home/zerlinshen/singlecell_factory/bridges/local_r_pipeline_macbook/`
- Recommended independent remote R workspace:
  - `/home/zerlinshen/multiomics_r_factory/`

Rationale:
- `singlecell_factory` should remain the main compute/analysis engine for remote single-cell workflows.
- The R layer is broader than scRNA-seq alone and may later cover spatial transcriptomics, polished publication graphics, multi-omics summaries, and other R-native plotting/report tasks.
- Therefore the long-term cleaner architecture is: analysis engine (`singlecell_factory`) + independent R workspace (`multiomics_r_factory`) + explicit bridge between them.

Intended use:
- Use the bridge mirror inside `singlecell_factory` when tight co-location with pipeline outputs is convenient.
- Use `/home/zerlinshen/multiomics_r_factory/` as the preferred long-term home for broader R analysis and figure workflows.
- Use compact manifest-backed bundles rather than forcing direct `.h5ad` conversion in R for NC2024-scale cohorts.
- Keep Mac-side work focused on reviewing and organizing the remote-generated figures/reports.

#### Bundle schema v2.1 (additive over v2)

The default emitter in `scripts/export_singlecell_r_bundle.py` writes
`schema_version = "singlecell_r_bundle_v2.1"`. v2.1 is fully backwards-compatible
with v2: the only addition is an `extensions: dict[str, dict]` manifest field.
Force legacy emission with `--schema-version v2`.

Known extension keys:

| Extension key | Modality | Producer (Python) | Reader (R) |
|---|---|---|---|
| `protein` | CITE-seq / ADT | `maybe_export_protein(...)` | `multiomics_r_factory/R/protein_module.R::load_protein_extension(bundle)` |
| `spatial` | spatial transcriptomics | `maybe_export_spatial(...)` | `multiomics_r_factory/R/spatial_module.R::load_spatial_extension(bundle)` |
| `multimodal_obsm` | WNN / MOFA embeddings (EXPERIMENTAL) | `maybe_export_multimodal_obsm(...)` | `multiomics_r_factory/R/integration_module.R::load_multimodal_extension(bundle)` |

Public Python API for registering an extension:

```python
from scripts.export_singlecell_r_bundle import add_extension
add_extension(manifest, "protein", version="1.0", files=[...], **fields)
```

`add_extension(...)` is idempotent for the same `name`. The R reader accepts
both `singlecell_r_bundle_v2` and `singlecell_r_bundle_v2.1`
(`ACCEPTED_V2_SCHEMAS`); unknown extension keys are skipped with
`[singlecell_r_bundle_v2.1] skipping unknown extension '<name>'` for
forward-compat. Per-modality loaders are fail-soft (`data=NULL` on load
failure, no exception). The bundle-level `claim_guard` string also applies at
the modality level until per-modality validation exists.

#### Bridge contract (Phase A hardening)

- **Atomic bundle export**: writer emits `<output>.tmp.<pid>` then `os.replace`,
  manifest written last; partial bundles are rolled back on failure.
- **Manifest sweep on the R side**: `validate_bundle_v2_manifest` checks
  manifest presence and per-file size; partial bundles are refused.
- **claim_guard surfacing**: the R reader emits the guard via `message()` to
  stderr and attaches it as `attr(expr_sparse, "claim_guard")` on the returned
  matrix.
- **h5ad backend order**: `zellkonverter` -> `SeuratDisk` ->
  `Seurat::ReadH5AD` (last; opt-in via `options(h5ad.allow_seurat_legacy = TRUE)`).
- **Test override**: `RSCRIPT_BIN` env var overrides the default Rscript binary
  in `tests/conftest.py` (default `r_multiomics_arrow` for production parity).

#### Bundle export CLI flags (new)

| Flag | Default | Purpose |
|---|---|---|
| `--schema-version {v1,v2,v2.1}` | `v2.1` | Force a specific bundle schema |
| `--include-protein` | off | Emit `protein` extension when present on AnnData |
| `--protein-obsm-key` | `protein_clr` | obsm key holding the protein matrix |
| `--protein-isotype-controls` | empty | Comma-separated isotype-control proteins |
| `--include-spatial` | off | Emit `spatial` extension when coords are present |
| `--spatial-obsm-key` | `spatial` | obsm key holding (x,y) coordinates |
| `--no-spatial-image-paths` | off | Suppress library image-path strings (paths only, never image bytes) |
| `--include-multimodal-obsm` | off | Emit `multimodal_obsm` (EXPERIMENTAL) |
| `--multimodal-obsm-keys` | `X_wnn X_mofa` | obsm keys to register as multimodal embeddings |

#### Cross-language parity sentinels

Two hard gates live in `tests/test_python_r_parity.py`: top-marker Jaccard
overlap >= 0.6 and PCA cosine similarity >= 0.95 between the Python AnnData
output and the R-side bundle reconstruction. The DAG-count assertion in
`tests/test_modular.py` is set to 29 modules (3 new bridge_ready modules
registered alongside the existing 25 + 1 paper_repro variants).

#### Known issues

- `tests/test_singlecell_r_bundle_export.py` carries 2 pre-existing v1-schema
  failures (asserts `schema_version == "singlecell_r_bundle_v1"` and reads
  `obs.csv.gz`). The default emitter is now v2.1/parquet; these are tracked
  separately and are not part of the v2.1 schema work.

### AI / Automation Checklist

Before launching a run, verify:

1. `--sample-root` exists and includes `outs/filtered_feature_bc_matrix`.
2. If `rna_velocity` is enabled, pass `--velocity-bam` explicitly.
3. If `rna_velocity` is enabled without `--velocity-gtf`, ensure `--transcriptome-dir` points to a reference containing `genes/genes.gtf(.gz)` or `genes.gtf(.gz)`.
4. In restricted-network environments, remove `validate_cbioportal` from `--optional-modules`.
5. Export `MPLCONFIGDIR` and `NUMBA_CACHE_DIR` to avoid startup/cache issues.

For Codex / Claude Code agents:
- project Codex skills live in `.codex/skills/` using standard Codex project
  skill management
- project Claude skills live in `.claude/skills/`
- reusable global skills live under `~/.codex/skills/` and
  `~/.claude/skills/`
- `codex_skills/` remains only as a legacy compatibility mirror for historical
  project-local skills

### Paper-Driven Continuous Optimization

To continuously improve the pipeline by learning from external papers/repos, enable `paper_repro` and provide a spec file:

```bash
python -m workflow.modular.cli \
  --project LUSC_paper_repro \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,paper_repro \
  --paper-spec-json docs/paper_repro/paper_spec.json
```

`paper_repro` outputs:
- provenance ledger (paper path, repo URL, commit, license)
- figure parity checks against current run outputs
- machine-readable reproduction report for iterative pipeline optimization

### Full Analysis (recommended)

```bash
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,cell_cycle,differential_expression,annotation,\
trajectory,pseudo_velocity,cnv_inference,pathway_analysis,cell_communication,\
gene_regulatory_network,immune_phenotyping,tumor_microenvironment,\
gene_signature_scoring
```

### With Checkpointing (crash recovery)

```bash
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression,annotation \
  --checkpoint
```

Resume from a failed module:

```bash
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression,annotation \
  --checkpoint --resume-from annotation
```

Resume behavior details:
- `--resume-from` requires an existing `.checkpoints/` directory from a prior run of the same `--project`.
- The pipeline reuses the latest run directory for that project that contains checkpoints.
- If checkpoint `after_<resume-from-1>` is missing, it automatically searches backward for the nearest earlier checkpoint in execution order.
- If no earlier checkpoint exists, the run exits with an explicit `FileNotFoundError`.

### Parallel Execution

```bash
python -m workflow.modular.cli \
  --project LUSC_3k_Analysis \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression,annotation,immune_phenotyping \
  --parallel-workers 4
```

### Immuno-Oncology Deep Analysis

```bash
python -m workflow.modular.cli \
  --project LUSC_immuno \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,differential_expression,annotation,\
immune_phenotyping,tumor_microenvironment,gene_signature_scoring,pathway_analysis
```

### RNA Velocity from BAM (no loom file)

```bash
python -m workflow.modular.cli \
  --project LUSC_velocity \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,rna_velocity \
  --velocity-bam data/raw/lung_carcinoma_3k_count/outs/possorted_genome_bam.bam \
  --transcriptome-dir /path/to/refdata-gex-GRCh38-2020-A
```

If `--transcriptome-dir` is set, `genes.gtf(.gz)` is auto-discovered from:
`<transcriptome-dir>/genes/genes.gtf(.gz)` (or `<transcriptome-dir>/genes.gtf(.gz)`).
You can still pass `--velocity-gtf` explicitly to override auto-discovery.

## Output Structure

Each run creates an independent timestamped folder under `--output-dir` (CLI default: `/home/zerlinshen/singlecell_factory/results`):

```
<output-dir>/<project>_<timestamp>/
├── final_adata.h5ad
├── run_manifest.json
├── module_status.csv
├── .checkpoints/                    # only when --checkpoint is enabled
└── <module_name>/                   # one subfolder per executed module
```

Current module output folders (as implemented):
- `annotation`: `cell_type_annotation.csv`, `cluster_majority_cell_type.csv`, `umap_cell_type.png`
- `batch_correction`: `umap_batch_before.png`, `umap_batch_after.png`
- `cell_communication`: `cell_communication_liana.csv`, `cell_communication_lr.csv`, `cell_communication_dotplot.png`, `cell_communication_heatmap.png`
- `cell_cycle`: `cell_cycle_scores.csv`, `cell_cycle_umap.png`
- `cell_fate`: `fate_probabilities.csv`, `terminal_states.csv`, `fate_heatmap.png`, `fate_umap_cellrank.png`
- `clustering`: `pca_variance_explained.png`, `umap_leiden.png`
- `cnv_inference`: `cnv_scores.csv`, `cnv_classification.json`, `cnv_score_umap.png`, `cnv_heatmap.png`
- `composition`: `composition_counts.csv`, `composition_proportions.csv`, `composition_test_results.csv`, `composition_barplot.png`, `composition_boxplot.png`
- `differential_expression`: `marker_genes.csv`, `marker_genes_all.csv`, `marker_top5_by_cluster.csv`, `de_dotplot_top5.png`, `de_heatmap_top5.png`, `de_volcano.png`
- `doublet_detection`: `doublet_scores.png`
- `evolution`: `evolution_clone_assignment.csv`, `evolution_clone_stats.csv`, `evolution_clone_markers.csv`, `evolution_clone_umap.png`
- `gene_regulatory_network`: `tf_activity_per_cell.csv`, `tf_activity_per_cluster.csv`, `tf_top_per_cluster.json`, `tf_activity_heatmap.png`, `tf_activity_umap.png`
- `gene_signature_scoring`: `gene_signature_scores.csv`, `gene_signature_per_cluster.csv`, `signature_heatmap.png`, `signature_umap.png`, `signature_correlation.png`
- `immune_phenotyping`: `immune_phenotyping.csv`, `immune_subtype_summary.csv`, `umap_immune_subtype.png`, `immune_signature_heatmap.png`
- `metacell`: `metacell_assignments.csv`, `metacell_summary.csv`, `metacells.h5ad`, `metacell_size_hist.png`, `metacell_umap.png`
- `paper_repro`: `paper_repro_registry.csv`, `paper_repro_figures.csv`, `paper_repro_report.json` (and `paper_repro_spec.template.json` when no spec is provided)
- `pathway_analysis`: `pathway_enrichment.csv`, `pathway_activity_per_cluster.csv`, `pathway_enrichment_bar.png`, `pathway_activity_heatmap.png`
- `pseudo_velocity`: `pseudo_velocity_speed.csv`, `pseudo_velocity_per_cluster.csv`, `pseudo_velocity_arrows*.png`, `pseudo_velocity_stream*.png`, `pseudo_velocity_speed_umap.png`, `pseudo_velocity_speed_boxplot.png`
- `pseudobulk_de`: `pseudobulk_counts.csv`, `pseudobulk_de_results.csv`, `pseudobulk_volcano.png`, `pseudobulk_heatmap.png`
- `qc`: `qc_violin_pre_filter.png`, `qc_violin_post_filter.png`, `qc_scatter_pre_filter.png`, `qc_scatter_post_filter.png`
- `rna_velocity`: `velocity_confidence.csv`, `velocity_top_genes.csv`, `velocity_stream_umap.png`, `velocity_grid_umap.png`, `velocity_length_distribution.png` (+ dynamical-mode plots)
- `trajectory`: `dpt_pseudotime.csv`, `pseudotime_per_cluster.csv`, `pseudotime_top_genes.csv`, `pseudotime_dpt_umap.png`, `paga_trajectory.png`
- `tumor_microenvironment`: `tme_scores_per_cell.csv`, `tme_scores_per_cluster.csv`, `checkpoint_expression.csv`, `tme_cyt_umap.png`, `tme_tis_umap.png`
- `validate_cbioportal`: `cbioportal_mutation_summary.csv`, `cbioportal_validation_report.json`

Each run is fully independent. Multiple runs accumulate under `results/`:

```
results/
├── LUSC_3k_Analysis_20260404_062118/    # Run 1
├── LUSC_immuno_20260405_091500/         # Run 2
└── LUSC_tme_20260405_140000/            # Run 3
```

## CLI Parameters

### General

| Parameter | Description |
|---|---|
| `--project` | Run name (used in output directory naming) |
| `--sample-root` | Dataset root directory |
| `--outs-dir` | Explicit path to Cell Ranger `filtered_feature_bc_matrix` (default: `<sample-root>/outs/filtered_feature_bc_matrix`) |
| `--output-dir` | Output root (default: `/home/zerlinshen/singlecell_factory/results`) |
| `--optional-modules` | Comma-separated module list (dependencies auto-included) |
| `--paper-spec-json` | JSON spec for paper-driven provenance + figure reproduction checks |
| `--paper-repro-strict` | Fail run when paper_repro has unresolved metadata/figure checks |
| `--markers-json` | Optional custom marker dictionary JSON for annotation |
| `--checkpoint` | Save checkpoints after each module for crash recovery |
| `--resume-from MODULE` | Resume from a specific module using saved checkpoints |
| `--parallel-workers N` | Number of parallel workers (default: 1 = sequential) |

For NC2024-scale full-cohort runs, prefer:

```bash
SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only python -m workflow.modular.cli ... --scale-mode massive --checkpoint
```

This preserves `.checkpoints/after_<module>.json` status/metadata sidecars while skipping full AnnData checkpoint writes in `massive` mode. It does not change numerical analysis results; it only avoids repeated 10GB+ checkpoint writes when disk headroom is more important than full-object resume after every optional module.

### Cell Ranger

| Parameter | Default | Description |
|---|---|---|
| `--fastq-dir` | empty | FASTQ directory for running `cellranger count` when needed |
| `--transcriptome-dir` | empty | Cell Ranger reference dir (also used for RNA-velocity GTF auto-discovery) |
| `--sample-id` | `lusc` | Sample ID passed to Cell Ranger |
| `--localcores` | 8 | CPU cores for Cell Ranger |
| `--localmem` | 64 | RAM (GB) for Cell Ranger |
| `--force-cellranger` | false | Force rerun of `cellranger count` |
| `--no-run-cellranger-if-missing` | false | Do not run Cell Ranger even if `outs-dir` is missing |

### QC

| Parameter | Default | Description |
|---|---|---|
| `--min-genes` | 200 | Minimum genes per cell |
| `--max-genes` | 7000 | Maximum genes per cell |
| `--min-counts` | 500 | Minimum UMI counts |
| `--max-counts` | 50000 | Maximum UMI counts |
| `--max-mito-pct` | 20 | Maximum mitochondrial % |
| `--max-ribo-pct` | 50 | Maximum ribosomal % |
| `--min-cells` | 3 | Minimum cells per gene after filtering |

### Doublet Detection

| Parameter | Default | Description |
|---|---|---|
| `--expected-doublet-rate` | 0.06 | Expected doublet fraction for Scrublet |
| `--no-remove-doublets` | false | Keep doublets (mark only, do not filter) |

### Clustering / DE / Annotation

| Parameter | Default | Description |
|---|---|---|
| `--n-top-genes` | 3000 | Number of HVGs |
| `--n-pcs` | 40 | PCA dimensions |
| `--n-neighbors` | 15 | k-NN neighbors |
| `--leiden-resolution` | 0.8 | Leiden clustering resolution |
| `--scale-data` | false | Apply `sc.pp.scale()` before PCA |
| `--de-method` | `wilcoxon` | DE method (`wilcoxon`, `t-test`, `t-test_overestim_var`, `logreg`) |
| `--de-n-genes` | 300 | Max genes ranked per cluster |
| `--de-pval-threshold` | 0.05 | Adjusted p-value cutoff for significant DE |
| `--de-logfc-threshold` | 0.25 | Minimum absolute log fold-change cutoff for DE |
| `--annotation-confidence-threshold` | 0.1 | Minimum annotation confidence score before assigning `Unknown` |
| `--reference-adata` | empty | Optional reference `h5ad` for annotation label transfer |
| `--reference-label-key` | `cell_type` | Label column in reference `obs` used for transfer |
| `--reference-k` | 15 | K neighbors for reference mapping |
| `--reference-min-confidence` | 0.6 | Min confidence needed to override marker label |
| `--reference-override-mode` | `conservative` | `conservative` (override Unknown/low-confidence only) or `all` |

### Batch Correction

| Parameter | Default | Description |
|---|---|---|
| `--batch-key` | sample | Batch column in adata.obs |
| `--batch-method` | harmony | Method: harmony/bbknn/combat/scanorama/scvi/mnn/fastmnn |
| `--scvi-max-epochs` | 200 | Max epochs for scVI training when `--batch-method scvi` |
| `--scvi-n-latent` | 30 | Latent dimension for scVI embedding |
| `--no-scvi-early-stopping` | false | Disable scVI early stopping (default behavior is enabled) |

### Trajectory / Cell Cycle / CNV / Signatures

| Parameter | Default | Description |
|---|---|---|
| `--regress-cell-cycle` | false | Regress out `S_score` and `G2M_score` after cell-cycle scoring |
| `--trajectory-root-cluster` | empty | Leiden cluster ID used as DPT root |
| `--cnv-reference-group` | empty | Reference group for CNV normalization |
| `--cnv-window-size` | 100 | Sliding-window size for CNV smoothing |
| `--signature-json` | empty | Custom signature JSON for `gene_signature_scoring` |

### cBioPortal Validation

| Parameter | Default | Description |
|---|---|---|
| `--cbioportal-genes` | empty | Comma-separated gene list for validation |
| `--cbioportal-study` | `lusc_tcga_pan_can_atlas_2018` | cBioPortal study ID |
| `--no-cbioportal-de-genes` | false | Disable auto-inclusion of top DE genes |
| `--cbioportal-top-n` | 20 | Number of top DE genes used for validation |

### RNA Velocity

| Parameter | Default | Description |
|---|---|---|
| `--velocity-loom` | - | Loom file with spliced/unspliced counts |
| `--velocity-bam` | - | Path to possorted_genome_bam.bam (Cell Ranger BAM output) |
| `--velocity-gtf` | - | Path to genes.gtf(.gz); optional if `--transcriptome-dir` is set (auto-discovery) |
| `--velocity-mode` | stochastic | scVelo mode: stochastic/dynamical |
| `--velocity-n-jobs` | 4 | Parallel workers for BAM extraction and scVelo dynamics |
| `--velocity-min-shared-counts` | 20 | Minimum shared counts for scVelo gene filtering |
| `--velocity-n-pcs` | 30 | PCA components for scVelo moments |
| `--velocity-n-neighbors` | 30 | Neighbors for scVelo moments |

If no loom file is provided, the module can extract spliced/unspliced counts directly from Cell Ranger BAM output using pysam (parallelized by chromosome). `genes.gtf(.gz)` is taken from `--velocity-gtf` or auto-discovered from `--transcriptome-dir`.

Important input rules:
- The pipeline does **not** auto-generate a loom file.
- The pipeline does **not** auto-generate a GTF file; it only resolves an existing one from `--velocity-gtf` or `--transcriptome-dir`.
- BAM extraction requires `--velocity-bam` (typically `<sample-root>/outs/possorted_genome_bam.bam`) plus a resolvable GTF.
- `--velocity-bam` is currently **not auto-inferred** from `--sample-root`; pass it explicitly when enabling `rna_velocity`.

The module runs scVelo on an internal adata copy and transfers only cell-level results back, preserving the shared gene index and enabling parallel execution with other modules.
RNA velocity BAM classification uses strict exon/intron rules only (no lossy fast-path).
Policy update (April 5, 2026): all lossy velocity shortcuts were removed from pipeline code.

Velocity BAM extraction cache:
- Default cache dir: `/tmp/singlecell_factory_velocity_cache`
- Override: `SCF_VELOCITY_CACHE_DIR=/path/to/cache`
- Cache key includes BAM/GTF identity + barcode/gene ordering + velocity extraction mode; repeated runs with the same inputs can skip BAM parsing.

On hosts with >=16 logical CPU cores, if `--velocity-n-jobs` is left at default (`4`), BAM extraction workers are auto-bumped to `8` and recorded in `run_manifest.json -> metadata.velocity_extract_n_jobs`.

**Dynamical mode** additionally outputs latent time UMAP and phase portraits for top velocity genes.

## Multi-Backend Support

| Module | Primary | Fallback |
|---|---|---|
| `pathway_analysis` | gseapy (MSigDB) | decoupler (PROGENy) -> built-in Hallmark (BH FDR) |
| `cell_communication` | LIANA (multi-method consensus) | Manual L-R scoring (18 TME pairs) |
| `gene_regulatory_network` | decoupler + DoRothEA | Manual TF-target scoring (12 key TFs) |
| `cnv_inference` | infercnvpy | Vectorized sliding-window smoothing |
| `batch_correction` | Harmony | BBKNN / Combat / Scanorama |

## Performance Optimizations (v5.0)

### GPU Acceleration
- **rapids-singlecell**: Auto-detected GPU backend for:
  - `clustering`: PCA/neighbors/UMAP/Leiden
  - `batch_correction`: post-correction neighbors/UMAP/Leiden
  - `differential_expression`: `rank_genes_groups` for supported methods
  - `evolution`: clone marker ranking
- Safe fallback behavior:
  - If GPU backend is unavailable, modules transparently use CPU scanpy.
  - If GPU clustering fails mid-run, the module reruns on a pristine CPU input object (no partial GPU-state reuse).
  - If `--batch-method scvi` is selected but `scvi-tools` is unavailable (or scVI training/input checks fail), `batch_correction` is marked `skipped` with explicit reason instead of crashing downstream optional workflow.
- Backend metadata tracking: `clustering_backend`, `de_backend`, `batch_post_backend` recorded in `run_manifest.json` for reproducibility auditing.

### Memory & I/O
- **Memory guard**: Parallel worker count is constrained by estimated AnnData copy size + available RAM to avoid OOM in branch execution.
- **Raw matrix policy**: Normalized/log1p matrix is stored in `adata.raw` before downstream analyses to improve biological interpretability for marker/score modules.
- **Zarr checkpoints**: 3-5x faster checkpoint I/O via `adata.write_zarr()` with automatic h5ad fallback for incompatible key names.
- **Massive metadata-only checkpoints**: `SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only` keeps checkpoint JSON sidecars but skips full AnnData checkpoint files in `--scale-mode massive`.
- **Async figure I/O**: Full render+write offloaded to background thread pool (PNG encoding no longer blocks the main thread).

### Computation
- **Clustering**: PCA on HVGs only, optional scaling via `--scale-data`, CPU/GPU branches now both produce UMAP outputs.
- **QC**: Optimized `sc.pp.calculate_qc_metrics` (`log1p=False`, `percent_top=None`).
- **Differential Expression**: Default switched to `wilcoxon` with configurable `--de-method` and `--de-n-genes`.
- **Trajectory**: Pseudotime-gene correlation refactored to sparse-friendly computation, avoiding full-matrix densification.
- **Cell communication**: Vectorized `np.where` scoring — eliminates O(n^2) Python loops.
- **CNV inference**: `scipy.ndimage.uniform_filter1d` vectorized sliding window.
- **Pseudo-velocity baseline**: Fully vectorized numpy broadcasting replaces per-cell Python loop.
- **RNA velocity**: BAM extraction parallelized by chromosome (~4-5x speedup).

### Pipeline Orchestration
- **Cost-aware tier scheduling**: Heavy modules (rna_velocity, cnv_inference) start first in parallel tiers.
- **Thread-safe status tracking**: `threading.Lock` on module status writes.
- **Copy-on-write branching**: Parallel modules get isolated AnnData copies with proper merge-back of obs/obsm/uns/metadata/module_dirs. Structural mutations (X, layers, varm, obsp) are detected and logged as warnings if an appending module makes changes that merge-back cannot capture.
- **Mutating module discovery**: Modules that modify adata structure declare `mutates_structure = True` as a class attribute. The pipeline automatically discovers these at registry build time and forces them into sequential execution. A static fallback set provides safety for modules that omit the declaration.
- **Module contract validation**: Modules declare `requires_keys` (e.g., `{"obs": ["leiden"]}`) and `provides_keys`. Before each `mod.run(ctx)`, the pipeline checks that required keys exist in adata. Optional modules with unmet requirements are skipped with a warning; mandatory modules raise immediately.
- **O(1) dependency resolution**: Topological sort uses `collections.deque` for O(1) queue operations; sequential execution logic is DRY-extracted into a single `_run_sequential()` helper shared by mutating, appending, and memory-fallback paths.
- **Runtime telemetry**: `metadata.module_runtime_sec` + `metadata.pipeline_wall_seconds` are persisted in each `run_manifest.json`.
- **Checkpoint/resume**: `--checkpoint` + `--resume-from` for crash recovery.

### Performance Benchmark Results (LUSC 3K Dataset)

Benchmark script output (`output/performance_benchmark_20260405_optfix/benchmark_compare.json`):

| Metric | Baseline | Optimized | Change |
|---|---|---|---|
| Wall-clock time | 10.04 s | 9.16 s | **-8.79%** |
| Peak RSS memory | 799 MB | 564 MB | **-29.45%** |
| Clusters found | 14 | 15 | ARI = 0.54 |

Latest end-to-end real-data run on **April 6, 2026** (`results/LUSC_GPU_FINAL_V2_20260406_002301`):
- Pipeline wall time: **81.907 s**
- Cells: `2588 -> 2382` after QC -> `2377` after doublet removal
- Clusters: `18`
- Module status: **24/24 ok** (3 mandatory + 21 optional)
- Key metadata: `clustering_backend=cpu`, `de_backend=cpu` (host has no detected CUDA backend)
- Heaviest modules by wall-time: `rna_velocity (31.377s)`, `validate_cbioportal (16.026s)`, `evolution (12.573s)`

GPU-change verification runs on **April 6, 2026**:
- `results/CODEX_GPU_SMOKE_POSTFIX_20260406_010537`: post-fix 11-module smoke (includes clustering/batch_correction/differential_expression/evolution), **11/11 ok**
- `results/CODEX_GPU_FULL_LOCAL_NOCKPT_20260406_010310`: local full run without external validation module, **22/23 ok**
- In that local full run, `rna_velocity` failed due missing spliced/unspliced inputs (no loom/BAM+GTF), which is expected and unrelated to GPU acceleration paths.

README profile verification run on **April 6, 2026**:
- `results/README_FULL_LOCAL_20260406_20260406_012249`: executed exactly with the documented `Full local analysis` command, **23/23 ok**, `pipeline_wall_seconds=40.962`.
- Detailed execution/audit report: `reports/README_FULL_LOCAL_20260406_运行与审查报告.md`.

RNA velocity bottleneck benchmark on real LUSC data (strict mode, same inputs, **April 5, 2026**):
- Cold run: `results/LUSC_VEL_STRICT_ONLY_20260405_173025`
- Warm run (same cache dir): `results/LUSC_VEL_STRICT_ONLY_WARM_20260405_173113`

| Metric | Cold | Warm | Change |
|---|---:|---:|---:|
| Pipeline wall time | 37.143 s | 15.116 s | **-59.3%** |
| `rna_velocity` module time | 26.422 s | 4.418 s | **-83.3%** |
| Velocity layer loading (`metadata.velocity_extract_seconds`) | 22.006 s | 0.090 s | **-99.6%** |
| Mean velocity confidence | 0.86117 | 0.86117 | identical |

Output reproducibility checks (cold vs warm): identical file hashes for `rna_velocity/velocity_confidence.csv`, `rna_velocity/velocity_top_genes.csv`, and `clustering/umap_leiden.png`.

RNA velocity extraction parallel scaling on the same dataset (strict mode, cold cache):
- `n_jobs=4`: `38.323 s`
- `n_jobs=8`: `26.148 s` (**31.8% faster**)

### Best Practice Configuration Recommendations

1. **Large datasets (>50k cells)**: start with `--parallel-workers 2-4`, then scale up only if RAM headroom is sufficient.
2. **DE for publication**: keep `--de-method wilcoxon`; use `t-test_overestim_var` only when speed is the top priority.
3. **Single-sample runs**: `pseudobulk_de` is expected to skip (requires >=2 biological samples).
4. **Batch correction**: only enable when a real batch column exists; otherwise keep skipped to avoid unnecessary recomputation.
5. **Runtime tuning loop**: inspect `run_manifest.json -> metadata.module_runtime_sec` and optimize the top 3 slowest modules first.
6. **RNA velocity repeated runs**: set `SCF_VELOCITY_CACHE_DIR` to a fast local SSD path and reuse it across reruns.
7. **RNA velocity policy**: keep strict exon/intron classification for all runs; optimize with cache and parallelism instead of lossy shortcuts.
8. **High-core hosts**: default velocity extraction auto-uses 8 workers on >=16 cores; set `--velocity-n-jobs` explicitly to override.

## Quality Verification Checklist

1. Check `module_status.csv` — mandatory modules should be `ok`; optional modules may be `ok` or `skipped` with a clear reason
2. Review `qc/qc_violin_*.png` — verify filter thresholds are reasonable
3. Review `clustering/umap_leiden.png` — clusters should be well-separated
4. Check `annotation/umap_cell_type.png` — cell types should be biologically coherent
5. Review `immune_phenotyping/immune_signature_heatmap.png` — CD8_exhausted should correlate with exhaustion score
6. Check `tumor_microenvironment/checkpoint_dotplot.png` — checkpoint expression patterns

## Methodology & References

Every analysis module uses publicly recognized, peer-reviewed methods. Below is the complete methodology audit with citations.

Module coverage check: **25 / 25 core modules documented and citation-aligned**. Phase C adds 3 bridge-ready modality modules — `protein_adt`, `spatial_ingest` + `spatial_neighborhoods`, `multimodal_integration` (EXPERIMENTAL) — registered in `workflow/modular/pipeline.py` and the catalog (DAG-count assertion = 29).

---

### 01. Cell Ranger Data Loading (`cellranger`)

| Item | Detail |
|---|---|
| **Method** | 10X Genomics Cell Ranger count matrix loading |
| **Implementation** | `scanpy.read_10x_mtx()` |
| **Input format** | Market Matrix (MTX) sparse format from Cell Ranger |
| **Reference** | Zheng et al., *Nature Communications*, 2017. DOI: [10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049) |

---

### 02. Quality Control (`qc`)

| Item | Detail |
|---|---|
| **Method** | Threshold-based cell and gene filtering |
| **Metrics** | n_genes_by_counts, total_counts, pct_counts_mt, pct_counts_ribo, pct_counts_hb |
| **Implementation** | `scanpy.pp.calculate_qc_metrics()`, `scanpy.pp.filter_genes()` |
| **Gene classes** | Mitochondrial (MT-), Ribosomal (RPS/RPL), Hemoglobin (HBA/HBB) |
| **Defaults** | 200-7000 genes, 500-50000 UMI, <=20% mito, <=50% ribo |
| **Reference** | Luecken & Theis, *Molecular Systems Biology*, 2019. DOI: [10.15252/msb.20188746](https://doi.org/10.15252/msb.20188746) |

---

### 03. Doublet Detection (`doublet_detection`)

| Item | Detail |
|---|---|
| **Method** | Scrublet — synthetic doublet simulation |
| **Algorithm** | Simulates doublets by averaging random cell pairs, scores real cells by k-NN similarity to synthetic doublets |
| **Implementation** | `scrublet.Scrublet.scrub_doublets(min_counts=2, min_cells=3, n_prin_comps=30)` |
| **Fallback** | For tiny/degenerate datasets, pipeline falls back to all-singlets (records `doublet_method=fallback_all_singlets`) instead of aborting |
| **Parameters** | expected_doublet_rate=0.06, min_gene_variability_pctl=85 |
| **Reference** | **Wolock et al., *Cell Systems*, 2019.** DOI: [10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005) |

---

### 04. Clustering (`clustering`)

| Item | Detail |
|---|---|
| **Normalization** | CPM (counts-per-million) + log1p, via `scanpy.pp.normalize_total()` |
| **HVG selection** | Seurat flavor HVG selection (`scanpy.pp.highly_variable_genes`) |
| **PCA** | Truncated SVD via ARPACK solver on HVGs only (memory-optimized) |
| **Neighbor graph** | k-nearest neighbors (k=15 default) in PCA space |
| **Dimensionality reduction** | UMAP (Uniform Manifold Approximation and Projection) |
| **Clustering** | Leiden community detection algorithm (igraph, undirected, resolution=0.8) |
| **References** | **McInnes et al., *JOSS*, 2018.** DOI: [10.21105/joss.00861](https://doi.org/10.21105/joss.00861) (UMAP); **Traag et al., *Scientific Reports*, 2019.** DOI: [10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z) (Leiden); **Stuart et al., *Cell*, 2019.** DOI: [10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031) (Seurat v3 HVG) |

---

### 05. Cell Type Annotation (`annotation`)

| Item | Detail |
|---|---|
| **Method** | Marker gene set scoring with strategy-aware label assignment |
| **Strategy (default)** | `cluster_voting` — aggregate raw mean expression per leiden cluster, score the cluster-aggregated expression with each cell-type marker panel, broadcast the argmax cell type label to every cell in that cluster. Robust at >100k cohorts. |
| **Strategy (fallback)** | `cell_argmax` — per-cell `sc.tl.score_genes` then per-cell argmax then per-cluster majority vote. Drifts on >100k cohorts (used in the v1 NC2024 run, replaced by `cluster_voting` in v2); retained for backward-compat only. Selected via `--annotation-strategy cell_argmax`. |
| **Implementation** | `workflow/modular/modules/annotation.py`; `score_genes` still runs per cell to produce confidence scores and to gate the optional reference-mapping override. |
| **Confidence metric** | max_score - second_max_score (vectorized via `np.partition`) |
| **Marker panels** | 10 cell types: Tumor epithelial, T cell, NK cell, B cell, Myeloid/Macro, Fibroblast, Endothelial, Plasma cell, Mast cell, Dendritic cell |
| **Unknown threshold** | Cells with confidence < 0.1 labeled "Unknown" |
| **References** | **Tirosh et al., *Science*, 2016.** DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) (gene set scoring); **Sanchez-Mejias et al., *Nature Communications*, 2024.** DOI: [10.1038/s41467-024-48700-8](https://doi.org/10.1038/s41467-024-48700-8) (NC2024 NSCLC atlas — paper-aligned annotation reference for the v2 cohort) |

---

### 06. Cell Cycle Scoring (`cell_cycle`)

| Item | Detail |
|---|---|
| **Method** | S-phase and G2M-phase gene set scoring |
| **Gene sets** | 47 S-phase genes + 49 G2M-phase genes (Tirosh/Regev lab) |
| **Implementation** | `scanpy.tl.score_genes_cell_cycle()` |
| **Optional** | `scanpy.pp.regress_out(["S_score", "G2M_score"])` to remove cell cycle effects |
| **Reference** | **Tirosh et al., *Science*, 2016.** DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |

---

### 07. CNV Inference (`cnv_inference`)

| Item | Detail |
|---|---|
| **Method (primary)** | Expression-based sliding-window CNV inference |
| **Algorithm** | Center gene expression by reference mean, smooth per chromosome using `scipy.ndimage.uniform_filter1d`, compute per-cell variance as CNV score |
| **Gene positions** | Auto-fetched from Ensembl BioMart via `pybiomart` if not pre-annotated |
| **Classification** | Percentile thresholding (default 75th) for malignant vs normal |
| **Method (fallback)** | `infercnvpy.tl.infercnv()` + `infercnvpy.tl.cnv_score()` |
| **References** | **Patel et al., *Science*, 2014.** DOI: [10.1126/science.1254257](https://doi.org/10.1126/science.1254257); **Tirosh et al., *Science*, 2016.** DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |

---

### 08. Differential Expression (`differential_expression`)

| Item | Detail |
|---|---|
| **Statistical test** | Wilcoxon rank-sum by default (configurable via `--de-method`) |
| **Implementation** | `scanpy.tl.rank_genes_groups(method=<de_method>, use_raw=False, n_genes=<de_n_genes>, pts=True)` |
| **Multiple testing** | Benjamini-Hochberg FDR correction (scanpy default) |
| **Significance filter** | adjusted p-value < 0.05 AND |log2FC| > 0.25 (configurable via `--de-pval-threshold`, `--de-logfc-threshold`) |
| **Output** | Ranked genes per cluster (`--de-n-genes`, default 300) with scores/p-values/fold changes |
| **References** | **Wilcoxon, *Biometrics Bulletin*, 1945.** DOI: [10.2307/3001968](https://doi.org/10.2307/3001968); **Benjamini & Hochberg, *JRSS-B*, 1995.** DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) |

---

### 09. Gene Regulatory Network (`gene_regulatory_network`)

| Item | Detail |
|---|---|
| **Method (primary)** | decoupler + DoRothEA regulon database |
| **Algorithm** | Univariate Linear Model (ULM) — infers TF activity from target gene expression |
| **Regulon confidence** | Levels A (high), B (medium), C (curated) from DoRothEA |
| **Implementation** | `decoupler.get_dorothea(organism="human", levels=["A","B","C"])`, `decoupler.run_ulm()` |
| **Fallback** | Gene set scoring for 12 curated TF-target sets (TP53, MYC, STAT1, NFKB1, HIF1A, SOX2, FOXP3, E2F1, NOTCH1, SMAD3, JUN, AP-1) |
| **References** | **Garcia-Alonso et al., *Genome Research*, 2019.** DOI: [10.1101/gr.240663.118](https://doi.org/10.1101/gr.240663.118) (DoRothEA); **Badia-i-Mompel et al., *Bioinformatics Advances*, 2022.** DOI: [10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016) (decoupler) |

---

### 10. Gene Signature Scoring (`gene_signature_scoring`)

| Item | Detail |
|---|---|
| **Method** | Scanpy gene set scoring (`sc.tl.score_genes`) |
| **Built-in signatures** | 10 cancer hallmark panels (see table below) |
| **Custom input** | User-provided JSON: `{"name": ["GENE1", "GENE2", ...]}` |
| **Minimum genes** | Requires >= 2 genes present per signature |

**Built-in Cancer Signatures:**

| Signature | Genes | Reference |
|---|---|---|
| Proliferation | MKI67, TOP2A, PCNA, MCM2, MCM6, CDK1, CCNB1, CCNB2 | Standard oncology panel |
| Apoptosis resistance | BCL2, BCL2L1, MCL1, BIRC5, XIAP, CFLAR | BCL2 family |
| Angiogenesis | VEGFA, VEGFB, FLT1, KDR, PECAM1, ANGPT2, NRP1 | VEGF/FLT pathway |
| Invasion/metastasis | MMP2, MMP9, MMP14, SNAI1, TWIST1, VIM, CDH2 | MMP/EMT markers |
| EMT mesenchymal | VIM, CDH2, FN1, SNAI2, ZEB1, ZEB2, TWIST1, MMP2 | Tan et al., *EMBO Mol Med*, 2014. DOI: [10.15252/emmm.201404208](https://doi.org/10.15252/emmm.201404208) |
| EMT epithelial | CDH1, EPCAM, KRT8, KRT18, KRT19, CLDN4, OCLN | Tan et al., *EMBO Mol Med*, 2014. DOI: [10.15252/emmm.201404208](https://doi.org/10.15252/emmm.201404208) |
| Stemness | POU5F1, NANOG, SOX2, KLF4, MYC, LIN28A, SALL4, BMI1 | Malta et al., *Cell*, 2018. DOI: [10.1016/j.cell.2018.03.034](https://doi.org/10.1016/j.cell.2018.03.034) |
| Hypoxia | VEGFA, SLC2A1, HK2, LDHA, PGK1, CA9, BNIP3, ENO1 | Buffa et al., *Br J Cancer*, 2010. DOI: [10.1038/sj.bjc.6605450](https://doi.org/10.1038/sj.bjc.6605450) |
| DNA damage response | BRCA1, BRCA2, ATM, ATR, RAD51, CHEK1, CHEK2, TP53 | DDR pathway |
| Glycolysis | HK2, PFKP, PKM, LDHA, ENO1, GAPDH, TPI1, ALDOA | Warburg effect |

---

### 11. Trajectory / Pseudotime (`trajectory`)

| Item | Detail |
|---|---|
| **Methods** | PAGA trajectory graph + Diffusion Pseudotime (DPT) + gene expression dynamics |
| **PAGA** | Partition-based graph abstraction — discovers trajectory topology between clusters |
| **DPT** | Diffusion map embedding + geodesic distance from root cell for temporal ordering |
| **Gene trends** | Top 30 genes correlated with pseudotime, smoothed heatmap visualization |
| **Implementation** | `scanpy.tl.paga()`, `scanpy.tl.diffmap()`, `scanpy.tl.dpt()` |
| **References** | **Wolf et al., *Genome Biology*, 2019.** DOI: [10.1186/s13059-019-1663-x](https://doi.org/10.1186/s13059-019-1663-x) (PAGA); **Haghverdi et al., *Nature Methods*, 2016.** DOI: [10.1038/nmeth.3971](https://doi.org/10.1038/nmeth.3971) (DPT) |

---

### 12. Cell Communication (`cell_communication`)

| Item | Detail |
|---|---|
| **Method (primary)** | LIANA multi-method consensus scoring |
| **Algorithm** | Aggregates 4+ L-R scoring methods (CellPhoneDB, NATMI, SingleCellSignalR, Connectome) into consensus rank |
| **Implementation** | `liana.mt.rank_aggregate(groupby="cell_type", resource_name="consensus")` |
| **Fallback** | Manual L-R scoring: pre-computed mean expression matrix, outer-product scoring for 18 curated TME L-R pairs (PD-L1/PD-1, VEGFA/KDR, TGFB1/TGFBR2, etc.) |
| **Reference** | **Dimitrov et al., *Nature Communications*, 2022.** DOI: [10.1038/s41467-022-30755-0](https://doi.org/10.1038/s41467-022-30755-0) |

---

### 13. Immune Phenotyping (`immune_phenotyping`)

| Item | Detail |
|---|---|
| **Method** | Gene set scoring for 15 immune subtypes + 3 functional signatures |
| **Implementation** | `scanpy.tl.score_genes()` per subtype panel |
| **Assignment** | Argmax scoring restricted to cells annotated as immune parent types |

**15 Immune Subtypes:**

| Subtype | Key Markers | Source |
|---|---|---|
| CD4 naive | CCR7, LEF1, SELL, TCF7, IL7R | Zheng et al., *Cell*, 2017 |
| CD4 memory | IL7R, S100A4, ANXA1, CCR6, AQP3 | Zheng et al., *Cell*, 2017 |
| Treg | FOXP3, IL2RA, CTLA4, IKZF2, TNFRSF18 | Zheng et al., *Cell*, 2017 |
| Th1 / Th2 / Th17 | TBX21 / GATA3 / RORC + cytokines | Zhang et al., *Nature*, 2018 |
| CD8 effector | CD8A, GZMB, PRF1, NKG7, GNLY, IFNG | Zheng et al., *Cell*, 2017 |
| CD8 memory | CD8A, IL7R, GPR183, SELL, TCF7 | Zhang et al., *Nature*, 2018 |
| CD8 exhausted | CD8A, PDCD1, LAG3, HAVCR2, TIGIT, TOX, CTLA4 | Zheng et al., *Cell*, 2017 |
| NK cytotoxic | NKG7, GNLY, KLRD1, KLRF1, NCAM1, FCGR3A | Tirosh et al., *Science*, 2016 |
| Macro M1 | CD68, NOS2, IL1B, TNF, CXCL10, CD80 | Zhang et al., *Nature*, 2018 |
| Macro M2 | CD68, CD163, MRC1, MSR1, TGFB1, IL10 | Zhang et al., *Nature*, 2018 |
| cDC1 / cDC2 / pDC | CLEC9A / CD1C / LILRA4 + subtype markers | Zhang et al., *Nature*, 2018 |

**Functional Scores:**
- **Exhaustion**: LAG3, HAVCR2, TIGIT, PDCD1, CTLA4, TOX, TOX2, ENTPD1
- **Cytotoxicity**: GZMA, GZMB, GZMK, GZMH, PRF1, NKG7, GNLY, FASLG
- **Activation**: CD69, CD38, HLA-DRA, ICOS, TNFRSF9, TNFRSF4

| Reference | DOI |
|---|---|
| **Zheng et al., *Cell*, 2017** — "Landscape of Infiltrating T Cells in Liver Cancer" | [10.1016/j.cell.2017.05.035](https://doi.org/10.1016/j.cell.2017.05.035) |
| **Zhang et al., *Nature*, 2018** — "Lineage tracking reveals dynamic relationships of T cells in colorectal cancer" | [10.1038/s41586-018-0694-x](https://doi.org/10.1038/s41586-018-0694-x) |
| **Tirosh et al., *Science*, 2016** — "Dissecting the multicellular ecosystem of metastatic melanoma" | [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |

---

### 14. Tumor Microenvironment (`tumor_microenvironment`)

| Signature | Genes | Method | Reference |
|---|---|---|---|
| **CYT** (Cytolytic Activity) | GZMA, PRF1 | Geometric mean: sqrt(GZMA x PRF1) | Rooney et al., *Cell*, 2015. DOI: [10.1016/j.cell.2014.12.033](https://doi.org/10.1016/j.cell.2014.12.033) |
| **TIS** (T-cell Inflamed, 18-gene) | CD27, CD274, CD276, CD8A, CMKLR1, CXCL9, CXCR6, HLA-DQA1, HLA-DRB1, HLA-E, IDO1, LAG3, NKG7, PDCD1LG2, PSMB10, STAT1, TIGIT, TGFB1 | Gene set score | Ayers et al., *JCI*, 2017. DOI: [10.1172/JCI91190](https://doi.org/10.1172/JCI91190) |
| **IFN-gamma** (10-gene) | IFNG, STAT1, CCR5, CXCL9, CXCL10, CXCL11, IDO1, PRF1, GZMA, HLA-DRA | Gene set score | Ayers et al., *JCI*, 2017. DOI: [10.1172/JCI91190](https://doi.org/10.1172/JCI91190) |
| **Immunosuppression** | IDO1, CD274, PDCD1LG2, HAVCR2, LAG3, TGFB1, IL10, VEGFA, ARG1 | Gene set score | Composite |
| **Immune ESTIMATE** | CD2, CD3D, CD3E, CD3G, CD8A, CD19, CD79A, MS4A1, LCK, GZMB, NKG7, PRF1, GNLY | Gene set score | Yoshihara et al., *Nat Commun*, 2013. DOI: [10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612) |
| **Stromal ESTIMATE** | DCN, LUM, COL1A1, COL3A1, COL5A1, FAP, ACTA2, VIM, FN1, SPARC, POSTN, THY1 | Gene set score | Yoshihara et al., *Nat Commun*, 2013. DOI: [10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612) |

**Checkpoint Molecules Profiled:** PD-1 (PDCD1), PD-L1 (CD274), PD-L2 (PDCD1LG2), CTLA-4, LAG-3, TIM-3 (HAVCR2), TIGIT, VISTA (VSIR), IDO1, B7-H3 (CD276)

---

### 15. Pathway Analysis (`pathway_analysis`)

| Backend | Method | Implementation | Reference |
|---|---|---|---|
| **gseapy** (primary) | Over-Representation Analysis with MSigDB Hallmark gene sets | `gseapy.enrich(gene_sets="MSigDB_Hallmark_2020")` | Subramanian et al., *PNAS*, 2005. DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102); Liberzon et al., *Cell Systems*, 2015. DOI: [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004) |
| **decoupler + PROGENy** | Multivariate Linear Model (MLM) pathway activity scoring | `decoupler.run_mlm()` with PROGENy top 300 footprint genes | Schubert et al., *Nature Communications*, 2018. DOI: [10.1038/s41467-017-02391-6](https://doi.org/10.1038/s41467-017-02391-6) |
| **Fallback** | Hypergeometric overlap test with built-in Hallmark gene sets | `scipy.stats.hypergeom.sf()` + Benjamini-Hochberg FDR | Standard hypergeometric test |

---

### 16. Pseudo-Velocity (`pseudo_velocity`)

| Item | Detail |
|---|---|
| **Method** | Local pseudotime gradient estimation in UMAP space with full visualization |
| **Algorithm** | For each cell, compute k-NN (k=15), calculate velocity = mean(spatial_direction x pseudotime_gradient) across neighbors. Produces quiver arrows, IDW-interpolated stream plots, and per-cluster speed statistics |
| **Implementation** | Vectorized NumPy broadcasting + `sklearn.neighbors.NearestNeighbors` + `scipy.spatial.cKDTree` for grid interpolation |
| **Visualizations** | Arrow field, stream plot, speed UMAP overlay, per-cluster speed boxplot |
| **Reference** | Inspired by La Manno et al., *Nature*, 2018. DOI: [10.1038/s41586-018-0414-6](https://doi.org/10.1038/s41586-018-0414-6) (RNA velocity concept adapted to pseudotime gradients) |

---

### 17. Batch Correction (`batch_correction`)

| Method | Algorithm | Reference |
|---|---|---|
| **Harmony** (default) | Iterative PCA-based soft clustering with batch diversity maximization | Korsunsky et al., *Nature Methods*, 2019. DOI: [10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0) |
| **BBKNN** | Batch-balanced k-nearest neighbors graph construction | Polanski et al., *Bioinformatics*, 2020. DOI: [10.1093/bioinformatics/btz625](https://doi.org/10.1093/bioinformatics/btz625) |
| **ComBat** | Empirical Bayes linear model batch adjustment | Johnson et al., *Biostatistics*, 2007. DOI: [10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037) |
| **Scanorama** | Mutual nearest neighbor panoramic stitching | Hie et al., *Nature Biotechnology*, 2019. DOI: [10.1038/s41587-019-0113-3](https://doi.org/10.1038/s41587-019-0113-3) |

---

### 18. RNA Velocity (`rna_velocity`)

| Item | Detail |
|---|---|
| **Method** | scVelo — RNA velocity from splicing kinetics |
| **Stochastic mode** | Moment-based estimation of splicing/degradation rates |
| **Dynamical mode** | Full transcriptome kinetics model with latent time inference |
| **Implementation** | `scvelo.tl.velocity()`, `scvelo.tl.velocity_graph()`, `scvelo.tl.velocity_confidence()` |
| **Architecture** | Copy-on-write execution (`adata.copy()`), then transfer only cell-level outputs back to main adata (non-mutating, parallel-safe) |
| **Requirements** | Spliced/unspliced count layers (from loom file, velocyto, or BAM extraction) |
| **Eligibility policy** | If no spliced/unspliced layers and no usable loom or BAM+GTF source are available, the module should be recorded as `skipped`, not treated as an unexpected failure. |
| **BAM extraction** | If no loom file is provided, spliced/unspliced counts are extracted directly from Cell Ranger BAM output (`possorted_genome_bam.bam`) using pysam + GTF-based exon/intron classification |
| **Outputs** | `velocity_stream_umap.png`, `velocity_grid_umap.png`, `velocity_length_distribution.png`, `velocity_confidence.csv`, `velocity_top_genes.csv` (+ `velocity_latent_time_umap.png`, `velocity_phase_portraits.png` in dynamical mode) |
| **Reference** | **Bergen et al., *Nature Biotechnology*, 2020.** DOI: [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3) |

---

### 19. cBioPortal Validation (`validate_cbioportal`)

| Item | Detail |
|---|---|
| **Method** | Cross-reference DE genes with TCGA large-cohort mutation data |
| **API** | cBioPortal REST API v2 (public, no authentication required) |
| **Endpoints** | Gene lookup, sample counts, mutation frequency profiling |
| **Default study** | LUSC TCGA Pan-Cancer Atlas 2018 |
| **Reference** | **Cerami et al., *Cancer Discovery*, 2012.** DOI: [10.1158/2159-8290.CD-12-0095](https://doi.org/10.1158/2159-8290.CD-12-0095); **Gao et al., *Science Signaling*, 2013.** DOI: [10.1126/scisignal.2004088](https://doi.org/10.1126/scisignal.2004088) |

---

### 20. Tumor Evolution (`evolution`)

| Item | Detail |
|---|---|
| **Method** | CNV-based clonal clustering + pseudotime-ordered evolution + phylogenetic reconstruction |
| **Clone identification** | Hierarchical clustering (Ward's method) on per-cell CNV profiles using correlation distance |
| **Phylogenetic tree** | Clone centroid dendrogram from `scipy.cluster.hierarchy` |
| **Evolution timeline** | Pseudotime distribution per clone reveals temporal ordering of clonal expansion |
| **Clone markers** | Wilcoxon rank-sum test for clone-specific DE genes with BH FDR correction |
| **Implementation** | `scipy.cluster.hierarchy.linkage()`, `scipy.spatial.distance.pdist()`, `scanpy.tl.rank_genes_groups()` |
| **References** | **Patel et al., *Science*, 2014.** DOI: [10.1126/science.1254257](https://doi.org/10.1126/science.1254257); **Gao et al., *Nature Biotechnology*, 2021.** DOI: [10.1038/s41587-020-00795-2](https://doi.org/10.1038/s41587-020-00795-2) (CopyKAT); **Tirosh et al., *Science*, 2016.** DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |

---

### 21. Pseudobulk Differential Expression (`pseudobulk_de`)

| Item | Detail |
|---|---|
| **Method** | Aggregate raw counts by biological sample plus grouping columns, then perform bulk-style DE |
| **Primary backend** | `pydeseq2` (`DeseqDataSet`, `DeseqStats`) |
| **Fallbacks** | Mann-Whitney U, then Wilcoxon rank-sum |
| **Multiple testing** | Benjamini-Hochberg FDR (`_bh_adjust`) |
| **Confirmatory mode** | Requires an explicit contrast contract: `--pseudobulk-contrast-col`, `--pseudobulk-contrast-a`, `--pseudobulk-contrast-b` (or `--pseudobulk-contrast-json`) |
| **Exploratory mode** | Optional `group_vs_rest` output only when `--pseudobulk-exploratory-group-vs-rest` is explicitly enabled |
| **Eligibility policy** | Missing counts, missing contrast contract, or too few biological replicates should be recorded as `skipped`, not as silent paper-grade output |
| **Implementation** | `workflow/modular/modules/pseudobulk_de.py` |
| **References** | **Love et al., *Genome Biology*, 2014.** DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8); **Wilcoxon, 1945** DOI: [10.2307/3001968](https://doi.org/10.2307/3001968) |

---

### 22. Cell Fate Mapping (`cell_fate`)

| Item | Detail |
|---|---|
| **Method (primary)** | CellRank Markov-state fate probability inference |
| **Fallback** | Diffusion-style transition probabilities on neighbor connectivities + pseudotime terminal-state heuristics |
| **Implementation** | `cellrank.kernels.PseudotimeKernel`, `cellrank.estimators.GPCCA`; fallback in `workflow/modular/modules/cell_fate.py` |
| **References** | **Lange et al., *Nature Methods*, 2022.** DOI: [10.1038/s41592-021-01346-6](https://doi.org/10.1038/s41592-021-01346-6); **Haghverdi et al., 2016** DOI: [10.1038/nmeth.3971](https://doi.org/10.1038/nmeth.3971) |

---

### 23. Composition Analysis (`composition`)

| Item | Detail |
|---|---|
| **Method (primary)** | Bayesian compositional DA with scCODA (via `pertpy`) |
| **Fallback** | Mann-Whitney U (`n_groups=2`) or Kruskal-Wallis (`n_groups>2`) + BH FDR |
| **Implementation** | `workflow/modular/modules/composition.py` |
| **References** | **Büttner et al., *Nature Communications*, 2021.** DOI: [10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6); **Benjamini & Hochberg, 1995** DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) |

---

### 24. Metacell Aggregation (`metacell`)

| Item | Detail |
|---|---|
| **Method (primary)** | SEACells archetype-based metacell construction |
| **Fallback** | MiniBatchKMeans clustering in PCA space |
| **Output contract** | Informational module; writes standalone `metacells.h5ad`, does not replace `adata.X` |
| **Implementation** | `workflow/modular/modules/metacell.py` |
| **References** | **Persad et al., *Nature Biotechnology*, 2023.** DOI: [10.1038/s41587-023-01716-9](https://doi.org/10.1038/s41587-023-01716-9); **Sculley, KDD 2010.** DOI: [10.1145/1772690.1772862](https://doi.org/10.1145/1772690.1772862) |

---

### 25. Paper-Driven Reproduction (`paper_repro`)

| Item | Detail |
|---|---|
| **Method** | Spec-driven provenance tracking + figure parity checks against pipeline outputs |
| **Primary checks** | Source completeness (`paper_path`, `repo_url`, `repo_commit`, `license`) and figure similarity (`mae` or `pearson`) |
| **Outputs** | `paper_repro_registry.csv`, `paper_repro_figures.csv`, `paper_repro_report.json` |
| **Implementation** | `workflow/modular/modules/paper_repro.py` |
| **References** | **Sandve et al., *PLOS Computational Biology*, 2013.** DOI: [10.1371/journal.pcbi.1003285](https://doi.org/10.1371/journal.pcbi.1003285) |

## Automated Reference Management

The pipeline includes an automated mechanism to keep the `Complete Citation List` up-to-date:

1. **Module-Level Citations**: New modules should follow [`workflow/modular/modules/module_template.py`](workflow/modular/modules/module_template.py) and define a `__references__` dictionary. The template also includes `requires_keys`, `provides_keys`, and `mutates_structure` declarations for pipeline contract validation.
2. **Auto Discovery**: [`scripts/update_references.py`](scripts/update_references.py) parses `__references__` and also scans module source for DOI patterns, then enriches metadata via Crossref API when available.
3. **Pre-commit Hook**: [`.githooks/pre-commit`](.githooks/pre-commit) enforces sync; if references changed, commit is blocked until README is staged.
4. **Enable Hook**: run `git config core.hooksPath .githooks`.

## Validation Protocol (LUSC Real Run)

Use [`scripts/validate_optimizations.py`](scripts/validate_optimizations.py) to validate acceleration without sacrificing reliability.

- **Baseline set**: LUSC 3K (`data/raw/lung_carcinoma_3k_count/outs/filtered_feature_bc_matrix`).
- **Control experiment**: baseline vs optimized with same hardware, same seed, repeated runs (`--repeats`).
- **Core metrics**: ARI, pairwise Precision/Recall/F1, marker-gene F1, runtime, data-quality audit (`data_quality_audit.csv`).
- **Stats**: 95% CI + significance test (`runtime_p_value`, `accuracy_p_value`) with threshold `p < 0.05`.
- **Tolerance**: max metric loss threshold configurable (`--loss-threshold-pct`, default `0.5`).
- **Robustness tests**: extreme dropout, low cell count, NaN injection sanitization.
- **Rollback policy**: if thresholds fail, report emits `ROLLBACK_REQUIRED` and automated rollback steps.

Example:

```bash
python scripts/validate_optimizations.py --repeats 5 --warmup-runs 1 --seed 42
```

Outputs are written to `reports/validation/`:
- `optimization_validation_report.json`
- `optimization_validation_report.md`
- `perf_runs.csv`
- `accuracy_runs.csv`
- `data_quality_audit.csv`
- `robustness_results.csv`
- `runtime_boxplot.png`
- `metric_distribution.png`
- `marker_f1_distribution.png`

## Standardized Module Docs

Use [`docs/MODULE_TECH_DOC_TEMPLATE.md`](docs/MODULE_TECH_DOC_TEMPLATE.md) for every new module technical document. The template standardizes:
- Functional description
- Implementation details and complexity
- Benchmark and validation evidence
- Citation and best-practice references

---

### Complete Citation List

| # | Reference | DOI | Used by |
|---|---|---|---|
| 1 | Patel et al., *Science*, 2014 | [10.1126/science.1254257](https://doi.org/10.1126/science.1254257) | `evolution` (auto-detected from module source) |
| 2 | Tirosh et al., *Science*, 2016 | [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) | `evolution` (auto-detected from module source) |
| 3 | Gao et al., *Nature Biotechnology*, 2021 | [10.1038/s41587-020-00795-2](https://doi.org/10.1038/s41587-020-00795-2) | `evolution` (auto-detected from module source) |
| 4 | Argelaguet et al. et al., *Genome Biology*, 2020 | [10.1186/s13059-020-02015-1](https://doi.org/10.1186/s13059-020-02015-1) | `multimodal_integration` (MOFA factor model for multi-omics integration.) |
| 5 | Argelaguet et al., *Genome Biology*, 2020 | [10.1186/s13059-020-02015-1](https://doi.org/10.1186/s13059-020-02015-1) | `multimodal_integration` (auto-detected from module source) |
| 6 | Hao et al. et al., *Cell*, 2021 | [10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048) | `multimodal_integration` (Seurat WNN: weighted nearest neighbours for multimodal joint embeddings.) |
| 7 | Hao et al., *Cell*, 2021 | [10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048) | `multimodal_integration` (auto-detected from module source) |
| 8 | Sandve et al. et al., *PLOS Computational Biology*, 2013 | [10.1371/journal.pcbi.1003285](https://doi.org/10.1371/journal.pcbi.1003285) | `paper_repro` (Guides provenance capture and reproducibility reporting.) |
| 9 | Sandve et al., *PLoS Computational Biology*, 2013 | [10.1371/journal.pcbi.1003285](https://doi.org/10.1371/journal.pcbi.1003285) | `paper_repro` (auto-detected from module source) |
| 10 | Stoeckius et al. et al., *Nature Methods*, 2017 | [10.1038/nmeth.4380](https://doi.org/10.1038/nmeth.4380) | `protein_adt` (CITE-seq protocol; CLR normalization for ADT counts.) |
| 11 | Stoeckius et al., *Nature Methods*, 2017 | [10.1038/nmeth.4380](https://doi.org/10.1038/nmeth.4380) | `protein_adt` (auto-detected from module source) |
| 12 | Mulè et al. et al., *Nature Communications*, 2022 | [10.1038/s41467-022-29356-8](https://doi.org/10.1038/s41467-022-29356-8) | `protein_adt` (DSB normalization (stretch-goal stub).) |
| 13 | Mulè et al., *Nature Communications*, 2022 | [10.1038/s41467-022-29356-8](https://doi.org/10.1038/s41467-022-29356-8) | `protein_adt` (auto-detected from module source) |
| 14 | Chen et al. et al., *Science*, 2015 | [10.1126/science.aaa6090](https://doi.org/10.1126/science.aaa6090) | `spatial_ingest` |
| 15 | Chen et al., *Science*, 2015 | [10.1126/science.aaa6090](https://doi.org/10.1126/science.aaa6090) | `spatial_ingest` (auto-detected from module source) |
| 16 | 10x Genomics et al., *Unknown*, 2020 | N/A | `spatial_ingest` (Spatial barcoded array; tissue_positions_list.csv schema.) |
| 17 | 10x Genomics et al., *Unknown*, 2023 | N/A | `spatial_ingest` (Subcellular in-situ assay; same (x,y) per-cell schema.) |
| 18 | Palla et al. et al., *Nature Methods*, 2022 | [10.1038/s41592-021-01358-2](https://doi.org/10.1038/s41592-021-01358-2) | `spatial_neighborhoods` |
| 19 | Palla et al., *Nature Methods*, 2022 | [10.1038/s41592-021-01358-2](https://doi.org/10.1038/s41592-021-01358-2) | `spatial_neighborhoods` (auto-detected from module source) |
| 20 | Haghverdi et al., *Nature Methods*, 2016 | [10.1038/nmeth.3971](https://doi.org/10.1038/nmeth.3971) | `trajectory` (auto-detected from module source) |
| 21 | Wolf et al., *Genome Biology*, 2019 | [10.1186/s13059-019-1663-x](https://doi.org/10.1186/s13059-019-1663-x) | `trajectory` (auto-detected from module source) |
---

## Results

Each pipeline run creates an independent folder under `results/`. Example verified run:

```
results/LUSC_REVIEW_OPT_FIX_20260405_20260405_142200/
├── cellranger/                   (data loading)
├── qc/                           (4 QC plots)
├── doublet_detection/            (1 doublet score plot)
├── clustering/                   (2 plots: PCA elbow, UMAP)
├── annotation/                   (3 plots + 2 CSV)
├── cell_cycle/                   (1 plot + 1 CSV)
├── cnv_inference/                (2 plots + 1 CSV + 1 JSON)
├── differential_expression/      (3 plots + 3 CSV)
├── gene_regulatory_network/      (1 plot + 1 CSV + 1 JSON)
├── gene_signature_scoring/       (3 plots + 2 CSV)
├── trajectory/                   (4 plots + 3 CSV)
├── cell_communication/           (1 plot + 1 CSV)
├── immune_phenotyping/           (5 plots + 2 CSV)
├── tumor_microenvironment/       (5 plots + 3 CSV)
├── pathway_analysis/             (1 plot + 1 CSV)
├── pseudo_velocity/              (4 plots + 2 CSV)
├── final_adata.h5ad
├── evolution/                    (4 plots + 3 CSV)
├── module_status.csv             (22/22 modules: all OK)
└── run_manifest.json
```

## Legacy Workflows

Old entry points are still available: `workflow/standard.py`, `workflow/velocity.py`, `workflow/benchmark.py`.

Recommended: use `python -m workflow.modular.cli` for all new analysis.
