# scRNA-seq Pipeline Audit Report

## Executive Summary
Overall audit judgment: **this project is not ready to be trusted as a final validation pipeline in its current default state**.

What is solid:
- Cell Ranger outputs are structurally valid and plausible for PBMC 1k-like input.
- Core Scanpy mechanics (normalization, HVG, PCA, neighbors, UMAP, Leiden, marker ranking) execute when launched correctly.
- `adata.raw` is preserved before HVG subsetting, so full-gene marker visualization is possible.

What is blocking:
- The advertised end-to-end entry point is functionally broken.
- Default QC settings overfilter this PBMC-like dataset (151/1221 cells retained, 12.4%).
- “Doublet detection enabled” is not operational in this environment (fails and silently degrades to all non-doublet).
- No actual cell-type annotation is performed despite repeated claims.
- Reproducibility and module-boundary guarantees in docs are overstated and partially false.

## Project Overview
Audited scope included:
- Code: `run_end_to_end.sh`, `downstream_scanpy/run_scanpy_pipeline.py`, `downstream_scanpy/scanpy_pipeline_optimized.py`
- Data and outputs: `pbmc_1k_v3_count/outs/*`, `cellranger_runs/pbmc1k_practice/outs/*`
- Documentation and claims: `README.md`, `QUICKSTART.md`, `PROJECT_STRUCTURE.md`, `PIPELINE_IMPROVEMENTS.md`, `CHANGES_SUMMARY.md`, `COMPLETED.md`, `downstream_scanpy/pipeline_comparison.md`

Practical verification runs executed during audit:
- Optimized default run: `output/audit_default/`
- Simple run: `output/audit_simple/`
- Optimized run with adaptive thresholds: `output/audit_suggested/`

## Reconstructed Pipeline
### True entry points found
1. `run_end_to_end.sh` (intended top-level orchestrator; currently broken)
2. `downstream_scanpy/scanpy_pipeline_optimized.py` (working when invoked with correct CWD/env)
3. `downstream_scanpy/run_scanpy_pipeline.py` (working when invoked with correct CWD/env)

### Actual execution order implemented (optimized script)
1. **Data loading**: `sc.read_10x_mtx(...)` from `filtered_feature_bc_matrix`
2. **QC metrics**: mito/ribo flags + `calculate_qc_metrics`
3. **QC suggestion step**: MAD/percentile-based advisory thresholds
4. **Doublet detection attempt**: Scrublet call
5. **QC filtering**: `min_genes`, `max_genes`, `max_mito_pct`, optional predicted_doublet exclusion
6. **Gene filtering**: `min_cells`
7. **Save intermediate**: `after_qc.h5ad`
8. **Normalization/log1p**
9. **HVG selection** (`n_top_genes=2000`)
10. **Scaling** (HVG subset)
11. **PCA**
12. **Neighbors graph**
13. **UMAP**
14. **Leiden clustering**
15. **Marker detection** (`rank_genes_groups`, Wilcoxon, `use_raw=True`)
16. **Marker plotting** (top markers + canonical PBMC panel)
17. **Save outputs**: processed h5ad, CSVs, figures, logs

### Cell Ranger upstream stage (verified from real run metadata)
- Command recorded in `pbmc_1k_v3_count/_cmdline`: Cell Ranger count over PBMC v3 FASTQs + GRCh38-2024-A reference.
- Matrix integrity verified:
  - Filtered matrix dims: **38,606 genes x 1,221 cells**
  - NNZ: **4,099,222**
- `metrics_summary.csv` is internally plausible (e.g., median genes/cell 3,290; fraction reads in cells 95.6%).

## Code Audit
### A. Entry point and orchestration
- `run_end_to_end.sh` is not in sync with current Scanpy CLI:
  - Uses obsolete `--outdir` and omits required `--project`.
  - `--skip-cellranger` and `--skip-scanpy` flags only print messages; they do not alter control flow.
  - Conda env check uses fragile grep pattern and falsely reports missing env.
- Net effect: advertised end-to-end route is non-functional and misleading.

### B. Data loading and path handling
- Both Python pipelines default to `base_output_dir='../output'` (relative path).
- This is CWD-dependent and can write outside the repo if launched from root or another directory.
- This directly weakens module boundary guarantees.

### C. QC implementation
- Correct: mitochondrial and ribosomal flags are computed (for `gene_symbols` mode).
- Critical logic gap: `max_ribo_pct` exists in config but is never applied in filtering.
- Default QC for PBMC-like data is overly strict in this dataset:
  - `pct_counts_mt < 5` retains only 157/1221 cells before other filters.
  - Full default mask retains 151/1221 cells (12.4%).

### D. Normalization / HVG / scaling
- Correct sequencing: normalize -> log1p -> set `adata.raw` -> HVG subset -> scaling.
- `adata.raw` preservation avoids the classic full-gene-loss bug.
- Minor inconsistency: comment says “Store raw counts” but `adata.raw` stores log-normalized matrix.

### E. PCA / neighbors / UMAP / clustering
- Core implementation is standard Scanpy practice.
- Reproducibility weakness: no fixed `random_state` for UMAP/Leiden.

### F. Marker detection and annotation
- Marker ranking is run with `use_raw=True` (good).
- Canonical marker plotting checks `adata.raw.var_names` (good in gene-symbol mode).
- **No cell-type annotation assignment is implemented**:
  - No `cell_type` column added to `adata.obs`
  - No annotation decision logic output
  - Claims of “cell type annotation” in docs/code header are not supported by implementation.

### G. Doublet handling
- Code attempts Scrublet, but in this environment it fails due missing `skimage` when threshold is `None`.
- Exception is swallowed; pipeline proceeds with `predicted_doublet=False` for all cells.
- `doublet_threshold` CLI/config is parsed but never used in Scrublet call.
- Result: doublet removal is effectively non-operational while giving a “feature present” impression.

### H. Robustness and reproducibility
- Import-time failure occurs without `NUMBA_CACHE_DIR` workaround in this environment.
- No `environment.yml`/`requirements.txt` provided despite scripts/docs referencing env creation.
- Documentation is inconsistent and partially stale (`--outdir` examples remain in multiple files).

## Scientific Audit
### PBMC suitability of current defaults
- For this PBMC-like 10x dataset, default QC is scientifically too stringent:
  - Median mito is ~6.76%; fixed cutoff at 5% is inappropriate as default here.
  - Resulting 12.4% retention is a severe distortion and collapses biological structure.

### Dimensional reduction and clustering
- Method choices are standard.
- Default output quality is poor because QC pre-filtering removes most biological diversity.
- Under adaptive thresholds (`min_genes=1640`, `max_mito_pct=15.8`), clustering becomes biologically plausible (8 clusters, 1097 cells).

### Marker-based annotation defensibility
- Current pipeline does not perform annotation decisions; only plotting/marker tables.
- With default filtering, cluster marker evidence is weak/mixed for PBMC lineage calling.
- With suggested thresholds, canonical PBMC structure (T/NK/B/monocyte axes) is much clearer and broadly defensible as first-pass annotation.

### Contaminants / artifacts
- No ambient RNA correction module.
- No robust doublet removal in observed environment (operational failure).
- No explicit RBC/platelet contaminant handling step beyond QC/marker interpretation.
- Therefore downstream biological claims should be conservative.

## Accuracy / Reliability Assessment
Rating scale: **High / Moderate / Low confidence**

| Dimension | Confidence | Basis |
|---|---|---|
| Code correctness confidence | **Moderate-Low** | Core Scanpy steps work; orchestration and several config behaviors are broken/inconsistent. |
| Reproducibility confidence | **Low** | Broken wrapper, missing env spec, hidden env-variable requirements, stale docs. |
| Annotation confidence | **Low** | No implemented annotation logic; only marker visualization/ranking. |
| Biological interpretation confidence (default run) | **Low** | Overfiltering leaves 151 cells and collapses expected PBMC structure. |
| Biological interpretation confidence (suggested thresholds) | **Moderate** | Lineages become plausible, but no robust doublet/ambient control and no formal annotation framework. |
| Risk of false annotation | **High** | Manual interpretation required; mixed clusters under default settings; no explicit confidence framework in code. |

Strong evidence:
- Cell Ranger matrix integrity and metrics are plausible.
- Adaptive-threshold run recovers expected PBMC-like structure.

Weak evidence:
- Publication-strength cell typing, doublet control, contaminant control, and claim calibration.

## Reproducibility Assessment
Verified issues:
1. `run_end_to_end.sh` cannot currently drive a valid full run.
2. Conda env check in wrapper is brittle and fails in this environment despite existing `sc10x`.
3. No reproducible environment lockfile/spec (`environment.yml`, requirements, lock) is provided.
4. Pipelines are sensitive to launch directory due relative output path.
5. Hidden runtime requirements (`NUMBA_CACHE_DIR`, writable Matplotlib config) are undocumented.

Positive:
- When directly run from `downstream_scanpy/` with explicit env setup, scripts complete and generate consistent output structure.

## Major Problems Found
1. **Broken end-to-end wrapper** (`run_end_to_end.sh`) due obsolete CLI usage and ineffective skip flags.
2. **Scientifically problematic default QC** for this PBMC-like dataset (12.4% retention).
3. **Doublet detection is non-operational in practice** under current environment and silently degrades.
4. **No real annotation module** despite repeated claims of marker-based annotation.
5. **Boundary/reproducibility guarantees are overstated** (CWD-dependent output path, stale docs, missing env spec).
6. **`gene_ids` mode breaks QC marker logic** (MT/ribo detection and canonical marker lookup rely on `var_names` gene symbols).

## Minor Problems Found
1. `max_ribo_pct` exists but is unused.
2. `doublet_threshold` exists but is unused.
3. No fixed random seeds for stochastic steps.
4. Comments/documentation contain stale statements and old output conventions.
5. Heavy marker plotting logs are noisy and can obscure key warnings.

## What Is Already Good
1. `adata.raw` is set before HVG subsetting.
2. Marker ranking is performed using raw layer in optimized pipeline.
3. Output structure (when run correctly) is organized and useful (`data/`, `figures/`, `tables/`, `logs/`, `cache/`).
4. Adaptive QC suggestion mechanism is practically useful and empirically improved this dataset.
5. Cell Ranger input/output integrity is solid.

## High-Priority Fixes
1. **Fix orchestration immediately**:
   - Update `run_end_to_end.sh` to pass `--project`.
   - Implement actual skip logic.
   - Fix conda env detection.
2. **Anchor output paths to script location**, not CWD.
3. **Make doublet detection truly functional**:
   - Use `sc.pp.scrublet` (current API) and wire `doublet_threshold`.
   - Fail hard or explicit status when doublet detection is requested but unavailable.
4. **Implement explicit annotation stage**:
   - Add cluster-level scoring/rules and `adata.obs['cell_type']` + confidence flags.
5. **Ship reproducible environment files** and document required runtime env vars.
6. **Change default QC for PBMC-like data** to avoid destructive overfiltering; current default is not fit-for-purpose for this dataset.

## Medium-Priority Improvements
1. Apply `max_ribo_pct` or remove it from config/docs.
2. Add random seeds for reproducible UMAP/Leiden.
3. Add automated QC report of filter contribution (which rule removed how many cells).
4. Add optional ambient RNA correction and explicit RBC/platelet contamination checks.
5. Clean outdated docs referencing `--outdir` and old `results_*` folders.

## Publication-Readiness Assessment
Current status: **Not publication-ready**.

Reasons:
- Execution/reproducibility defects in official wrapper.
- No robust operational doublet handling in current environment.
- No actual annotation module despite claims.
- Default run quality is not biologically reliable for PBMC interpretation.
- Insufficient provenance hardening (env lock, seed control, validated parameter presets).

Acceptable current use:
- Learning/demo: yes (with caution and manual supervision)
- Internal exploratory analysis: conditional yes (after QC parameter correction and explicit caveats)
- Publication-grade claims: no

## Comparison to Standard 10x Expected Output
### Benchmark against standard 10x-to-Scanpy expectations

| Item | Observed in project | Judgment |
|---|---|---|
| Cell Ranger filtered matrix structure | Valid MTX + barcodes + features; dims 38,606 x 1,221 | **Matches standard expectation** |
| Cell Ranger QC plausibility | 1,221 cells, median genes 3,290, high fraction reads in cells | **Matches standard expectation** |
| AnnData dims after load | 1,221 x 38,606 | **Matches standard expectation** |
| QC metadata fields after QC calc | Standard Scanpy QC fields present | **Matches standard expectation** |
| Normalization/log1p behavior | Standard normalize_total + log1p | **Matches standard expectation** |
| HVG selection behavior | 2,000 HVGs selected as configured | **Matches standard expectation** |
| PCA/neighbors/UMAP/Leiden outputs | Present and internally consistent in processed h5ad | **Matches standard expectation** |
| Default filtering outcome | 151/1221 cells retained (12.4%), 3 clusters | **Suspicious deviation** |
| Suggested-threshold outcome | 1097/1221 cells retained (89.8%), 8 clusters | **Acceptable deviation / closer to standard** |
| Major PBMC lineage recovery (default) | Weak/partially collapsed | **Suspicious deviation** |
| Major PBMC lineage recovery (suggested) | T/NK/B/mono axes broadly recoverable | **Matches/acceptable** |
| Doublet handling | Claimed but not operational in tested env | **Likely workflow artifact / likely error** |
| Marker-based annotation deliverable | No explicit annotation assignments produced | **Likely error vs stated goals** |

### Explicit classification summary
- **Matches standard expectation**: input integrity, load dimensions, core Scanpy transforms, existence of embeddings/clusters/DE tables.
- **Acceptable deviation**: adaptive-threshold run yielding 8 clusters from PBMC-like data.
- **Suspicious deviation**: default QC collapse to 12.4% cells and under-resolved biology.
- **Likely error/workflow artifact**: non-functional wrapper, non-operational doublet mode, absent annotation module despite claims.

## Module Boundary and Dependency Isolation Assessment
### Classification
**Mostly separated but with risk** (not fully separated).

### Is the scRNA-seq module fully independent from other project-dependent results?
**No.**

### Boundary problems identified
1. Output root is relative (`../output`) and depends on CWD; this can route outputs outside intended module boundary.
2. `run_end_to_end.sh` routes Scanpy output intent to `pbmc_1k_v3_count/analysis` (mixing downstream artifacts with Cell Ranger run folder) and is currently CLI-incompatible.
3. Documentation still contains obsolete `--outdir` and old folder conventions, increasing risk of users writing/reading from wrong locations.
4. Duplicate Cell Ranger run directories (`pbmc_1k_v3_count/` and `cellranger_runs/pbmc1k_practice/`) are not inherently wrong but create provenance ambiguity if not explicitly controlled.

### What should be reorganized
1. Hard-anchor output root to repository path derived from `__file__`.
2. Add explicit `--base_output_dir` CLI option with absolute-path normalization.
3. Keep Cell Ranger raw outputs and Scanpy outputs in separate top-level namespaces consistently.
4. Remove or archive stale docs and legacy path examples.
5. Add a single source-of-truth run manifest recording input path, script version, parameters, and output folder.

## Final Verdict
**One-sentence verdict:** The project contains a workable Scanpy core but fails strict final-audit criteria due broken orchestration, non-operational doublet handling, absent true annotation logic, and scientifically unsafe default filtering for this PBMC-like dataset.

Top 5 reasons:
1. End-to-end runner is currently not valid for the current CLI.
2. Default QC overfilters to 12.4% cells and distorts biology.
3. Doublet module is effectively disabled in tested environment.
4. Claimed annotation step is not implemented as actual label assignment.
5. Reproducibility and boundary guarantees are undermined by path/env/documentation issues.

Acceptability by use case:
- **Learning/demo**: **Acceptable with caution**
- **Internal lab exploration**: **Conditionally acceptable after high-priority fixes**
- **Figure generation**: **Conditionally acceptable only after QC and annotation fixes**
- **Publication-grade analysis**: **Not acceptable**

Explicit mandatory answer:
**Does this project’s real output look consistent with a standard and valid 10x-to-Scanpy analysis pipeline?**  
**Partially**: upstream 10x outputs and Scanpy mechanics are consistent, but the current default downstream settings and execution layer are not reliably valid for final biological interpretation.

## Real Run Reproduction and Official Result Comparison
### Real execution scope performed
Real reruns were executed from the actual official 10x Cell Ranger filtered matrix input:
- `pbmc_1k_v3_count/outs/filtered_feature_bc_matrix`

Executed reruns:
1. Optimized default: `output/rerun_opt_default`
2. Simple default: `output/rerun_simple_default`
3. Optimized suggested thresholds (`min_genes=1640`, `max_mito_pct=15.8`): `output/rerun_opt_suggested`

Note on upstream rerun boundary:
- A full FASTQ->Cell Ranger rerun was not repeated in this final check due runtime cost.
- Upstream input integrity was still validated using real files:
  - Two independent Cell Ranger outputs exist in project (`pbmc_1k_v3_count` and `cellranger_runs/pbmc1k_practice`)
  - Decompressed hashes of `matrix.mtx`, `features.tsv`, and `barcodes.tsv` are identical between them.

### Rerun vs existing project outputs (execution-based)
Compared pairs:
- `audit_default` vs `rerun_opt_default`
- `audit_simple` vs `rerun_simple_default`
- `audit_suggested` vs `rerun_opt_suggested`

Observed:
- `processed.h5ad` files: **exact md5 match** for all three pairs.
- (optimized) `after_qc.h5ad`: **exact md5 match** for both optimized pairs.
- `cluster_counts.csv`, `filtering_stats.csv` (where present), `marker_genes.csv`: **exact md5 match**.
- AnnData internals: exact match on dimensions, raw dimensions, HVG count (2000), obsm/obsp availability, `rank_genes_groups` presence.
- Figure files: same file set. Most PNGs are byte-identical; violin plots differ by md5 in optimized runs (consistent with plotting stochastic jitter), while quantitative outputs remain identical.

### Rerun vs official/standard expected PBMC-like behavior
Official baseline from real Cell Ranger outputs:
- Filtered matrix: **38,606 genes x 1,221 cells**
- Cell Ranger graph clustering: **7 clusters** over all 1,221 cells

Rerun outputs:
- Optimized default:
  - 151 cells retained (12.4%), 3 clusters
  - Mechanically reproducible, biologically overfiltered
- Optimized suggested:
  - 1,097 cells retained (89.8%), 8 clusters
  - Broad PBMC lineages evident in canonical markers:
    - CD14/S100A8/S100A9/Lyz monocytes
    - MS4A1/CD79A B-lineage clusters
    - IL7R/CD3D T-lineage clusters
    - NKG7/GNLY NK cluster
    - FCGR3A+ myeloid-like cluster

Classification:
- **Exact matches**:
  - Rerun quantitative outputs vs existing project outputs (h5ad/tables) for corresponding configurations.
  - Core object availability (PCA/neighbors/UMAP/Leiden/marker outputs).
- **Acceptable differences**:
  - Non-identical violin PNG bytes with unchanged quantitative content.
  - 8 Scanpy clusters vs 7 Cell Ranger clusters in permissive run (expected method-level difference).
- **Suspicious differences**:
  - Default rerun retains only 12.4% cells and under-represents PBMC diversity.
- **Likely workflow errors**:
  - Persistent doublet module non-operation (Scrublet failure fallback).
  - Default QC preset unsuitable for this PBMC practice dataset.

Direct conclusion:
- **Can this project be rerun successfully from the real input data?**  
  **Yes for the core Scanpy workflow from real Cell Ranger matrix input; no for the advertised end-to-end wrapper in its current state.**
- **Are the rerun outputs consistent with the existing project outputs?**  
  **Yes. Quantitative outputs are reproducible and hash-identical for matched configurations.**
- **Are the rerun outputs broadly consistent with the official/standard expected results?**  
  **Conditionally yes**: with adaptive/permissive QC settings they are broadly consistent with expected PBMC structure; with current default settings they are not.
