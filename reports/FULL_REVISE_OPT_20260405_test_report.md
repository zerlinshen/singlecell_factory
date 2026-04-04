# Full scRNA Pipeline Revision & Test Report

Date: 2026-04-05
Workspace: /home/zerlinshen/singlecell_factory

## 1) Code Revisions Applied

### A. Pipeline robustness (mandatory doublet stage)
- File: workflow/modular/modules/doublet_detection.py
- Change:
  - Added robust fallback when Scrublet is unstable on tiny/degenerate inputs.
  - New behavior: mark all cells as singlets instead of aborting pipeline.
  - Added metadata key `doublet_method` with values:
    - `scrublet`
    - `fallback_all_singlets`
  - Threshold plotting is now conditional when Scrublet threshold exists.

### B. Timezone-safe manifest timestamp
- File: workflow/modular/pipeline.py
- Change:
  - Replaced deprecated `datetime.utcnow()` with timezone-aware UTC timestamp generation.

### C. RNA velocity module hardening (already integrated in this workspace)
- File: workflow/modular/modules/rna_velocity.py
- Highlights:
  - Numeric version guard for NumPy2/scVelo compatibility patch.
  - Non-mutating copy-on-write behavior retained.
  - Improved BAM counting worker architecture and warning logs.

### D. Documentation sync
- README.md:
  - Added RNA velocity output directory details.
  - Added doublet fallback behavior note.
- PROTOCOL.md:
  - Added explicit RNA velocity module section and output files.
  - Added doublet fallback behavior note.

## 2) Automated Test Matrix

### Command A
```bash
NUMBA_CACHE_DIR=/tmp/numba_cache PYTHONPATH=. pytest -q --no-cov
```
Result: **40 passed, 2 warnings**

### Command B (with coverage gates)
```bash
NUMBA_CACHE_DIR=/tmp/numba_cache PYTHONPATH=. pytest -q
```
Result: **40 passed, 2 warnings**
Coverage summary:
- TOTAL: **90.87%** (threshold: 90%)
- benchmark.py: 87%
- standard.py: 92%
- velocity.py: 96%

### Command C (timing profile)
```bash
NUMBA_CACHE_DIR=/tmp/numba_cache PYTHONPATH=. pytest --no-cov --durations=10
```
Result: **40 passed, 2 warnings**
Slowest test:
- tests/test_modular.py::test_modular_pipeline_minimal (~0.76s)

## 3) Real Data End-to-End Pipeline Run

### Command
```bash
NUMBA_CACHE_DIR=/tmp/numba_cache PYTHONPATH=. python -m workflow.modular.cli \
  --project FULL_REVISE_OPT_20260405 \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,cell_cycle,differential_expression,annotation,trajectory,pseudo_velocity,cnv_inference,pathway_analysis,cell_communication,gene_regulatory_network,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,evolution \
  --parallel-workers 4
```

### Output
- run dir: results/FULL_REVISE_OPT_20260405_20260405_014053
- final_adata.h5ad: generated
- run_manifest.json: generated
- module_status.csv: generated
- files under run dir (<=2 depth): 72

### Module Status
- Mandatory:
  - cellranger: ok
  - qc: ok
  - doublet_detection: ok
- Optional (requested): all **14/14 ok**
  - clustering, cell_cycle, differential_expression, annotation,
    trajectory, pseudo_velocity, cnv_inference, pathway_analysis,
    cell_communication, gene_regulatory_network, immune_phenotyping,
    tumor_microenvironment, gene_signature_scoring, evolution

### Key Metadata Snapshot
- raw_cells: 2588
- cells_after_qc: 2382
- doublet_method: scrublet
- doublets_detected: 5
- cells_after_doublet_removal: 2377
- n_clusters: 17
- de_significant_genes: 6551
- n_clones: 3

## 4) Warnings Observed (non-fatal)

1. Matplotlib tight_layout warning in QC plotting.
2. FutureWarning from scanpy PCA argument (`use_highly_variable`).
3. Colorbar warnings in some plotting functions.
4. Pandas `groupby(observed=...)` future behavior warning in evolution.

All warnings were non-blocking; no module failed.

## 5) Residual Gap

- RNA velocity real-data run was not included in the end-to-end command above because required velocity inputs are incomplete in current data tree:
  - BAM exists.
  - GTF/loom not found under `data/`.

To run RNA velocity end-to-end, provide either:
- `--velocity-loom <file.loom>`
or
- `--velocity-bam <bam> --velocity-gtf <genes.gtf.gz>`

