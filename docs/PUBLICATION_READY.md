# Publication-Ready Methods Template — NC2024 NSCLC Reproduction

This file provides copy-paste-ready Methods section text and citation patterns for manuscripts using this pipeline's NC2024 outputs.

---

## Methods Section Template

### Single-cell RNA-seq preprocessing and analysis

Single-cell preprocessing and analysis followed Sanchez-Mejias et al. (2024) [1] with the modifications described below. Raw count matrices were loaded from ArrayExpress E-MTAB-13526 (CellRanger v3.1.0 output). Doublets were identified per sample using Scrublet [2] (expected_doublet_rate=0.06, n_prin_comps=30). Cells failing QC thresholds (minimum counts, minimum genes, maximum mitochondrial fraction) were excluded. The cohort was split into a tumor sub-dataset and a background+healthy sub-dataset prior to integration, matching the original study design.

Highly variable genes were selected (Scanpy, n_top_genes=2000) and principal components computed using sparse-aware truncated SVD (ARPACK) to avoid in-memory densification of the ~884k-cell matrix — mathematically equivalent to standard PCA on sparse input. Batch effects were corrected using Harmony [3] (batch_key=sample, n_pcs=15), run independently on each sub-dataset. Cells were clustered using the Leiden algorithm [4] (resolution=1.0) on the 15-dimensional harmonised embedding and visualised by UMAP [5].

Cluster marker genes were identified by Wilcoxon signed-rank test [6] with Bonferroni correction, requiring expression in ≥30% of cells per cluster (min_pct=0.30). Cell type identities were assigned by manual review of canonical lineage markers.

Full environment, argument hashes, and package versions for each cohort run are recorded in the run ledger files:
- `ops/run_ledger/nc2024_tumor_20260426T041139.json`
- `ops/run_ledger/nc2024_bh_20260426T052950.json`

Parameter rationale and paper alignment details are documented in `ops/nc2024_methodology_audit/AUDIT_2026-04-26_v2.md`.

### Cell type proportions and alignment

Cell type proportions and cluster-to-cell-type alignment are reported in `ops/nc2024_methodology_audit/ALIGNMENT_REPORT_v2_2026-04-26.md`.

---

## Citation Patterns

Use these inline patterns when drafting the manuscript:

```
[paper alignment]
"Single-cell preprocessing followed Sanchez-Mejias et al. (2024) [1] with
sparse-exact engine modifications (see ops/run_ledger/<run>.json for full
environment and arguments)."

[cell type proportions]
"Cell type proportions were derived from Harmony-integrated Leiden clusters
annotated as described in ALIGNMENT_REPORT_v2_2026-04-26.md."

[parameter rationale]
"Parameter choices (n_pcs=15, leiden_resolution=1.0, de_correction=bonferroni,
de_min_pct=0.30) are justified against the source paper in
AUDIT_2026-04-26_v2.md."
```

> **These are NOT the pipeline defaults.** As of 2026-07-28 the canonical
> profile is `n_pcs=40`, `leiden_resolution=0.8` on both the CLI and the
> `config.py` dataclasses. Reproduce the clustering geometry above with
> `--scientific-profile paper-15pc --acknowledge-scientific-non-equivalence`;
> `de_correction=bonferroni` / `de_min_pct=0.30` remain launcher-level flags
> (`--de-correction`, `--de-min-pct`). Confirm against the run's
> `resolved_scientific_parameter_diff` before quoting any of these in a
> manuscript.

---

## References

1. Sanchez-Mejias et al., *Nature Communications* 15:4388 (2024). DOI: 10.1038/s41467-024-48700-8
2. Wolock et al., *Cell Systems* (2019). DOI: 10.1016/j.cels.2018.11.005
3. Korsunsky et al., *Nature Methods* (2019). DOI: 10.1038/s41592-019-0619-0
4. Traag et al., *Scientific Reports* (2019). DOI: 10.1038/s41598-019-41695-z
5. McInnes et al., *JOSS* (2018). DOI: 10.21105/joss.00861
6. Wilcoxon, *Biometrics Bulletin* (1945). DOI: 10.2307/3001968

---

## Key Paths for Supplementary Material

| Item | Path |
|---|---|
| Methodology audit | `ops/nc2024_methodology_audit/AUDIT_2026-04-26_v2.md` |
| Alignment report | `ops/nc2024_methodology_audit/ALIGNMENT_REPORT_v2_2026-04-26.md` |
| Tumor run ledger | `ops/run_ledger/nc2024_tumor_20260426T041139.json` |
| B/H run ledger | `ops/run_ledger/nc2024_bh_20260426T052950.json` |
| Git commit | `fdc48efb1becd940472388d5a41c4c281b6b158a` |
| Conda env | `sc_gpu` |
