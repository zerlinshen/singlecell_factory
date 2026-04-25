# 2026-04-24 - NC2024 CSS fidelity audit

## Objective

- Quantify the memory boundary behind the repeated `large` capacity failures.
- Create a rerunnable clean-vs-CSS benchmark lane for the fresh NC2024 full-cohort `massive` result.
- Measure exact cluster fidelity, broad label fidelity, rare-group behavior, and `ELF3` distribution without launching another full-cohort stage-1 run.

## Starting Context

- Fresh direct `massive` rerun:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329`
- Fresh final object:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/final_adata.h5ad`
- Carry-forward lessons:
  - treat `large` as a capacity probe
  - do not rerun prepare
  - do not rerun full-cohort stage-1 without a new blocker

## What Was Run

- Host:
  - `ubuntu-tail`
- Working directory:
  - `/home/zerlinshen/singlecell_factory`
- Execution mode:
  - `benchmark`
- New rerunnable script:
  - `/home/zerlinshen/singlecell_factory/scripts/nc2024_css_fidelity_audit.py`
- New test:
  - `/home/zerlinshen/singlecell_factory/tests/test_nc2024_css_fidelity_audit.py`
- Canonical 100k benchmark:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_CSS_FIDELITY_AUDIT_100K_AUTO_20260423_203921`
- Evidence-only 50k benchmark with the same script:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_CSS_FIDELITY_AUDIT_AUTO_20260423_203821`
- Superseded first-pass 50k benchmark before full-cohort ELF3 summary was added:
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_CSS_FIDELITY_AUDIT_AUTO_20260423_203603`

## Outcome

- `success`
- No prepare rerun was launched.
- No new `NC2024_NSCLC_FULL_COHORT_STAGE1_*` result directory was created.
- The benchmark used the existing fresh `final_adata.h5ad`.

## Evidence

- Script checks:
  - `python -m py_compile scripts/nc2024_css_fidelity_audit.py`
  - `pytest -q -o addopts="" tests/test_nc2024_css_fidelity_audit.py`
  - result: `6 passed`
- 100k benchmark metrics:
  - subset size: `100000`
  - HVGs used: `1000`
  - clean clusters: `19`
  - CSS clusters: `17`
  - clean-vs-CSS ARI: `0.4190380411518998`
  - clean-vs-CSS NMI: `0.6539602054106193`
  - CSS-vs-full-massive Leiden ARI: `0.5478705324241598`
  - CSS-vs-full-massive Leiden NMI: `0.7219989456991266`
  - peak RSS: `12779.69140625 MB`
- Memory math confirmed from HDF5:
  - final object `X` is CSR `810218 x 30374`
  - CSR core: `11.169010004 GB`
  - dense `float32`: `98.438246128 GB`
  - dense `float64`: `196.876492256 GB`
- `ELF3` full-object distribution:
  - positive cells: `56809 / 810218` (`7.01%`)
  - tumor epithelial: `36158 / 54158` positive (`66.76%`)
  - non-immune bucket: `36915 / 70716` positive (`52.20%`)
  - Myeloid/Macro: `9975 / 283360` positive (`3.52%`), `98` high cells

## Interpretation

- The 96GB memory issue is explained by dense materialization risk, not by the sparse final object itself.
- CSS preserves a useful broad biological representation but should not be described as exact Leiden-cluster equivalence.
- The benchmark supports the existing caution that rare groups need targeted follow-up; exact cluster identity is the weak point.
- `ELF3` is primarily a tumor epithelial / non-immune signal in the full object, with a smaller Myeloid/Macro-positive subset suitable for later targeted biological audit.

## What Changed

- Added remote script:
  - `/home/zerlinshen/singlecell_factory/scripts/nc2024_css_fidelity_audit.py`
- Added remote test:
  - `/home/zerlinshen/singlecell_factory/tests/test_nc2024_css_fidelity_audit.py`
- Added local phase-3 reproduction package:
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/reproduction_packages/NC2024_phase3_css_fidelity_audit_20260424`
- Updated local reproduction README to point to the phase-3 package and canonical audit path.

## Cautions for the Next Run

- Do not treat the clean-vs-CSS ARI/NMI as a full 810k all-at-once gold-standard comparison; the 810k all-at-once lane remains capacity-limited.
- Use the 100k benchmark as the current best controlled estimate for CSS fidelity.
- Use `elf3_full_cohort_group_summary.csv` for unbiased full-object `ELF3` distribution; the benchmark subset intentionally enriches ELF3-high cells for audit sensitivity.
- If deeper rare-cluster precision is needed, run targeted rare-group audits rather than another full-cohort stage-1 run.

## Improvement Ideas

- Add a focused `ELF3` Myeloid/Macro follow-up comparing marker context, patient distribution, and tumor-epithelial ambient/doublet signals.
- Add per-cell-type stratified clean-vs-CSS benchmarks for specific rare groups if paper-facing claims require cluster-level precision.

## Classification

- `canonical`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_CSS_FIDELITY_AUDIT_100K_AUTO_20260423_203921`
  - `/home/zerlinshen/singlecell_factory/scripts/nc2024_css_fidelity_audit.py`
  - `/home/zerlinshen/singlecell_factory/tests/test_nc2024_css_fidelity_audit.py`
- `evidence-only`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_CSS_FIDELITY_AUDIT_AUTO_20260423_203821`
- `superseded`
  - `/home/zerlinshen/singlecell_factory/results/NC2024_NSCLC_CSS_FIDELITY_AUDIT_AUTO_20260423_203603`
- `failed exploratory`
  - none
