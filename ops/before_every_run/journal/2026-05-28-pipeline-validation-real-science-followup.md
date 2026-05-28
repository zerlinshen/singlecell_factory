# Pipeline Validation Real-Science Follow-Up - 2026-05-28

## Objective

Finish the remaining real-data scientific checks for the integration-biology and
multiomics validation round without weakening claim boundaries, then clean
reproducible bulky intermediates.

## Evidence Root

- `/home/zerlinshen/projects/pipeline-validation-20260528/`
- Plan: `RUN_PLAN.md`
- Cleanup record: `cleanup_records/cleanup_20260528.md`

## What Was Attempted

- Targeted population sensitivity for LUSC AT1, LUSC DC mature, LUSC cDC2,
  Trevino inhibitory interneuron, and Trevino intermediate progenitor.
- LUSC true-count-compatible pseudobulk DE with dataset-aware/matched design.
- Unbiased/cell-matched 3D-genome validation using GM12878 Rao DpnII Hi-C
  compartments and B35T1NC Micro-C author TADs.
- Storage cleanup of reproducible checkpoint `.h5ad` files and generated contact
  TSVs after compact evidence was retained.

## What Succeeded

- Targeted sensitivity:
  - PASS: LUSC AT1, LUSC cDC2, Trevino inhibitory interneuron, Trevino
    intermediate progenitor.
  - REVIEW: LUSC DC mature, because batch mixing improves only with an excessive
    local-purity penalty.
- LUSC origin DE:
  - Uses 71/87 true-count-compatible samples; excludes 16 fractional-count
    samples.
  - Matched/dataset-aware tumor-primary vs normal-adjacent contrast has 43
    samples across 5 datasets, permutation q<0.1 = 0, and 12/12 sentinel
    directions correct.
- 3D genome:
  - GM12878 unbiased Hi-C chr19 A/B compartments pass published-subcompartment
    concordance after sign fixing (`0.814`, 543 bins compared).
  - `hic_tad.py` now computes cross-boundary insulation and excludes sparse
    non-variable tails from compartment eigendecomposition instead of downgrading
    the whole chromosome.

## What Failed Or Remains Risky

- LUSC tumor-stage DE remains conditional single-dataset evidence (6 samples, 1
  dataset), not a final biological result.
- The current lightweight factory TAD boundary caller remains below null on
  B35T1NC Micro-C author TADs (`F1=0.049`, `recall_over_null=0.40` best tested
  setting). Keep factory TAD boundaries exploratory/visual-QC only until a
  stronger caller is integrated.
- Independent two-lane code review was not available in this Codex App session;
  the final review is therefore a conservative self-audited quality gate, not a
  merge-ready independent approval.

## What Changed

- Code:
  - `workflow/modular/modules/hic_tad.py`
  - `tests/test_wave2b_hic_smoke.py`
- Docs/run memory:
  - `README.md`
  - `docs/MULTIOMICS_MODULE_RATIONALE.md`
  - `AGENTS.md`
  - this journal entry and `LATEST.md`
- Validation scripts:
  - `scripts/run_targeted_population_sensitivity.py`
  - `scripts/run_lusc_matched_de_sanity.py`
  - `scripts/run_unbiased_hic_validation.py`

## Artifact Status

- Canonical evidence-only validation root:
  `/home/zerlinshen/projects/pipeline-validation-20260528/`
- Canonical compact summaries:
  - `targeted_sensitivity/targeted_population_sensitivity_summary.json`
  - `matched_de/lusc_matched_de_summary.json`
  - `genome3d/unbiased_hic_validation_summary.json`
- Superseded/deleted:
  - generated `genome3d/*.contacts.tsv.gz`
  - Trevino linked-pipeline `.checkpoints/*.h5ad`
  - validation-round `__pycache__`
- Preserved:
  - raw/public inputs
  - `.cool` and BED ground-truth fixtures
  - Trevino final h5ad, bundle/R outputs, manifests, reports, and ledgers

## Verification

- `python3 -m py_compile` on all three validation scripts and
  `workflow/modular/modules/hic_tad.py`: pass.
- `jq empty` on all three validation summary JSON files: pass.
- `rg '\bNaN\b|Infinity'` across validation summary JSON files: no matches.
- `pytest --no-cov -q tests/test_wave2b_hic_smoke.py`: 14 passed in 0.49s.
- `git diff --check`: pass.
- Cleanup verification:
  - no Trevino checkpoint `.h5ad` remains under `.checkpoints/`;
  - no generated contact TSV remains under the validation `genome3d/`;
  - validation root is 1.8M;
  - `/home/zerlinshen` has 168G available.

## Next Operator Memory

- Do not call this round a clean global pass.
- Supported claim: bounded real-data validation improved materially and several
  prior blockers are now resolved.
- Non-final claim boundaries: LUSC DC mature, tumor-stage DE, and factory TAD
  biology remain review/conditional.
- Next real scientific work should integrate or wrap a stronger TAD caller and
  validate it against the retained B35T1NC/cooltools-style ground truth before
  using TAD boundaries as biological evidence.
