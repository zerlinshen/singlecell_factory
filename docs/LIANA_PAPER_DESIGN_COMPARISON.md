# LIANA lane vs De Zuani 2024 paper design — comparison table

**Purpose:** pin down exactly how the factory's cell–cell communication
evidence lanes (Wave-6 LUAD-vs-LUSC LUCA LIANA, Wave-8/9 patient-stratified
recurrence) differ from the baseline paper's design, and why the claim class
stays `partial`.

Baseline paper: De Zuani et al., *Nat Commun* 15:4388 (2024),
DOI `10.1038/s41467-024-48700-8` (PMC `PMC11116453`).
Paper scRNA data: E-MTAB-13526. Paper spatial data: E-MTAB-13530.

## Paper design (verbatim from Methods / Results)

| Aspect | Paper (CellPhoneDB) |
|---|---|
| Cohort | Own cohort: 25 treatment-naive NSCLC patients (LUAD + LUSC) + background + 2 healthy donors |
| Cells in analysis | LUAD n=105,749 cells; LUSC n=230,066 cells |
| Conditions | tumour / background / healthy as **independent datasets**, CellPhoneDB run separately per condition |
| Engine | CellPhoneDB (statistical permutation framework), ref 29 |
| Scaling | Did **not** run complete datasets; **stratified 50% subsample** preserving cell-type / patient / sample proportions (permutation tests do not scale to >10⁶ cells) |
| Pair filter 1 | LR pair expressed in **≥30% of cells** in the cell-type cluster of interest |
| Pair filter 2 | mean **log(1+expression) > 1.0** for the LR pair |
| Pair filter 3 | **Bonferroni-adjusted p < 0.01** |
| Contrast logic | tumour-exclusive LR lists (pairs **not detected** in background or healthy) → LUAD-vs-LUSC differential wiring (Fig. 2B–D); ICI/checkpoint co-expression differences |
| Multi-condition claim basis | Tumour-exclusive lists compared across histologies at cohort level |

## Factory lanes (this workspace)

| Aspect | Wave-6 group lane | Wave-8 patient lane | Wave-9 expansion (this lane) |
|---|---|---|---|
| Cohort | LUCA core (Salcher 2022), tumour_primary lung | same | same |
| Cells | ~3k per histology (cap, seed 23) | ~1k per donor, 3+3 donors | ~1k per donor, **8+8 donors** from eligible pool (47 LUAD / 22 LUSC with ≥800 cells & ≥6 types) |
| Conditions | tumour only, LUAD vs LUSC | tumour only, per donor | tumour only, per donor |
| Engine | LIANA `rank_aggregate` consensus (multi-method) | same | same |
| Pair statistics | LIANA consensus ranks (`magnitude_rank`, `specificity_rank`) | same + donor recurrence fractions | same + donor recurrence fractions |
| Contrast logic | panel row counts + Jaccard | panel recurrence across 3+3 donors | panel recurrence across 8+8 donors |
| Exclusivity logic | none (no background/healthy arm) | none | none |

## Consequence for claims

| Claim | Status | Why not higher |
|---|---|---|
| NC2024-F2-02 | `partial` | Method lane on external cohort, tumour-only; no paper exclusivity design, no figure parity |
| NC2024-F2-03 | `partial` | VEGF/EGFR panel contrast exists but is subsample + descriptive, not the paper's Bonferroni-filtered cohort statistics |
| NC2024-F2-04 | `partial` (strengthened by Wave-8/9 donor recurrence) | Recurrence across 8+8 donors is descriptive; paper's claim rests on tumour-vs-background exclusivity + multi-condition cohort design we have not replicated |

## What would be required to move toward `direct`

1. **Background/healthy arm**: run the same engine on LUCA `normal_adjacent` +
   healthy donors of matched histology, and compute tumour-exclusivity
   analogous to the paper (pairs in tumour not detected in background/healthy).
2. **Engine parity run**: one bounded CellPhoneDB execution on a LUCA
   stratified subsample to quantify LIANA-consensus vs CellPhoneDB panel
   overlap (engine sensitivity), so the engine difference is measured rather
   than asserted.
3. **Paper-cohort lane**: E-MTAB-13526 (16 GB raw, deleted locally per
   retention policy 2026-05-25; re-fetchable) with the paper's own filters
   (≥30% expression, mean log1p > 1.0, Bonferroni p < 0.01) on a stratified
   50% subsample — the only route to paper-design parity.
4. Multi-condition statistical model across donors (e.g. pseudo-bulk
   interaction scores with donor as unit), not just recurrence fractions.

Until at least (1) and (2) land, all F2 claims remain `partial` and any
tumour-exclusivity or figure-parity language is **out of bounds**.
