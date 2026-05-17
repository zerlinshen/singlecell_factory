# STATUS: one-off, promote-on-reuse

---
prompt_id: PLOT-3
plot: peak-gene-linkage-volcano
plan_revision: v4.2
authored_at: 2026-05-16
---

# Visual-Verdict Prompt — PLOT-3 (Peak-Gene Linkage Volcano)

## Target Artifact

The peak-gene linkage volcano plot produced by the Wave-5 v4.2 canonical
run. Axes: **x = |Pearson|**, **y = -log10(FDR)**, with the **top-1000
peak-gene links** highlighted. Path resolved at verdict time from the
canonical run ledger; do **not** hard-code.

## Operator Instructions

You are performing a **biological-substance** visual verdict on the
peak-gene linkage volcano. The check is whether the v3-method peak-gene
linkage output exhibits the expected confidence concentration and panel-
member enrichment among the top-ranked links.

Inspect the plot and answer each numbered checklist item with one sentence
of evidence drawn from what you actually see on the rendered figure. Do
not speculate beyond what the figure shows.

## Biological-Substance Checklist

1. **|Pearson| concentration of top-1000.** Is the **|Pearson|
   distribution of the top-1000 highlighted links in the 0.8–1.0 range**
   (high-confidence linkages) consistent with the v3 method? Note the
   approximate min/max of the highlighted points along the x-axis.

2. **Panel-member enrichment among top-50.** Are **panel-member genes
   among the top-50 visually highlighted** (e.g., **OLIG2 and SOX10**
   expected near the top by |Pearson|, with additional panel members such
   as PAX6 / SOX2 / NEUROD2 visible in the top region)? Identify any
   panel-member gene labels you can read.

3. **FDR separation top-1000 vs background.** Does the **FDR distribution
   (y-axis)** show a clear **separation between the top-1000 highlighted
   points and the background point cloud** — i.e., are the top-1000 at
   substantially higher -log10(FDR) than the bulk?

4. **Artifact freedom.** Is the plot free of obvious technical artifacts —
   specifically: vertical bands at |Pearson| = 1.0 driven by sparse
   detection, horizontal banding from FDR ties dominating the top, or a
   highlight set that does not visually correspond to a peak in the
   |Pearson|/FDR joint distribution?

## Verdict Emission

Emit a single verdict value from `{pass, partial, fail}` plus one sentence
of evidence per checklist item.

- `pass` — all four checklist items show clear biological-substance
  evidence.
- `partial` — at least one item is ambiguous or marginal but no item is
  outright failing.
- `fail` — any item shows a clear failure (e.g., top-1000 |Pearson| is
  centered well below 0.8, no panel members in top-50, no FDR separation,
  obvious sparse-detection vertical band).

Write the verdict and per-item evidence into:

    .omc/research/wave5/visual_verdicts/PLOT-3_verdict.md

The verdict file is read by the biology validation gate (AC-VAL-4 visual
verdicts). Do not edit this prompt during the verdict pass.
