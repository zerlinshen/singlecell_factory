# STATUS: one-off, promote-on-reuse

---
prompt_id: PLOT-1
plot: UMAP-RNA
plan_revision: v4.2
authored_at: 2026-05-16
---

# Visual-Verdict Prompt — PLOT-1 (UMAP-RNA Cluster Separability)

## Target Artifact

The RNA UMAP plot produced by the Wave-5 v4.2 canonical run (path resolved
at verdict time from the canonical run ledger; do **not** hard-code).

## Operator Instructions

You are performing a **biological-substance** visual verdict on the RNA
UMAP. This is *not* a rendering-quality check — it is a check on whether
the embedding recovers PCW21 cortical-lineage structure.

Inspect the plot and answer each numbered checklist item with one sentence
of evidence drawn from what you actually see on the rendered figure. Do
not speculate beyond what the figure shows.

## Biological-Substance Checklist

1. **Lineage cluster count.** Does the UMAP show **≥ 5 visually-separable
   clusters** consistent with PCW21 cortical lineages (radial glia / IPC /
   excitatory neuron / inhibitory neuron / OPC / microglia)? Note whether
   the clusters are distinct point-clouds or only density modes.

2. **Radial-glia marker concentration.** Are the radial-glia
   marker-expression overlays (**PAX6, VIM, SOX2**) concentrated in **one
   or two contiguous clusters** rather than scattered across the entire
   embedding? Identify the cluster(s) that show the strongest expression.

3. **OPC-lineage separation.** Are the OPC markers (**OLIG2, SOX10**)
   concentrated in a **distinct cluster** that does **not** overlap with
   neuronal-lineage clusters (ExN / IN)? Note overlap, if any.

4. **Artifact freedom.** Is the layout free of obvious artifacts —
   specifically: single-point clusters, ring artifacts, density-only blobs
   with no internal lineage structure, or a single dominant blob
   containing all cells?

## Verdict Emission

Emit a single verdict value from `{pass, partial, fail}` plus one sentence
of evidence per checklist item.

- `pass` — all four checklist items show clear biological-substance
  evidence.
- `partial` — at least one item is ambiguous or marginal but no item is
  outright failing.
- `fail` — any item shows a clear failure (e.g., no separable clusters,
  RG markers scattered everywhere, OPC overlapping ExN, obvious
  artifact).

Write the verdict and per-item evidence into:

    .omc/research/wave5/visual_verdicts/PLOT-1_verdict.md

The verdict file is read by the biology validation gate (AC-VAL-4 visual
verdicts). Do not edit this prompt during the verdict pass.
