# STATUS: one-off, promote-on-reuse

---
prompt_id: PLOT-2
plot: UMAP-ATAC
plan_revision: v4.2
authored_at: 2026-05-16
---

# Visual-Verdict Prompt — PLOT-2 (UMAP-ATAC Cluster Separability)

## Target Artifact

The ATAC UMAP plot produced by the Wave-5 v4.2 canonical run (path resolved
at verdict time from the canonical run ledger; do **not** hard-code).

## Operator Instructions

You are performing a **biological-substance** visual verdict on the ATAC
UMAP. The check is whether the chromatin-accessibility embedding recovers
the same major PCW21 cortical lineages as the RNA embedding (PLOT-1) and
whether joint RNA/ATAC structure is preserved.

Inspect the plot and answer each numbered checklist item with one sentence
of evidence drawn from what you actually see on the rendered figure. Do
not speculate beyond what the figure shows.

## Biological-Substance Checklist

1. **Lineage recovery vs RNA.** Does the ATAC UMAP recover the **same
   major lineages** as the RNA UMAP — i.e., **≥ 4 visually-separable
   clusters** corresponding to the principal PCW21 lineages (radial glia
   / IPC+ExN / IN / OPC, and microglia where present)?

2. **Glial vs neuronal peak-set differences.** Are peak-set differences
   between **glial vs neuronal clusters** visually evident — e.g., via
   cell-type colorization, marker-peak overlays, or accessibility-score
   gradients separating OPC/glia from ExN/IN?

3. **Joint structure with RNA.** Is there evidence of **joint structure**:
   does the cluster topology in ATAC roughly align with the RNA
   embedding's topology (in particular, is the
   **RG → IPC → ExN developmental axis** preserved)?

4. **Artifact freedom.** Is the layout free of obvious technical artifacts
   — specifically: single-point clusters, ring artifacts,
   doublet-driven smears, or a single dominant blob containing all
   cells?

## Verdict Emission

Emit a single verdict value from `{pass, partial, fail}` plus one sentence
of evidence per checklist item.

- `pass` — all four checklist items show clear biological-substance
  evidence.
- `partial` — at least one item is ambiguous or marginal but no item is
  outright failing.
- `fail` — any item shows a clear failure (e.g., ATAC does not recover
  major lineages, no glial/neuronal peak-set separation, RG→IPC→ExN axis
  absent, obvious artifact).

Write the verdict and per-item evidence into:

    .omc/research/wave5/visual_verdicts/PLOT-2_verdict.md

The verdict file is read by the biology validation gate (AC-VAL-4 visual
verdicts). Do not edit this prompt during the verdict pass.
