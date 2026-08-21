---
name: plot-qc-review
description: Use after a pipeline module emits QC/diagnostic figures, before accepting run artifacts, or whenever the user asks to check plots. Vision-based rubric review of single-cell figures (UMAP, QC violins, marker dotplots/heatmaps, elbows, distribution comparisons) with a PASS/NOTES/FAIL verdict per plot.
---

# plot-qc-review — vision QC gate for singlecell_factory figures

Use the model's image-reading ability (`ReadMediaFile`) to actually LOOK at every emitted
figure before accepting it. Text logs cannot catch a broken embedding, a batch-split UMAP,
or an empty panel — this gate exists for exactly that.

## Procedure

1. Collect the figures: the run's figure/ROI output dirs (see the run contract /
   `ops/briefs`), or the paths the user gave. If more than ~6 figures, review in batches
   and aggregate verdicts at the end.
2. For EACH figure: `ReadMediaFile` it. If the returned image is downsampled and any
   text (axis labels, legends, cluster ids) is unreadable, re-read with
   `full_resolution: true` or a `region` crop. **Never pass a figure you could not read.**
3. Apply the rubric below. Judge only what is visible — do not invent biology.
4. Emit the verdict contract (below). Any FAIL means the artifacts are NOT accepted;
   route the finding to `execute-and-recover-pipeline` with the concrete fix.

## Rubric

### Universal technical checks (every figure)
- Empty panels, missing axes/labels/legends, clipped text, legend overflow.
- Rendering adequacy for the cell count: a 40k/900k-cell embedding rendered as a few
  hundred visible points indicates a downsampling or alpha bug, not biology.
- Color scales: continuous scales actually varying; categorical palettes not collapsing
  distinct groups into one color.

### Per-type biological red flags
- **UMAP/t-SNE embeddings**: one structureless blob (integration/embedding failure);
  clusters separated strictly by batch/library (uncorrected batch effect); long streaks
  or grid stripes (artifact); single cells isolated in huge voids (outlier handling bug).
- **QC violins/histograms (n_genes, n_counts, percent_mt)**: degenerate (empty or
  single-spike) distributions; filtering thresholds drawn far from any density mass;
  percent_mt near zero everywhere for fresh tissue (suspicious, likely wrong MT prefix).
- **Marker dotplots/heatmaps**: markers expressed everywhere (no specificity — wrong
  genes or normalized-away signal); labels swapped vs. known biology; a whole row/column
  constant (dead gene or dead cluster).
- **Elbow/scree**: no elbow in the plotted range — chosen PC count is unjustified by
  the figure itself.
- **Doublet/ambient-RNA score distributions**: perfectly uniform or single-value scores
  (scorer silently failed).
- **Comparison panels (metric distributions, clustering comparison)**: identical
  distributions across methods (comparison not actually computed); axis ranges hiding
  the differences being claimed.

### Degenerate spikes: determinism vs silent scorer failure

A histogram collapsed to a single spike at a perfect value (e.g. F1 = 1.0 in all
repeats) sits between two very different truths: a deterministic pipeline legitimately
reproducing itself, or a scorer that silently stopped varying (constant broadcast,
cached metric, repeats not actually rerun). The plot alone cannot separate them — so
do not FAIL from the plot alone, and do not PASS either. This is exactly the
verify-from-data case:

1. Verdict **NOTES**, and state the suspicion plainly.
2. Before the figure is cited as evidence, check the underlying per-repeat metrics
   table (CSV/JSON next to the figure): FAIL if the metric was constant-broadcast
   rather than recomputed per repeat; upgrade to PASS if the repeats genuinely
   recomputed to the same value.
3. Cosmetic tells worth recording: an axis extending past the metric's possible
   maximum (e.g. F1 axis drawn to 1.4) is auto-scaling, not biology — note it, don't
   fail it.

## Verdict contract

Per figure:
```
<file>: PASS | NOTES | FAIL — <what you see, with panel/region evidence> — <fix if FAIL>
```

Overall:
- Any FAIL → overall FAIL: list the failing evidence, name the module/params most
  likely responsible, propose the concrete re-run or re-plot fix. Do NOT accept the run.
- Only NOTES → acceptable with caveats recorded in the run brief.
- All PASS → say so, with one line of evidence per figure (not just "looks fine").
