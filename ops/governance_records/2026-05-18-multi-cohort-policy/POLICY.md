# Multi-Cohort Processing Policy

Control-plane governance record only; not a scientific output.

- Timestamp UTC: `2026-05-18T09:14:26Z`
- Decision: per-sample filtering is allowed and encouraged, but article-level
  annotation must happen on a cohort-level merged/sketch/metacell atlas with
  explicit batch-aware integration and label transfer back to all cells.
- Rationale: single-patient annotation loses cross-patient comparability and
  can confound donor/batch effects with cell-type or state calls.

## Default Execution Pattern

1. per-sample/per-patient QC and doublet filtering
2. merge clean cells or construct balanced sketch/metacell atlas
3. cohort-level normalization/HVG/PCA/integration with batch keys
4. cohort-level clustering and annotation
5. label transfer to all cells
6. patient-level composition, pseudobulk, DE, figures, and conclusions

## Scope

This is a policy/documentation change. It does not change runtime defaults or
rerun data by itself.
