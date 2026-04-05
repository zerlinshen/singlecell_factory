# LUSC Optimization Validation Report

- Dataset: `LUSC_3K_real_run`
- Repeats: `10`
- Acceptance loss threshold: `0.1%`

## Performance
- Baseline mean time: `2.585s` (95% CI `2.540` - `2.630`)
- Optimized mean time: `2.516s` (95% CI `2.483` - `2.549`)
- Mean speedup: `2.63%`
- Runtime significance p-value: `0.024414`

## Accuracy & Reliability
- Cluster ARI mean: `1.0000`
- Pairwise Precision/Recall/F1 mean: `1.0000` / `1.0000` / `1.0000`
- Marker F1 mean: `1.0000`
- Max relative loss across key metrics: `0.0000%` (threshold `0.1%`)
- Accuracy significance p-value: `0.000000`

## Algorithm Shift (Wilcoxon -> t-test_overestim_var)
- Marker F1 (same clustering, different DE test): `0.8282`
- Interpretation: >=0.85 means algorithm substitution is close; <0.85 suggests keeping wilcoxon in production.

## Robustness
- `extreme_dropout_90pct`: **PASS**
- `low_cell_count_80`: **PASS**
- `nan_injection_sanitized`: **PASS**

## Rollback Decision
- Final decision: **ROLLBACK_REQUIRED**
- Reason: algorithm-shift marker consistency too low