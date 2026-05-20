# 2026-05-19 human conclusion log for clustering/final-run strategy

## What changed

Added a human-readable conclusion log for the NC2024/Cell controller-validation discussion. The log captures:

- why the current NC `scale_mode=massive -> CSS` run is evidence-only and not final scientific clustering;
- why the final NC path should move to sparse-exact or validated hybrid clustering;
- why Cell/Trevino low ARI is likely a parity-boundary/information gap rather than a simple pipeline failure;
- why GPU DE fallback is a RAPIDS version/API issue, separate from GPU clustering;
- the approved next plan: preserve current NC DE evidence, produce Cell parity report, benchmark NC clustering routes, use sparse-exact/hybrid for final real-data claims, and validate GPU DE in a separate RAPIDS environment.

## Human conclusion artifacts

- NC ledger:
  `/home/zerlinshen/projects/nc-reproduction/ledger/human_review/2026-05-19-clustering-final-strategy-human-conclusion.md`
- NC run evidence:
  `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-19T1010Z-82f1964/evidence/human_conclusion_clustering_final_strategy_20260519.md`
- Cell/Trevino ledger:
  `/home/zerlinshen/projects/wave5-trevino/ledger/human_review/2026-05-19-clustering-final-strategy-human-conclusion.md`
- Cell/Trevino run evidence:
  `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-19T1005Z-82f1964/evidence/human_conclusions/human_conclusion_clustering_final_strategy_20260519.md`

## Classification

This is a human-facing conclusion and decision log. It does not promote the NC partial run to source-of-truth and does not change the Cell/Trevino source-of-truth. It records the strategy change for future final reruns.

