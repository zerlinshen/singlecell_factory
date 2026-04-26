---
name: bioinformatics-autosteer
description: Low-friction workflow steering for bioinformatics projects. Use when working on single-cell, spatial transcriptomics, WGS, Hi-C/3C, proteomics, large dataset execution, remote analysis orchestration, or local artifact handoff, especially when deciding which existing skills to apply first.
---

Use this skill to bias toward the right helper skills early instead of treating every task as generic coding.

1. Prefer these existing skills automatically when relevant.
- `before-every-run` first for any remote execution or resume on `ubuntu-tail`.
- `singlecell-remote-workflow` for Ubuntu-to-Mac analysis handoff, remote run status, and local artifact bridges.
- `execute-and-recover-pipeline` for failed jobs, checkpoint/resume, and evidence-based reruns.
- `singlecell-factory-module-delivery` when the task is to introduce or evolve a module in the remote `singlecell_factory` pipeline.
- `develop-and-integrate-module` for pipeline/module changes in `singlecell_factory`.
- `review-and-validate-quality` for final risk review, scientific output checks, and doc synchronization.
- Follow [before-every-run/references/shenxin-routing-contract.md](/Users/zerlinshen/.codex/skills/before-every-run/references/shenxin-routing-contract.md) for lane order.

2. For large datasets, assume staging first.
- Start with lighter module sets.
- Prefer clustering-first runs for very large objects.
- Treat swap as a safety buffer, not a substitute for RAM.

3. For multi-omics projects, preserve expensive assets.
- Protect raw data and references by default.
- Treat outputs and caches as reproducible unless explicitly marked important.

4. Keep the workflow practical.
- Remote side for heavy compute.
- Local side for organization, interpretation, comparison, and handoff artifacts.
- When in doubt, optimize for successful staged execution over ambitious all-at-once runs.
