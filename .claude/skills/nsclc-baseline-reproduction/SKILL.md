---
name: nsclc-baseline-reproduction
description: Use when reproducing or validating the current high-IF NSCLC baseline article. Produces a claim inventory, reproducible task list, and claim-validation matrix grounded in the stable 40k reference, the 40k massive/CSS benchmark, and the 900k scale-validity lane.
---

Use this skill for the current baseline reproduction workflow.

## Inputs it should rely on
- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run/LATEST.md`
- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/first reproduction - Nature Communications 2024 NSCLC single-cell/00_START_HERE_AGENT_HANDOFF.md`
- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/NC2024_SINGLE_CELL_CLAIM_MAPPING.md`
- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/NC2024_SINGLE_CELL_REPRODUCTION_STATUS.md`
- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/NC2024_SINGLE_CELL_GAP_ANALYSIS.md`
- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/NC2024_SINGLE_CELL_COMPARISON_REPORT.md`

## Operating rules

1. Use the 40k clean vs 40k CSS benchmark as the primary controlled comparator for overlapping biological modules.
2. Use the 900k CSS proxy lane only as scale-validity / downstream-coherence evidence, not as strict biological fidelity proof.
3. Judge support by artifact class, not by Python-versus-R identity.
4. Keep unsupported spatial and cross-modal claims explicit.
5. Route compute execution to existing workflow skills; do not turn this skill into a second execution engine.
6. Before requesting fresh remote evidence, invoke `before-every-run` first.

## Expected outputs
- updated claim inventory when needed
- updated reproducible task list
- updated claim-validation matrix
- concise summary of what is direct, partial, and unsupported now

## Default interpretation vocabulary
- `direct`
- `partial`
- `unsupported`
- `full_object`
- `bundle`
- `proxy`
- `40k_fidelity`
- `900k_survivability`
