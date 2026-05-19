# Ambient RNA Correction Policy

**Status:** ✅ DECIDED 2026-05-19 — option (a) **Integrate ambient correction as a permanent factory module**
**Canonical plan:** `/home/zerlinshen/.omc/plans/nc-cell-clustering-final-strategy-plan.md` (APPROVED v4, 2026-05-19)
**Decision gate:** Plan Component G, G-G1
**Downstream gate:** G-G2 — feeds Component D's G-D3 (iii) claim-coverage subsection
**Author:** zerlinshen + Claude (Opus 4.7)

## Decision

**Option (a) chosen: integrate ambient RNA correction as a permanent factory module.**

Rationale: `singlecell_factory` is a **general-purpose single-cell production tool** intended to serve all future single-cell analysis pipelines — not a one-shot reproduction harness for NC2024 NSCLC + Cell/Trevino specifically. Tumor scRNA-seq is the dominant intended use case and is the modality most vulnerable to ambient RNA contamination (free-floating mRNA from lysed cells produces malignant↔immune marker bleed, particularly in NSCLC and other epithelial cancers). A permanent module amortizes the engineering cost across every future tumor project, eliminates the per-claim sensitivity-analysis bookkeeping that option (d) would have required, and is the only option compatible with Principle 4 ("Parity before claims") at industrial scale.

Options (b), (c), (d) were rejected:
- **(b) require pre-corrected input** — rejected because NC and Cell/Trevino public data are not pre-corrected, and asking each future project to run SoupX externally without factory audit trail violates the project's reproducibility contract (factory CLAUDE.md: "Always prioritize reproducibility, deterministic outputs, and stable pipeline contracts").
- **(c) exclude ambient-driven claims** — rejected because it would systematically narrow every future project's claim envelope on the most scientifically interesting axis (malignant↔immune crosstalk), not just NC.
- **(d) sensitivity-analysis-per-claim** — rejected because the bookkeeping cost is paid per claim per project forever, whereas (a)'s engineering cost is paid once.

## Scope of this policy

This policy locks the **what** (factory must have a permanent ambient correction module) but **defers the how** (tool choice, integration location, dependency surface) to a follow-up plan.

## Tool selection — DEFERRED to follow-up plan

Three candidate tools, each with tradeoffs to evaluate in the follow-up plan:

| Tool | Language | Method | Inputs needed | Compatibility |
|---|---|---|---|---|
| **SoupX** | R | Empirical ambient profile from cell-free droplets | `raw_feature_bc_matrix` (the unfiltered Cell Ranger output) | Existing `r_multiomics_factory` R bridge — already in the factory ecosystem. Most peer-reviewed tumor scRNA-seq usage. |
| **DecontX** (celda) | R | EM on cluster structure to estimate per-cell contamination fraction | Filtered matrix only (no raw required) | Same R bridge. Doesn't require raw matrix — survives data sources that lost it. |
| **cellbender remove-background** | Python | Variational autoencoder over droplets | `raw_feature_bc_matrix` | Pure Python, GPU-friendly, fits the rapids-singlecell stack. Newest method; state-of-the-art per recent benchmarks. |

Tool selection criteria (the follow-up plan will weigh these):
1. **Compatibility with NC + Cell/Trevino actual input availability** — does `raw_feature_bc_matrix` exist for these datasets? If not, only DecontX is usable for them. (Check at follow-up plan time.)
2. **Compatibility with all future expected datasets** — tumor scRNA-seq generally has the raw matrix; rarer single-cell modalities may not.
3. **Engineering cost** — Python module is cheapest (no new bridge); R modules reuse existing R bridge but add cross-language complexity.
4. **GPU compatibility** — cellbender is GPU-native and aligns with the factory's rapids-singlecell direction.
5. **Peer-review weight** — SoupX has the longest production track record.

## Effect on the canonical plan's Component D

Per the plan's G option-(a) decision-trigger language: **"Open a follow-up plan for ambient module integration (out of scope for this plan). Block D's commit on integration milestone."**

This means:
- **Component D's methodology-binding subsections** are now blocked on the ambient module integration milestone, in addition to G-C3 (C winner) and G-F3 (Wilcoxon-mandatory, ✅ DONE).
- **Component D's non-methodology subsections** (already drafted in the skeleton at `projects/nc-reproduction/ledger/designs/nc-final-science-rerun-design.md`) remain unblocked and can be finalized.
- **Component D's G-D3 (iii) "claim families explicitly cannot support" subsection** becomes minimal once the ambient module integrates — most ambient-driven claims become supportable. Pre-integration, this subsection should enumerate "all ambient-driven claims" as currently-blocked-pending-module.
- **NC + Cell/Trevino reproductions** wait for the ambient module before final claims are issued. Controller-validation runs (already complete) are NOT retroactively invalidated; the new module gates final-claim runs only.

## Interim policy (between this decision and module integration)

Per Plan Principle 1 ("Plan, do not execute"), this policy decision does NOT itself launch any rerun. While the follow-up plan is being designed and the module is being built:
- No final-claim NC or Cell rerun launches.
- Existing 2026-05-19 controller-validation evidence stays as-is, labeled "controller-validation" not "final-science".
- Documentation in this repo and in `nc-reproduction` + `wave5-trevino` ledgers refers to claims as **provisional pending ambient correction**.

## Required follow-up plan

A new strategic plan needs to be drafted at `~/.omc/plans/factory-ambient-correction-module-integration.md`. It must cover:
- Tool selection (per criteria above)
- Module surface: name, CLI flag, config schema, dependency surface (env updates if Python; R bridge usage if R)
- Where the stage sits in the pipeline (between QC and doublet detection? between doublet detection and clustering?)
- Provenance: how the ambient profile and per-cell contamination fractions get recorded in `adata.uns` + manifest
- Tests: parity / regression / contract gates
- Compatibility with the existing `scale_mode` presets (no silent activation; explicit opt-in or default-on with disclosure)
- Effect on `r_multiomics_factory` if R tool chosen (cross-factory contract update)
- Audit gates and approval signoff

## ADR cross-reference

This decision updates Plan v4's ADR (Decision section) and Component G's status. Plan should be updated to reflect:
- G-G1: ✅ DECIDED (this doc)
- G-G2: requires the follow-up plan's integration milestone to land before D's claim-coverage subsection finalizes

## Changelog

- **2026-05-19:** Decision recorded. Option (a) chosen after user clarified that the factory is a general-purpose single-cell production tool, not a one-shot reproduction harness. Tool selection deferred to follow-up plan.
