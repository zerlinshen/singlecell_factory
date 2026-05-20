# Ambient RNA Correction Policy

**Status:**
- ✅ Phase 1 — WHAT — DECIDED 2026-05-19: integrate ambient correction as a permanent factory module
- ✅ Phase 2 — HOW — DECIDED 2026-05-20: **tool = DecontX**, **activation = conditional (QC-triggered, not default-on)**, **install scope = lightweight (`bioconductor-celda` 1.26.0 in `r_multiomics` env)**

**Canonical plan:** `/home/zerlinshen/.omc/plans/nc-cell-clustering-final-strategy-plan.md` (APPROVED v4, 2026-05-19; nc-reproduction project itself ABORTED 2026-05-20 — policy survives because it is factory-level, not project-level)
**Decision gate:** Plan Component G, G-G1 (Phase 1) + this doc Phase 2 (Phase 2)
**Downstream effect:** Component D / NC reproduction is moot (project aborted); policy now applies to the next single-cell project (Nature Genetics 2025 pan-cancer 3D genome — DOI 10.1038/s41588-025-02188-0)
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

---

# Phase 2 (HOW) — DECIDED 2026-05-20

## Tool selection: DecontX (locked)

**Chosen: DecontX (Bioconductor `celda` package, R-based).**

Rationale per Phase 1 criteria:

| Criterion | DecontX evaluation |
|---|---|
| 1. Compatibility with `raw_feature_bc_matrix` availability | **DecontX does NOT require raw matrix** — runs on the filtered cellranger matrix alone. Survives data sources that lost the raw matrix. SoupX and cellbender both require raw, which makes them brittle for re-analysis pipelines. |
| 2. Compatibility with all future datasets | Filtered-matrix-only requirement = best portability. Tumor scRNA, spatial-derived single-cell, and re-analyzed public datasets all retain filtered matrices but may not retain raw. |
| 3. Engineering cost | Reuses the existing R bridge pattern (`_r_scripts/run_wnn.R` → subprocess via Python). No new env, no new bridge type. `bioconductor-celda` 1.26.0 already installed into the `r_multiomics` env (R 4.5.3) at policy-decision time — verified `decontX` function loadable. |
| 4. GPU compatibility | DecontX is CPU-only. Acceptable because (a) it runs once per sample at QC stage, not on hot path; (b) decontamination is not a tight loop, run-time is minutes at <50k cells; (c) avoids the cellbender heavyweight CUDA dependency. |
| 5. Peer-review weight | Yang S et al. 2020 *Genome Biology* (DOI 10.1186/s13059-020-1950-6) — peer-recognized; broadly cited in tumor scRNA literature. Comparable validation to SoupX; newer than SoupX in 10x v3 chemistry contexts. |

**Rejected alternatives:**
- **SoupX** — needs raw matrix; raw matrix availability is dataset-dependent; harder to enforce factory-wide.
- **cellbender remove-background** — heavyweight (CUDA + variational autoencoder; per-sample 30-60 min on GPU); kept on the table for a possible Phase 3 (advanced) lane if a future dataset has known severe contamination AND raw matrix.

## Activation policy: **CONDITIONAL, NOT DEFAULT-ON**

**Core principle:** Ambient correction is *not* a free lunch. DecontX models contamination as a per-cell-type mixture and can over-correct when:
- Cell-type clusters are poorly resolved (early QC stage),
- Sample is already clean (`pct_top_50 < 0.5`, low mt%, low predicted doublets),
- Cohort batch effects dominate over ambient effects.

Therefore the factory **MUST NOT** run DecontX unconditionally. It is a conditional QC remediation step, activated **per-sample** based on initial QC metrics.

### Trigger rules (factory default thresholds)

Activation decision is made **per-sample** at the QC stage, immediately after `qc.py` computes the standard metrics but **before** doublet detection. Sample triggers DecontX when **any** of the following conditions fire:

| Trigger ID | Condition | Default threshold | Rationale |
|---|---|---|---|
| **T1** `top50_high` | Median `pct_counts_in_top_50` per cell > threshold | `50.0` (interpreted as **percent**, matching scanpy's 0–100 convention) | High top-50 fraction implies cells dominated by few highly-expressed genes — characteristic ambient signature where free-floating housekeeping/lineage markers flood all cells. |
| **T2** `mt_excess` | Median `pct_counts_mt` per cell > threshold | `15.0` (interpreted as **percent**, matching scanpy's 0–100 convention) | Elevated mitochondrial fraction indicates cell lysis; lysed cells release mRNA into the droplet supernatant → ambient contamination. |
| ~~**T3** `doublet_excess`~~ | ~~Scrublet `predicted_doublet_rate` > expected × multiplier~~ | ~~multiplier = `2.0`~~ | **REMOVED 2026-05-20** — DAG ordering (`qc → ambient_correction → doublet_detection`) makes `obs["predicted_doublet"]` always absent when ambient triggers are evaluated. The original idea was sound (doublet over-prediction can reflect ambient signal forming pseudo-doublet expression profiles), but it is unreachable in the current pipeline shape. Doublet-driven ambient re-trigger is deferred to a potential future Phase 3 (re-invoke ambient after doublet detection). CLI flag `--ambient-trigger-doublet-mult` and config field `trigger_doublet_multiplier` were dropped at the same time. |
| **T4** `low_cell_count_correlation` | Spearman ρ between `total_counts` and `n_genes_by_counts` across cells < threshold | `0.85` | Healthy droplets show tight count↔gene correlation. Ambient-tainted samples show many low-count droplets with elevated gene counts (ambient genes bleed in). Inverse: low ρ → ambient suspect. |
| **T5** `cross_sample_variance` | Coefficient of variation of housekeeping gene expression across samples within same cohort > threshold | `0.50` | Healthy cohort housekeeping should be stable; ambient contamination is per-sample → high CV signals sample-specific contamination. (Cohort-level trigger, not per-sample.) |

**Trigger combination:** any single trigger among T1, T2, T4 fires → DecontX runs on that sample. T5 fires → DecontX runs on **all samples in the cohort**. (T3 was removed 2026-05-20 — see strikethrough row above.)

**No trigger fires** → DecontX is **skipped**. The sample proceeds directly to doublet detection. The decision is recorded in `adata.uns["ambient_correction"]` and `module_status.csv` as `skipped_no_trigger` with the QC metrics that justified skipping.

### Threshold override (per-project)

Defaults above are factory-level. Projects can override via:
- CLI flag: `--ambient-trigger-top50 0.40` (project tuning)
- Config: `AmbientConfig.trigger_thresholds = {...}` (deeper project pinning)
- Env var: `SC_AMBIENT_TRIGGERS_DISABLE=1` (force-skip; e.g., paper-faithful reproduction where original paper did NOT run ambient correction — this is the De Zuani 2024 mode and would have been the right setting if that project had survived)

**Recorded in manifest:** all threshold values used, regardless of source (default / CLI / config / env), are written to `run_manifest.json` under `ambient.trigger_thresholds`.

### Re-QC after DecontX

If DecontX runs on a sample, the post-correction flow is:
1. Apply `decontXcounts` as the new `.X` (count matrix).
2. **Re-run QC** (`qc.py`) on the corrected counts to recompute `pct_counts_in_top_50`, `pct_counts_mt`, `total_counts`, `n_genes_by_counts`.
3. **Re-run doublet detection** (`doublet_detection.py`) on the corrected counts (doublet predictions on raw vs. corrected can differ).
4. Continue normal pipeline (clustering, batch correction, ...).

**Audit trail:** both pre- and post-correction QC metrics are stored in `adata.uns["ambient_correction"]["qc_pre"]` and `["qc_post"]`. The pipeline manifest also records DecontX's per-cell contamination fraction summary (median, q25, q75, max) under `ambient.contamination_summary`.

### Provenance fields written to `adata.uns["ambient_correction"]`

```python
{
  "engine": "decontx",                          # or "skipped"
  "tool_version": "celda 1.26.0",
  "decision": "triggered" | "skipped_no_trigger",
  "triggers_fired": ["T1_top50_high", "T4_low_count_correlation"],  # empty when skipped
  "trigger_thresholds": {...},                  # the resolved threshold dict
  "qc_pre": {"top50_median": ..., "mt_median": ..., ...},
  "qc_post": {"top50_median": ..., "mt_median": ..., ...},   # null when skipped
  "contamination_summary": {"median": ..., "q25": ..., "q75": ..., "max": ...},  # null when skipped
  "decontx_runtime_seconds": 123.4,             # null when skipped
}
```

## Pipeline placement

DecontX runs as a **new optional module** named `ambient_correction`, positioned in the DAG between `qc` and `doublet_detection`. Stage order:

```
cellranger → qc → ambient_correction (conditional) → doublet_detection → ... (existing chain)
```

The module is `optional` (not in `MANDATORY_STAGES`) because:
- The conditional policy itself may decide to skip (no-trigger samples bypass).
- Paper-faithful reproductions may want to opt out entirely (`SC_AMBIENT_TRIGGERS_DISABLE=1`).
- The factory should not force ambient correction on lightweight tests or sandbox runs.

But the module is **enabled by default** in production project runs — the conditional logic decides per-sample whether the *correction itself* runs, distinct from whether the *module is invoked*. The module is always invoked; it self-skips when no triggers fire.

## Implementation artifacts (this delivery)

- **R driver script:** `workflow/modular/modules/_r_scripts/decontx_run.R` ✅ DELIVERED 2026-05-20 — accepts h5ad input via zellkonverter, runs `celda::decontX`, writes decontaminated h5ad + JSON summary.
- **Python module:** `workflow/modular/modules/ambient_correction.py` ✅ DELIVERED 2026-05-20 — `AmbientCorrectionModule` class with 5-trigger evaluator, conda subprocess to R (`conda run -n r_multiomics Rscript ...`), full provenance to `adata.uns["ambient_correction"]`, cross-language debugging design (env isolation, stderr pass-through, keep-temp-on-failure, repro-cmd logging, dry-run mode, no-silent-skip on R missing).
- **Pipeline wiring:** `workflow/modular/pipeline.py` `_build_registry()` ✅ DONE 2026-05-20 — registered as `AmbientCorrectionModule()`. `module_catalog.py` `MODULE_SPECS` ✅ DONE 2026-05-20 — `ambient_correction.depends_on=("qc",)`, `doublet_detection.depends_on=("qc","ambient_correction")` (DAG forces ambient before doublet). Added to `MANDATORY_MODULES` so it always runs (self-skips internally when no trigger fires).
- **CLI flags:** `cli.py` ✅ DONE 2026-05-20 — 8 flags: `--ambient-disable-triggers`, `--ambient-dry-run-triggers-only`, `--ambient-trigger-top50/-mt/-count-corr/-cohort-cv`, `--ambient-decontx-max-iter`, `--ambient-batch-column`. `AmbientCorrectionConfig` lives in `config.py` and is constructed from these args at `parse_args()`-end. (Originally 9 flags; `--ambient-trigger-doublet-mult` was dropped 2026-05-20 with the T3 removal.)
- **Tests:** parity test for "no triggers → skipped"; contract test for "triggered → post-QC re-runs" — TODO for next pass, written after first NG2025 run exercises the trigger evaluator on real data.

## Out of scope (Phase 3 / future)

- cellbender lane as advanced opt-in for severe-contamination + raw-matrix-available cases.
- Empty-droplet ambient profile estimator (Lun et al. emptyDrops) as a separate sanity check.
- Cross-sample contamination diffusion modeling (currently each sample is independent in DecontX; cohort effects only via T5 trigger).

## Changelog

- **2026-05-19:** Phase 1 decision recorded. Option (a) chosen after user clarified that the factory is a general-purpose single-cell production tool, not a one-shot reproduction harness. Tool selection deferred to follow-up plan.
- **2026-05-20:** Phase 2 decision recorded. Tool = DecontX; activation = conditional (5-trigger rule, threshold override allowed); installed `bioconductor-celda` 1.26.0 in `r_multiomics` env; R + Python artifacts to be written; pipeline wiring DEFERRED pending user policy review. **Note:** the De Zuani 2024 NSCLC reproduction project (nc-reproduction) that originally motivated Phase 1 was ABORTED 2026-05-20 due to cell-calling divergence; the policy survives because it is factory-level and now serves the next single-cell project (NG 2025 pan-cancer 3D genome, DOI 10.1038/s41588-025-02188-0).
- **2026-05-20 (post-review):** Phase 2 hardening landed (Team B). Five behavioral fixes:
  1. **T3 removed** — DAG ordering makes `obs["predicted_doublet"]` unreachable at ambient-trigger time. The doublet-driven re-trigger idea is deferred to potential Phase 3. CLI flag `--ambient-trigger-doublet-mult` and `AmbientCorrectionConfig.trigger_doublet_multiplier` were dropped to match. Trigger ID set is now T1, T2, T4, T5.
  2. **Missing-QC fail-loud** — when both `pct_counts_in_top_50` AND `pct_counts_mt` are absent, the module now RAISES `RuntimeError` with `decision="qc_metrics_unavailable"` recorded in `adata.uns`, instead of silently returning `decision="skipped_no_trigger"`. Principle 9: a "no triggers fired" label for an "no triggers evaluable" run is a silent algorithmic compromise.
  3. **R script relocated** — `decontx_run.R` moved from `singlecell_factory/workflow/modular/modules/_r_scripts/` to `r_multiomics_factory/scripts/` (sibling of `run_wnn.R`). Python resolves the path via env var `R_MULTIOMICS_SCRIPT_DIR` if set, else via suite-root sibling discovery (`<suite>/r_multiomics_factory/scripts/decontx_run.R`).
  4. **Conda env declared** — `r_multiomics_factory/envs/r_multiomics.yml` now declares the minimal env (`r-base=4.5.3`, `bioconductor-celda=1.26.0`, `bioconductor-singlecellexperiment`, `bioconductor-zellkonverter`, `r-matrix`).
  5. **Subprocess hardening** — switched from `subprocess.run(timeout=...)` to `Popen` + `communicate(timeout=...)` + `os.killpg(SIGKILL)` with `start_new_session=True`, so an R grandchild wedged in `celda::decontX` is killed together with its conda parent. `shutil.which("conda")` preflight gives a clear error when conda is missing. `_raise_with_artifacts` annotated `NoReturn` to match its raise contract.
