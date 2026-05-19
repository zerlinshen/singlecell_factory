# Component E Executor Brief — GPU DE Environment Validation

**Status:** ⏸ DRAFT — pending explicit "build E env" approval
**Canonical plan:** `/home/zerlinshen/.omc/plans/nc-cell-clustering-final-strategy-plan.md` (APPROVED 2026-05-19)
**Plan section:** Component E (GPU DE environment validation)
**Hard upstream dependencies:** none
**Hard non-mutation invariant:** **do NOT modify the existing `sc_gpu_stable` environment** under any circumstances. F-3 + Principle 6 enforce DE rigor at the code level; this brief only verifies whether GPU DE is admissible on any final-claim path.
**Created:** 2026-05-19

## 1. Goal

Establish whether GPU DE is admissible on any final-claim path. The current `sc_gpu_stable` env ships `rapids-singlecell 0.13.4` which does not expose `rapids_singlecell.tl.rank_genes_groups`; DE silently falls back to CPU. F-3 already banned the silent CPU Welch fallback, so any GPU-DE workflow today either uses `sc.get.rank_genes_groups_df` (if available) or raises. Component E builds a separate test environment to evaluate whether a newer `rapids-singlecell` version closes this gap with sufficient CPU parity to be admitted into final-claim runs.

## 2. Target Environment

- **Name:** `sc_gpu_de_test` (clearly distinct from `sc_gpu_stable`; the name itself prevents `conda activate sc_gpu_stable` muscle memory).
- **Manifest:** `singlecell_factory/ops/envs/sc_gpu_de_test.yaml` — full conda lock with versions pinned. The pinned `rapids-singlecell` version must expose `tl.rank_genes_groups`. As of the audit, the target version is whatever the next rapids-singlecell release ≥ 0.14 ships; the exact version is selected by the executor session.
- **CUDA / RAPIDS base:** matching the `sc_gpu_stable` toolchain (same CUDA major version, same NumPy/scipy major versions) so the parity tests below isolate the `rapids_singlecell.tl.rank_genes_groups` change rather than infrastructure drift.

## 3. Non-Mutation Invariant (G-E1)

- **Forbidden:** any change to `sc_gpu_stable` (conda env path, packages, environment hooks, shell scripts that activate it). The verifier check below MUST pass before AND after the new env is built.
- **Verifier check:** record SHA256 of `conda env export -n sc_gpu_stable --no-builds` before and after the new env build. Equal hash = no mutation.
- **Rollback procedure:** if `sc_gpu_stable` mutation is detected at any point, the new env is destroyed immediately, the root cause is documented in `singlecell_factory/ops/env_validation/rollback_<timestamp>.md`, and the next attempt starts from a clean conda clone.

## 4. G-E2 API-Presence Gate

After the new env is built:

1. Activate `sc_gpu_de_test`.
2. Run a smoke check:
   ```python
   import rapids_singlecell as rsc
   assert callable(rsc.tl.rank_genes_groups), "G-E2 fail: API absent"
   ```
3. Run a toy DE on a 100-cell synthetic AnnData with 2 clusters and confirm it returns a markers DataFrame without errors.

**Pass criterion:** smoke + toy DE both succeed; the env passes G-E2. **Fail:** GPU DE remains inadmissible — env is destroyed, finding logged.

## 5. G-E3 Parity Gate (all three tiers must pass)

CPU vs GPU DE concordance must be reported on a defensible metric, at **all three** scale tiers independently. A partial pass (e.g., small + medium pass, full fails) does NOT admit GPU DE — it logs as "GPU DE bounded by scale" and surfaces a follow-up investigation but does not unlock D's GPU DE branch.

| Tier | Cell count | Source | Concordance metric (proposed; executor session can refine) |
|---|---|---|---|
| Small | ~5,000 | NC tumor 5k slice or matched synthetic | top-200 markers per cluster: Jaccard overlap ≥ 0.95 AND Spearman rank correlation ≥ 0.95 |
| Medium | ~100,000 | NC tumor 100k slice | top-200 markers per cluster: Jaccard overlap ≥ 0.92 AND Spearman ≥ 0.92 |
| Full | NC tumor full cohort (877k) | the canonical input | top-200 markers per cluster: Jaccard overlap ≥ 0.90 AND Spearman ≥ 0.90 |

Both CPU and GPU runs at each tier use **Wilcoxon** (per F-3 / Principle 6). Welch is not benchmarked; it is banned.

Executor session may refine the exact thresholds with citation (e.g., references to published rapids-singlecell vs scanpy DE parity studies); thresholds must be declared BEFORE running parity tests, never after.

## 6. Decision Triggers

| Outcome | Action |
|---|---|
| G-E1 violated (`sc_gpu_stable` mutated) | Immediate rollback; root-cause note in ledger; retry from clean clone. |
| G-E2 fails (API absent in target version) | GPU DE inadmissible; D uses CPU Wilcoxon. Try newer `rapids-singlecell` version or wait. |
| G-E2 passes, G-E3 fails any tier | GPU DE inadmissible for final claims. Result logged as "GPU DE bounded by scale [TIER]" or "GPU DE concordance gap"; D uses CPU Wilcoxon. |
| All three G-E3 tiers pass independently | GPU DE admitted into D's available DE engines. D's methodology subsection may select GPU DE. |

## 7. Output Directory Layout

```
singlecell_factory/ops/envs/
└── sc_gpu_de_test.yaml                     (conda lock, F-3-compliant)

singlecell_factory/ops/env_validation/
├── G-E1-isolation-evidence.md              (before/after sc_gpu_stable SHA256)
├── G-E2-api-presence-evidence.md           (smoke + toy run output)
├── G-E3-parity-small-tier.md
├── G-E3-parity-medium-tier.md
├── G-E3-parity-full-tier.md
└── gpu_de_parity_report.md                 (final admission status; one of {admitted, bounded, inadmissible})
```

The final `gpu_de_parity_report.md` updates the canonical plan's ADR follow-up #11 (`gpu_de_admission_status` one-liner).

## 8. Launch Pre-conditions

Before this brief can be acted on:
- [ ] Explicit "build E env" approval from human.
- [ ] Confirmation that the chosen `rapids-singlecell` target version exposes `tl.rank_genes_groups` (verified via release notes / source inspection before build).
- [ ] Disk and conda-cache budget verified on the host that will run the parity tests.

**Building the env IS a real-system action** (writes to conda envs dir) — per `user_collab_style.md` real-run gates, this requires explicit go.

## 9. Owner / Agent Routing

- Env build + manifest: `oh-my-claudecode:executor`.
- API-presence smoke + parity methodology: `oh-my-claudecode:scientist`.
- G-E1 isolation verification + rollback: `oh-my-claudecode:verifier`.
- Final admission report: `oh-my-claudecode:critic` reviews G-E3 thresholds for honesty.
