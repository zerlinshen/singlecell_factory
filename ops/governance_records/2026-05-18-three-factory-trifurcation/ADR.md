# ADR: Three-Factory Trifurcation — Phase 1

- **Date:** 2026-05-18
- **Status:** Phase 1 complete; smoke gate green; Phase 2 + Phase 3 + C1 codex review + C2 reproduction SOP upgrade outstanding
- **Plan:** `/home/zerlinshen/.omc/plans/factory-trifurcation-consensus-plan-2026-05-18.md` (v2)
- **Spec:** `/home/zerlinshen/.omc/specs/deep-interview-factory-trifurcation-2026-05-18.md`

## Decision

Split `/home/zerlinshen/` into a three-factory sibling architecture:

| Repo | Role | Language |
|---|---|---|
| `singlecell_factory` | Python compute pipeline | Python |
| `r_multiomics_factory` | R analysis (renamed from `multiomics_r_factory`) | R |
| `plotting_factory` (new) | Cross-language plotting helpers | Python + R subtrees |

Each factory is an independent git repository with its own dependency
management, contracts vendoring, and governance lifecycle. Cross-factory
edits travel via bridge symlinks (three-dot relative paths) and vendored
byte-identical files validated by a cross-repo parity script.

Project data continues to live under `/home/zerlinshen/projects/<project-id>/`
and is not affected by the trifurcation.

## Drivers

1. **Separation of concerns.** Plotting was leaking through R analysis modules
   and Python pipeline modules. Cross-cutting plot code was duplicated and
   theme-config drift had landed in three places. A dedicated `plotting_factory`
   collapses that surface.
2. **Future omics extensibility.** ATAC, Hi-C, VDJ, ribo, and spatial lanes
   are all converging in `r_multiomics_factory`. Keeping pure R analysis in a
   `r_`-prefixed sibling clarifies the import direction and makes future
   per-language factories (e.g., a Julia or Rust sibling) a low-friction add.
3. **Time-budget reality.** User signaled a ~10-day window before real analysis
   must resume. M1 minimum-viable approach (Phase 1 in 1-2 days) defers the
   higher-risk modality plot extraction + figure_bundle schema work to
   Phases 2-3 without blocking analysis resumption.

## Alternatives Considered

| Option | Why rejected |
|---|---|
| **Op1: Keep monorepo** (no trifurcation) | Forecloses the separation-of-concerns goal; codex's recent consolidation was tooling-level only, not directory-level. |
| **Op2: Two factories** (Py + R only, no plotting) | Leaves theme-config + plot helpers duplicated; doesn't address the cross-cutting plotting surface that triggered this work. |
| **Op3: Plotting as a subdirectory** of `r_multiomics_factory` | Couples plot lifecycle to R repo's release cadence; Python plotting would still need a separate home. |
| **Op4: Single plotting factory without language subdirs** | Mixes language-specific tooling (renv vs pip) at top level; harder for agents to reason about which lane owns a file. |
| **Op5: Plotting as a sub-package inside `singlecell_factory`** | Same coupling problem as Op3 from the Python side; doesn't address R plot helpers. |

Per Round-2 deep-interview deliberation, the user explicitly authorized the
most aggressive (Op2-style physical separation) approach: "plot is plot,
analysis is analysis." That confirmation drove the dual-language
`plotting_factory` with `python/` and `r/` top-level subtrees.

## Why Chosen

- **Matches user vision** of clean separation with explicit language lanes
- **Fits time budget**: Phase 1 file moves + renames + path updates land in
  1-2 days; analysis resumption unblocked at Phase 1 acceptance
- **Future-proof for omics expansion**: each new modality adds R-side analysis
  in `r_multiomics_factory` and plotting helpers in `plotting_factory/r/`,
  not in `singlecell_factory`
- **Preserves codex work**: 11 codex commits in SC (boundary `e44dc5c..e782eaf`)
  + 5 codex commits in R (boundary `e3b73d2..cba6162`) are preserved in Phase 0
  rollback floor; no codex content was discarded
- **Atomic rollback floor**: Phase 0 SHA anchors (SC `13c2c885`, R `cba6162`)
  let any phase reset cleanly

## Consequences

### Adopted
- Bridge symlinks now span TWO target sets: `local_r_pipeline_macbook` →
  `r_multiomics_factory` and `local_plot_pipeline` → `plotting_factory`. The
  symlink validator (`scripts/ci/check_bridge_symlink.sh`) was hardened to
  resolve targets with `readlink -f`, assert expected sibling repo prefix,
  reject dangling links, reject non-symlink files under bridge paths, and
  reject self-bridges.
- `renv.lock` lifecycle now spans TWO R repos (`r_multiomics_factory` +
  `plotting_factory`). Phase 2 will prune `plotting_factory/renv.lock` to
  the minimum plotting-only dependency set.
- Cross-factory CI gates extend: contracts parity checked across all three
  repos via vendored byte-identical hash files. `.expected_sha256` must be
  refreshed whenever `bundle_schema.yaml` content changes.
- Vendored `tools/check_contracts_cross_repo.sh` accepts both
  `multiomics_r_factory` (legacy) and `r_multiomics_factory` (current) repo
  names during the transition window. Phase 1.5 reference sweep completed;
  the legacy fallback can be removed at Phase 2 entry per the plan §1.4.5
  step 5 follow-up.

### Surfaced for follow-up
- **C1 codex review** flagged 4 pre-existing test failures in
  `tests/test_modular.py::test_clustering_gpu_fallback*` (codex wave6 work
  changed `clustering.py` fallback policy default to raise without updating
  these tests). NOT caused by Phase 1; deferred to C1.
- `tests/test_modular.py::test_new_modules_in_dag` assertion was 31 (stale);
  pipeline currently registers 42 modules (wave6 multi-omics expansion).
  Phase 1 lead commit refreshed the assertion to 42.
- Architecture validator warns about missing `results/nc2024_tumor_20260426_v2`
  run dir — historical, exit 0, acceptable per plan §1.6.

## Follow-ups

1. **Phase 2 (~3-5 days):** Extract plot halves from 6 modality modules in
   `r_multiomics_factory`; migrate Python plots into `plotting_factory/python/`;
   activate cross-factory CI gate; extend Wave5 v5.1 evidence gate boundary.
   **DONE.**
2. **Phase 3 (~2-3 days):** Define `figure_bundle` schema in
   `plotting_factory/schema/`; activate per-plot YAML config schemas; ship
   visual-regression smoke; flip cross-factory CI gate to blocking;
   finalize governance update. **DONE (2026-05-18).**
   - §3.1: `figure_bundle_schema.yaml` + JSON validator + sync-contracts extension — DONE
   - §3.2: 5 per-plot YAML config schemas — DONE
   - §3.3: Visual regression smoke (5 plot types, plotting_factory `8216331`) — DONE
   - §3.4: `cross_factory_contract_gate.sh` promoted to blocking; extended to verify
     `figure_bundle_schema.yaml` parity across all 3 repos (sha256 `a98ed120...`) — DONE
   - §3.5: Governance finalization across all 3 repos (SC `4eafc8e`, R `d31c3c8`) — DONE
3. **C1 codex review (parallel with Phases 1-2):** Five-cluster review of
   the 11 codex SC commits + 5 codex R commits.
4. **C2 reproduction-methodology-upgrade (parallel with Phase 2):** Doc
   (`docs/PAPER_REPRODUCTION_SOP.md`) + skill
   (`~/.claude/skills/paper-reproduction-from-upstream/SKILL.md`) +
   contract validator (`contracts/project_run_contract.yaml`).
5. **C3 real-run-3D-genome project:** Deferred per user; revisit after
   Phase 3 governance finalization.

## Phase 1 commit ledger

| Repo | SHA | Subject |
|---|---|---|
| singlecell_factory | `0c1c542` | docs: add drift-detector required references (Phase 1 pre-flight drift fix) |
| r_multiomics_factory | `02acb02` | tools: support transition name in check_contracts_cross_repo.sh (Phase 1.4.5) |
| r_multiomics_factory | `f78edf1` | rename: multiomics_r_factory → r_multiomics_factory (Phase 1.2) |
| plotting_factory | `8357b41` | Bootstrap plotting_factory skeleton (Phase 1.1) |
| r_multiomics_factory | `cfbfc5e` | Migrate pure plot files to plotting_factory (Phase 1.3) |
| plotting_factory | `a5d0a8a` | Import pure plot files from r_multiomics_factory (Phase 1.3) |
| singlecell_factory | `10c6bcf` | bridges: re-point + add local_plot_pipeline + harden check_bridge_symlink.sh (Phase 1.4) |
| singlecell_factory | `295b8b7` | rename: SC string references sweep (Phase 1.5 A-E) |
| singlecell_factory | `0ca023f` | smoke gate: refresh contracts hash + update modules-in-dag assertion (Phase 1.6 lead fix) |
| r_multiomics_factory | `db9f8c8` | contracts: refresh .expected_sha256 (Phase 1.6 lead fix) |

**Phase 0 rollback floor:**
- SC `wave6-trevino-v5.1` pre-Phase-0: `13c2c885c981262fc634125cc46f4abc9c08298f`
- R `master` pre-Phase-0: `cba6162583510f84df40bc8f86047fd53ba94ce6`

**Phase 0 post-commit anchors (referenced as rollback target if Phase 1 reverted):**
- SC: `e782eaf2dbc52219725d14450c3c2891206f32aa`
- R: `478dc9ce7b1e3619b7bf12dc40f0b8e916e36adc`
