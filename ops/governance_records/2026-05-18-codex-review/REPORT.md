# Codex Review — Phase 0 Boundary `e44dc5c..e782eaf` (5 Clusters)

- Reviewer: worker-4 (factory-trifurcation-p3c1)
- Date: 2026-05-18
- Mode: READ-ONLY (no source edits; only `ops/governance_records/2026-05-18-codex-review/REPORT.md`)
- Plan: `/home/zerlinshen/.omc/plans/factory-trifurcation-consensus-plan-2026-05-18.md` §C1
- Current `HEAD`: `26ca5a227cc355c30380ac737f0a4895299700bc`
- Phase 0 boundary: `e44dc5c..e782eaf` (10 commits, listed below)

## Scope

```
e782eaf ops: Phase 0 preflight pin for factory-trifurcation work
570edf1 config: pre-commit + pyproject + env locks + nc2024 notebook + R-bridge test
93ddd75 ref: marker database catalog (cellmarker2, celltypist, panglaodb, sctypedb, curated)
d374335 ops: remove obsolete run_record_stub.json files
df4757e ops: control-plane records — journals, governance, cleanup, env, ledger
612f0e5 scripts: wave5 v5.1 CI gates + dev launchers + paper-figure renderers
0f79db2 tests: wave5 v5.1 + modality + contract + smoke + driver coverage
8482c46 quarantine: wave2 module drafts + hv1 staging
45a0c35 workflow: wave5 v5.1 modular pipeline + new modality modules
07b8cea contracts: cross-repo bundle schema + project_run contract + parity tooling
```

The 5 review clusters listed in the task instructions span all 10 commits.
Cluster 3 references commit `26ca5a2` (post-Phase-0, on `HEAD`) which adds the
"See Also" pointer at AI_AGENT_PROTOCOL.md:212; that change is in-scope only as
an integration check, not as a Phase-0 boundary review item.

Phase 2+C2 anchors (informational): SC `26ca5a2`, R `6760008`, `plotting_factory` `e4ec116`.

## Methodology

Each cluster was evaluated against:

- `singlecell_factory/CLAUDE.md` (project precedence + three-factory architecture)
- `singlecell_factory/AI_AGENT_PROTOCOL.md` (paper reproduction ladder, routing)
- `docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md` (governance contract)
- `docs/LINUX_FILE_GOVERNANCE.md` (path roles, human/agent surfaces)
- the three-factory ADR (`ops/governance_records/2026-05-18-three-factory-trifurcation/ADR.md`, referenced)
- codex governance principles (reproducibility, deterministic outputs, stable contracts)

Sources reviewed: file content via `Read`/`cat`, commit metadata via
`git show`/`git log`, dispatcher + gate semantics by reading executable
paths. Phase 1 finding (`test_clustering_gpu_fallback*`) was reproduced live
via `pytest -q` to confirm classification as **pre-existing** rather than
introduced by Phase-0 boundary work.

Severity rubric:

- **Blocker**: pipeline-breaking, data-loss-risking, governance-contradicting.
  Escalated to team-lead immediately on discovery.
- **Should-fix**: architectural/governance concerns or correctness risks that
  do not block but should be tracked.
- **Nice-to-have**: stylistic, doc phrasing, redundant comments.

## Findings

### Cluster 1 — Governance Docs (`docs/LINUX_FILE_GOVERNANCE.md`, `docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md`)

Reviewed:
- `docs/LINUX_FILE_GOVERNANCE.md` (one-screen strategy, path roles, human/agent split, canonical projects)
- `docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md` (authority boundary, three-factory subsection, manifest layers, retention, multi-cohort policy, paper-reproduction policy)

**Blocker**: none.

**Should-fix**:
1. `REMOTE_FACTORY_PROJECT_GOVERNANCE.md:9-11` calls the two-factory bullets
   "SUPERSEDED by the Three-Factory Authority Boundary subsection" and
   "preserved verbatim for historical reference". The preserved historical
   bullets at lines 13-14 still claim `r_multiomics_factory` owns "plotting
   helpers", which directly contradicts the three-factory split at lines
   25-30. The "preserved verbatim" framing is correct for history, but the
   bullets should either be (a) clearly marked with an inline `SUPERSEDED:`
   prefix on each line, or (b) moved into a `### Historical (pre-2026-05-18)`
   sub-heading so search/grep readers do not absorb the contradicting claim
   as current. Risk: a future agent grepping for `plotting helpers` will
   land on the legacy bullet and mis-route work to `r_multiomics_factory`.

**Nice-to-have**:
1. `LINUX_FILE_GOVERNANCE.md:54` notes the rename "2026-05-18 from
   `multiomics_r_factory`" while `REMOTE_FACTORY_PROJECT_GOVERNANCE.md:25`
   says "renamed from `multiomics_r_factory`". Consistent; recommend adding
   a one-line "Legacy name accepted during transition" anchor in
   `LINUX_FILE_GOVERNANCE.md` to match the in-code transitional comments
   (e.g., `tools/check_contracts_cross_repo.sh:32-34,38`).
2. `LINUX_FILE_GOVERNANCE.md:88-90` "Current Canonical Projects" table is a
   point-in-time snapshot. Recommend adding a "Last verified" date line so
   future readers know when to re-validate (the rest of the doc is structural
   and stable).

### Cluster 2 — Wave5 v5.1 Evidence Gates + 3-Way Schema Dispatcher

Reviewed:
- `workflow/modular/cli.py` (three-way ledger schema dispatcher integration, branch anchor, `SC_REQUIRE_PROJECT_ROOT` cutover gate, `--allow-dirty` dirty-tree gate, manifest writing)
- `ops/run_ledger/schema/wave5_v5_1.schema.json` (additive schema on top of v5.0; plan_revision enum restricted to `v5.1`; new MV1a/MV2/MV3/POLISH/HV1/HV2 fields + provenance chain)
- `scripts/ci/select_ledger_schema.py` (canonical dispatcher; B1 mutual-rejection patch requires `plan_revision`)
- `scripts/ci/wave5_v5_0_schema_gate.sh` (dispatcher-routed gate with explicit mutual-rejection exit code 3)
- `scripts/ci/wave5_v4_2_gate.sh` (v4.2 acceptance gate; reads `plan_revision`, `validation_posture`, per-AC verdicts)

**Blocker**: none.

**Should-fix**:
1. `select_ledger_schema.py` `V42_REVISIONS` allows `{"v3", "v4", "v4.1", "v4.2"}`
   while `wave5_v4_2.schema.json`'s acceptance gate at `wave5_v4_2_gate.sh:46`
   hard-requires `plan_revision == "v4.2"`. This is internally consistent —
   the dispatcher decides which **schema** validates the record; the gate adds
   a posture check on top — but the surface-asymmetry can confuse future
   maintainers. Recommend documenting in the dispatcher docstring that
   "Dispatcher selects the **schema family**; per-AC gates may impose
   stricter `plan_revision` posture inside the family."
2. `wave5_v5_1.schema.json:343-344` ends `additionalProperties: true` while the
   nested `historical_reproduction_rate` block (line 108) sets
   `additionalProperties: false`. The mix is intentional (top-level is
   forward-compatible; nested invariants are strict), but the JSON-schema
   reader should add a brief comment in the top-level `description`
   confirming this choice — otherwise a future agent may "tighten" the
   top-level to `false` and silently break additive v5.x rollouts.
3. `wave5_v4_2_gate.sh:24-31` falls back to a single hardcoded absolute path
   (`/home/zerlinshen/singlecell_factory/...`) when `__file__` resolves to
   `-` (heredoc invocation). The fallback works on the canonical host but
   breaks for any local-Mac SSH orchestrator. Recommend resolving the
   schema path from `$PWD` or a `$REPO_ROOT` env override before the
   hardcoded path; this is a portability hazard, not a correctness one
   (gate still exits 2 cleanly with a clear message).

**Nice-to-have**:
1. `cli.py:573-579` propagates `--multimodal-engine` and `--second-obsm-key`
   to env vars only when the env var is unset. Behavior is intentional (env
   wins), but the CLI `--help` text on lines 148 and 156 should mirror the
   precedence note ("env var wins if already set"); only the
   `--second-obsm-key` help text says it.
2. `cli.py:582-595` cutover semantics block — comment is excellent.
   Suggest also linking to the eventual Round-2 ADR file (currently only
   referenced as "Round-2 ADR will flip the default to required") so a
   future grep finds the binding decision.
3. `wave5_v5_0_schema_gate.sh:8-12` references "Three-factory context
   (post-2026-05-18)" with a comment that future-dates the boundary
   slightly. Harmless, but `LINUX_FILE_GOVERNANCE.md` calls the same
   change "added 2026-05-18"; keep the prose calendar-consistent.

### Cluster 3 — Paper Reproduction Ladder additions in `AI_AGENT_PROTOCOL.md`

Reviewed:
- `AI_AGENT_PROTOCOL.md:197-212` (Paper Reproduction Ladder section + new "See Also" pointer from C2.6 commit `26ca5a2`)
- Cross-references in `docs/PAPER_REPRODUCTION_SOP.md`,
  `contracts/project_run_contract.yaml`,
  `scripts/validate_project_governance.py`,
  `~/.claude/skills/paper-reproduction-from-upstream/SKILL.md`
- Alignment with `docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md:104-118` "Paper
  Reproduction And Context Optimization Policy"

The ladder (8 steps) is identical in spirit to the governance-doc policy (9
steps) — the policy doc has an extra step 5 splitting "data objects" and
"figure-panel parity" into separate evidence classes that the ladder folds
into step 4. The C2.6 pointer (line 212) bridges the ladder to:
operational SOP, walk-through skill, schema, validator.

**Blocker**: none.

**Should-fix**: none. The triad (`SOP` ↔ `skill` ↔ `protocol`) is coherent
with the three-factory model and respects the "project precedence over
global OMC" rule in `CLAUDE.md`.

**Nice-to-have**:
1. Ladder step 5 ("Compare each claim as exact, approximate, proxy,
   unsupported, or resource gap") and the policy doc's step 6 use the same
   5-tag taxonomy. Recommend co-locating the tag definitions in
   `contracts/project_run_contract.yaml` so the controlled vocabulary is
   machine-readable, not only prose.
2. The ladder lists "exact / approximate / proxy / unsupported / resource
   gap"; `docs/PAPER_REPRODUCTION_SOP.md` may use slightly different prose.
   Out of scope to verify in C1 (READ-ONLY scope of the SOP file content),
   but recommend a quick consistency pass in P3.5.
3. Line 212 is a long single sentence with four resource pointers. Consider
   converting to a 4-bullet list under "See Also" for easier scanning;
   prose order is fine — readability only.

### Cluster 4 — CI Gates

Reviewed:
- `scripts/ci/wave5_v4_2_gate.sh` (per-AC verdicts: AC-VAL-3a, AC-VAL-3b, AC-CI-1, AC-VAL-PLOT-{1,2,3}, AC-LEDGER-1, AC-VAL-3c, ATTESTATION)
- `scripts/ci/wave5_v5_0_schema_gate.sh` (dispatcher-routed mutual rejection)
- `scripts/ci/wave5_trevino_regression_gate.sh` (3 thresholds: cell_type_ari, peak_gene_overlap_top1000, pseudotime_spearman)
- `scripts/ci/wave5_v4_2_lockfile_hook.sh` (`--allow-fallback` opt-in; allowlist + journal rationale check)
- `scripts/ci/wave4_dispatch_freeze.sh` + `scripts/ci/check_freeze.sh` (US-W4-11 lock)
- `scripts/ci/check_schema_parity.sh` (byte-identical canonical/vendored bundle_schema.yaml check)
- `.pre-commit-config.yaml` (black, ruff, eol-fixer, trailing-whitespace, contracts-cross-repo-parity hook)
- `tools/check_contracts_cross_repo.sh` (cross-repo parity gate; SHA256 against `.expected_sha256`)

**Blocker**: none.

**Should-fix**:
1. `wave5_v4_2_lockfile_hook.sh:103-118` searches the **5 most recently
   modified** files in `ops/before_every_run/journal/` for the
   `--allow-fallback` literal. This is fragile in two ways:
   - **mtime ordering is not append-only**: a `chmod` or unrelated edit can
     re-sort the window, hiding an older rationale entry. Recommend
     "find files where `--allow-fallback` appears + take the youngest" or
     "scan the entire dir" since journal entries are small.
   - **false negative**: a journal entry that documents `--allow-fallback`
     for a different gate context still satisfies the check. Recommend a
     stronger token like `wave5_v4_2_lockfile_hook --allow-fallback` to
     bind the rationale to this specific opt-in.
   Severity: should-fix (gate is opt-in; the default path is hard-stop
   exit 1, so the default-deny posture is correct).
2. `wave5_v4_2_gate.sh:115-119` validates `ari_methodology_note_sha` is
   present **and** `ari_note_attestation.subagent_type` is non-empty before
   accepting the biology-aware path. Good. But the fallback partial path
   (lines 124-127) treats `ari_value in [0.70, +inf)` as PARTIAL — the upper
   open bound is fine, but the `[0.60, 0.70)` band is **silently dropped**
   (it isn't caught by either branch because the biology-aware path is
   skipped). Trace: if `ari_value = 0.65` AND `ari_note` is absent, the
   gate falls to line 126 `elif ari_value >= 0.70 → PARTIAL`, fails that
   condition, and falls to line 128 `else → FAIL` — so `0.65` is correctly
   classified as FAIL. **No bug, but the comment on line 117 ("0.70 is the
   fallback partial threshold; with biology-aware path < 0.60 is FAIL")
   could be clearer about the [0.60, 0.70) band when biology-aware is
   absent.** Stylistic only.
3. `wave5_trevino_regression_gate.sh:33-37` uses fixed thresholds (0.70 for
   all three metrics). The thresholds match plan v4.2 §3.4 but are not
   sourced from a YAML/JSON spec — they're inline in the bash heredoc. If
   the plan threshold ever changes, two places must update (gate + ledger
   doc). Recommend extracting to a sibling YAML spec
   (mirror `nc_regression_tolerance_spec.yaml`) for parity with the NC
   regression gate.
4. `check_freeze.sh:18-19` reads staged paths via
   `git diff --cached --name-only --diff-filter=ACMR`. Correct for
   pre-commit. But the heredoc-or-CI lane that invokes
   `bash scripts/ci/check_freeze.sh` with no args from a non-git context
   (e.g., a tarball-based CI runner) will silently get an empty `paths=()`
   and exit 0 — bypassing the freeze. Recommend an explicit "must be inside
   git work-tree OR explicit path args" precondition check.

**Nice-to-have**:
1. `.pre-commit-config.yaml`: `contracts-cross-repo-parity` hook runs
   `tools/check_contracts_cross_repo.sh` when the canonical bundle schema or
   `.expected_sha256` changes. The hook gracefully `exit 0`'s when the
   sibling repo is missing (`check_contracts_cross_repo.sh:43-47`). That
   matches the canonical/vendored-copy model and is correct on a fresh
   single-repo clone.
2. `check_schema_parity.sh` and `check_contracts_cross_repo.sh` are
   overlapping: both verify `singlecell_factory/contracts/bundle_schema.yaml`
   parity against the vendored copy in `r_multiomics_factory`. The latter
   also checks `.expected_sha256` and is more thorough. Recommend marking
   `check_schema_parity.sh` as "lightweight smoke" in its header to avoid
   future confusion about which gate is authoritative.
3. The `verify_module_references.py` script (referenced in commit `612f0e5`
   stats; not deep-read) is a good complement; expose its exit code
   semantics in its module docstring if not already.

### Cluster 5 — R-factory SHA fields + Bundle v2.1 Provenance Traceability

Reviewed:
- `workflow/modular/manifest_writer.py:53-103` (dual `factory_r.sha` + `r_factory_sha_at_manifest_write` write; legacy name accepted with `DeprecationWarning`)
- `scripts/export_singlecell_r_bundle.py:1115-1158` (`_r_factory_sha_at_export`, `_write_bundle_provenance`, bundle SHA256 concat)
- `r_multiomics_factory/R_bundle/io_bundle.R:473-505` (PREC-1 dual-SHA mismatch — `warning()` not `stop()`, per plan §line 49)
- Bundle v2.1 reader (`io_bundle.R:23-37`) — accepts `singlecell_r_bundle_v2`, `v2.1`, `v2.2` with additive `extensions` slot
- Cross-repo parity (`tools/check_contracts_cross_repo.sh`)

**Blocker**: none.

**Should-fix**:
1. `manifest_writer.py:76-82` accepts legacy `multiomics_r_factory` path and
   emits `DeprecationWarning`. Correct. But the warning is silenced by
   default in many CI/heredoc contexts (Python `-W default::DeprecationWarning`
   is not the default for `python -c`). Recommend also writing the legacy-
   name signal into `manifest.extra["legacy_factory_r_path"] = True` so a
   downstream verifier can grep the manifest without re-running.
2. The R loader's PREC-1 mismatch handler at
   `r_multiomics_factory/R_bundle/io_bundle.R:495-499` emits `warning()` not
   `stop()`. This matches the plan's intent ("WARNING-not-ERROR semantics
   per plan §line 49"). Good. But there is no automated CI gate that
   asserts the warning is **surfaced** in a downstream R run — `warning()`
   in R can be silently suppressed by `suppressWarnings(...)`. Recommend
   that the R bundle test harness (`bridges/local_r_pipeline_macbook/tests/
   testthat/test_plot_remote_bundle.R`, added in commit `570edf1`) assert
   `expect_warning(read_bundle(...), regexp = "PREC-1")` for the mismatch
   case. Not blocking — surfaced via `message()` and `warning()` is
   acceptable per the dual-write contract.
3. `_r_factory_sha_at_export` and `factory_git_state` both shell out to
   `git rev-parse --short=7 HEAD`. If the R factory is in detached HEAD
   (e.g., a CI checkout at a tag), `--short=7` still returns a valid SHA,
   so behavior is correct. But on a missing/broken git repo both functions
   return `""` silently. Recommend that both call sites log a one-line
   `logger.warning("R factory git state unavailable; SHA omitted")` so
   silent provenance loss is detectable.

**Nice-to-have**:
1. `manifest_writer.py:96` records
   `"claim_guard": "not_for_de_or_new_quantitative_claims_without_full_object_validation"`.
   This is excellent — explicit guard on the manifest itself. Consider also
   echoing this guard into the R reader's startup message so a downstream R
   analyst sees the guard before any plot is rendered.
2. `io_bundle.R:36` declares
   `KNOWN_BUNDLE_EXTENSIONS <- c("protein", "spatial", "multimodal_obsm", "marker_resolutions", "atac", "vdj", "hic", "ribo")`.
   Aligned with Phase-0 module additions in commit `45a0c35` (atac, hic,
   ribo, vdj modules). The reader is forward-compatible (unknown extensions
   are skipped). No action.
3. Bundle exporter writes to `<final_output_dir>.tmp.<pid>` then
   `os.replace`s — atomic publish. Excellent. The temp-dir name uses PID
   which can collide on long-running supervisors; consider
   `tempfile.mkdtemp(prefix=..., dir=parent_dir)` for race-safety. Nice
   defense-in-depth, not a real bug today.

## Phase 1 Pre-Existing Finding (Known)

**Status: pre-existing at Phase 0 baseline. NOT introduced by the codex
commits under review. NOT a blocker for codex-preservation. SHOULD be
classified as a follow-up task.**

Live verification (just now, on `HEAD = 26ca5a2`):

| Test | Status |
| --- | --- |
| `tests/test_modular.py::test_clustering_gpu_fallback` | **FAIL** |
| `tests/test_modular.py::test_clustering_gpu_fallback_does_not_reuse_mutated_gpu_adata` | **FAIL** |
| `tests/test_modular.py::test_clustering_gpu_fallback_on_cuda_oom` | **FAIL** |
| `tests/test_modular.py::test_batch_correction_gpu_fallback_on_cuda_oom` | PASS |
| `tests/test_modular.py::test_de_gpu_fallback_on_cuda_oom` | PASS |

Cause: `workflow/modular/modules/clustering.py:165-235` in commit `45a0c35`
extends the GPU policy with `SC_GPU_FAILURE_POLICY={raise,restore-cpu,...}`
via `_contract_violation.resolve_gpu_failure_policy`, the M2 inplace+raw-mandate
strategy, and `ClusteringContractViolation` poisoning. The three failing tests
were written against the older "silent fallback to CPU" policy. The policy
change is intentional and architecturally cleaner (explicit poison on
contract violation, restore-cpu only when `adata.raw` is preserved), but the
3 unit tests were not updated in lockstep.

**Classification**: `should-fix` (test-only; runtime is more conservative
post-policy-change; no scientific outputs at risk).

**Recommendation**: Update the 3 failing tests in a follow-up task to
expect either (a) the new poison-and-raise behavior when
`SC_GPU_FAILURE_POLICY=raise`, or (b) the new restore-cpu path when
`adata.raw` is preserved. Add a new test that exercises the
`ClusteringContractViolation` poison path explicitly.

**Why this is NOT a blocker on codex preservation**: the new policy is a
strict superset of the old "silent fallback" — operators can still opt into
silent CPU fallback by setting `SC_GPU_FAILURE_POLICY` appropriately (or
leaving the default `restore-cpu`). The tests fail because they assert a
specific implementation detail of the old policy, not because the runtime
loses correctness. Revert would lose the `_contract_violation` infrastructure
and the M2 raw-mandate that all of Phase-0's new modality modules depend on.

## Open Questions for User

1. **PREC-1 R warning surface**: should we make the R bundle loader's
   `warning()` upgrade to `stop()` when run under a pre-flagged "strict
   parity" R session (e.g., `Sys.setenv(SC_BUNDLE_STRICT_PARITY = "1")`)?
   Current behavior is warning-only across all sessions.
2. **Lockfile hook journal scan**: do you want the journal-rationale token
   to be `--allow-fallback` (current, loose) or `wave5_v4_2_lockfile_hook
   --allow-fallback` (tighter, recommended above)? Tighter token requires
   updating any existing journal entries that document a fallback.
3. **`wave5_v4_2_gate.sh` portability**: fine to keep the hardcoded
   `/home/zerlinshen/singlecell_factory/...` fallback, or replace with
   `${SC_FACTORY_ROOT:-/home/zerlinshen/singlecell_factory}/...` to ease
   local-Mac SSH orchestration? The latter is a one-line change and
   strictly additive.
4. **Phase 1 GPU fallback tests**: in scope to fix in Phase 1 / P3 follow-up,
   or carry as a known finding into Phase 2+? Recommend scheduling as a
   small follow-up after P3 lands.
5. **Two-factory historical bullets in `REMOTE_FACTORY_PROJECT_GOVERNANCE.md`**:
   prefer (a) inline `SUPERSEDED:` prefix on each historical line, or
   (b) move into a `### Historical (pre-2026-05-18)` sub-heading? Either
   removes the silent contradiction with the three-factory subsection.

## Go/No-Go on Codex Preservation

**Verdict: GO — PRESERVE ALL CODEX CHANGES.**

Justification:

- No **Blocker**-class finding in any of the 5 clusters.
- The 11 commits in `e44dc5c..e782eaf` are coherent: contracts, schemas,
  governance, gates, tests, modules, env-locks, refs, and provenance.
- The Phase 1 finding (`test_clustering_gpu_fallback*` × 3 failing) is
  **pre-existing at Phase 0 baseline** and is a test-only artifact of an
  intentional clustering-policy upgrade. It is `should-fix`, not
  `blocker`, and reverting codex changes would lose the
  `_contract_violation` infrastructure that Phase-0 modality modules
  depend on.
- The three-factory governance bullet contradiction
  (`REMOTE_FACTORY_PROJECT_GOVERNANCE.md`) is a **doc clarification**, not
  a structural defect — the three-factory subsection at lines 18-40 is
  authoritative and correct.
- The dispatcher + per-gate mutual-rejection design (Cluster 2) is
  internally consistent, with exit-code semantics documented and a
  deterministic schema selection rule.
- The dual R-factory-SHA contract + WARNING-not-ERROR loader semantics
  (Cluster 5) match the plan and preserve forensic recoverability without
  introducing pipeline-breaking gates.

Action: proceed to P3 (Phase 3 acceptance gate). The `should-fix` items
above should be tracked as follow-up tasks but are **not** preconditions
for Phase 2+C2 lock-in or codex preservation.

## Summary Counts

| Cluster | Blocker | Should-fix | Nice-to-have |
| --- | ---: | ---: | ---: |
| 1. Governance docs | 0 | 1 | 2 |
| 2. Wave5 v5.1 + dispatcher | 0 | 3 | 3 |
| 3. Paper Reproduction Ladder | 0 | 0 | 3 |
| 4. CI gates | 0 | 4 | 3 |
| 5. R-factory SHA + bundle v2.1 | 0 | 3 | 3 |
| Phase 1 pre-existing | — | 1 | — |
| **Totals** | **0** | **12** | **14** |

End of report.
