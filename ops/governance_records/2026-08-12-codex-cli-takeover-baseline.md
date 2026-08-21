# Codex CLI takeover baseline — 2026-08-12

## Purpose and scope

This record anchors the Codex Desktop takeover of the mature Codex CLI project
without copying the repository or importing Claude Code state wholesale.

- Canonical primary root:
  `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory`
- Retired path:
  `/home/zerlinshen/singlecell_factory` (absent; do not recreate or symlink)
- Linked read-only scope:
  `r_multiomics_factory`, `plotting_factory`, the suite root, and the top-level
  project directories under `/home/zerlinshen/projects`
- Explicitly excluded from synchronization:
  authentication files, credential stores, local-only settings, session logs,
  caches, and transient runtime locks
- This takeover step launches no pipeline, changes no source file, resets no
  worktree, creates no compatibility link, and makes no commit.

## Authority read before takeover

The takeover read and applied these sources in order:

1. `AI_AGENT_PROTOCOL.md`
2. `AGENTS.md`
3. `CODEX.md`
4. `ops/before_every_run/LATEST.md`
5. `ops/before_every_run/journal/2026-08-07-post-audit-triage-and-literature-grounding.md`
6. `.omx/plans/wave5-v51-codex-takeover-consensus-20260517.md`
7. `.omc/state/mission-state.json` and `.omx/ultragoal/goals.json`

## Takeover decision

The code and Codex project assets are already on the same filesystem used by
Codex Desktop. The safe operation is an in-place context takeover, not a file
copy. Project-local `AGENTS.md`, `CODEX.md`, `.codex/skills/`, history, and OMX
evidence remain authoritative at the canonical root.

The May 2026 OMX state is retained as historical evidence and must not be
auto-resumed:

- all missions in `.omc/state/mission-state.json` are `done`; its last update
  is `2026-05-28T03:50:18.808Z`;
- Ultragoal `G001` is marked `in_progress` but explicitly `superseded` by G002;
- G002 is `review_blocked`, while the later G003 records the independent review
  as `complete`;
- the latest operational journal dated 2026-08-07 reports `triage_only`, no
  incomplete run in the prior 14-day window, and no pipeline run launched.

Therefore, current working-tree evidence and the August journal take precedence
over stale workflow flags. No old mode is resumed by this takeover.

## Primary repository Git baseline

Captured before this record was created.

- branch: `wave6-trevino-v5.1`
- HEAD: `d2fc4a34e6a1fc7a4a941e7866340a4640e02670`
- upstream: `origin/wave6-trevino-v5.1`
- ahead/behind: `+0/-0`
- origin: `git@github.com:zerlinshen/singlecell_factory.git`
- staged changes: none
- unstaged binary diff SHA-256:
  `d82b2b622ff9356bba351af473bff7dd7869145f9100a8b7cfcf98f194aec283`
- staged binary diff SHA-256 (empty input):
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`
- diff stat: 9 tracked files changed, 295 insertions, 20 deletions

### Modified tracked files

The hashes below are working-tree content hashes, not Git blob IDs.

```text
10dcaf62cf26f2f5c3960d3de2c42f36fb029b9d1f6d2d7e0d914cfc047ac660  AGENTS.md
9fec7a9b31e7422619404ce4e918cfe4e5425077276d6b96901980828b32aabb  README.md
45a615e2ac101b6fbf196c0c75cb50becec344c1121a8deb1c21f38cfb54a3b4  contracts/.expected_sha256
8ed3c592fc8b4daa8f4133694c6fa868b559f6aba2a26e8ea439165d66458b28  contracts/bundle_schema.yaml
4d800c4248234942dca39a2191d2e36227ced033d95d9e23ce15c1fba4d29e9b  docs/MULTIOMICS_MODULE_RATIONALE.md
51f00b52f52bd90eb81b5fb1272d5a141a977900a03cc64cc7b694330a4aa02f  ops/before_every_run/LATEST.md
2208ec37c8562982a4bdb9c1c45896a739063d1dd2614c9d25a3b96662a0ae19  scripts/export_singlecell_r_bundle.py
2e35822c9db93772b7e74ec5f2e2a38f57ba9d7eaf29bbb3bd9e433a88af2964  tests/test_singlecell_r_bundle_export.py
9b22fa5c5f36fddd1bbf17ae2634a358930e5a38d875c2dede165d7bdc9bb5f4  tests/test_wave2b_ribo_smoke.py
```

### Pre-existing untracked files

These hashes were captured before this takeover record was added.

```text
f522f24708f16b6daab8615f7a17f7dbdffc6e23fb9200326c0a7192e79958f0  .kimi-code/skills/h5ad-sanity-check/SKILL.md
f6c74e3d8244e18845b65eeeba916669393168457604a0c0ae9730cda609542b  .kimi-code/skills/h5ad-sanity-check/h5ad_sanity.py
45bf87dc8cbb2411dfd7f5cba4c18a9d01b9e437b4e09c945fbcd320648ca1e8  .kimi-code/skills/plot-qc-review/SKILL.md
6accab4f5fecddf9ecbf3f85b4b2a8d3f0ce24fed682f7c9bd8cd4aaab672a77  ops/before_every_run/journal/2026-08-06-gate-chain-smoke-project-guards.md
f66efd2fe87a1a553425683d5dbe7818b6b26d30ee4921a5c2335ed6ae403d6c  ops/before_every_run/journal/2026-08-07-post-audit-triage-and-literature-grounding.md
```

### Worktrees

```text
/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory
  branch: wave6-trevino-v5.1
  HEAD: d2fc4a34e6a1fc7a4a941e7866340a4640e02670

/home/zerlinshen/.worktrees/pipeline-remediation-20260802-verify/singlecell_factory
  detached HEAD: f67cfad32660a9c301f401890f7bdb7b2cb3a60d

/home/zerlinshen/.worktrees/publication-baseline-cleanup-20260803/suite/singlecell_factory
  branch: codex/publication-baseline-cleanup-20260803
  HEAD: dcd1172a72b5bc7f45cee5ec002cb431d4df1125
```

## Linked suite status (read-only)

### Suite root

- path: `/home/zerlinshen/Bioinformatics Research Pipeline`
- branch: `master`
- HEAD: `a60d15d8c7e23cea65f49b63344d9b57eec58af4`
- upstream divergence: `+0/-0`
- dirty: 6 modified tracked files and 5 untracked files
- takeover rule: observe only; do not fold these changes into the primary repo

### r_multiomics_factory

- branch: `master`
- HEAD: `4707744df82b02e1cd41c2653170520f79cd03aa`
- upstream divergence: `+0/-0`
- dirty: 8 modified tracked files and 4 untracked paths
- takeover rule: read-only dependency and contract peer

### plotting_factory

- branch: `master`
- HEAD: `9bb839d3c40e6b971f351f1cfb100b2b7493fc19`
- upstream divergence: `+0/-0`
- clean at capture time
- takeover rule: read-only presentation dependency

## Linked project directory inventory (read-only)

There were 27 non-hidden top-level project directories at capture time:

```text
a549-crossstudy-compatibility-v1
claim-integrity-validation-20260802
figure-quality-fix
gpu-probe-20260527
h2170-fosl2-ko-validation
h2170-serpine2-loop-verify
h2170-serpine2-loops
hgmm-smoke
integration-bench-20260526
ledger
lusc-gt-concordance-20260728
lusc-integration-gate-20260527
nc-reproduction
pipeline-optimization-20260529
pipeline-optimization-20260610
pipeline-publication-baselines
pipeline-scientific-audit-20260805
pipeline-validation-20260527
pipeline-validation-20260528
pipeline-validation-20260719
pipeline-verify-20260526
production-projects
raw-fasta-smoke
reproductions
round9-singlecell-comparison
skill-integration-assessment
wave5-trevino
```

No project directory was entered for mutation and no project output was
created, moved, or deleted.

## Current operational truth carried forward

- Yost / NG2025 / LKB1 is retired and must not be reopened or presented as live
  work.
- The August 2026 canonical evidence is the Wave-9 LUCA/NSCLC work recorded in
  `ops/before_every_run/LATEST.md` and its linked project evidence.
- The 2026-08-07 triage found no failed or incomplete recent run requiring
  automatic recovery.
- The current dirty tree contains substantive bundle-schema/export/test work.
  Ownership and intended outcome must be identified before any edit, commit,
  reset, run, or cross-repository contract synchronization.
- The dirty-tree baseline is evidence to preserve, not cleanup material.

## Safe next execution entry

Before changing code, the next operator must:

1. Recompute `git status --porcelain=v2 --branch` and the unstaged binary diff
   SHA-256.
2. Compare them with this record and explain any drift.
3. Review the current diff by concern: bundle schema/export, tests, docs, and
   project-guard additions.
4. Identify which concern the user wants resumed.
5. Run only the focused validation for that selected concern; do not launch a
   broad pipeline merely to prove takeover.

## Stop condition for this takeover phase

This phase is complete when this record exists, its captured hashes reproduce
for every pre-existing modified/untracked file, and a post-write status check
shows that the only new takeover artifact is this governance record.
