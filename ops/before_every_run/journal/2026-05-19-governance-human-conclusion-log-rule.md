# 2026-05-19 governance update: human conclusion log requirement

## Objective

Make human-readable conclusion logs a durable governance requirement rather
than a one-off convention from the NC2024/Cell clustering discussion.

## Confirmation

Existing rule already required remote run journaling through:

- `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/LATEST.md`
- `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/journal/`
- `before-every-run` skill preflight/post-run journaling

Gap found: the rule did not explicitly require human-facing scientific decision
logs when a discussion changes interpretation, final-run strategy, benchmark
lane choice, claim support, or source-of-truth status.

## Changes Made

- Updated suite root README:
  `/home/zerlinshen/Bioinformatics Research Pipeline/README.md`
- Added suite governance README:
  `/home/zerlinshen/Bioinformatics Research Pipeline/governance/README.md`
- Updated singlecell governance/readme surfaces:
  `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/README.md`
  `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md`
- Updated cross-factory repo READMEs:
  `/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory/README.md`
  `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/README.md`
- Updated skill surfaces:
  `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/.codex/skills/before-every-run/SKILL.md`
  `/Users/zerlinshen/.codex/skills/before-every-run/SKILL.md`

## New Durable Rule

When a run discussion changes scientific interpretation, final-run strategy,
claim support, benchmark lane choice, source-of-truth status, or human-facing
next actions, write a human-readable conclusion log in addition to operational
run memory.

Required locations:

- project ledger:
  `/home/zerlinshen/projects/<project-id>/ledger/human_review/<date>-<topic>.md`
- run evidence when tied to a concrete run:
  `/home/zerlinshen/projects/<project-id>/runs/<run-id>/evidence/<topic>.md`
- factory run memory index:
  `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/LATEST.md`
  and matching journal entry under `ops/before_every_run/journal/`

The log must separate observed facts, interpretation, decisions, rejected or
downgraded routes, next actions, and artifact classification (`canonical`,
`evidence-only`, `superseded`, or `failed exploratory`).

## Verification

- `git diff --check` passed for changed markdown/skill files in
  `singlecell_factory`, `r_multiomics_factory`, and `plotting_factory`.
- Confirmed rule text is present in:
  suite README, suite governance README, singlecell README, singlecell remote
  governance doc, r_multiomics README, plotting README, and remote
  `before-every-run` skill.

## Artifact Classification

- This governance update is `canonical` operational policy.
- It does not change scientific source-of-truth runs.
- It formalizes the already-created NC2024/Cell clustering human conclusion log
  as the first example of the policy.

