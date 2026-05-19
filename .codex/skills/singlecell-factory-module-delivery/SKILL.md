---
name: singlecell-factory-module-delivery
description: Add or evolve a module in the remote `singlecell_factory` pipeline safely. Use when introducing a new module, integrating a paper method as a module, extending module wiring, or changing module contracts/outputs in `workflow/modular`. Before implementation, read the current run memory and route through the Shenxin execution chain so module work stays aligned with remote validation, docs, and handoff.
---

# Singlecell Factory Module Delivery

Use this skill when the task is to **introduce a new module** or materially evolve an existing one in the remote `singlecell_factory` pipeline.

This is the project-specific module-introduction surface. It is narrower and more operational than the generic `develop-and-integrate-module` skill.

## Read First

Before doing any module work, read:

- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run/LATEST.md`
- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/first reproduction - Nature Communications 2024 NSCLC single-cell/00_START_HERE_AGENT_HANDOFF.md`
- [before-every-run/references/shenxin-routing-contract.md](/Users/zerlinshen/.codex/skills/before-every-run/references/shenxin-routing-contract.md)

Then choose the next skill lane:

1. `before-every-run`
2. `bioinformatics-autosteer`
3. `singlecell-factory-module-delivery`
4. `execute-and-recover-pipeline`
5. `review-and-validate-quality`
6. `readme-sync-enforcer`

## When To Use This Skill

Use it when you need to:

- add a new module under `workflow/modular/modules/`
- add new module dependencies or registry wiring
- add new module CLI/config parameters
- integrate a paper method into the modular pipeline
- define new `requires_keys` / `provides_keys`
- introduce a new output contract that must still fit project conventions

Do not use it for:

- ordinary run monitoring
- rerunning a failed pipeline without code change
- local report packaging only

## Delivery Workflow

## 1. Define the module contract first

Before any implementation, write down:

- module name
- purpose
- expected inputs
- expected outputs
- required upstream modules
- output files or figures
- whether the module is optional or mandatory

If the module is paper-driven, capture:

- paper title
- method summary
- critical assumptions
- known failure modes

## 2. Wire it in the remote pipeline, not as a one-off script

Treat the remote repo as canonical:

- `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory`

Required integration surfaces usually include:

- module file in `workflow/modular/modules/`
- pipeline registry and dependencies
- config/CLI plumbing if the module adds parameters
- module status / output behavior

## 3. Follow the current large-cohort safety pattern

When the module touches large objects:

- preserve sparse/lazy structures where possible
- materialize only at narrow third-party boundaries
- avoid hidden densification on full-cohort paths

This is especially important for NC2024-style full-cohort work.

## 4. Validate on the correct execution lane

Default validation choice:

- `execution_mode = debug_massive`

Use `controller_validation` only when the module is already stable enough that the orchestration path itself needs proof.

Validation expectations:

- identify first failing module
- prove that the new module runs or that the frontier moved
- keep exact run directory and evidence

## 5. Close the loop

After the module change:

- update run memory through `before-every-run`
- update docs through `readme-sync-enforcer`
- update the skill system if the new module changes workflow expectations

If the new module materially changes how skills should be used, also update:

- `/Users/zerlinshen/Downloads/1. Codex/SKILL_PROTOCOL.md`
- `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/SKILL_PLAYBOOK.md`

## Verification Checklist

- module contract is explicit
- registry/dependency wiring is complete
- outputs are categorized and documented
- remote validation run exists
- run memory was updated
- README/PROTOCOL changes were made if behavior changed
- `SKILL_PROTOCOL.md` was updated if this new module changes workflow usage
