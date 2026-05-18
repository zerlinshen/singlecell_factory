# Remote Factory/Project Governance

This document is the phase-1 remote governance contract for `singlecell_factory` and `/home/zerlinshen/projects`.

It is a control-plane document. It is not a scientific result and must not be cited as biological evidence.

## Authority Boundary

- `singlecell_factory` is a factory: it owns code, contracts, validators, schemas, docs, tests, and execution behavior.
- `r_multiomics_factory` is a factory: it owns R-side code, reports, plotting helpers, and vendored bundle-contract consumers.
- `/home/zerlinshen/projects/<project-id>` is the project/data layer: it owns inputs, runs, figures, reports, manifests, evidence, and project-specific conclusions.
- The local Mac Reproduction Trail is downstream review/sync. It is not the implementation target for phase 1 remote governance.

## Human vs Agent Surfaces

Use two complementary entrypoints for every governed project/run.

Human-facing surfaces should answer: what figure matters, what conclusion is supported, and what decision it informs. Prefer concise biological names such as `figures/`, `reports/`, `conclusions.md`, and `human_summary.md`. These files should link back to run evidence but should not require reading executor logs first.

Agent-facing surfaces should answer: how was it run, can it be resumed, and can another executor audit it. Prefer explicit machine-readable names such as `manifest.json`, `python/run_manifest.json`, `python/module_status.csv`, `r/**/*manifest*.json`, `logs/`, and `.checkpoints/`.

Phase 1 defines this split in the factory contract and validator. It does not rename legacy runs or move existing project artifacts.

## Manifest Layers

Do not collapse manifest layers into one file.

1. Run-root `manifest.json` is the cross-factory run envelope at:
   `/home/zerlinshen/projects/<project-id>/runs/<run-id>/manifest.json`
2. Producer-native Python provenance remains beside Python outputs, for example:
   `python/run_manifest.json`, `python/**/run_manifest.json`, and `python/**/module_status.csv`.
3. Producer-native R provenance remains under `r/` when present.
4. Bundle provenance remains under `python/bundle/` or the producer-native bundle directory.

Legacy runs may lack one or more layers. Record that as validation `warning` or `info` when substitute provenance exists. Do not backfill legacy run files during phase 1.

## Governance Reports

Factory-side validation reports belong in control-plane locations such as:

- `ops/governance_records/<governance-round>/`
- `docs/`

They are not scientific outputs. Scientific claims still require project run evidence.

## Dirty-State Pinning

Before governance edits or validation, record:

- branch
- full commit SHA
- `git status --short`
- tracked diff hash from `git diff | sha256sum`
- untracked inventory and inventory hash from `git ls-files --others --exclude-standard | sort`

The untracked inventory must be persisted in a durable control-plane record, not only `/tmp`.

## Purpose-Driven Module Selection

Before selecting optional modules, record the scientific question, modalities present, expected human deliverables, required upstream artifacts, candidate modules, safe skip conditions, and wet-lab/manuscript decision. The machine-readable version of this checklist lives in `contracts/project_run_contract.yaml` under `purpose_driven_module_selection`.

Phase 1 adds metadata, documentation, and tests only. It must not change runtime defaults unless a later execution plan explicitly scopes that behavior change.

Use module catalog metadata to choose modules by:

- scientific question
- modalities present
- required upstream artifacts
- emitted artifacts
- manifest impact
- safe skip conditions
- wet-lab or manuscript decision needs

Running every optional module is not the governance default.

## Paper Reproduction And Context Optimization Policy

Paper reproduction is faithful-first, then context-optimized.

1. Clone or stage the upstream paper repository, scripts, and supplementary methods when available.
2. Pin repository URL, commit/tag, DOI, data accessions, license, environment files, and raw/processed input boundaries.
3. Reproduce from raw data first when public raw data exists and host capacity allows it.
4. If raw data is missing or impractical, use the earliest public computable input and label the boundary explicitly.
5. Reproduce both data objects and figures; data-object parity and figure-panel parity are separate evidence classes.
6. Compare each claim as exact, approximate, proxy, unsupported, or resource gap.
7. Map paper methods to existing modules, parameter changes, new reusable modules, or paper-specific scripts.
8. Add/update factory modules only when the method should become reusable across projects.
9. After faithful reproduction, optimize for our context: wet-lab decision, cohort scale, memory limits, available modalities, downstream hypotheses, and project-specific best parameters.

Reproduction project policies must record `upstream_repository`, `raw_data_reproduction`, `data_object_reproduction`, `figure_reproduction`, `module_gap_decisions`, and `context_optimization_decisions`.

## Project Run Retention And Final Backup Policy

Project run outputs are not a permanent warehouse. The default retention rule is
to keep the latest validated run result as the active source of truth and delete
older bulky run directories after provenance is captured. This applies to
reproduction projects and real analysis projects unless a project-specific rule
requires stricter retention.

Every project must define a project-level retention/final-backup rule that
records:

- project purpose and analysis type;
- method family and selected/best parameters;
- module list, filtering thresholds, batch/integration keys, sketch/metacell
  settings, random seed, and any paper-specific target choices;
- latest retained run id/path and superseded run ids/paths;
- final backup contents, including root `manifest.json`, producer-native
  manifests, `module_status.csv`, launch command/log, parameter file,
  environment pins, final labeled object or compact atlas, figures/reports, and
  conclusion summary;
- never-delete list: raw data, prepared canonical inputs, source code,
  environments, launch scripts, governance records, and currently cited report
  assets.

Factory governance cleanup records must capture what was deleted, why deletion
was safe, what replaced it, and where the final backup lives.

## Multi-Cohort Processing Policy

The default policy for multi-patient or multi-cohort datasets is:

1. Filter and QC per sample/patient first.
2. Merge cleaned cells, or build a balanced sketch/metacell atlas when full
   memory pressure is too high.
3. Perform cohort-level integration, clustering, and annotation with explicit
   batch keys.
4. Transfer labels back to all cells.
5. Report patient-level composition, pseudobulk, DE, and figures from the
   full labeled cohort.

This policy preserves memory headroom without sacrificing cross-cohort
annotation precision. Per-patient runs are smoke tests or QC evidence, not final
article-level annotation evidence unless confirmed by a cohort-level atlas.

## Project Validation

Use `scripts/validate_project_governance.py` for read-only validation of real projects. Phase-1 pilot projects:

- `/home/zerlinshen/projects/wave5-trevino`
- `/home/zerlinshen/projects/nc-reproduction`

The validator writes reports only when `--output-dir` is provided. Use a factory control-plane output directory, not a project run directory.
