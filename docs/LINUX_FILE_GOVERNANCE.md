# Linux File Governance

This is the remote Linux file-governance map for the single-cell and
multiomics analysis system.

It is a control-plane document. It is not biological evidence, and it does not
modify pipeline behavior. File-governance work must not change `workflow/`,
module logic, CLI defaults, or run semantics unless a later task explicitly
requests pipeline implementation work.

## One-Screen Strategy

```text
/home/zerlinshen/
├── singlecell_factory/              # Python factory: code/contracts/docs/control plane
│   ├── workflow/                    # pipeline implementation
│   ├── contracts/                   # canonical schemas and run/project contracts
│   ├── docs/                        # factory docs and governance maps
│   ├── ops/
│   │   ├── before_every_run/        # remote run memory
│   │   ├── governance_records/      # policy, validation, cleanup records
│   │   └── run_ledger/              # factory-level execution ledger
│   ├── results/                     # legacy/smoke only for governed projects
│   └── output/                      # legacy/smoke only
├── r_multiomics_factory/            # R factory: R-native analysis modules
├── plotting_factory/                # Plotting factory: cross-language plot helpers (python/ + r/ subtrees)
├── projects-bootstrap/              # project creation helper
└── projects/                        # project/data layer
    ├── README.md                    # projects root index
    ├── FILE_GOVERNANCE.md           # project file policy
    └── <project-id>/
        ├── README.md                # project human entrypoint
        ├── project.yaml             # project metadata and factory pins
        ├── ledger/
        │   └── project_retention_policy.yaml
        ├── inputs/                  # protected input handles
        ├── configs/                 # project-specific configs
        ├── notebooks/               # exploratory notebooks
        └── runs/
            ├── README.md            # run layer rules
            └── <run-id>/
                ├── manifest.json    # cross-factory run envelope
                ├── python/          # singlecell_factory outputs
                ├── r/               # r_multiomics_factory outputs
                ├── logs/            # launch/runtime logs
                └── evidence/        # optional human reports/figures/reviews
```

## Path Roles

| Path | Role | Rule |
| --- | --- | --- |
| `/home/zerlinshen/singlecell_factory` | Python factory | Code, contracts, validators, docs, tests, control-plane records |
| `/home/zerlinshen/r_multiomics_factory` | R factory | R-native analysis modules (renamed 2026-05-18 from `multiomics_r_factory`) |
| `/home/zerlinshen/plotting_factory` | Plotting factory | Cross-language plot helpers; `python/` + `r/` subtrees; theme/schema/contracts (added 2026-05-18) |
| `/home/zerlinshen/projects-bootstrap` | Bootstrap helper | Creates governed project roots |
| `/home/zerlinshen/projects` | Scientific project layer | Inputs, configs, ledgers, runs, reports, figures, conclusions |
| `/home/zerlinshen/singlecell_factory/ops/governance_records` | Control plane | Policy, cleanup, validation, and structure decisions |
| `/home/zerlinshen/singlecell_factory/results` | legacy/smoke | No new governed project science |

## Human-Facing And Agent-Facing Split

Use `human-facing` surfaces for:

- current source-of-truth pointer
- concise reports
- final figures
- conclusion summaries
- caveats for biological or manuscript decisions

Use `agent-facing` surfaces for:

- `manifest.json`
- `python/run_manifest.json`
- `module_status.csv`
- launch logs
- command records
- validation JSON
- cleanup manifests
- run ledgers

Both are necessary. Human readers should not need logs first. Agents should not
need to infer provenance from prose.

## Current Canonical Projects

| Project | Current source of truth | Notes |
| --- | --- | --- |
| `nc-reproduction` | `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88` | Structure validation on filtered P15_T1; future article-scale work requires per-sample filtering then cohort merge |
| `wave5-trevino` | `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88` | Public-resource Cell figure reproduction; not FASTQ/fragments/BPNet exact parity |

## What This Round Executed

- Added `/home/zerlinshen/projects/README.md`.
- Added `/home/zerlinshen/projects/FILE_GOVERNANCE.md`.
- Added run-layer README files for current project roots.
- Added this factory-side governance map.
- Added a dated governance report under
  `ops/governance_records/2026-05-18-linux-file-governance-architecture/`.

## Three-Factory Round 2026-05-18

Completed Phase 0 through Phase 3 of the three-factory restructure.

**Factories established:**
- `singlecell_factory` — Python upstream analysis (this repo)
- `r_multiomics_factory` — R-native analysis (renamed from `multiomics_r_factory`)
- `plotting_factory` — dual-language visualization library (`python/` + `r/` subtrees)

**Phase summary:**
- Phase 0: repo bootstrap and skeleton for `plotting_factory`; `r_multiomics_factory` rename
- Phase 1: R plotting modules migrated to `plotting_factory/r/`
- Phase 2: Python plotting modules extracted to `plotting_factory/python/`; bridge symlinks wired
- Phase 3: cross-language theme tokens, per-plot YAML config schemas, `figure_bundle_schema.yaml`
  contracts, visual regression smoke tests, CI gate promotion to blocking, and governance finalization

**CI gate**: `scripts/ci/cross_factory_contract_gate.sh` verifies `figure_bundle_schema.yaml`
parity across all three factories; gate is blocking (non-advisory) as of Phase 3.

**ADR**: `docs/adr/` records factory-trifurcation decisions.

This round is documentation and governance only: it did not launch runs, migrate
scientific artifacts, delete data, or change source-of-truth run IDs.

## Hard Boundary

This round is documentation and governance only: do not modify pipeline. It did
not launch runs, migrate artifacts, delete data, or change source-of-truth run
IDs.
