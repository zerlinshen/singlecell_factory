# Paper reproduction ladder policy and skill sync — 2026-05-18

- Date/time UTC: `2026-05-18T09:48:31Z`
- Objective: Persist the new reproduction operating rule in skills, repository-level docs, project policies, and run memory.
- What changed: Reproduction now follows `clone/stage upstream repo -> pin commit/data/license -> raw-data reproduction when possible -> earliest public computable input when raw data is unavailable -> reproduce data objects and figures -> classify claim parity/resource gaps -> module-gap analysis -> add/update reusable modules -> context optimization`.
- Canonical governance record: `/home/zerlinshen/singlecell_factory/ops/governance_records/2026-05-18-paper-reproduction-ladder-policy/`.
- Updated repository surfaces: `README.md`, `AI_AGENT_PROTOCOL.md`, `PROTOCOL.md`, `docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md`, `contracts/project_run_contract.yaml`.
- Updated project policies: `/home/zerlinshen/projects/wave5-trevino/ledger/project_retention_policy.yaml` and `/home/zerlinshen/projects/nc-reproduction/ledger/project_retention_policy.yaml` now include upstream repo/raw-data/data-object/figure/module-gap/context-optimization fields.
- Updated skills locally and remotely: `reproduction-workspace-governance`, `nsclc-baseline-reproduction`, `develop-and-integrate-module`, `singlecell-factory-module-delivery`, and `remote-run-report-bridge`.
- Artifact classification: policy/docs are `canonical`; no pipeline run was launched; no data was deleted.
- Next operator reminder: do not start future paper reproductions directly from our preferred pipeline if an upstream repo/raw-data path exists. Faithful reproduction comes first, context optimization second.
