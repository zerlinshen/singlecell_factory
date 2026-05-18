# Paper Reproduction Ladder Policy

Created: `2026-05-18T09:48:10Z`

This governance round updates skills and repository-level operating docs so paper reproduction follows a faithful-first ladder before context optimization.

## Required Ladder

1. Clone/stage upstream paper repo, scripts, and supplementary methods when available.
2. Pin repo URL, commit/tag, DOI, data accessions, license, environment, and input boundary.
3. Run raw-data reproduction first when public raw data and host capacity allow it.
4. If raw data is missing or infeasible, use the earliest public computable input and state that boundary.
5. Reproduce both data objects and figure panels.
6. Classify claim support as exact, approximate, proxy, unsupported, or resource gap.
7. Map paper methods to existing module, parameter change, new reusable module, or paper-specific script.
8. Add/update factory modules only when the method should be reusable.
9. Optimize for our context after faithful reproduction: biological question, wet-lab decision, cohort scale, memory, modalities, and downstream hypotheses.

## Updated Surfaces

- `README.md`
- `AI_AGENT_PROTOCOL.md`
- `PROTOCOL.md`
- `docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md`
- `contracts/project_run_contract.yaml`
- `/home/zerlinshen/projects/wave5-trevino/ledger/project_retention_policy.yaml`
- `/home/zerlinshen/projects/nc-reproduction/ledger/project_retention_policy.yaml`
- local and remote Codex skills for reproduction governance, NSCLC reproduction, module delivery, and report bridging.
