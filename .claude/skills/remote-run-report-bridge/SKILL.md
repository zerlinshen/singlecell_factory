---
name: remote-run-report-bridge
description: Use for bioinformatics runs that finish on a remote server and need to be pulled back into the local Mac workflow, summarized, organized into local figure/report artifacts, and packaged into a human-readable report with conclusions and suggested next analyses.
---

Use this skill when a remote bioinformatics run has produced usable outputs and the next job is to turn them into a local report.

1. Prefer the remote server as the analysis source of truth.
- Read `run_manifest.json`, `module_status.csv`, and key module outputs first.
- Distinguish clearly between clean completion and partial completion with usable artifacts.

2. Choose the bridge by payload size.
- Prefer a compact CSV bundle or selected outputs when plotting is the goal.
- Only pull full `h5ad`/large objects when object-level local analysis is truly needed.

3. Use the local MacBook for reporting and organization.
- Bridge remote outputs into the local coordination workspace.
- Prefer lightweight imported figure/report artifacts over maintaining a local plotting codebase here.
- Produce a concise summary with figures, short conclusions, and suggested next analysis directions.

4. Report with evidence.
- Cite the exact remote run directory and local output directory.
- Mention module failures or caveats that affect interpretation.

5. Keep the workflow practical.
- Optimize for quick interpretation, not maximal data transfer.
- Prefer stable, reproducible report generation over fragile all-in-one conversions.

6. Close the loop.
- If the report changes the stable operating guidance:
  - `readme-sync-enforcer`
- If the report materially strengthens or weakens paper support:
  - `nsclc-baseline-reproduction`
