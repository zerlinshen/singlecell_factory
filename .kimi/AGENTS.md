# Kimi Entry Point: singlecell_factory

Kimi agents entering `/home/zerlinshen/singlecell_factory` must start here:

1. Read `AI_AGENT_PROTOCOL.md`.
2. Read `.kimi/skills/project-onboarding/SKILL.md`.
3. Read `ops/before_every_run/LATEST.md` before any run, rerun, recovery, or
   result-truth decision.
4. Use the `kimi-bioinformatics-agent` user skill and the read-only
   `bioinformatics-workflow` plugin when available.

This is the upstream execution and run-truth workspace. Do not claim success
without concrete artifacts such as `run_manifest.json`, `module_status.csv`,
expected AnnData outputs, logs, and downstream handoff evidence when relevant.

Mac-led SSH operation and direct remote operation are both valid. The same
protocol and verification requirements apply in either mode.

