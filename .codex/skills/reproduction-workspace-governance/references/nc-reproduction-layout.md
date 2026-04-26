# NC-Style Reproduction Layout

Use this layout for active paper reproduction workspaces where humans need fast
access to final evidence and agents need durable run memory.

## Human-Facing Layer

- `00_HUMAN_START_HERE.md`: one-page current state, most important PDFs, figure
  comparisons, fidelity/gap summaries, and active ROI topics.
- `human_review/README.md`: human navigation map.
- `human_review/figures/README.md`: final figures, paper-vs-ours boards, and
  report assets.
- `human_review/comparisons/README.md`: claim, data, fidelity, and gap
  comparisons.
- `human_review/roi/README.md`: index of regions or genes of interest.
- `human_review/roi/<TOPIC>.md`: evidence and open questions for one user-facing
  biological topic, such as `ELF3`.

## Agent-Facing Layer

- `agent_runs/README.md`: how future agents should read run history.
- `agent_runs/TEMPLATE_RUN_README.md`: required fields for each run.
- `agent_runs/YYYY-MM-DD-<slug>/RUN.md`: one meaningful run or workflow round.

Each run record should contain:
- objective
- local and remote roots
- data source
- execution mode
- commands or wrappers
- outputs
- validation evidence
- issues and fixes
- current artifact classification
- handoff notes

## Remote R Plotting Contract

For the NC2024 single-cell factory workstream, R plotting/reporting is
remote-side by default. The local Mac workspace is for review, organization,
and packaging. Do not describe local R plotting as the active bridge unless a
new verified contract explicitly restores it.
