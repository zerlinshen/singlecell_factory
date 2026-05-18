# Paper Reproduction SOP

Control-plane document. Not a scientific output.

**Version:** 1.0 (2026-05-18, C2.1)
**Cross-references:**
- Factory routing: `ops/governance_records/2026-05-18-three-factory-trifurcation/ADR.md`
- Skill: `~/.claude/skills/paper-reproduction-from-upstream/SKILL.md` (C2.2)
- Contract schema: `contracts/project_run_contract.yaml` (C2.3)
- Agent protocol: `AI_AGENT_PROTOCOL.md` § Paper Reproduction Ladder

---

## 1. Scope

This SOP covers the end-to-end workflow for reproducing published single-cell or multiome papers using the three-factory architecture. It applies from the moment a paper is selected for reproduction through the final population of `ledger/project_retention_policy.yaml` with evidence and parameter decisions.

**In scope:**
- Forking or cloning the upstream paper repository
- Staging raw data and recording provenance
- Evaluating data object and figure parity
- Routing method gaps to the correct factory
- Recording context-tuning decisions

**Out of scope:**
- End-to-end pipeline execution (see `PROTOCOL.md` and `$execute-and-recover-pipeline`)
- Module contract changes (see `$develop-and-integrate-module`)
- R-side bundle reporting (see `$singlecell-remote-workflow`)
- Visual regression testing of plot outputs (Phase 3 of trifurcation plan)

---

## 2. When to Use

Trigger this SOP when:

- A new paper has been selected for reproduction and no `projects/<project-id>/upstream/` directory exists yet.
- Parity comparison is needed between an existing run and the upstream paper's claimed outputs.
- A reproduction project's `ledger/project_retention_policy.yaml` is missing any of the six C2 fields: `upstream_repository`, `raw_data_reproduction`, `data_object_reproduction`, `figure_reproduction`, `module_gap_decisions`, `context_optimization_decisions`.
- An agent needs to decide whether a paper method should become a reusable factory module or remain a paper-specific script.

Do **not** trigger this SOP for:
- Projects that are not manuscript reproductions (analysis_type ≠ manuscript_reproduction).
- Incremental re-runs on a project whose upstream is already pinned.

---

## 3. Fork and Clone Steps

Target path convention: `projects/<project-id>/upstream/<repo-name>/`

### 3.1 Identify the upstream repository

1. Locate the paper's code repository (GitHub, GitLab, Zenodo, or supplementary link).
2. Record the canonical URL and the paper DOI before cloning.
3. Check the license. If the license prohibits reproduction or redistribution, escalate to the user before proceeding.

### 3.2 Fork (recommended) or clone

**Option A — GitHub fork (preferred for public repos):**
```bash
# Fork via gh CLI first, then clone your fork
gh repo fork <upstream-url> --clone=false
git clone git@github.com:<your-org>/<forked-repo>.git \
    /home/zerlinshen/projects/<project-id>/upstream/<repo-name>/
cd /home/zerlinshen/projects/<project-id>/upstream/<repo-name>/
git remote add upstream <upstream-url>
git fetch upstream
```

**Option B — Direct clone (when forking is impractical):**
```bash
git clone <upstream-url> \
    /home/zerlinshen/projects/<project-id>/upstream/<repo-name>/
```

### 3.3 Pin the commit

```bash
cd /home/zerlinshen/projects/<project-id>/upstream/<repo-name>/
# Pin the exact state used for reproduction
FORK_SHA=$(git rev-parse HEAD)
FORK_DATE=$(date -u +%Y-%m-%d)
echo "Fork SHA: $FORK_SHA"
echo "Fork date: $FORK_DATE"
```

Record `$FORK_SHA` and `$FORK_DATE` immediately in the project retention policy (§6).

### 3.4 Environment staging

```bash
# Inspect environment files in the upstream repo
ls environment.yml requirements.txt renv.lock setup.py pyproject.toml 2>/dev/null
# Stage a conda/pip/renv environment separately — do not install globally
```

Do not modify the upstream clone after pinning. Treat it as read-only reference material.

---

## 4. DOI and License Recording

All provenance fields are **mandatory** before any reproduction run is launched. Record them in `ledger/project_retention_policy.yaml` under the `upstream_repository` key (see §6 for schema).

| Field | Where to find it |
|---|---|
| `url` | GitHub/GitLab URL, or Zenodo deposit URL |
| `commit_or_tag` | `git rev-parse HEAD` or release tag in the upstream repo |
| `doi` | Paper DOI from journal page or supplementary; Zenodo DOI if code-only deposit |
| `license` | `LICENSE` file in the upstream repo; use SPDX identifier (e.g., `MIT`, `CC-BY-4.0`, `GPL-3.0`) |
| `forked_at` | ISO-8601 date: `date -u +%Y-%m-%d` |
| `fork_path` | Absolute path: `/home/zerlinshen/projects/<project-id>/upstream/<repo-name>/` |

If the upstream repo has no LICENSE file, record `license: unlicensed_no_file_present` and escalate to the user before using any code.

---

## 5. Raw Data Ingest

### 5.1 Determine the earliest public computable input

Work backwards from the paper's figure panels:
1. FASTQ / fragment files (rawest; compute-intensive)
2. Per-cell-barcode count matrices or peak matrices (e.g., GEO processed matrices)
3. Author-supplied RDS / H5AD / processed tables

Prefer the earliest tier that is feasible within resource constraints. Record the chosen boundary in `raw_data_reproduction.boundary`.

### 5.2 Stage raw data

Destination: `projects/<project-id>/inputs/upstream/`

```bash
mkdir -p /home/zerlinshen/projects/<project-id>/inputs/upstream/
# Example: download GEO processed matrices
# Use $download-large-file skill for large files (>1 GB)
```

**Storage policies:**
- Do not stage FASTQ files unless the project explicitly requires raw-data reproduction and storage capacity is confirmed.
- Never place staged raw data inside the factory tree (`singlecell_factory/`, `r_multiomics_factory/`, `plotting_factory/`).
- Record the accession ID, download date, and checksum in `inputs/upstream/provenance.txt`.

### 5.3 Record raw data reproduction outcome

After staging, populate `raw_data_reproduction` in the project retention policy:
- `attempted`: `true` if any reproduction run was launched from raw data.
- `success`: `true` if the run produced outputs matching the paper's reported metrics.
- `boundary`: free-text description of what stopped reproduction if `success: false` (e.g., `"GPU OOM at module 7"`, `"missing reference genome hg38"`, `"FASTQ not publicly available"`).

---

## 6. Parity Evaluation Taxonomy

Use exactly five classes when comparing reproduced outputs to the paper's claimed outputs. Apply these consistently across both data objects and figure panels.

| Class | Definition | Example |
|---|---|---|
| `exact` | Output matches the paper's reported value within floating-point precision or pixel-identical rendering | Cell count identical; DEG p-value matches to 4 decimal places |
| `approximate` | Output is directionally correct and within an acceptable tolerance (method-defined) | Cluster proportions within ±5%; UMAP topology preserved with different random seed |
| `proxy` | A related but not identical metric was reproduced due to unavailable inputs or method differences | Pseudobulk DE reproduced instead of single-cell DE due to missing raw counts |
| `unsupported` | The paper claims a result but the evidence to reproduce it is absent from public inputs | Figure relies on internal dataset not deposited |
| `resource_gap` | Reproduction is technically feasible but blocked by compute or storage constraints | BPNet model training requires 8× A100 GPUs not available |

Assign one class per data object and per figure panel. Never leave a parity field blank — if unknown, use `unsupported` with a rationale note.

---

## 7. Module Gap Decision Routing

After parity evaluation, classify each method or analysis step that is missing or different from the factory's current capabilities. Route each gap to exactly one target.

### 7.1 Gap taxonomy

| Gap class | Definition |
|---|---|
| `exact` | Factory already produces this output with no changes |
| `approximate` | Factory output is within parity tolerance with a parameter change only |
| `proxy` | Factory produces a related output; the gap is methodological |
| `unsupported` | Gap is not addressable without new code |
| `resource_gap` | Gap is addressable in principle but blocked by compute or storage |

### 7.2 Factory routing targets

For each `unsupported` gap, select one routing target:

| Target | When to use |
|---|---|
| `singlecell_factory` | The method is a reusable Python analysis step (module contract required; use `$develop-and-integrate-module`) |
| `r_multiomics_factory` | The method is a reusable R analysis step (language-native; must not mix plotting logic) |
| `plotting_factory` | The gap is purely a visualization or theme-config issue (Python or R plot helper) |
| `paper_specific_script` | The method is one-off and not worth generalizing; keep under `projects/<id>/` |
| `parameter_change_only` | No new code needed; existing module behavior changes with config tuning |

**Decision rule (three-factory ADR):**
- If a function body returns a ggplot object or composes one → `plotting_factory`
- If a function is called from `workflow/modular/modules/` → `singlecell_factory`
- If a function is R-native analysis (not plotting) and reusable → `r_multiomics_factory`
- If a function is called only from `projects/` scripts → `paper_specific_script`
- If no new function is needed → `parameter_change_only`

See `ops/governance_records/2026-05-18-three-factory-trifurcation/ADR.md` for the full factory authority boundary description.

### 7.3 Record each gap decision

Populate `module_gap_decisions` in the project retention policy as a list:

```yaml
module_gap_decisions:
  - module: <method-name>
    gap_class: unsupported          # taxonomy class from §6
    target: singlecell_factory      # routing target from §7.2
    rationale: >-
      <Free text: why this target; what makes it reusable or not>
```

---

## 8. Context-Tuning Record Schema

After faithful reproduction, record every parameter that was tuned away from the upstream paper's settings. This schema captures the justification for divergence and the parity impact.

Each entry in `context_optimization_decisions` uses this template:

```yaml
context_optimization_decisions:
  - parameter: <parameter-name>           # e.g., n_neighbors, resolution, min_cells
    upstream_value: <upstream-default>    # value from the paper's code/methods section
    our_value: <our-chosen-value>         # value actually used in the run
    rationale: >-
      <Free text: why we diverged — e.g., dataset size difference, resource constraint,
       library version change, batch structure difference>
    parity_class: approximate             # parity class from §6 for this parameter's effect
```

**Rules:**
- Every parameter that differs from the upstream paper's stated settings must have an entry.
- `parity_class` must be drawn from the taxonomy in §6.
- Do not leave `rationale` blank — a brief sentence is sufficient.
- After parameter tuning is complete, run the parity gate again and update `data_object_reproduction` and `figure_reproduction` fields accordingly.

---

## Appendix A: Quick Reference Checklist

```
[ ] 1. Upstream repo URL, DOI, license identified
[ ] 2. Fork/clone to projects/<id>/upstream/<repo>/ and pin commit SHA + date
[ ] 3. DOI + license recorded in project_retention_policy.yaml (upstream_repository)
[ ] 4. Raw data staged to projects/<id>/inputs/upstream/ (or boundary recorded)
[ ] 5. raw_data_reproduction fields populated (attempted, success, boundary)
[ ] 6. Each data object evaluated with parity taxonomy → data_object_reproduction
[ ] 7. Each figure panel evaluated with parity taxonomy → figure_reproduction
[ ] 8. Each method gap routed to correct factory target → module_gap_decisions
[ ] 9. Each parameter divergence documented → context_optimization_decisions
[ ] 10. project_retention_policy.yaml committed to project ledger
```

## Appendix B: Related Skills and Tools

- `$paper-reproduction-from-upstream` (C2.2) — agent skill that walks through this SOP interactively
- `$develop-and-integrate-module` — use when a gap routes to `singlecell_factory` or `r_multiomics_factory`
- `$before-every-run` — run preflight journal before any reproduction pipeline launch
- `$singlecell-factory-module-delivery` — final delivery checklist when a paper method becomes a factory module
- `$nsclc-baseline-reproduction` — project-specific skill for NC2024 reproduction
- `$download-large-file` — robust download for large GEO/Zenodo files
- `$reproduce-run-retention` — cleanup and retention after a reproduction run completes
