# AGENTS.md

This project is a bioinformatics analysis workspace centered on single-cell and multi-omics workflows.

## Project Scope

This workspace should support:
- scRNA-seq analysis
- spatial transcriptomics (`空间转录组 / 空转`)
- WGS analysis
- Hi-C / 3C style chromatin interaction analysis
- proteomics and cross-omics integration

The current active use case is:
- remote analysis on Ubuntu
- local plotting and interpretation on the MacBook

## Current Working Defaults

- Remote server alias: `ubuntu-tail`
- Remote analysis repo: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory`
- Remote conda execution pattern: `/home/zerlinshen/conda/bin/conda run -n sc_gpu ...`
- Local plotting project: `/Users/zerlinshen/Downloads/1. Antigravity/R pipeline for scRNA analysis`

Default architecture:
- heavy compute runs on the Ubuntu server
- local Mac work focuses on plotting, figure generation, interpretation, and light bridge scripts

## Primary Rules

1. Protect raw assets.
- Never delete or overwrite reference genomes, annotation references, or raw downloaded datasets unless explicitly requested.
- `data/raw/` and reference folders are considered high-value assets.

2. Treat results as reproducible unless told otherwise.
- Old run outputs under `results/`, `output/`, caches, and temporary exports can be cleaned when they are clearly redundant or test-only.
- Before deleting result folders, preserve any run the user explicitly identifies as important.

3. Prefer resumable and restart-safe workflows.
- Use checkpointed runs when available.
- Prefer recovery and resume over blind reruns.
- For large downloads, always prefer resumable transfer.

4. Keep remote and local roles cleanly separated.
- Remote side owns analysis, preprocessing, batch jobs, and large downloads.
- Local side owns figure production, report polishing, and bridge scripts.

5. Preserve scientific traceability.
- When summarizing results, cite concrete files such as `run_manifest.json`, `module_status.csv`, checkpoints, marker tables, and output artifacts.
- Distinguish clearly between:
  - fully completed runs
  - partially completed runs with usable outputs
  - failed runs

## Domain-Specific Guidance

### scRNA-seq

- Prefer `singlecell_factory` on the remote server for major analysis runs.
- Prefer `--checkpoint` for non-trivial jobs.
- Prefer `--gpu-mode auto`; fall back to CPU-compatible execution when GPU libraries are unstable.
- For local visualization, prefer slim handoff bundles over moving giant `.h5ad` files when plotting is the main goal.

### Spatial Transcriptomics

- Treat spatial data as larger and more visualization-heavy than standard scRNA.
- Preserve coordinate files, image assets, and spot metadata.
- Avoid destructive cleanup of image-linked assets.
- Prefer workflows that keep spatial metadata aligned with expression matrices.

### WGS

- Preserve FASTQ/BAM/CRAM/VCF and reference resources by default.
- Treat intermediate files as potentially expensive to regenerate.
- Be especially cautious with deletion because storage is large but reprocessing cost is also high.

### Hi-C / 3C

- Preserve contact matrices, valid pairs, and reference/index assets unless told otherwise.
- Keep normalization method and resolution metadata traceable.
- Avoid silent replacement of matrix outputs from different parameter settings.

### Proteomics

- Preserve raw vendor or converted spectrum files by default.
- Keep sample sheets, contrast tables, and normalization metadata paired with results.
- Distinguish exploratory outputs from publication-ready summaries.

## Operational Preferences

- Prefer concise, direct progress updates.
- Make reasonable assumptions and keep moving unless a decision is truly risky.
- Avoid asking for permission for routine operational work.
- Do not use `git` unless explicitly requested.
- Do not run destructive commands against source code or raw data.

## Cleanup Policy

Safe-to-clean targets when redundant:
- `results/`
- `output/`
- temporary handoff bundles
- cache directories
- stale failed run directories

Do-not-clean targets by default:
- `data/raw/`
- reference folders
- pipeline source code
- environment definitions

## Reporting Style

When reporting bioinformatics progress:
- state what was run
- state where outputs live
- state whether the run finished cleanly
- state the next scientifically meaningful step

When reporting storage cleanup:
- state what was preserved
- state what was deleted
- state reclaimed disk space

## Current Priority

Near-term priority for this workspace:
- maintain a robust remote-to-local single-cell workflow
- support large official dataset downloads on the remote server
- preserve raw/reference assets while aggressively trimming redundant test outputs
- keep the project ready to expand into spatial transcriptomics, WGS, Hi-C/3C, and proteomics

## Skill Routing Defaults

For bioinformatics work in this project, proactively prefer the most relevant installed skills instead of waiting for explicit invocation:
- `singlecell-remote-workflow` for remote analysis plus local plotting handoff
- `execute-and-recover-pipeline` for job recovery, resume, and large-run triage
- `develop-and-integrate-module` for pipeline code changes
- `review-and-validate-quality` before declaring workflow changes done
- `bioinformatics-autosteer` as the low-friction reminder layer for choosing the right workflow skill early
- `remote-run-report-bridge` for remote result pullback, local R plotting, and PDF-style summary packaging
- `readme-sync-enforcer` whenever pipeline behavior, fallback rules, or execution modes change
- `test-run-cleanup-policy` after official 10x test runs, benchmark reruns, or superseded experiment outputs
- `skill-usage-playbook` when choosing which installed skill should drive the next step in this project
