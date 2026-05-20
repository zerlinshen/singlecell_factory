# 2026-05-19 - raw FASTQ/FASTA-reference local-to-remote smoke

## Objective

- Start a remote run from the local Mac Codex session over `ssh ubuntu-tail`.
- Validate the remote runtime coupling from raw 10x FASTQ plus a Cell Ranger reference containing `fasta/genome.fa`.
- Preserve per-run parameters/logs while enforcing the project policy: keep the latest/best validated run as canonical and archive/delete superseded bulky outputs only after preserving provenance.

## Starting Context

- Active factory repo: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory`.
- Conda env required: `/home/zerlinshen/conda/envs/sc_gpu_stable`.
- Cell Ranger: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/envs/cellranger-10.0.0/bin/cellranger`, version 10.0.0.
- FASTQ source: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/envs/cellranger-10.0.0/external/cellranger_tiny_fastq`.
- Reference: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ref/reference/refdata-gex-GRCh38-2024-A`.
- FASTA evidence: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ref/reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa`.
- Project root: `/home/zerlinshen/projects/raw-fasta-smoke`.

## What Was Run

Two local-launched remote runs were attempted:

1. `2026-05-19T0928Z-82f1964`
   - optional modules: `clustering,annotation`
   - result: process exit 0; Cell Ranger, QC, doublet detection, and clustering passed; optional annotation failed because the tiny matrix retained only 377 genes and no valid default marker genes.
   - classification: `superseded/evidence-only`.

2. `2026-05-19T0931Z-82f1964`
   - optional modules: `clustering`
   - result: clean pass for Cell Ranger, QC, doublet detection, and clustering.
   - classification: `canonical` latest/best raw-reference smoke.

Both launches recorded `logs/launch_command.txt`, `logs/params.json`, `logs/pipeline.log`, `manifest.json`, producer `run_manifest.json`, `module_status.csv`, and run ledger JSON.

## Outcome

- Local-to-remote execution path works.
- Cell Ranger preflight, FASTQ folder check, reference check, align/count, and matrix output completed.
- Clean canonical run: `/home/zerlinshen/projects/raw-fasta-smoke/runs/2026-05-19T0931Z-82f1964`.
- Clean run module status: `cellranger`, `qc`, `doublet_detection`, `clustering` all `ok`.
- Shape: raw cells `1142`, raw genes `38606`; after QC `1142` cells and `377` genes; clusters `9`.

## Retention Action

- Preserved first-run records under `/home/zerlinshen/projects/raw-fasta-smoke/ledger/run_records/2026-05-19T0928Z-82f1964`.
- Deleted superseded bulky first-run outputs:
  - `/home/zerlinshen/projects/raw-fasta-smoke/runs/2026-05-19T0928Z-82f1964`
  - `/home/zerlinshen/projects/raw-fasta-smoke/inputs/upstream/tinygex_cellranger`
- Retained latest/best generated outputs:
  - `/home/zerlinshen/projects/raw-fasta-smoke/runs/2026-05-19T0931Z-82f1964`
  - `/home/zerlinshen/projects/raw-fasta-smoke/inputs/upstream/tinygex_cellranger_clean`
- Cleanup record: `/home/zerlinshen/projects/raw-fasta-smoke/ledger/cleanup_records/2026-05-19T0933Z-latest-best-retention`.
- Protected raw/reference inputs were not deleted.

## Cautions for the Next Run

- System `python3` is not enough for pipeline execution; use `sc_gpu_stable`.
- Cell Ranger is not globally on PATH; prepend the repo-bundled `envs/cellranger-10.0.0/bin`.
- Tiny FASTQ is appropriate for operational raw/reference smoke only, not article-level biology.
- Do not force annotation on tiny smoke matrices without marker-gene coverage; record it as a resource/data-boundary gap rather than fabricating labels.
- The project policy now encodes the rule: record parameters/logs for every launched run, but keep only latest/best validated generated outputs as canonical.

## Classification

- `canonical`: `/home/zerlinshen/projects/raw-fasta-smoke/runs/2026-05-19T0931Z-82f1964`.
- `evidence-only`: archived first-run records under `ledger/run_records/2026-05-19T0928Z-82f1964`.
- `superseded`: deleted first-run bulky outputs listed above.
- `failed exploratory`: none retained as bulky run output.
