---
name: singlecell-remote-workflow
description: Run and bridge the single-cell workflow between the Ubuntu server and the local Mac coordination workspace. Use when working on singlecell_factory runs, large 10x/Chen Lab downloads, remote result triage, or handoff from AnnData outputs to local organized figure/report artifacts. Before any remote execution on `ubuntu-tail`, first use `before-every-run` to read the latest run memory and carry forward the latest lessons.
---

Use this skill for the end-to-end single-cell workflow on this machine.

0. Before any remote execution, invoke `before-every-run`.
- Read:
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run/LATEST.md`
  - `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/LATEST.md`
- State the execution mode before running anything:
  - `execution_mode = debug_massive`
  - `execution_mode = controller_validation`
  - `execution_mode = benchmark`

1. Treat the Ubuntu server as the analysis engine.
- SSH alias: `ubuntu-tail`
- Remote repo: `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory`
- Conda env: `sc_gpu`
- Prefer `python -m workflow.modular.cli ... --checkpoint --gpu-mode auto`
- If GPU compatibility is unstable, keep results and rerun fragile modules with CPU fallback.

2. Treat the Mac as the organization and interpretation layer.
- Local coordination workspace: `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail`
- Prefer small handoff bundles over copying giant `.h5ad` files when the goal is organizing already-processed outputs.
- Keep local work focused on artifact organization, comparison boards, handoff docs, and result packaging.

3. For large official downloads.
- Prefer resumable downloads and keep partial files for resume.
- Store remote raw data under `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/data/raw/`
- For the current project, monitor the WCH run plus the 10x NSCLC download queues before starting duplicate work.

4. For WCH / LUSC work.
- Prepared WCH input lives at `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/data/raw/wch_lung_cancer_atlas/prepared_input.h5ad`
- Best stable run so far: `WCH_LUSC_STAGE_IIII_CPU_STABLE_20260420_213646`
- Read `run_manifest.json`, `module_status.csv`, and module subdirectories before deciding whether to rerun.

5. For NC2024 full-cohort work, use these defaults.
- `execution_mode = debug_massive`
  - run direct `massive`
  - do not rerun prepare if canonical prepared input is still valid
  - stop at the first failing module
  - patch only that module
  - rerun direct `massive`
- `execution_mode = controller_validation`
  - use the official `large -> massive` orchestration path
  - treat `large` as a capacity probe, not the main completion lane
- `execution_mode = benchmark`
  - use only for controlled comparisons of fidelity/performance lanes

6. For local handoff.
- Prefer CSV bundle handoff when full `h5ad` transfer is too large.
- Current compact-bundle source pattern remains useful: `obs.csv`, `X_pca.csv`, `X_umap.csv`, `marker_expr.csv`
- Keep the local side focused on organizing processed outputs rather than maintaining a local plotting pipeline.

7. After a meaningful remote run, write back the outcome.
- Update:
  - `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/journal/`
  - `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/ops/before_every_run/LATEST.md`
- Refresh the local mirror under:
  - `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/before-every-run/`
- Classify artifacts explicitly:
  - `canonical`
  - `evidence-only`
  - `superseded`
  - `failed exploratory`

8. When reporting status.
- Separate `finished with usable outputs` from `cleanly completed`.
- Always cite concrete files, checkpoint names, or module statuses.

9. Hand off to the right next skill.
- After successful remote execution that needs local packaging:
  - `remote-run-report-bridge`
- After a workflow/result should be treated as stable:
  - `review-and-validate-quality`
- After repeated reruns or superseded experiments:
  - `test-run-cleanup-policy`
