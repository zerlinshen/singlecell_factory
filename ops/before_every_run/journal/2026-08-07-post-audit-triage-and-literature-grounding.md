# Post-audit triage + literature grounding — 2026-08-07

## Scope

- Task: apply `execute-and-recover-pipeline` (triage lane) and `literature-grounding`
  after the 2026-08-05 scientific audit; verify Wave-1 fixes and documented
  residuals against primary literature via the academic-search MCP.
- Execution mode: **triage_only** — no pipeline run launched. Justification: no
  failed/incomplete run exists in the 14-day window, and all deferred items are
  explicit owner decisions (see below).
- Companion deliverable:
  `/home/zerlinshen/projects/pipeline-scientific-audit-20260805/reports/03_remediation/LITERATURE_VALIDATION_2026-08-07.md`

## Triage findings (evidence-based)

- `pipeline-scientific-audit-20260805`: 8/8 runs `complete` (latest
  `2026-08-05T1459Z-4e5392e`, lane `W9_nsclc_visium_spatial_liana`, Visium
  E-MTAB-13530).
- `lusc-gt-concordance-20260728`: 3 failed runs on 07-28 morning
  (cellranger→differential_expression chain) were superseded the same day by 7
  consecutive `complete` runs from `2026-07-28T1342Z` onward — nothing to recover.
- `claim-integrity-validation-20260802/r5`: targeted engine-marking probe, not a
  production failure — honesty stamps behaved as designed
  (`grn_claimable=false` exploratory fallback; cell_communication loud error on
  Ensembl-indexed axis).
- `pipeline-publication-baselines`: idle (only retired run dir).
- Latest run memory reviewed: journal `2026-05-27-q22-lusc-integration-gate-gpu.md`
  (Q22 PASS, harmony gate on 92,430-cell LUSC cohort) and
  `ops/run_ledger/wave5_trevino_20260517T1005Z-7b539a5.v5.1.json` (WNN fix APPROVED).

## Literature validation outcome

5/5 claims **SUPPORTED**, no Wave-1 decision needs revision:

1. ORA≠GSEA prerank boundary — Subramanian 2005 PNAS `10.1073/pnas.0506580102`; GSEApy Fang 2023 `10.1093/bioinformatics/btac757`.
2. WNN neighbors rebuild — Hao 2021 **Cell** `10.1016/j.cell.2021.04.048` (bibliographic correction: not Nat Biotechnol).
3. Ambient absolute-% residual — SoupX `10.1093/gigascience/giaa151`; CellBender Nat Methods 2023 `10.1038/s41592-023-01943-7`; FastCAR PMID 38030970 ("highly sample-specific").
4. Doublet recall residual — Xi & Li Cell Systems 2021 `10.1016/j.cels.2020.11.008`; DoubletFinder `10.1016/j.cels.2019.03.003`; scDblFinder `10.12688/f1000research.73600.1`.
5. NSCLC starter pack off-lung fail-honest — CellTypist Science 2022 `10.1126/science.abl5197`; SingleR Nat Immunol 2019 `10.1038/s41590-018-0276-y`.

Primary anchors verified via `academic-search:get_paper_by_id`; no hand-written
citation fields.

## Deferred — owner decisions (not acted on)

- NC2024 LIANA communication lane; full spatial env install; 900k re-runs.
- Wave-2: P1-6 doctor spatial/squidpy preflight; P2 metacell engine stamp;
  P2 cell_cycle DAG order.

## Housekeeping

- No factory code changed; no run outputs touched; no cleanup candidates beyond
  what prior journals already recorded.
- MCP note: `biomcp`/`playwright`/`context7` installed 2026-08-06 in
  `~/.kimi-code/mcp.json`; they load after session restart and can corroborate
  this validation via PubTator3/Europe PMC federation on a future pass.
