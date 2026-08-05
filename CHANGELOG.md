# Changelog

Notable changes to `singlecell_factory`. Most recent first.

## 2026-08-05 — Wave-3 cell-cycle pre-layout order + spatial/comm evidence

- **cell_cycle option A**: depends on doublet_detection; clustering
  `runs_after` cell_cycle so regress_out can affect embeddings when both
  modules are requested. Post-layout regress still uses Wave-2 exploratory stamps.
- **C4 Visium**: squidpy installed; Moran's I 8/8 markers I>0.2 on 10x mouse brain.
- **Communication**: bounded LUSC LIANA lane (claimable engine) under audit project;
  NC2024-F2 rows remain partial (method lane, not full paper multi-condition).

## 2026-08-05 — Wave-2 operator safety + secondary claim stamps

Follow-up to Wave-1 claim honesty (`806a384`). Plan:
`projects/pipeline-scientific-audit-20260805/reports/03_remediation/WAVE2_RALPLAN.md`.

- **scfactory**: recipe/module preflight for required packages (squidpy);
  `visium_neighborhoods` declares `required_python_packages`; doctor
  `readiness.spatial_analytics`.
- **metacell**: SEACells vs MiniBatchKMeans engine/claim stamps; exploratory
  UMAP title on fallback.
- **cell_cycle**: design C — `cell_cycle_regressed_after_layout` + force
  clustering claim exploratory when regress runs after layout.
- **config/cli**: default `--tissue` / `tissue` is `unspecified` (not lung).
- **marker_db_loader**: does not invent lung when tissue unspecified.
- Tests: metacell/cell_cycle/pathway depth + visium preflight.

## 2026-08-05 — Wave-1 scientific claim honesty (audit remediation)

Governed scientific audit project:
`/home/zerlinshen/projects/pipeline-scientific-audit-20260805/`
(master report under `reports/00_master/`; consensus plan under
`reports/03_remediation/`).

- **pathway_analysis**: stop labeling `gseapy.enrich` (ORA) as confirmatory
  rank-based GSEA; prefer `prerank`; stamp `gseapy_ora` / `gseapy_mixed`
  distinctly (`pathway_claimable` false when mixed).
- **annotation**: `DEFAULT_MARKERS` tumor/immune starter pack is
  non-claimable; tissue-mismatch warnings; exploratory UMAP titles.
- **trajectory**: rebuild neighbor graph when preferring `X_wnn` so PAGA/DPT
  geometry matches the advertised embedding.
- **integration_select**: add `__references__` (integration benchmark +
  Harmony / scVI).
- **scfactory run**: forward `--allow-dirty` to the modular CLI.
- Tests: `test_pathway_claim_honesty.py`, `test_annotation_claim_honesty.py`,
  trajectory contract updates.

## 2026-08-05 — First full sync to GitHub (main: `82f1f6f..d192836`)

This repository's full working history was pushed to
`github.com:zerlinshen/singlecell_factory` on 2026-08-05 after a proxy-layer
network fault was resolved (github.com traffic routes through the local Clash
Verge `MESL` proxy group; see `~/.grok/GITHUB_PUSH_GUIDE.md`).

### main — wave6 modular factory hardening (28 commits, 329 files, +45k/-1.2k)

Sparse / memory-safe execution engines (opt-in via env flags):

- `SC_DE_ENGINE=sparse` — sparse Welch differential expression with float64 parity (`ce634d0`)
- `SC_CELLCOMM_ENGINE=sparse` — one-pass groupby-mean cell communication (`16a5459`)
- `SC_CNV_ENGINE=chunked` — sparse-aware chunked CNV smoothing (`d5068b6`)
- `SC_CLUSTERING_ENGINE=sparse_exact` — disables CSS routing for exact clustering (`5dc2a65`)

Memory governance:

- Observational memory guard (`SC_MEM_GUARD=on`) across DE / CNV / cell_comm /
  clustering / pseudobulk_de / metacell (`32f73a4`)
- Pre-flight budget check + watchdog cooperative abort, preventing kernel
  OOM-kill (`acac170`)
- Global densify policy with CI grep-ban and annotated allowlist (`cb2a5cd`);
  `plan_densify` guards on legacy dense paths in DE / CNV / evolution modules
  (`92fa664`, `0406652`, `6102127`)

Bundle v2 (Python <-> R interchange):

- v2 schema (parquet + mtx) with Mac-side asset aggregation pipeline
  (`b602abc`, `e1a4e37`)
- Round-trip subprocess tests preventing silent contract breakage (`357c34c`)

NC2024 paper alignment:

- Paper-faithful parameters + tumor/B+H cohort subset flag + new launch
  script (`b560c36`); methodology audit (`66cc909`); sparse-exact launcher
  (`b602abc`)

Testing / ops infrastructure:

- 4-tier staircase test fixtures (nano / small_real / medium_real / full_real)
  with NC2024 subset generator (`c341431`)
- Opt-in performance regression gate with 5k synthetic baseline (`b91649b`)
- Run audit-trail ledger capturing CLI / env / git / RSS / parity per run
  (`d64489c`)
- Scale mode decomposed into 4 independent capability flags (`b9ef446`)

### Branch: `codex/publication-baseline-cleanup-20260803`

Publication-baseline cleanup line on top of main: scientific-inference input
contract enforcement (`a5d2abd`), batch-confounded clustering can no longer be
recorded as claimable (`5b18f88`), memory-ratchet gate fails fast (`5cc0b9c`),
NG2025 live pointers retired (`dcd1172`), plus validation/test isolation fixes.

### Branches: `wave5-v4.2-execution`, `wave5-v5.0-figures`

Wave5 canonical-run ledger work: v3 + v4.2 ledger schemas, CI gates, lockfile
pre-commit hook, execution journal (`wave5-v4.2-execution`); figure-parity
rendering with 3-layer comparison PDF and schema isolation (`wave5-v5.0-figures`).
