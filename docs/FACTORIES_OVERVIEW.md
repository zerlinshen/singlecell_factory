# Factories Overview

_singlecell_factory + r_multiomics_factory + plotting_factory technical reference_

generated_at: 2026-05-02T13:11:04+00:00 (Three-factory trifurcation: 2026-05-18; R repo renamed under previous_name policy; new `plotting_factory` introduced)

## 1. Executive summary

`singlecell_factory` is a Python/AnnData heavy-compute pipeline organized as 33 concrete module files across the layers `ingest`, `quality_control`, `latent_structure`, `markers`, `annotation`, `state_dynamics`, `biology`, `aggregation`, `reproduction`, `spatial`, `multimodal`, `genomic_optional`, `external_validation`, `covariates`, and `reporting`. The catalog declares 29 active module specs (3 mandatory + 26 optional) with explicit dependencies, modality tags, and `bridge_ready` flags. The CLI entrypoint is `scfactory run`.

`r_multiomics_factory` is the R-side analysis workspace (trifurcation rename applied 2026-05-18; see RENAME_NOTE.md for previous_name). It consumes the bundle v2.1 contract emitted by `scripts/export_singlecell_r_bundle.py`, reads it via `R_bundle/io_bundle.R`, and routes Seurat objects to per-modality analysis modules (`R/*_module.R`: annotation, batch_integration, integration, preprocessing, protein, spatial, atac, hic, ribo, vdj, etc.). As of Phase 2 (2026-05-19), all plot helpers have been migrated out of this repo into `plotting_factory`.

`plotting_factory` (introduced 2026-05-18 / Phase 1-3) is the dual-language visualization library at `/home/zerlinshen/Bioinformatics Research Pipeline/plotting_factory/`. It owns `python/` + `r/` subtrees, shared `theme/` tokens, per-plot YAML config schemas in `schema/`, and the canonical `figure_bundle_schema.yaml` contract (vendored byte-identically into the other two repos). It contains: 5 general R plot helpers (composition_plots, dim_plots, expression_plots, qc_plots, theme_config), 6 R modality plot modules (atac_plots, hic_plots, protein_plots, ribo_plots, spatial_plots, vdj_plots), 20 Python wave5 figure renderers, and 5 visual-regression smoke tests.

Together: 97+ collected tests on the named contract surface, 4 starter recipes, 3 modality extensions (protein/spatial/multimodal_obsm), an NC2024 paper-faithful reproduction lane, plus the new C2 paper-reproduction-from-upstream pathway anchored by `docs/PAPER_REPRODUCTION_SOP.md`.

Use `scfactory run` for routine pipeline runs, the documented bundle v2.1 contract for agent-mediated handoff to R, and `--recipe nc2024_paper` for paper-faithful reproduction.

## 2. Architecture

Two side-by-side stacks bridged by the bundle v2.1 contract. Python writes a temp dir, validates manifest + SHA256 sweep, then `os.replace` swaps it into place; R reads via `io_bundle.R::read_bundle_v2`, surfaces `CLAIM_GUARD` as both a stderr message and a matrix attribute, and routes to module plotters.

### Python (singlecell_factory) layers

- **ingest** (2 modules, 1 bridge_ready): cellranger, protein_adt
- **quality_control** (2 modules, 0 bridge_ready): qc, doublet_detection
- **latent_structure** (2 modules, 1 bridge_ready): clustering, batch_correction
- **covariates** (1 modules, 0 bridge_ready): cell_cycle
- **markers** (2 modules, 2 bridge_ready): differential_expression, pseudobulk_de
- **annotation** (1 modules, 1 bridge_ready): annotation
- **state_dynamics** (4 modules, 0 bridge_ready): trajectory, pseudo_velocity, rna_velocity, cell_fate
- **genomic_optional** (2 modules, 0 bridge_ready): cnv_inference, evolution
- **biology** (6 modules, 3 bridge_ready): pathway_analysis, cell_communication, gene_regulatory_network, immune_phenotyping, tumor_microenvironment, gene_signature_scoring
- **external_validation** (1 modules, 0 bridge_ready): validate_cbioportal
- **reporting** (1 modules, 1 bridge_ready): composition
- **aggregation** (1 modules, 1 bridge_ready): metacell
- **reproduction** (1 modules, 1 bridge_ready): paper_repro
- **spatial** (2 modules, 2 bridge_ready): spatial_ingest, spatial_neighborhoods
- **multimodal** (1 modules, 1 bridge_ready): multimodal_integration

### R (r_multiomics_factory) modules

annotation_module, atac_module, batch_integration_module, cli_utils, hic_module, integration_module, io_bridge, marker_module, marker_db_module, pipeline_steps, preprocessing_module, protein_module, ribo_module, spatial_module, vdj_module. As of Phase 2 (2026-05-19), the plot halves of all six modality modules and the five general plot helpers have moved to `plotting_factory/r/` (see below).

### plotting_factory (dual-language plotting)

`r/` subtree:
- `r/core/`: theme_config (palettes, fonts, provenance footer)
- `r/general/`: composition_plots, dim_plots, expression_plots, qc_plots (formerly in `r_multiomics_factory/R/`)
- `r/modality/`: atac_plots, hic_plots, protein_plots, ribo_plots, spatial_plots, vdj_plots (extracted from the 6 modality modules in Phase 2.1)

`python/` subtree:
- `python/wave5/`: 20 figure renderers migrated from `singlecell_factory/scripts/dev/` (Phase 2.2). Standalone scripts only; nothing called from `workflow/modular/modules/`.

`schema/` subtree (Phase 3.2):
- `pca_scatter.yaml`, `umap_scatter.yaml`, `marker_dotplot.yaml`, `marker_heatmap.yaml`, `qc_violin.yaml` — each declares `required_input_columns`, `optional_grouping_keys`, `color_palette_tokens`, `figure_size_defaults`, `language_renderer_hints`.

`contracts/` subtree:
- `figure_bundle_schema.yaml` — canonical lives in `singlecell_factory/contracts/`, vendored here + into `r_multiomics_factory/contracts/`. Parity validated by `scripts/ci/cross_factory_contract_gate.sh` (blocking since Phase 3.4).

`tests/`:
- `test_visual_regression_smoke.py` — structural-hash comparison (image dimensions + non-zero-pixel ratio + dominant channel, 5% tolerance) against 5 reference PNGs in `tests/fixtures/`.

## 3. Capability matrix

| module | layer | modality | depends_on | bridge_ready | status |
|---|---|---|---|---|---|
| cellranger | ingest | singlecell_rna | - | no | stable |
| qc | quality_control | singlecell_rna | cellranger | no | stable |
| doublet_detection | quality_control | singlecell_rna | qc | no | stable |
| clustering | latent_structure | singlecell_rna | doublet_detection | yes | stable |
| cell_cycle | covariates | singlecell_rna | clustering | no | stable |
| batch_correction | latent_structure | singlecell_rna | clustering | no | stable |
| differential_expression | markers | singlecell_rna | clustering | yes | stable |
| annotation | annotation | singlecell_rna | clustering | yes | stable |
| trajectory | state_dynamics | singlecell_rna | clustering | no | stable |
| pseudo_velocity | state_dynamics | singlecell_rna | trajectory | no | stable |
| rna_velocity | state_dynamics | singlecell_rna_splicing | clustering | no | stable |
| cnv_inference | genomic_optional | copy_number | clustering | no | stable |
| pathway_analysis | biology | singlecell_rna | differential_expression | no | stable |
| cell_communication | biology | singlecell_rna | annotation | no | stable |
| gene_regulatory_network | biology | singlecell_rna | clustering | no | stable |
| validate_cbioportal | external_validation | external_cancer_genomics | differential_expression | no | stable |
| immune_phenotyping | biology | singlecell_rna | annotation | yes | stable |
| tumor_microenvironment | biology | singlecell_rna | annotation | yes | stable |
| gene_signature_scoring | biology | singlecell_rna | clustering | yes | stable |
| evolution | genomic_optional | copy_number_state_dynamics | cnv_inference, trajectory | no | stable |
| pseudobulk_de | markers | singlecell_rna | differential_expression | yes | stable |
| cell_fate | state_dynamics | singlecell_rna | trajectory | no | stable |
| composition | reporting | singlecell_rna | annotation | yes | stable |
| metacell | aggregation | singlecell_rna | clustering | yes | stable |
| paper_repro | reproduction | singlecell_rna | clustering | yes | stable |
| protein_adt | ingest | protein_adt | qc | yes | stable |
| spatial_ingest | spatial | spatial_transcriptomics | qc | yes | stable |
| spatial_neighborhoods | spatial | spatial_transcriptomics | spatial_ingest | yes | stable |
| multimodal_integration | multimodal | multimodal_joint_embedding | clustering | yes | experimental |

## 4. Bundle v2.1 contract

- Schema version string: `singlecell_r_bundle_v2.1`. Back-compat literal: `singlecell_r_bundle_v2`. The exporter's `ExportConfig.schema_version` accepts `v1`, `v2`, or `v2.1`.
- Top-level manifest fields: `schema_version`, `source`, `bundle`, `cell_alignment`, `expression`, `files`, `extensions`.
- Known extension keys (`KNOWN_EXTENSION_KEYS`): `protein`, `spatial`, `multimodal_obsm`.
    - `protein`: CLR-normalized ADT matrix (`protein.parquet`); fields include `normalization`, optional `isotype_controls`.
    - `spatial`: (x,y) coordinates + library metadata (`spatial.parquet`); image data is STRING-pointer only, never serialized.
    - `multimodal_obsm` (EXPERIMENTAL): one parquet per joint embedding (`X_wnn`, `X_mofa`); same plotting-only `claim_guard`.

- `CLAIM_GUARD` literal: `not_for_de_or_new_quantitative_claims_without_full_object_validation`. Surfaced in R as a `message()` call AND as `attr(expr_sparse, "claim_guard")`. Meaning: bundle data is plotting-only; do not draw new quantitative claims from it without full-object validation.
- Atomicity: writers stage every file under a `TemporaryDirectory`; SHA256 is computed for each, recorded in the manifest, and cross-checked against on-disk byte sizes. The manifest is written last, then `os.replace(temp_dir, final_output_dir)` performs the atomic swap; any pre-existing target directory is moved to a sibling backup and restored on failure.

## 5. Workflow walkthroughs

### nc2024_paper

_NC2024 NSCLC paper-faithful baseline reproduction_  

**Command**: `python scripts/scfactory.py run --recipe nc2024_paper <input>`

Optional modules:
- clustering
- differential_expression
- annotation
- trajectory
- paper_repro

Env:
- `SC_CLUSTERING_ENGINE`=sparse_exact

Bundle config:
- `enabled`: False
- `include_protein`: False
- `include_spatial`: False
- `include_multimodal_obsm`: False
- `schema_version`: v2.1

When to use:

> Honors NC2024 paper-aligned parameters.
> Uses sparse-exact clustering (Phase 7+ default for ~300k+ cohorts).
> Bundle is OFF by default — paper_repro emits its own parity artifacts;
> pass --bundle explicitly if you also want the v2.1 R bundle.

### quick_explore

_Fast first-look exploration: clustering + DE only, bundle on_  

**Command**: `python scripts/scfactory.py run --recipe quick_explore <input>`

Optional modules:
- clustering
- differential_expression

Bundle config:
- `enabled`: True
- `include_protein`: False
- `include_spatial`: False
- `include_multimodal_obsm`: False
- `schema_version`: v2.1

When to use:

> Most common use case for new users on a single sample.
> Skips trajectory/annotation to keep wall-time low; rerun with another
> recipe (e.g. nc2024_paper) once you've inspected the bundle.

### cite_seq_full

_CITE-seq full analysis: RNA + ADT (CLR) + bundle with protein extension_  

**Command**: `python scripts/scfactory.py run --recipe cite_seq_full <input>`

Optional modules:
- clustering
- differential_expression
- annotation
- trajectory
- protein_adt

Bundle config:
- `enabled`: True
- `include_protein`: True
- `include_spatial`: False
- `include_multimodal_obsm`: False
- `schema_version`: v2.1

When to use:

> Pairs CITE-seq RNA + protein (ADT) with CLR-normalized protein outputs
> in the v2.1 bundle (protein extension on).
> Inputs may be a CITE-seq .h5ad with obsm['protein_counts'] or a
> cellranger sample-root with an Antibody Capture modality in the H5.

### visium_neighborhoods

_Visium spatial transcriptomics: ingest + neighborhood analytics + spatial bundle_  

**Command**: `python scripts/scfactory.py run --recipe visium_neighborhoods <input>`

Optional modules:
- clustering
- differential_expression
- annotation
- spatial_ingest
- spatial_neighborhoods

Bundle config:
- `enabled`: True
- `include_protein`: False
- `include_spatial`: True
- `include_multimodal_obsm`: False
- `schema_version`: v2.1

When to use:

> Visium / generic spatial workflow. spatial_ingest attaches (x,y) coords +
> library metadata; spatial_neighborhoods computes Moran's I, neighborhood
> enrichment, and co-occurrence.
> Requires squidpy in the active environment for full neighborhood
> analytics; without it, spatial_neighborhoods will degrade to a no-op
> with a warning. The v2.1 bundle's spatial extension is enabled.

## 6. scfactory CLI reference

### scfactory

```
usage: scfactory [-h] {run,doctor,report} ...

User-friendly wrapper around the modular pipeline.

positional arguments:
  {run,doctor,report}
    run                Auto-detect modality and dispatch to
                       workflow.modular.cli
    doctor             Read-only environment health check
    report             Render a single self-contained HTML report from a run
                       directory

options:
  -h, --help           show this help message and exit
```

### scfactory run

```
usage: scfactory run [-h] [--project PROJECT] [--out OUT] [--bundle]
                     [--bundle-out BUNDLE_OUT]
                     [--optional-modules OPTIONAL_MODULES] [--recipe RECIPE]
                     [--list-recipes] [--dry-run]
                     [input]

Run the modular pipeline on an .h5ad or sample-root.

Module precedence (highest wins):
  1. --optional-modules <list>   (escape hatch; overrides all)
  2. --recipe <name>             (preset modules + env + bundle)
  3. auto-detected modality      (fills any remaining gap)

--bundle and recipe.bundle.enabled are OR-ed: passing --bundle always forces bundle on even if the recipe disables it.

positional arguments:
  input                 .h5ad file or sample-root directory (omit only with
                        --list-recipes)

options:
  -h, --help            show this help message and exit
  --project PROJECT     Run name (default: scfactory_<modality>)
  --out OUT             Output directory (default: <repo>/results)
  --bundle              After run, export an R bundle via
                        scripts/export_singlecell_r_bundle.py (forces bundle
                        on even if recipe disables it)
  --bundle-out BUNDLE_OUT
                        Bundle output path (default: <out>/<project>/r_bundle)
  --optional-modules OPTIONAL_MODULES
                        Override auto-detected optional module list (comma-
                        separated; passthrough escape hatch — wins over
                        --recipe)
  --recipe RECIPE       Apply a preset recipe from recipes/<name>.yaml
                        (modules + scale-mode + env + bundle config)
  --list-recipes        List available recipes (one per line) and exit 0
  --dry-run             Print what would run, do not execute
```

### scfactory doctor

```
usage: scfactory doctor [-h] [--json]

options:
  -h, --help  show this help message and exit
  --json      Emit machine-readable JSON
```

### scfactory report

```
usage: scfactory report [-h] [--out OUT] [--title TITLE] run_dir

positional arguments:
  run_dir        Path to a run directory containing run_manifest.json

options:
  -h, --help     show this help message and exit
  --out OUT      Output HTML path (default: <run_dir>/report.html)
  --title TITLE  Report title (default: project name from manifest)
```

## 7. NC2024 claim coverage

Source: `r_multiomics_factory/docs/NSCLC_CLAIM_VALIDATION_MATRIX.md`. 15 claims parsed: 4 direct, 7 partial, 4 unsupported.

| claim_id | family | module_family | lane | status | evidence |
|---|---|---|---|---|---|
| NC2024-F1-01 | broad cell-type atlas | core scRNA analysis | 40k_fidelity | direct | 40k clean run + annotation outputs |
| NC2024-F1-02 | tumour/background composition shift | composition | 40k_fidelity | direct | 40k clean run composition outputs |
| NC2024-F1-03 | immune composition shift | immune/TME | 40k_fidelity | direct | 40k clean immune outputs |
| NC2024-F1-04 | Treg/NK/exhaustion pattern | immune/TME | 40k_fidelity | direct | 40k clean immune outputs |
| NC2024-F1-05 | CAML-like hybrid state | scRNA-core / tumour-state | 40k_fidelity | partial | 40k marker/annotation outputs |
| NC2024-F2-01 | abundance correlation | composition | 40k_fidelity | partial | 40k composition outputs |
| NC2024-F2-02 | communication differences | cell-cell communication | none | unsupported | no current communication run |
| NC2024-F2-03 | VEGF/EGFR interaction pattern | cell-cell communication | none | unsupported | no current communication run |
| NC2024-F2-04 | LUSC-specific checkpoint wiring | cell-cell communication | none | unsupported | no current communication run |
| NC2024-F3-01 | macrophage state diversity | immune/TME | 900k_survivability | partial | 900k proxy v3 immune + TME outputs |
| NC2024-F3-02 | foetal-like macrophage programme | pathway/regulatory | none | unsupported | current runs lack dedicated pathway/regulatory analysis for this claim |
| NC2024-F4-01 | spatial localization | spatial | bundle_v2.1 | partial | AnnData -> bundle v2.1 `spatial` extension -> `R/spatial_module.R::plot_spatial... |
| NC2024-F4-02 | spatial interaction context | spatial + communication | bundle_v2.1 | partial | bundle v2.1 `spatial` extension + upstream `spatial_neighborhoods.py` (squidpy-... |
| NC2024-ABS-02 | macrophage vs cytotoxic inverse relationship | immune/TME | 900k_survivability | partial | 900k immune subtype summary + composition + TME |
| NC2024-ABS-03 | similar composition but subtype-specific divergence | subtype comparison | 40k_fidelity | partial | 40k clean + 40k benchmark + docs |

## 8. Test coverage and quality gates

Total collected on the named contract surface: **97** tests.

| file | tests | protects against |
|---|---:|---|
| `tests/test_r_bundle_contract.py` | 16 | Python <-> R bundle v2/v2.1 round-trip subprocess contract; guards CLAIM_GUARD surfacing in R. |
| `tests/test_python_r_parity.py` | 0 | Cross-language sentinels (scanpy vs Seurat) — fires when defaults silently diverge. |
| `tests/test_modular.py` | 54 | Modular pipeline contracts: dependency resolution, ctx.set_module_dir, manifest/status outputs. |
| `tests/test_scfactory.py` | 16 | scfactory CLI wrapper: auto-detection planner, recipe precedence, doctor JSON shape. |
| `tests/test_singlecell_r_bundle_export.py` | 5 | Bundle export internals: schema versions, atomic publish, extension manifest entries. |
| `tests/test_densify_audit.py` | 2 | Grep-ban: every .toarray()/.todense() in modules must carry a densify-allowed comment. |
| `tests/test_module_catalog.py` | 4 | Module catalog stability: mandatory chain, default optional list, layer coverage. |

## 9. Known limitations and roadmap

- `multimodal_integration` is **EXPERIMENTAL**: off by default, gated by `SC_MULTIMODAL_ENGINE`. Emits the bundle v2.1 `multimodal_obsm` extension only when explicitly enabled.
- Tissue images in the spatial extension are **STRING-pointer only**. Image bytes are never serialized into the bundle; readers must resolve paths themselves.
- Two bundle-v1 schema tests are pinned to the v1 literal for back-compat coverage and are intentionally NOT removed.
- No cell-cell communication / pathway scoring / RNA-velocity / trajectory parity checks against R yet — listed as roadmap.
- `--include-multimodal-obsm` is opt-in only and is never auto-passed even on the multimodal modality.
- A pre-existing repo-doc-sync hook on `singlecell_factory` blocks commits until citation/doc references are staged.
- Follow-up: a dedicated "How to add a new module" extension guide is on the roadmap but not part of this overview.

## Appendix: file index

| path | summary |
|---|---|
| `r_multiomics_factory/AGENTS.md` | Start with `AI_AGENT_PROTOCOL.md` for onboarding, architecture, and task |
| `r_multiomics_factory/AI_AGENT_PROTOCOL.md` | Canonical onboarding index for AI agents entering |
| `r_multiomics_factory/CLAUDE.md` | Claude Code must start with `AI_AGENT_PROTOCOL.md` before editing this |
| `r_multiomics_factory/CODEX_PROFILE.md` | Codex must start with `AI_AGENT_PROTOCOL.md` for onboarding, architecture, and |
| `r_multiomics_factory/MODULE_UPDATE_SKILL.md` | Use this checklist every time a new module is added or an existing module is changed. |
| `r_multiomics_factory/README.md` | AI agents must start with [AI_AGENT_PROTOCOL.md](AI_AGENT_PROTOCOL.md). That |
| `r_multiomics_factory/SKILL_PLAYBOOK.md` | Start with `AI_AGENT_PROTOCOL.md` before choosing a skill. This playbook routes |
| `r_multiomics_factory/docs/ARTICLE_CLAIM_INVENTORY_NSCLC_NATCOMM_2024.md` | Baseline article: |
| `r_multiomics_factory/docs/MODULE_CAPABILITY_MAP.md` | `singlecell_factory` is the compute and analysis engine. `r_multiomics_factory` is the downstream R/report wo... |
| `r_multiomics_factory/docs/NC2024_NSCLC_BASELINE_REPRODUCTION_MEMO.md` | - Journal: `Nature Communications` |
| `r_multiomics_factory/docs/NSCLC_BASELINE_REPRO_TASK_LIST.md` | This task list converts the baseline article into reproducible work items using the approved evidence schema. |
| `r_multiomics_factory/docs/NSCLC_CLAIM_VALIDATION_MATRIX.md` | This matrix is the claim authority table for the current baseline article. |
| `r_multiomics_factory/docs/NSCLC_REPRODUCTION_COMPARISON_MEMO.md` | This memo compares the three current NSCLC lanes and states what each lane proves. |
| `r_multiomics_factory/docs/PYTHON_R_PROVENANCE_AND_FIDELITY.md` | This document defines how to interpret Python-generated figures and R-generated figures in this workspace, an... |
| `singlecell_factory/AGENTS.md` | This file applies to `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory` and all subdirectories. |
| `singlecell_factory/AI_AGENT_PROTOCOL.md` | Canonical onboarding index for AI agents entering |
| `singlecell_factory/BEST_PRACTICES.md` | - Codex 项目规则：`AGENTS.md`（本项目根目录） |
| `singlecell_factory/CLAUDE.md` | Start with `AI_AGENT_PROTOCOL.md` for onboarding, read order, and task routing. |
| `singlecell_factory/CODEX.md` | Codex must start with `AI_AGENT_PROTOCOL.md` for onboarding, read order, and |
| `singlecell_factory/PROTOCOL.md` | AI agents should start with `AI_AGENT_PROTOCOL.md` before using this file. This |
| `singlecell_factory/README.md` | AI agents must start with [AI_AGENT_PROTOCOL.md](AI_AGENT_PROTOCOL.md). That |
| `singlecell_factory/docs/FACTORIES_OVERVIEW.md` | _singlecell_factory + r_multiomics_factory technical reference_ |
| `singlecell_factory/docs/MODULE_TECH_DOC_TEMPLATE.md` | - Module name: |
| `singlecell_factory/docs/PUBLICATION_READY.md` | This file provides copy-paste-ready Methods section text and citation patterns for manuscripts using this pip... |

