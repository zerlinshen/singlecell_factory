# Component G — Ambient RNA Correction Module — Executor Brief

**Parent plan:** `~/.omc/plans/factory-ambient-correction-module-integration.md` (APPROVED 2026-05-20)
**Parent strategic plan:** `~/.omc/plans/nc-cell-clustering-final-strategy-plan.md` (Component G, APPROVED v4 2026-05-19)
**Created:** 2026-05-20 by Claude (Opus 4.7) under autopilot authorization
**Status:** ⏸ DRAFT — pending executor session launch (multi-day implementation work)
**Default tool:** `cellbender` (pure Python, GPU-native PyTorch); fallbacks: DecontX (R bridge, clustering-circularity), SoupX (R bridge, slower)

This brief converts the strategic ambient plan into a tactical task list an executor agent (or a fresh-context Claude session) can run. It picks the concrete integration site, the env/module file skeletons, the test fixtures, and the dry-run sequence. It does NOT execute — actual implementation requires its own session.

---

## Critical preconditions (G-AMB-0 precondition gate)

**Before any of the G-AMB-1..G-AMB-5 gates run, the executor MUST complete the G-AMB-0 raw-matrix availability audit.** This is a hard precondition.

### G-AMB-0 audit checklist (anchor datasets: NC2024 NSCLC + Cell/Trevino)

| Audit step | Target | Expected | Action if failed |
|---|---|---|---|
| 0.1 | NC2024 — locate cellranger raw output for empty-droplet estimation | `raw_feature_bc_matrix.h5` per sample OR equivalent unfiltered .mtx with all barcodes | If only filtered .mtx.gz exists (as is the case for `data/raw/nc2024_nsclc_emtab13526/` per 2026-05-20 audit), check whether EMTAB-13526 published raw_feature_bc_matrix.h5 alongside filtered. If NOT, NC2024 is **cellbender-ineligible**; consider (a) re-running cellranger on FASTQs if available, (b) using DecontX (filtered-matrix-tolerant), or (c) deferring NC2024 ambient correction with documented limitation. |
| 0.2 | Cell/Trevino — locate cellranger raw output | Same as 0.1 | Same options. |
| 0.3 | Anchor coverage threshold | 100% (both anchor datasets must have raw matrices) | If <100% anchor coverage: **invalidate cellbender as the default**, escalate to G-AMB-1 with the audit conclusion → fallback to DecontX (and document the circularity downgrade). |
| 0.4 | Anticipated future cohort sampling (Architect-v2 synthesis) | 80% of likely-near-term cohorts have raw matrices (sample 3-5 candidate datasets) | If <80%: surface as a "long-term cellbender adoption risk" but does NOT block G-AMB-1. |

**As of 2026-05-20 (Claude's pre-brief check):** NC2024 raw input at `data/raw/nc2024_nsclc_emtab13526/` is **post-cellranger filtered .mtx.gz only** (D1_1, D1_2, ..., D2_4 — barcodes/features/matrix triples). EMTAB-13526 metadata is at `.idf.txt` + `.sdrf.txt`. The executor MUST verify whether EMTAB-13526 published raw_feature_bc_matrix.h5 elsewhere (likely not — that's the common case for ArrayExpress submissions; raw outputs are often FASTQs only, which would require re-running cellranger). **This is the single most consequential finding from the pre-brief audit and likely forces a fallback decision at G-AMB-1.**

---

## Module file plan

### New files

```
singlecell_factory/
├── workflow/modular/modules/ambient_correction.py    # new module class
├── ops/envs/sc_ambient_cellbender.yaml                # new conda env spec
├── ops/policy/ambient_tool_selection_signoff.md       # G-AMB-1 signoff record
├── tests/test_ambient_correction.py                   # unit + parity tests
└── tests/parity/test_ambient_module_contract.py       # MUTATING_MODULES + AnnData round-trip
```

### Module class skeleton (`workflow/modular/modules/ambient_correction.py`)

```python
from __future__ import annotations
import logging
import subprocess
from pathlib import Path
import numpy as np
from scipy import sparse
from ..context import PipelineContext

logger = logging.getLogger(__name__)

__references__ = {
    "Fleming_cellbender_2023": {
        "title": "Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender",
        "authors": "Fleming, Marioni, Babadi",
        "journal": "Nature Methods",
        "year": "2023",
        "doi": "10.1038/s41592-023-01943-7",
        "description": "cellbender remove-background: Bayesian deep-generative ambient RNA correction.",
    },
}


class AmbientCorrectionModule:
    """Ambient RNA correction. Runs after qc, before doublet_detection.

    Per-tool dispatch: cellbender (default), DecontX, SoupX (R bridges for the latter two).
    The module receives the QC-filtered AnnData and rewrites adata.X in-place with the
    ambient-corrected matrix. Upstream obs/var/uns annotations MUST be preserved.
    """

    name = "ambient_correction"
    mutates_structure = True  # Discovered by pipeline.module_dependencies() into MUTATING_MODULES
    provides_keys = {"obs": [], "var": [], "obsm": [], "uns": ["ambient_correction"]}

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Ambient correction requires AnnData")
        cfg = ctx.cfg.ambient_correction
        tool = (cfg.tool or "cellbender").lower()
        if not cfg.enabled:
            ctx.metadata["ambient_correction"] = {
                "applied": False, "skip_reason": cfg.disable_reason or "ambient_correction.enabled=false"
            }
            logger.warning("Ambient correction SKIPPED — reason=%s", cfg.disable_reason)
            return
        # ... dispatcher: tool == "cellbender" -> self._run_cellbender(adata, cfg, ctx)
        # ...               tool == "decontx" -> self._run_decontx(adata, cfg, ctx)
        # ...               tool == "soupx"   -> self._run_soupx(adata, cfg, ctx)

    # ----------- cellbender path -----------
    def _run_cellbender(self, adata, cfg, ctx):
        # Pre: raw_matrix_path must point to raw_feature_bc_matrix.h5 (NOT filtered)
        # Step 1: serialize sample-stratified launch (cellbender runs per-sample;
        #         downstream concatenate). On 877k 4-sample cohort: ~2-4h GPU wall time.
        # Step 2: subprocess into sc_ambient_cellbender env with cellbender remove-background.
        # Step 3: read corrected .h5 ← assert shape matches adata.shape ← assert barcode order matches
        # Step 4: adata.X = corrected_X  (in-place, preserves obs/var/uns)
        # Step 5: populate adata.uns["ambient_correction"] = {full v2 manifest dict}
        # Step 6: write ctx.metadata["ambient_correction"] = {applied=True, tool="cellbender", ...}
        raise NotImplementedError("cellbender path — implement in executor session")

    # ----------- decontx path (R bridge) -----------
    def _run_decontx(self, adata, cfg, ctx):
        raise NotImplementedError("DecontX path — implement in executor session")

    # ----------- soupx path (R bridge) -----------
    def _run_soupx(self, adata, cfg, ctx):
        raise NotImplementedError("SoupX path — implement in executor session")
```

### Config dataclass addition (`workflow/modular/config.py`)

```python
@dataclass
class AmbientCorrectionConfig:
    enabled: bool = False             # Default OFF — must be explicitly enabled per project
    tool: str = "cellbender"           # cellbender | decontx | soupx
    disable_reason: str = ""           # Free-text reason when enabled=False
    raw_matrix_path: str = ""          # Path to cellranger raw_feature_bc_matrix (cellbender/SoupX)
    expected_cells: int = 0            # cellbender --expected-cells (0 = auto-estimate)
    fpr: float = 0.01                  # cellbender --fpr (false-positive rate)
    epochs: int = 150                  # cellbender --epochs
    total_droplets_included: int = 25000  # cellbender --total-droplets-included
    learning_rate: float = 1e-4        # cellbender --learning-rate
    random_seed: int = -1              # -1 means use cfg.random_state
```

Add to `PipelineConfig`:
```python
ambient_correction: AmbientCorrectionConfig = field(default_factory=AmbientCorrectionConfig)
```

### CLI flag additions (`workflow/modular/cli.py`)

```
--enable-ambient-correction       (sets cfg.ambient_correction.enabled=True)
--ambient-tool {cellbender,decontx,soupx}   (default: cellbender)
--ambient-disable-reason TEXT     (companion to disabled state; mandatory if --no-ambient-correction)
--ambient-raw-matrix-path PATH    (required for cellbender + soupx; auto-derived from sample paths if unset)
--ambient-expected-cells INT      (cellbender knob)
--ambient-fpr FLOAT
--ambient-epochs INT
--ambient-total-droplets-included INT
--ambient-learning-rate FLOAT
--ambient-random-seed INT
```

### Pipeline wiring (`workflow/modular/pipeline.py`)

Add dependency edge:
```python
MODULE_DEPENDENCIES = {
    ...,
    "qc": [...],
    "ambient_correction": ["qc"],          # NEW: runs after qc
    "doublet_detection": ["ambient_correction"] if "ambient_correction" enabled else ["qc"],
}
```

`mutates_structure = True` on the class is auto-discovered into `MUTATING_MODULES` via `module_dependencies()`; do NOT manually edit `_MUTATING_MODULES_FALLBACK`.

---

## Env build plan (G-AMB-2)

### `sc_ambient_cellbender` env (separate from `sc_gpu_stable`)

```yaml
# ops/envs/sc_ambient_cellbender.yaml
name: sc_ambient_cellbender
channels:
  - conda-forge
  - pytorch
  - nvidia
dependencies:
  - python=3.11
  - pip
  - pytorch>=2.0
  - pytorch-cuda=12.1
  - cellbender>=0.3.0
  - anndata
  - scipy
  - numpy
  - pip:
      - tables  # for cellbender HDF5 IO
```

**Why separate env:** cellbender pins PyTorch + CUDA 12.x, which may conflict with `sc_gpu_stable`'s rapids-singlecell 0.13.4 stack (RAPIDS 24.x + CUDA 12.4). Isolation prevents cross-contamination. **G-E1 invariant applies:** `sc_gpu_stable` SHA256 MUST NOT change during ambient env build. Capture before/after hashes.

**Env build verifier (G-AMB-2 smoke test):**
```bash
mamba env create -f ops/envs/sc_ambient_cellbender.yaml
/home/zerlinshen/conda/envs/sc_ambient_cellbender/bin/cellbender remove-background --help | head -10
# Expected: cellbender CLI prints usage
# Then run cellbender on a tiny synthetic 1k-cell raw matrix fixture; confirm exit 0.
```

---

## Test matrix

### Unit tests (`tests/test_ambient_correction.py`)
- `test_module_skipped_when_disabled` — `cfg.ambient_correction.enabled=False` → run() returns early, metadata records `applied: False, skip_reason`.
- `test_module_requires_raw_matrix_path_for_cellbender` — `enabled=True, tool="cellbender", raw_matrix_path=""` → run() raises `ConfigError`.
- `test_disable_reason_required_when_disabled` — CLI gate: if `--no-ambient-correction`, must pass `--ambient-disable-reason`.
- `test_random_seed_propagates` — `random_seed=42` ends up in the cellbender subprocess env.

### Contract tests (`tests/parity/test_ambient_module_contract.py`)
- `test_mutates_structure_discovered` — after `module_dependencies()`, `"ambient_correction"` in `pipeline.MUTATING_MODULES`.
- `test_dependency_edge` — `"ambient_correction"` depends on `"qc"`; `"doublet_detection"` depends on `"ambient_correction"` when enabled.
- `test_anndata_roundtrip_preserves_obs` — fixture: upstream sets `adata.obs["__sentinel__"] = "preserved"` AND `adata.var["__sentinel_var__"] = 1`. After `AmbientCorrectionModule.run()`, both sentinels survive AND `adata.X` differs from input AND `adata.shape` unchanged.
- `test_manifest_provenance_full` — `adata.uns["ambient_correction"]` contains all v2 fields: `tool, applied, raw_matrix_sha256, input_matrix_sha256, output_matrix_sha256, ambient_profile_summary, parameters{}, runtime_seconds, peak_rss_gb, random_seed, cuda_version, torch_version`.

### Integration / dry-run (G-AMB-4 + G-AMB-5)

- **G-AMB-4 (small slice, ~10-50k cells from 1-2 samples)**
  - Arm A: assert mean ambient fraction estimated >1% AND a known ambient-prone gene set's counts decrease post-correction.
  - Arm B: if Arm A fails because the library is clean, assert `ambient_profile_summary.background_fraction < 0.5%` AND full v2 parameter propagation, AND audit-log line "clean-library detected".
  - At least ONE arm must pass.
- **G-AMB-5 (full 877k tumor cohort dry run, correction only — no clustering)**
  - Record runtime, peak RSS, ambient profile summary.
  - Determinism check: rerun on a 50k subset with same `random_seed`; assert `ambient_profile_sha256` matches.
  - If mismatch: record `ambient_seed_is_deterministic: false`, surface as follow-up, does NOT block G-AMB-5 pass.

---

## Approval gate flow

```
G-AMB-0 (precondition: raw-matrix audit) — executor + human review
         │
         ▼  (if anchor coverage == 100% → cellbender admissible; else fallback to DecontX)
G-AMB-1 (tool selection sign-off) — human + Architect + Critic
         │
         ▼
G-AMB-2 (env + dependency surface) — env builds; sc_gpu_stable SHA256 unchanged; smoke passes
         │
         ▼
G-AMB-3 (module contract) — mutates_structure discovery + AnnData round-trip + dependency edges
         │
         ▼
G-AMB-4 (small-slice parity, Arm A or Arm B passes)
         │
         ▼
G-AMB-5 (full-cohort dry run + determinism check) — UNLOCKS Component D's methodology-binding subsection
         │
         ▼ (CONTINGENT only on G-AMB-4 both arms fail OR G-AMB-5 errors)
G-AMB-6 (fallback trigger: re-route to DecontX or SoupX with documented downgrade)
```

---

## Estimated effort

| Phase | Wall clock | Notes |
|---|---|---|
| G-AMB-0 audit | 1-2h | File-system + EMTAB metadata inspection + possibly contact ArrayExpress |
| G-AMB-1 sign-off | 0.5-1h | Architect + Critic review of audit conclusion |
| G-AMB-2 env build | 2-4h | cellbender + PyTorch + CUDA env build; potential dep conflict triage |
| G-AMB-3 module code + tests | 1 day | Module class + config + CLI + pipeline wiring + tests |
| G-AMB-4 small-slice run | 2-3h | 10-50k cells × 1-2 samples on GPU |
| G-AMB-5 full-cohort dry run | 4-8h | 877k cells, 4 samples, GPU; plus determinism check rerun |
| Total | **~3-4 days minimum** | Assumes G-AMB-0 audit doesn't force a fallback |

If G-AMB-0 forces fallback to DecontX (R bridge): add 1 day for R env + R bridge wiring + parity testing. SoupX adds 1 more day (slower runtime).

---

## What this brief does NOT do

- **Does NOT pick a value for `ambient_expected_cells`** — that's a per-dataset call made at executor time once the raw matrix is loaded.
- **Does NOT lock cellbender as the chosen tool** — G-AMB-0 may invalidate cellbender for NC2024 if raw matrix unavailable.
- **Does NOT begin implementation** — executor must claim this brief, build the env, and proceed gate-by-gate with human checkpoints at G-AMB-1 and (if triggered) G-AMB-6.

---

## Pre-brief audit conclusions (2026-05-20, Claude under autopilot)

1. **NC2024 raw matrix availability:** As of repo state on 2026-05-20, only `<sample>-matrix.mtx.gz` filtered triples exist (no `raw_feature_bc_matrix.h5`). The executor's first task in G-AMB-0 is to check whether EMTAB-13526 published raw matrices on ArrayExpress / supplementary archives, or whether FASTQs+cellranger re-run is needed. **High likelihood that cellbender will NOT be admissible for NC2024 without re-running cellranger from FASTQs (potentially weeks of compute).** DecontX (filtered-matrix-tolerant) is the realistic fallback for NC2024.
2. **Cell/Trevino raw matrix availability:** Not audited in this brief — executor must check.
3. **rapids-singlecell env status:** `sc_gpu_stable` has rapids-singlecell 0.13.4; unaffected by ambient module (separate env).
4. **`sc_gpu_de_test` env status:** Built but cuvs-blocked at G-E2 PARTIAL FAIL (see `ops/env_validation/gpu_de_parity_report.md`). Unrelated to ambient module.

---

**Created by:** Claude (Opus 4.7), 2026-05-20, under autopilot authorization to draft this brief. Implementation is OUT OF SCOPE for this brief's authoring session and requires a separate executor session.
