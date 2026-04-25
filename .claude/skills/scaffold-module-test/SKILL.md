---
name: scaffold-module-test
description: Generate synthetic AnnData and isolated single-module tests to validate module contracts and outputs
triggers:
  - test module
  - scaffold test
  - module test
  - test isolation
argument-hint: "<module_name>"
---

# Scaffold Module Test

Generate synthetic test data and isolated test cases for a single pipeline module.

## When to Activate

- User creates a new module and needs tests
- User says "test module X" or "scaffold test for X"
- After `/develop-and-integrate-module` creates a new module

## Workflow

### Step 1: Analyze Module Contract

Read `workflow/modular/modules/<module_name>.py` and extract:
- `requires_keys`: what adata keys the module needs
- `provides_keys`: what the module adds
- `mutates_structure`: whether it modifies X/layers
- Optional dependencies (try/except imports)
- Constructor parameters from config

### Step 2: Generate Synthetic AnnData

Create a minimal synthetic AnnData that satisfies `requires_keys`:

```python
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

def make_synthetic_adata(n_obs=200, n_vars=500, seed=42):
    """Minimal synthetic AnnData for module testing."""
    rng = np.random.default_rng(seed)
    X = csr_matrix(rng.poisson(2, size=(n_obs, n_vars)).astype(np.float32))
    adata = AnnData(X)
    adata.var_names = [f"Gene_{i}" for i in range(n_vars)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_obs)]
    return adata
```

Then add required keys based on module's `requires_keys`:
- `obs["leiden"]` → assign random cluster labels
- `obsm["X_pca"]` → generate random PCA embeddings
- `obsm["X_umap"]` → generate random 2D UMAP
- `obs["dpt_pseudotime"]` → generate sorted random pseudotime
- `layers["counts"]` → copy of X as raw counts
- `var["mt"]`, `var["ribo"]` → random boolean masks

### Step 3: Generate Test Function

Write a test in `tests/test_<module_name>.py`:

```python
def test_<module_name>_runs_on_synthetic_data():
    """<ModuleName> module runs without error on synthetic data."""
    adata = make_synthetic_adata()
    # ... add required keys ...
    
    cfg = PipelineConfig(...)  # minimal config
    ctx = PipelineContext(cfg=cfg, run_dir=tmp_path, ...)
    ctx.adata = adata
    ctx.set_module_dir("<module_name>")
    
    mod = <ModuleClass>()
    mod.run(ctx)
    
    # Verify provides_keys are present
    for key_type, keys in mod.provides_keys.items():
        for key in keys:
            assert key in getattr(adata, key_type), f"Missing {key_type}.{key}"

def test_<module_name>_skips_on_missing_deps():
    """<ModuleName> gracefully skips when optional deps unavailable."""
    # ... test with missing optional package mock ...

def test_<module_name>_output_isolation():
    """<ModuleName> writes only to its module directory."""
    # ... verify no files written outside ctx.figure_dir / ctx.table_dir ...
```

### Step 4: Verify

Run the generated test:
```bash
pytest -q tests/test_<module_name>.py -v
```

### Step 5: Report

```
SCAFFOLD TEST REPORT
====================
Module: <module_name>
Tests generated: 3
Tests passing: X/3
Synthetic data: n_obs=200, n_vars=500
Required keys satisfied: [list]
Provides keys verified: [list]
```

## Examples

```
/scaffold-module-test clustering
/scaffold-module-test evolution
/scaffold-module-test cell_communication
```

## Notes

- Synthetic data should be minimal (200 cells, 500 genes) for fast tests
- Always use `tmp_path` fixture for output directories
- Mock optional dependencies to test fallback paths
- Add the generated test file to `tests/` and verify with `pytest`
