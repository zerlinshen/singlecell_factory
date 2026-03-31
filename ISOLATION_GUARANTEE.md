# Complete Isolation Guarantee

## Principle: Code and Results Must Never Mix

This pipeline enforces **strict separation** between:
1. **Analysis Code** (`downstream_scanpy/`) - Immutable
2. **Analysis Results** (`output/<project>/`) - Isolated per project

---

## Isolation Levels

### Level 1: Directory Structure Isolation

```
singlecell_factory/
├── downstream_scanpy/     ← CODE (never touched by pipeline runs)
│   ├── run_scanpy_pipeline.py
│   └── scanpy_pipeline_optimized.py
│
└── output/                ← RESULTS (all outputs here)
    ├── project_A/        ← Completely isolated
    ├── project_B/        ← Completely isolated
    └── project_C/        ← Completely isolated
```

**Guarantee**: Running the pipeline will NEVER create/modify files in `downstream_scanpy/`

---

### Level 2: Per-Project Cache Isolation

Each project has its own cache directory:

```
output/project_A/
└── cache/
    └── cache_scanpy_10x_mtx_*.h5ad   ← Only Project A uses this

output/project_B/
└── cache/
    └── cache_scanpy_10x_mtx_*.h5ad   ← Only Project B uses this
```

**Guarantee**: Cache files cannot leak between projects. Deleting Project A's cache does not affect Project B.

---

### Level 3: Configuration Isolation

Each project stores its complete configuration:

```
output/project_A/
├── logs/pipeline_config.json      ← Project A's exact settings
└── tables/pipeline_config.csv     ← Human-readable format

output/project_B/
├── logs/pipeline_config.json      ← Project B's exact settings
└── tables/pipeline_config.csv     ← Different from Project A
```

**Guarantee**: You can always reproduce a project by examining its stored config.

---

### Level 4: Data Isolation

Processed data is strictly separated:

```
output/project_A/
├── data/after_qc.h5ad      ← Project A's QC-filtered data
└── data/processed.h5ad     ← Project A's final results

output/project_B/
├── data/after_qc.h5ad      ← Project B's QC-filtered data (different!)
└── data/processed.h5ad     ← Project B's final results (different!)
```

**Guarantee**: Project A cannot accidentally read Project B's data.

---

## Proof of Isolation

### Test 1: Code Directory Remains Clean

```bash
# Before running
ls downstream_scanpy/
# Output: run_scanpy_pipeline.py, scanpy_pipeline_optimized.py

# After running any number of projects
ls downstream_scanpy/
# Output: run_scanpy_pipeline.py, scanpy_pipeline_optimized.py (unchanged!)
```

### Test 2: Parallel Projects Don't Interfere

```bash
# Terminal 1 - Start long-running analysis
python scanpy_pipeline_optimized.py --tenx_dir=... --project=big_analysis

# Terminal 2 - Run quick analysis simultaneously
python run_scanpy_pipeline.py --tenx_dir=... --project=quick_test

# Result: Both complete successfully with no conflicts
```

### Test 3: Deletion Safety

```bash
# Delete one project
rm -rf output/old_project/

# Other projects are completely unaffected
ls output/other_project/
# All files present and valid
```

---

## Common Scenarios

### Scenario 1: Sharing Results with Collaborators

```bash
# Simply zip the entire project folder
zip -r collaborator_A.zip output/my_project/

# Collaborator unzips and gets everything:
# - Data
# - Figures
# - Tables
# - Logs (provenance)
# - Cache (optional, can be excluded)
```

### Scenario 2: Parameter Exploration

```bash
# Run 1: Default parameters
python scanpy_pipeline_optimized.py --tenx_dir=... --project=exp_default

# Run 2: More permissive QC
python scanpy_pipeline_optimized.py --tenx_dir=... --project=exp_permissive \
  --min_genes=100 --max_mito_pct=15

# Run 3: Higher resolution clustering
python scanpy_pipeline_optimized.py --tenx_dir=... --project=exp_highres \
  --leiden_resolution=1.0

# All three exist independently
ls output/
# exp_default/  exp_permissive/  exp_highres/
```

### Scenario 3: Multiple Datasets

```bash
# Dataset A
python scanpy_pipeline_optimized.py --tenx_dir=/data/pbmc1k/ --project=pbmc1k

# Dataset B (completely different)
python scanpy_pipeline_optimized.py --tenx_dir=/data/tumor/ --project=tumor_sample

# No chance of cross-contamination
```

---

## What Gets Stored Where

| Data Type | Location | Per-Project? |
|-----------|----------|--------------|
| Source code | `downstream_scanpy/` | ❌ N/A |
| Processed .h5ad files | `output/<project>/data/` | ✅ Yes |
| Figures/plots | `output/<project>/figures/` | ✅ Yes |
| CSV tables | `output/<project>/tables/` | ✅ Yes |
| Execution logs | `output/<project>/logs/` | ✅ Yes |
| Cache files | `output/<project>/cache/` | ✅ Yes |
| Config files | `output/<project>/logs/` | ✅ Yes |

---

## Best Practices

### DO ✅

- Use descriptive project names: `pbmc1k_v2_stringent`
- Archive old projects: `tar czf old.tar.gz output/old_project/`
- Share entire project folders for reproducibility
- Delete cache/ to save space (will be regenerated if needed)

### DON'T ❌

- Don't modify code in `downstream_scanpy/` during analysis
- Don't manually move files between project directories
- Don't rely on global cache directories
- Don't run multiple projects with the same name simultaneously

---

## Troubleshooting Isolation Issues

### Issue: "Cache file locked"

**Cause**: Same project name used simultaneously
**Solution**: Use unique project names for parallel runs

### Issue: "Permission denied in downstream_scanpy/"

**Cause**: Attempting to write to code directory
**Solution**: Check that `--project` parameter is provided

### Issue: Results appear in wrong project

**Cause**: Hard-coded paths in custom modifications
**Solution**: Always use relative paths and the `dirs` dictionary

---

## Summary

| Guarantee | Implementation |
|-----------|---------------|
| Code isolation | All code in `downstream_scanpy/`, never written to |
| Result isolation | Each project in `output/<project>/` |
| Cache isolation | Project-specific `cache/` subdirectories |
| Config isolation | Stored per-project in `logs/` |
| Data isolation | All outputs use project's `data/` directory |

**This pipeline is designed to be completely safe for:**
- Multiple users on shared systems
- Parallel execution of different projects
- Systematic parameter exploration
- Long-term reproducibility

---

*Last updated: 2026-03-30*
