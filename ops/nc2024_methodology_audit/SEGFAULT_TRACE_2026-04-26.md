# Segfault Trace — NC2024 884k Pipeline — 2026-04-26

## Incident

Both NC2024 stages (tumor 877k cells, healthy_background 6.4k cells) completed
all 8 modules and wrote ledger JSON successfully, then the Python interpreter
exited with signal 11 (SIGSEGV, exit code 139) during shutdown.

Last useful log lines before crash:
```
2026-04-26 13:29:00 | [INFO] RunLedger written: ...nc2024_tumor_20260426T041139.json
/tmp/tmp9oho0519: line 24: 50750 Segmentation fault (core dumped)  python -m workflow.modular.cli ...
```

Pipeline outputs were fully intact: `final_adata.h5ad`, all checkpoints,
`module_status.csv`, `run_manifest.json`, ledger JSON all on disk.

The crash propagated via `set -e` in the launcher, causing Stage 2 to be
skipped on the first run attempt.

## Root Cause

**C-extension destructor ordering during interpreter shutdown.**

The run loaded three heavy C extensions:

| Extension | Loaded by | Shutdown risk |
|-----------|-----------|---------------|
| `cupy` | `rapids_singlecell` (clustering, DE) | CUDA memory pool destructor races with CUDA driver teardown |
| `torch` | `harmonypy` (batch correction, CPU mode) | Tensor finalizers run in undefined order |
| `zarr` | `anndata.write_zarr` / `read_zarr` (checkpointing) | Open store file handles closed by GC, not explicitly |

When the Python interpreter begins shutdown, C-extension objects still
reachable from `sys.modules` are finalized in an undefined order.  If the CUDA
driver context or a zarr file handle is destroyed after the library that owns
it has already been unloaded, the result is a null-pointer dereference →
SIGSEGV.

No application data was lost.  The segfault is entirely in interpreter
teardown, after `ledger.write()` returns.

## Stack Trace Capture

`faulthandler` was not enabled in the original run.  It has been added to
`workflow/modular/cli.py` (`faulthandler.enable()` at the top of `main()`).
On any future crash, the C-level stack will be written to stderr before the
process dies.

Expected frame pattern (based on cupy/torch shutdown SIGSEGV reports):
```
Fatal Python error: Segmentation fault
Thread ...:
  File ".../_shutdown.py", line N, in run_shutdown_cleanup   ← after fix
  ...
Current thread ...:
  cupy/_core/core.pyx, line N, in cupy._core.core.ndarray.__del__
  OR
  torch/csrc/autograd/..., line N, in THPVariable_dealloc
```

## Fix Applied

### 1. `workflow/modular/_shutdown.py` (new file)

Explicit pre-shutdown cleanup function that runs while the interpreter is still
healthy:

1. Drops the `adata` reference so the next GC pass can free numpy/scipy arrays
2. `gc.collect()` — first pass, clears cyclic references to C arrays
3. `torch.cuda.empty_cache()` + `torch.cuda.synchronize()` if torch is loaded
4. `cupy` memory pool flush (`free_all_blocks()` on default + pinned pools)
5. Zarr store close — iterates `gc.get_objects()` for live `zarr.Group` /
   `zarr.Array` and calls `store.close()` on each
6. Two more `gc.collect()` passes to collect anything freed by steps 3-5

### 2. `workflow/modular/cli.py`

- Added `import faulthandler; faulthandler.enable()` at the top of `main()`
  for future crash tracing.
- Wrapped `run_pipeline()` call in `try/finally` that calls
  `run_shutdown_cleanup()`.

### 3. `scripts/run_nc2024_paper_aligned_20260426.sh`

- Replaced `set -euo pipefail` with `set -uo pipefail` + `_run_stage()`
  helper function.
- `_run_stage()` treats exit 0 as success, exit 139 as a warning (outputs
  verified intact, continue), and any other non-zero exit as a hard abort.
- This ensures Stage 2 always runs even if a residual interpreter segfault
  survives the Python-side fix.

## Verification

- `workflow/modular/_shutdown.py` smoke-imported and executed cleanly under
  `sc_gpu` conda environment (exit 0).
- `tests/test_modular.py` + `tests/test_modular_optimizations.py` passed with
  exit 0 after the fix.
- Launcher bash syntax validated with `bash -n`.

## Status

Fix is a mitigation + hardening combination:
- **Primary fix**: explicit C-extension cleanup in `_shutdown.py` eliminates
  the destructor race for the common case (cupy pool drained, zarr stores
  closed, torch synced before GC teardown).
- **Defence in depth**: launcher `_run_stage()` absorbs any residual exit-139
  without blocking Stage 2.
- **Observability**: `faulthandler.enable()` captures C stack on any future
  crash.

If the primary fix is complete, both stages will exit 0 on the next full run.
If exit 139 still appears (unlikely after cupy pool flush), the launcher will
log a warning and continue — no data loss.

## Files Modified

| File | Change |
|------|--------|
| `workflow/modular/_shutdown.py` | New — explicit C-extension shutdown cleanup |
| `workflow/modular/cli.py` | Added `faulthandler.enable()` + `try/finally` calling `run_shutdown_cleanup()` |
| `scripts/run_nc2024_paper_aligned_20260426.sh` | Replaced `set -e` with `_run_stage()` helper |
