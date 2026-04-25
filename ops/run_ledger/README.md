# ops/run_ledger

Each file in this directory is a JSON audit record written by `RunLedger` at the end of every modular pipeline run. Records use schema version `run_ledger_v1` and capture CLI arguments, filtered environment variables (SC_*, SCF_*, CUDA_*, CONDA_* prefixes plus a small standard allowlist), git SHA/branch/dirty flag, conda environment, Python version, hostname, sample root, optional modules requested, per-module wall time and RSS peak, total pipeline wall time, process-level RSS peak sampled at 1-second intervals via psutil, and a SHA-256 checksum of `final_adata.h5ad`. Example (abbreviated):

```json
{
  "schema_version": "run_ledger_v1",
  "timestamp_utc": "2026-04-25T10:00:00Z",
  "project": "NC2024_NSCLC",
  "git_sha": "abc1234",
  "peak_rss_bytes": 42949672960,
  "final_adata_sha256": "e3b0c44298fc1c149afb...",
  "module_results": [
    {"name": "qc", "status": "ok", "wall_seconds": 12.3, "rss_peak_bytes": 0}
  ]
}
```

**Retention policy**: keep all ledger files indefinitely. Records are append-only and are never auto-pruned. They serve as the primary audit trail for production runs (including NC2024 sparse-exact 900k) and are tracked by git.
