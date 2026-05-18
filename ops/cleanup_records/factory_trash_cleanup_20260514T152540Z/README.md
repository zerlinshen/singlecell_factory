# Factory Trash Cleanup 2026-05-14

Created UTC: `2026-05-14T15:25:40Z`

Purpose:
- Restore free disk after `/` reached 100% usage and sandbox/tool startup failed.
- Remove factory-local trash while preserving raw data, references, source code,
  canonical v2 NC2024 outputs, and provenance for deleted v1 outputs.

Disk:
- Before cleanup: `/dev/nvme0n1p5 563G used 535G avail 0 use 100%`
- After cleanup: `/dev/nvme0n1p5 563G used 340G avail 195G use 64%`
- Reclaimed: about `195G`

Preserved source of truth:
- `/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2`
- `/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2`
- `/home/zerlinshen/singlecell_factory/ops/run_ledger/`
- raw data under `/home/zerlinshen/singlecell_factory/data`
- references under `/home/zerlinshen/singlecell_factory/ref`
- all source code and environment definitions

Deleted:
- v1 NC2024 output directories marked in `AGENTS.md` as annotation-bugged and
  not citation-safe.
- checkpoint-only sparse-exact probe directories that lacked final run evidence.
- tiny 2026-05-14 test result directories.
- factory-local `multiomics_r_factory/output` smoke artifacts.
- Python/test/matplotlib caches and generated `coverage.json`.
- `/tmp` pytest/numba/codex bubblewrap scratch directories.
- user/package caches: pip, uv, numba, matplotlib, Chrome cache, and conda
  package/index/tarball cache via `conda clean --all --yes`.

Provenance retained:
- v1 tumor and B/H `run_manifest.json`
- v1 tumor and B/H `module_status.csv`
- v1 tumor and B/H run-ledger JSON copies
- deletion target file inventory before removal:
  `deletion_targets.file_inventory.tsv`

Intentionally kept:
- Older full-cohort result directories still referenced by README and run-memory
  files, including the stage1/extended/rerun-all-eligible runs. These should be
  removed only with a doc/run-memory update that changes the cited source of
  truth.
