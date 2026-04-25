# Staircase Test Fixtures

A 4-tier scaffold for catching scale-related bugs at the smallest tier
that would actually trigger them, instead of waiting until 884k production
runs to discover OOMs / threshold-gated routing bugs.

## Tiers

| Tier | n_cells | Source | Purpose | Default? |
|---|---|---|---|---|
| `nano` | ~5k | synthetic CSR | Unit + parity numerics; covers algorithm correctness | ✅ runs in default `pytest` |
| `small_real` | ~100k | NC2024 prepared zarr (~10 samples stratified) | Memory-affecting changes; grouped-Scrublet trigger threshold | ⛔ opt-in via `pytest -m small_real` |
| `medium_real` | ~500k | NC2024 prepared zarr (~50 samples stratified) | Pre-publication runs; full module suite at near-production scale | ⛔ opt-in via `pytest -m medium_real` |
| `full_real` | 884,050 | original NC2024 prepared zarr (path pointer, not copied) | Publication-grade gate; one-shot before paper submission | ⛔ opt-in via `pytest -m full_real` |

## Why staircase

Today's OOM happened because:
- nano (5k) does not exceed the 100k threshold, so grouped Scrublet auto-routing was never exercised.
- production launch (884k) blew up at module 3 of 21.

The 100k tier catches the auto-threshold bug; the 500k tier catches density-dependent memory effects (HVG sparsity differs at scale); the 884k tier is the final publication-grade gate.

## Generate

```bash
# All tiers (consumes ~3-5 GB disk):
python scripts/build_staircase_fixtures.py

# One tier at a time:
python scripts/build_staircase_fixtures.py --tiers small_real
python scripts/build_staircase_fixtures.py --tiers medium_real

# full_real does not copy data — writes a path pointer to tests/data/staircase/full_real.path:
python scripts/build_staircase_fixtures.py --tiers full_real
```

Re-run after changes to the underlying prepared zarr. Seed is `42` by default for reproducible cell selection; override via `--seed`.

## Run

```bash
# nano runs in default suite (no marker needed)
pytest tests/

# Individual real-data tiers
pytest -m small_real tests/test_staircase_smoke.py
pytest -m medium_real tests/test_staircase_smoke.py

# full_real — only before publication; combine with paper-aligned launcher
bash scripts/run_nc2024_paper_aligned_20260425.sh
```

## What's tracked here

- This README
- `.gitignore` excludes `*.h5ad` / `*.zarr` (large) but allows `*.summary.json` (provenance)
- `build_staircase_fixtures.py` writes `<tier>.summary.json` next to each fixture with cell counts, sample provenance, seed, source SHA — keep these tracked for reproducibility audits

## Discipline

| Change category | Required gate |
|---|---|
| Pure-numerics module change | nano (default) |
| Memory / density / chunking change | nano + small_real |
| Paper-method-aligning change | nano + small_real + medium_real |
| Pre-submission publication run | nano + small_real + medium_real + full_real |
