# Two Real-Dataset Final Factory Validation

Control-plane validation report; scientific outputs remain under project roots.

- Generated UTC: `2026-05-19T05:17:53.771280+00:00`
- Verdict: `PASS_TWO_REAL_DATASET_FACTORY_BRIDGE_VALIDATED`
- Datasets: NC2024 NSCLC P15_T1 and Cell/Trevino public-resource reproduction.

## NC2024

- Run: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`
- Shape: `5281 x 19504`
- Bundle: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/python/bundle`
- R outputs: `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88/r` (8 files)
- Boundary: project/module architecture validation, not full article-scale cohort annotation.

## Cell/Trevino

- Pipeline run: `2026-05-17T2004Z-13c2c88`
- Evidence run: `20260517T1436Z-13c2c88`
- Shape: `55653 x 25519`
- Bundle: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/python/bundle`
- R outputs: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/r` (8 files)
- Quality gate: `conditional` with expected resource gaps documented.

## Three-Factory Result

- `singlecell_factory`: upstream run evidence, bundle export, governance validator.
- `r_multiomics_factory`: R bundle reading and plotting entrypoint.
- `plotting_factory`: source modules used by R plotting and visual smoke tests.

## Ownership Governance

- Rule: `owner-by-primary-output`.
- `singlecell_factory`: default global control plane for Python-heavy upstream and cross-factory validation.
- `r_multiomics_factory`: local owner for R-heavy/spatial primary scientific truth when R produces primary objects and biological interpretation.
- `plotting_factory`: presentation-only plotting surface; it must not own biological conclusions.
- Canonical policy: `docs/OWNER_BY_PRIMARY_OUTPUT_GOVERNANCE.md`.

## Hardening Gates

- Clean-room minimal real-data gate: `cleanroom_minimal_realdata.json`.
- Figure parity gate: `figure_parity_gate.json` (conditional until curated paper/process-data references are registered).
- CI governance gate: `ci_governance_gate.json`.

## Remaining Boundaries

- NC2024 retained P15_T1 is not full multi-cohort paper-scale truth.
- Cell/Trevino remains conditional public-resource reproduction, not raw FASTQ/fragments/BPNet parity.
