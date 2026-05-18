# NC2024 Architecture Validation

Control-plane validation report only; not a scientific output.

- Generated UTC: `2026-05-18T08:57:09Z`
- Verdict: `PROJECT_ARCHITECTURE_PASS_WITH_FACTORY_LEGACY_DEBT`
- Project root: `/home/zerlinshen/projects/nc-reproduction`
- Factory repo: `/home/zerlinshen/singlecell_factory`
- Branch/head: `wave6-trevino-v5.1` / `13c2c885c981262fc634125cc46f4abc9c08298f`

## Conclusion

The new architecture is valid for the NC2024 project layer: project metadata, run roots, root manifests where available, Python provenance where available, figures, and logs live under `/home/zerlinshen/projects/nc-reproduction`. Governance validation is control-plane only and writes under `ops/governance_records/`.

The factory is not fully pure yet because pre-project-architecture NC2024 outputs remain under `singlecell_factory/results/`. They are classified as legacy debt and were not moved or deleted during this validation pass.

## Governance Validator

- Overall status: `pass`
- Severity counts: `{'info': 2, 'warning': 1}`
- Control-plane only: `True`
- Report: `ops/governance_records/2026-05-18-nc2024-architecture-validation/project_validations/nc-reproduction_all-runs.governance_validation.json`

## NC2024 Project Runs

- `2026-05-14T0542Z-d192836`: root_manifest=True, python=True, r=True, logs=True, python_manifests=0, module_status=0, figures=12, evidence_json=0
- `2026-05-15T1815Z-wave4-200k-preflight`: root_manifest=False, python=True, r=True, logs=True, python_manifests=0, module_status=0, figures=11, evidence_json=0
- `2026-05-15T1833Z-wave4-800k-clustering`: root_manifest=True, python=True, r=True, logs=True, python_manifests=1, module_status=1, figures=11, evidence_json=0

## Factory Purity Check

- `ops/run_records` exists: `False`
- This validation record: `ops/governance_records/2026-05-18-nc2024-architecture-validation`
- Governance stray files under project root: `0`
- Legacy NC2024 result paths in factory `results/`: `48`
- Legacy NC2024 result total size: `30667532081` bytes
- Legacy path inventory: `ops/governance_records/2026-05-18-nc2024-architecture-validation/factory_legacy_nc2024_results.txt`

## Recommendation

Treat the architecture as accepted for new work. Do not write new NC2024 scientific outputs into the factory. Run a separate retention cleanup for legacy `results/NC2024*` after confirming which files, if any, still need local/report references.
