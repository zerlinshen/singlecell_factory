# Project Governance Validation: wave5-trevino

Control-plane validation report only; not a scientific output.

- Generated: `2026-05-18T08:42:52Z`
- Project root: `/home/zerlinshen/projects/wave5-trevino`
- Requested run: `all discovered runs`
- Overall status: `pass`
- Severity counts: `{'warning': 10, 'info': 1}`

## Runs

### `2026-05-17T2004Z-13c2c88`
- Root manifest: `/home/zerlinshen/projects/wave5-trevino/runs/2026-05-17T2004Z-13c2c88/manifest.json`
- Figures: `15`
- Evidence JSON: `0`
- Python run manifests: `1`
- Module status files: `1`
- Bundle provenance files: `0`
- R manifests/provenance: `0`

### `20260516T0931Z-d192836f1bb0`
- Root manifest: `missing_expected_manifest`
- Figures: `54`
- Evidence JSON: `0`
- Python run manifests: `1`
- Module status files: `1`
- Bundle provenance files: `0`
- R manifests/provenance: `0`
- Findings:
  - `warning` `missing_root_manifest`: Missing cross-factory run envelope; acceptable as legacy/partial state when substitute provenance exists.
  - `warning` `missing_r_dir`: Missing canonical run subdirectory; warning for legacy/partial runs.
  - `warning` `missing_logs_dir`: Missing canonical run subdirectory; warning for legacy/partial runs.

### `20260517T1005Z-7b539a5`
- Root manifest: `missing_expected_manifest`
- Figures: `6`
- Evidence JSON: `0`
- Python run manifests: `1`
- Module status files: `1`
- Bundle provenance files: `1`
- R manifests/provenance: `0`
- Findings:
  - `warning` `missing_root_manifest`: Missing cross-factory run envelope; acceptable as legacy/partial state when substitute provenance exists.
  - `warning` `missing_r_dir`: Missing canonical run subdirectory; warning for legacy/partial runs.
  - `warning` `missing_logs_dir`: Missing canonical run subdirectory; warning for legacy/partial runs.

### `20260517T1436Z-13c2c88`
- Root manifest: `missing_expected_manifest`
- Figures: `173`
- Evidence JSON: `4`
- Python run manifests: `0`
- Module status files: `0`
- Bundle provenance files: `0`
- R manifests/provenance: `0`
- Findings:
  - `warning` `missing_root_manifest`: Missing cross-factory run envelope; acceptable as legacy/partial state when substitute provenance exists.
  - `warning` `missing_r_dir`: Missing canonical run subdirectory; warning for legacy/partial runs.
  - `warning` `missing_logs_dir`: Missing canonical run subdirectory; warning for legacy/partial runs.
  - `info` `no_producer_native_manifest_found`: No producer-native manifest/module-status files found at expected paths; inspect legacy evidence manually if this run is cited.

## Project Findings
- `warning` `missing_project_yaml`: Missing project.yaml; warning for legacy projects, required for new-governance projects.
