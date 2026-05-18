# Remote Factory/Project Governance Run

Control-plane governance record only; not a scientific output.

## Objective

Implement phase-1 remote-first governance in `/home/zerlinshen/singlecell_factory` and validate `/home/zerlinshen/projects/wave5-trevino` plus `/home/zerlinshen/projects/nc-reproduction` without moving or deleting project run artifacts.

## Scope

- Remote factory docs/contracts/scripts/tests
- Read-only project validation
- No local Reproduction Trail implementation
- No project artifact migration, deletion, or consolidation

## Baseline

See `BASELINE_PIN.md`, `git_status_short.txt`, `tracked_diff.patch`, and `pre_governance_untracked_inventory.txt`.

## Final Verification

timestamp_utc: 2026-05-18T08:02:50Z
branch: wave6-trevino-v5.1
head: 13c2c885c981262fc634125cc46f4abc9c08298f
status_count_final: 960
tracked_diff_sha256_final: 9205c21de21cc38e064d978fa91dd965fcabdb2674da95e534c7bfe0391c3486
untracked_inventory_sha256_final: 43a74b326892a1c4c32fd40f95ff8b113ba2576a5ae396a6c195e94c128becad

Commands passed:
- `python3 -m py_compile scripts/validate_project_governance.py`
- `python3 -m unittest tests.test_project_governance_validator` (3 tests)
- `git diff --check -- README.md AI_AGENT_PROTOCOL.md PROTOCOL.md docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md contracts/project_run_contract.yaml scripts/validate_project_governance.py tests/test_project_governance_validator.py`
- `python3 scripts/validate_project_governance.py /home/zerlinshen/projects/wave5-trevino --output-dir ops/governance_records/2026-05-18-remote-factory-project-governance/project_validations`
- `python3 scripts/validate_project_governance.py /home/zerlinshen/projects/nc-reproduction --output-dir ops/governance_records/2026-05-18-remote-factory-project-governance/project_validations`
- `python3 scripts/validate_project_governance.py /home/zerlinshen/projects/wave5-trevino --run-id 20260517T1436Z-13c2c88 --output-dir ops/governance_records/2026-05-18-remote-factory-project-governance/project_validations`
- JSON assertion: every governance validation report has `overall_status == pass` and `control_plane_only == true`
- PyYAML parse check for `contracts/project_run_contract.yaml` when PyYAML is available

Validation report summary:
- `nc-reproduction_all-runs.governance_validation.json`: pass; severity={'info': 2, 'warning': 1}; runs=3
- `wave5-trevino_20260517T1436Z-13c2c88.governance_validation.json`: pass; severity={'warning': 4, 'info': 1}; runs=1
- `wave5-trevino_all-runs.governance_validation.json`: pass; severity={'warning': 10, 'info': 1}; runs=4

AI slop cleanup report:
- Scope: Ralph-owned governance files only.
- Behavior lock: targeted unittest, py_compile, diff check, YAML parse, real-project validator reports.
- Simplification: tightened excessive Markdown/YAML blank lines after adding human-vs-agent and purpose-driven module-selection contract text.
- No runtime default changes, no new dependencies, no project artifact migration.

Changed governance surfaces:
- `README.md`
- `AI_AGENT_PROTOCOL.md`
- `PROTOCOL.md`
- `contracts/project_run_contract.yaml`
- `docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md`
- `scripts/validate_project_governance.py`
- `tests/test_project_governance_validator.py`
- `ops/governance_records/2026-05-18-remote-factory-project-governance/`

Scope boundary:
- Project roots were used as read-only validation inputs.
- Validation reports were written only to the factory control-plane run record.
- Local Mac Reproduction Trail implementation remains out of scope for this phase.

## Post-Question Structure Decision

timestamp_utc: 2026-05-18T08:43:56Z

Decision implemented:
- New factory governance/control-plane records use `ops/governance_records/`.
- Project analysis run records should live with projects under `/home/zerlinshen/projects/<project-id>/runs/` and project ledgers.
- Existing historical `ops/run_records/` entries were not migrated in this pass; that is a separate retention/migration task.
- The factory remains a pure analysis environment plus contracts, validators, docs, and factory governance state.

Fresh verification:
- `python3 -m py_compile scripts/validate_project_governance.py`: pass
- `python3 -m unittest tests.test_project_governance_validator`: pass, 3 tests
- `git diff --check -- README.md AI_AGENT_PROTOCOL.md PROTOCOL.md docs/REMOTE_FACTORY_PROJECT_GOVERNANCE.md contracts/project_run_contract.yaml scripts/validate_project_governance.py tests/test_project_governance_validator.py`: pass
- Real project validations rerun with `--output-dir ops/governance_records/2026-05-18-remote-factory-project-governance/project_validations`: pass

Validation reports:
- `nc-reproduction_all-runs.governance_validation.json`: pass; severity={'info': 2, 'warning': 1}; runs=3
- `wave5-trevino_20260517T1436Z-13c2c88.governance_validation.json`: pass; severity={'warning': 4, 'info': 1}; runs=1
- `wave5-trevino_all-runs.governance_validation.json`: pass; severity={'warning': 10, 'info': 1}; runs=4

Fresh pins:
- status_count: 960
- tracked_diff_sha256: 80af5af257e438179d47e7e7c40f83822789af6e507486882bca260c41ecbd75
- untracked_inventory_sha256: f4c4695bac44fa4afd88b46ed0b0165123ad31433e8b1e1615cf871665686093
