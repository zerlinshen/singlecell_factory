# Shenxin Routing Contract

Use this routing order for the dominant Shenxin workflow:

1. `before-every-run`
2. `bioinformatics-autosteer`
3. `singlecell-remote-workflow`
4. `execute-and-recover-pipeline`
5. `remote-run-report-bridge`
6. `review-and-validate-quality`
7. `test-run-cleanup-policy`
8. `post-run-failure-cleanup`
9. `nsclc-baseline-reproduction`
10. `readme-sync-enforcer`
11. `skill-usage-playbook`
12. `develop-and-integrate-module`
13. `singlecell-factory-module-delivery`

## Default interpretation

- remote = compute and source of truth
- local = organization, review, comparison, handoff, report packaging
- `large` = capacity probe for NC2024 full cohort
- `massive` = debug/recovery lane for NC2024 full cohort

## Execution mode defaults

- `debug_massive`
  - direct `massive`
  - stop at first failing module
  - patch only that module
- `controller_validation`
  - validate the official `large -> massive` orchestration path
- `benchmark`
  - deliberate comparison only
