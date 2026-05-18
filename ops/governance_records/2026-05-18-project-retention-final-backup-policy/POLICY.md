# Project Retention And Final Backup Policy

Control-plane governance record only; not a scientific output.

- Timestamp UTC: `2026-05-18T09:23:19Z`
- Decision: latest-run-only is the default project storage rule, but each project
  must also keep a project-specific retention/final-backup policy recording the
  analysis type, method family, selected/best parameters, final backup contents,
  and never-delete assets.
- Scope: documentation and governance contract update only; no data deletion or
  pipeline execution is performed by this policy record.

## Required Project Policy Fields

- project purpose and analysis type
- method family and selected/best parameters
- module list and safe skip conditions
- filtering/QC thresholds
- batch/integration keys
- sketch/metacell settings when used
- random seed and paper-specific target choices
- latest source-of-truth run id/path
- superseded run ids/paths
- final backup contents
- never-delete assets

## Final Backup Contents

- root `manifest.json`
- producer-native manifests and `module_status.csv`
- launch command/wrapper and launch log
- parameter file and environment pins
- final labeled object or compact atlas
- human-facing figures/reports
- conclusion summary
