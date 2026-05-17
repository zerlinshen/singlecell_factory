# Wave 5.1 HV2 WNN Sweep Pre-Registration

This directory holds the SHA-pinned search space and aggregation rule that must
exist before any v5.1 WNN sweep result is generated.

Rules:

- Do not edit `search_space.json` or `aggregation_rule.json` after first sweep
  execution.
- Every sweep result row must record the SHA-256 of both files.
- If no configuration passes the guard, record
  `INCOMPLETE-NO-VIABLE-FIX`; do not relax the guard post hoc.

