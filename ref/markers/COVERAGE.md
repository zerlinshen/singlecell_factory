# Marker DB Coverage Matrix

This file tracks which `(tissue, condition)` pairs are populated in each vendored DB version.
AC-3 verification scripts MUST consult this file before evaluating whether a new tissue/condition
combination is supported. If the pair is absent, AC-3 verification SKIPs with status
`INSUFFICIENT_DB_COVERAGE`.

## Coverage table

| tissue | condition | cellmarker2/v2.0 | panglaodb/v2021-03-18 | celltypist/v2.6.0 | sctypedb/v1.0 |
|---|---|---|---|---|---|
| lung | NSCLC | populated | populated | populated | populated |

## Notes

- "populated" = `markers.tsv` contains at least one row for `(tissue, condition)`.
- "placeholder" = directory exists but `markers.tsv` has no rows for this pair.
- Initial commit covers only `(lung, NSCLC)` matching the NC fixture.
- New pairs are added per DB-snapshot refresh (see `LICENSE_AUDIT.md` refresh policy).
