# Marker DB License Audit

All vendored marker databases are offline snapshots. No runtime network calls are made.

| DB | Version | License | Source URL | Snapshot Date | Commercial Use |
|---|---|---|---|---|---|
| CellMarker 2.0 | v2.0 | CC BY 4.0 | http://bio-bigdata.hrbmu.edu.cn/CellMarker/ | 2024-01-15 | Yes (with attribution) |
| PanglaoDB | v2021-03-18 | CC BY 4.0 | https://panglaodb.se/ | 2021-03-18 | Yes (with attribution) |
| CellTypist | v2.6.0 | MIT | https://www.celltypist.org/ | 2024-03-01 | Yes |
| scTypeDB | v1.0 | MIT | https://github.com/IanevskiAleksandr/sc-type | 2024-01-15 | Yes |

## Usage constraints

- CellMarker 2.0 and PanglaoDB require attribution in publications.
- CellTypist and scTypeDB (MIT) have no additional publication requirements.
- All snapshots are read-only factory-side static assets under `ref/markers/<db>/<version>/`.
- Snapshots must NOT be modified after SHA256 is recorded in `MANIFEST.json`.

## Refresh policy

When a new DB version is vendored:
1. Add a new versioned subdirectory: `ref/markers/<db>/<new_version>/`.
2. Populate `MANIFEST.json` with the new `sha256` from the downloaded file.
3. Update this file with the new row.
4. Update `ref/markers/COVERAGE.md` with any new `(tissue, condition)` pairs.
5. The old version directory is retained for reproducibility.
