# NC2024 Project Latest-Run Retention Cleanup

Control-plane cleanup record only; not a scientific output.

- Timestamp UTC: `2026-05-18T09:19:39Z`
- Project: `/home/zerlinshen/projects/nc-reproduction`
- Kept run: `2026-05-18T0900Z-13c2c88`
- Deleted run count: `3`
- Deleted total bytes: `142437966218`

## Deleted Runs

- `2026-05-14T0542Z-d192836`: 220545479 bytes, root_manifest=True, python_manifests=0, module_status=0, files=25
- `2026-05-15T1815Z-wave4-200k-preflight`: 22732430819 bytes, root_manifest=False, python_manifests=0, module_status=0, files=13058
- `2026-05-15T1833Z-wave4-800k-clustering`: 119484989920 bytes, root_manifest=True, python_manifests=1, module_status=1, files=35813

## Preserved

- `/home/zerlinshen/projects/nc-reproduction/runs/2026-05-18T0900Z-13c2c88`
- raw/prepared inputs, configs, notebooks, launch scripts, source code, environments, and governance records

## Disk

Before:
```text
Filesystem      Size  Used Avail Use% Mounted on
/dev/nvme0n1p5  563G  300G  234G  57% /
```

After:
```text
Filesystem      Size  Used Avail Use% Mounted on
/dev/nvme0n1p5  563G  168G  367G  32% /
```
