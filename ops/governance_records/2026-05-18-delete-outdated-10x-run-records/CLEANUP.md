# Delete Outdated 10x/Test Run Records

Control-plane cleanup record only; not a scientific output.

- Timestamp UTC: `2026-05-18T08:53:21Z`
- Branch: `wave6-trevino-v5.1`
- Head: `13c2c885c981262fc634125cc46f4abc9c08298f`
- Decision: delete outdated factory `ops/run_records` stubs because current article validation lives in project-owned run/evidence structures.
- Deleted file count: `6`
- Deleted bytes: `3545`

## Deleted Files
- `ops/run_records/nc2024/20260426_201800/run_record_stub.json` (555 bytes, sha256 `1e25541e3a2d3e1a1fe69f50cf9214e8f05d37715180a704ffdf3615cc88db64`)
- `ops/run_records/nc2024/20260426_201801/run_record_stub.json` (552 bytes, sha256 `975d7fa6c75f54b6548fd874dcfdfc9a2719dacb1ce9068351eb49369480735c`)
- `ops/run_records/nc2024/20260426_202319/run_record_stub.json` (555 bytes, sha256 `fb2047fc092a9265ebaf7f6f98cc6aad474ba893b5c4f402e3627f9b513ca83e`)
- `ops/run_records/nc2024/20260426_202320/run_record_stub.json` (552 bytes, sha256 `3ad221f5115c98356ccb7139a48964e164d4f2852474a6348940e070a4182130`)
- `ops/run_records/wave2/wave2_us018_disk_preflight_20260515/run_record_stub.json` (517 bytes, sha256 `151ea55ed3e4eccbcbc81a66cfb17662681c6490798e028897622dfc82e21932`)
- `ops/run_records/wave5-trevino/20260516T0931Z-d192836f1bb0/run_record_stub.json` (814 bytes, sha256 `7a603ae0762072d26290a0fe634629bdcfa69f0ec81a78391547596fca470098`)

## Preserved

- `/home/zerlinshen/projects/*` project data and article-validation runs
- raw downloads and references
- factory code, contracts, validators, docs, tests, and environment definitions
- `ops/governance_records/*` factory governance records
