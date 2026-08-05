# Wave-9 LUCA patient LIANA OOM recovery — 2026-08-05

## Objective

Complete the Wave-9 LUCA tumour-primary patient LIANA expansion from 3+3 to
8+8 donors while preserving descriptive claim boundaries and recoverability.

## What was attempted

- Original Kimi launch under run `2026-08-05T1250Z-w9tmp`.
- First takeover relaunch under `2026-08-05T1336Z-4e5392e`.
- Root-cause fix and canonical relaunch under `2026-08-05T1345Z-4e5392e`.

## Observed failure

**Regression of previous issue.** The first process stopped without a Python
traceback. Kernel evidence shows OOM kill of PID 8855 at 91,747,284 kB anonymous
RSS. `scanpy.read_h5ad(..., backed="r")` backed `X`, but loaded atlas-wide
892,296×892,296 `obsp` connectivities/distances. The driver retained a main
atlas object and opened another per donor, reproducing the old dense/whole-slot
materialization failure family documented in April and May 2026.

## Fix

- Read only required HDF5 `obs` columns.
- Read only selected rows from CSR-encoded `X`.
- Do not materialize `obsp`, `obsm`, layers, raw, or unrelated metadata.
- Process donors sequentially.
- Write `donor_summary.json` and incremental `patient_liana_meta.json`.
- Add `--resume` with summary+CSV validation.
- Document the OOM trap in `README.md` and `PROTOCOL.md`.

## Successful result

- Canonical run:
  `/home/zerlinshen/projects/pipeline-scientific-audit-20260805/runs/2026-08-05T1345Z-4e5392e/`.
- 8 LUAD + 8 LUSC donors; 16/16 successful; 0 errors.
- 69.33 s wall; 1,788,184 kB maximum RSS; no swap activity during the run.
- Resume replay with CSV row/schema validation: 1.51 s, 227,936 kB maximum
  RSS, exit 0.
- Focused tests: 6 passed; `py_compile` and `git diff --check` passed.
- Full suite validation: 34/34 gates passed at 2026-08-05T14:14:38Z with
  `CONDA_PREFIX=/home/zerlinshen/conda`; raw-data governance passed after the
  blocked E-MTAB-13530 staging root received an explicit manifest.

## Scientific conclusion

F2-04 remains `partial`. Expanded recurrence is descriptive and shows mixed or
modest histology differences, not a clean LUSC-specific programme. It does not
reproduce the paper's background/healthy arms, CellPhoneDB multi-condition
statistics, causality, or figure parity.

## Artifact classification

- `canonical`: `2026-08-05T1345Z-4e5392e` W9.2 outputs and manifest.
- `failed exploratory`: `2026-08-05T1250Z-w9tmp` and
  `2026-08-05T1336Z-4e5392e`.
- `evidence-only`: kernel OOM record and failed-run logs.
- `superseded`: original whole-AnnData donor-loading route.

## Remaining risk and next operator memory

- Preserve sparse HDF5 row reads and one-donor-at-a-time execution for LUCA.
- Do not substitute `read_h5ad(backed="r")` without proving non-`X` slots stay lazy.
- W9.1 remains blocked: P11_T1/T2 targets were zero-filled placeholders despite
  matching advertised lengths. They were quarantined; no section is staged.
  Do not resume the same route without HDF5/tar/non-zero payload validation.
- `project.yaml`, project README, and retention policy source-of-truth did not
  change; the project README and claim records were updated for Wave-9 status.
- Full 34-gate validation is complete. No commit or push was performed.
