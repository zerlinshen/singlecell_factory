# 2026-04-30 - NC2024 Remote-To-Mac Handoff Protocol

## Objective

- Close the loop after the Rarrow remote validation by pulling review-ready
  figure artifacts back to the Mac workspace.
- Create a durable protocol so future remote figure/report outputs are not left
  remote-only when they are intended for local review.

## Starting Context

- Rarrow wrapper validation succeeded remotely on `ubuntu-tail`.
- B/H and tumor wrapper smoke outputs existed under
  `/home/zerlinshen/multiomics_r_factory/output/`.
- The local workspace already separated human-facing evidence under
  `human_review/` from agent run memory under `agent_runs/`, but it did not yet
  have an explicit remote-to-Mac handoff protocol.

## What Was Run

- Host/source: `ubuntu-tail`
- Local destination:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/human_review/figures/2026-04-30-r-arrow-wrapper-smoke`
- Transfer method: `rsync -az`

Remote sources:

```text
/home/zerlinshen/multiomics_r_factory/output/nc2024_bh_arrow_wrapper_smoke_20260430T100416+0800
/home/zerlinshen/multiomics_r_factory/output/nc2024_tumor_arrow_wrapper_smoke_20260430T100432+0800
/home/zerlinshen/singlecell_factory/ops/env_records/r_arrow_20260430T095459+0800
```

## Outcome

- `success`: B/H and tumor PNG/PDF smoke outputs were transferred to the Mac.
- `success`: wrapper logs and Rarrow runtime probe were transferred under
  `provenance/`.
- `success`: local `TRANSFER_MANIFEST.tsv` and `checksums.sha256` were created.
- `success`: durable protocol created at:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/REMOTE_TO_MAC_HANDOFF_PROTOCOL.md`

## Evidence

- Local handoff directory:
  `/Users/zerlinshen/Downloads/1. Codex/2.Reproduction Trail/2026-04-23 - Nature Communications 2024 NSCLC single-cell reproduction/human_review/figures/2026-04-30-r-arrow-wrapper-smoke`
- Local manifest:
  `TRANSFER_MANIFEST.tsv`
- Local checksums:
  `checksums.sha256`
- Local provenance:
  `provenance/nc2024_bh_wrapper_validation_v2.log`
  `provenance/nc2024_tumor_wrapper_validation.log`
  `provenance/r_multiomics_arrow_post_r_probe.txt`

## Problems Encountered

- The missing norm was real: prior docs said the Mac is for review/packaging but
  did not specify when to pull back, where to place files, or how to verify the
  transfer.

## What Changed

- Added `REMOTE_TO_MAC_HANDOFF_PROTOCOL.md`.
- Added the first protocol-compliant handoff:
  `human_review/figures/2026-04-30-r-arrow-wrapper-smoke/`.
- Updated local `human_review`, `agent_runs`, root README, and Codex profile
  surfaces to point to the new protocol and handoff.

## Cautions For The Next Run

- Remote artifacts remain canonical.
- Pull rendered figures/reports and small provenance files back to the Mac when
  human review is expected.
- Do not pull full `final_adata.h5ad`, prepared Zarr, raw data, or bulky bundle
  archives unless local object-level analysis is explicitly required.

## Classification

- canonical:
  - remote output directories listed above
  - `REMOTE_TO_MAC_HANDOFF_PROTOCOL.md`
- local-review-copy:
  - `human_review/figures/2026-04-30-r-arrow-wrapper-smoke/`
- evidence-only:
  - Rarrow wrapper smoke figures and logs
- superseded:
  - remote-only figure/report handoffs without a local review copy when human
    review is expected
- failed exploratory:
  - none
