# Wave-9 NSCLC Visium (W9.1) unblocked and completed — 2026-08-05

## Objective

Stage the bounded E-MTAB-13530 subset (16 sections, 4 NSCLC patients ×
T1/T2/B1/B2) and run the W9.1 paper-tissue spatial×communication lane
(`scripts/nsclc_visium_spatial_liana.py`) after the same-day zero-payload
block.

## What was attempted

- Content-first re-probe of the same EBI HTTPS payload path that had delivered
  zero-filled placeholders ~1 h earlier: HTTP 206 range GET on
  `P11_T1-filtered_feature_bc_matrix.h5` returned the HDF5 magic at offset 0;
  range GET on `P11_T1-spatial.tar` returned a real tar header. A full
  single-stream GET of the 2.9 MB matrix completed in 6.4 s (~459 kB/s).
- Full staging: all 32 payloads via plain HTTPS GET (curl, resume+retries,
  8 parallel workers), then per-file content validation: size match + HDF5
  magic + ≥1 tar member + non-zero payload byte, SHA-256 recorded in
  `DOWNLOAD_VALIDATION_20260805_HTTPS_STAGED.json`.
- Extraction into `sections/<SECTION>/{filtered_feature_bc_matrix.h5, spatial/}`;
  `sections_ready=16`.
- Lane run: `runs/2026-08-05T1459Z-4e5392e` (execution_mode = debug_massive).

## Observed behaviour worth preserving

- The earlier block was real (four zero-filled placeholders at advertised
  lengths, quarantined). The route became live again within about an hour.
  Treat EBI FIRE delivery as flaky and **always validate content, never size**.
- Two staging corrections vs the initial plan:
  1. The BioStudies-JSON-derived target list duplicated the two P25_B1 files
     and missed P25_B2 entirely; both P25_B2 files were fetched against the
     live directory listing and validated.
  2. `spatial.tar` members extract flat (no `spatial/` prefix); they were
     rearranged into `sections/<SECTION>/spatial/` for `scanpy.read_visium`.
     `/tmp/e13530_extract.sh` was fixed accordingly.

## Successful result

- Canonical run:
  `/home/zerlinshen/projects/pipeline-scientific-audit-20260805/runs/2026-08-05T1459Z-4e5392e/`.
- 16/16 sections successful; 0 errors; wall 75.81 s; maximum RSS
  2,300,244 kB (~2.3 GB); seed 41; 2,500 spots/section cap.
- Per-section LIANA consensus 3,937–32,408 rows; Leiden clusters 4–8.
- Paired tumour-minus-background contrast (4 patients × 2T+2B):
  checkpoint_any LIANA deltas P11 −26 / P19 +830 / P24 +337.5 / P25 +261
  (3/4 positive, heterogeneous); foetal-mac panel Moran's I deltas mostly
  positive in tumour sections.

## Scientific conclusion

F4-01/F4-02 remain `partial`. The lane is genuine NSCLC paper-cohort tissue
evidence (De Zuani cohort Visium, no longer mouse-brain demo), but descriptive:
Leiden spot clusters as proxy groups (no cell2location, no pathologist
annotation), LIANA consensus engine (not paper CellphoneDB multi-condition),
`figure_parity=false`. `checkpoint_immune` stays 0-by-construction on cluster
labels; use `checkpoint_any`. No claim was upgraded.

## Artifact classification

- `canonical`: `2026-08-05T1459Z-4e5392e` W9.1 outputs and manifest;
  staged raw pool `data/raw/nc2024_nsclc_visium_emtab13530/` (state
  `staged_valid`, `valid_sections_staged=16`).
- `evidence-only`: quarantined zero placeholders
  `downloads/invalid_zero_placeholders_20260805/` and
  `DOWNLOAD_VALIDATION_20260805.json` (the blocked-state finding stands as the
  record that size match is not success evidence).
- `superseded`: blocked-state MANIFEST/PROVENANCE wording (replaced by
  staged_valid state the same day).

## Remaining risk and next operator memory

- EBI payload route is unstable: it served zero-filled placeholders and real
  bytes within the same hour on 2026-08-05. Any re-fetch must re-validate
  content (HDF5 magic, tar members, non-zero payload) before use.
- `/tmp/e13530_targets.json` is now the corrected 32-file list (dupes removed,
  P25_B2 added); recreate from the live directory listing if lost.
- Wave-9 was validated (34/34) and then pushed the same day: sc `d9d4ea2`
  (`wave6-trevino-v5.1`), r `4707744` (`master`), suite `a60d15d` (`master`).
  The earlier entries' "no commit or push" notes describe the pre-push state.
