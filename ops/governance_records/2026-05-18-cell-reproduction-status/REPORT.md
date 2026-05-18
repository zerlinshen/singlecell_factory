# Cell/Trevino Reproduction Status

Control-plane status record only; not a scientific output.

- Generated UTC: `2026-05-18T09:20:32Z`
- Run dir: `/home/zerlinshen/projects/wave5-trevino/runs/20260517T1436Z-13c2c88`
- Decision: `conditional`
- Reason: Code review is clean for current scope; scientific gate is conditional because expected false checks/resource gaps remain explicitly documented rather than forced.
- Evidence JSON count: `17`
- Direct evidence panel count: `31`
- Progress snapshot panel count: `53`
- Unexpected false checks: `[]`
- Missing referenced paths: `[]`
- Hash mismatches: `[]`

## Expected Resource Gaps

- `figure2cd_public_author_cre_gene_resources.json:figure2c_public_link_count_matches_paper_64878`
- `figure3b_public_matrix_atac_pseudotime_transfer.json:all_marker_direction_checks_pass`
- `figure7_public_asd_bpnet_evidence.json:s6d_caption_like_case_matches_262`

## Resource Boundaries

- Current computation starts from GEO public matrices, public author RDS/tables, and Cell supplementary tables, not public FASTQ/SRA raw reads; no FASTQ-to-CellRanger parity is claimed.
- Figure 7 exact BPNet denominator universe/model weights/ref-alt tracks are absent from staged public resources; S6-derived summaries are generated with explicit gap panels.
- Non-computational/IHC/source-image panels are classified as not reproducible from single-cell matrices unless raw microscopy images are obtained.

## Code Review

- Recommendation: `APPROVE`
- Clean: `True`

## Retention

No Cell/Trevino run deletion performed in this pass; latest evidence run remains 20260517T1436Z-13c2c88.
