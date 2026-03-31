# Recommended Actions (Prioritized)

## Must fix now
1. Repair `run_end_to_end.sh` to match current Scanpy CLI (`--project`), implement real skip flags, and fix conda-env detection.
2. Make output location deterministic and module-safe (absolute path anchored to repo, not CWD-relative `../output`).
3. Make doublet detection operational:
   - Wire `doublet_threshold` into Scrublet call.
   - Add required dependency (`skimage`) or fail explicitly when unavailable.
4. Replace PBMC default QC settings that currently overfilter this dataset (151/1221 retained).
5. Implement explicit annotation output (`obs['cell_type']` + confidence) instead of marker plots only.

## Should fix next
1. Add reproducible environment spec (`environment.yml` or equivalent lockfile).
2. Fix `gene_ids` mode so mitochondrial/ribosomal QC and marker matching use gene symbols reliably.
3. Add fixed random seeds for UMAP/Leiden and log them.
4. Apply or remove `max_ribo_pct` consistently.
5. Update stale docs (`--outdir`, old `results_*` paths, overstated annotation/isolation claims).

## Nice to improve later
1. Add ambient RNA contamination checks/correction option.
2. Add post-clustering QC summaries (doublet risk, marker coherence, contamination flags).
3. Add regression tests for:
   - argument compatibility
   - output path safety
   - required output artifacts and columns in saved h5ad
4. Add a run manifest file containing input matrix path, software versions, parameters, and hashes.
