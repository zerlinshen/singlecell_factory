# Project: Single-cell RNA-seq analysis (Scanpy)

## Context
This project performs single-cell RNA-seq analysis using Scanpy.
Data source: 10x Genomics (Cell Ranger output).
Main object: AnnData (.h5ad)

## Goals
- QC filtering
- normalization
- HVG selection
- PCA / UMAP
- clustering (Leiden)
- cell type annotation
- downstream analysis (DEG, pathway)

## Preferences
- Use Scanpy (Python), not Seurat (R)
- Keep code reproducible and modular
- Avoid unnecessary complexity
- Prefer clear step-by-step pipelines

## Style
- concise
- practical
- no generic explanations
- focus on actionable code

## Commands
- when editing code: show full updated block
- when debugging: identify root cause first
