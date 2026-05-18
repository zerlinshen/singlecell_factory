#!/usr/bin/env python3
"""Build a pipeline-ready AnnData from public GSE162170 scRNA matrices.

This prepares the earliest-public-input RNA lane for Wave-5 raw/public Trevino
reproduction.  It streams gene-by-cell TSVs in chunks, sparsifies immediately,
and writes ``prepared_input.h5ad`` under a project input directory.  No
scientific output is written inside the factory tree.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


RNA_COUNTS = "GSE162170_rna_counts.tsv.gz"
RNA_SPLICED = "GSE162170_rna_spliced_counts.tsv.gz"
RNA_UNSPLICED = "GSE162170_rna_unspliced_counts.tsv.gz"
RNA_AMBIGUOUS = "GSE162170_rna_ambiguous_counts.tsv.gz"
RNA_METADATA = "GSE162170_rna_cell_metadata.txt.gz"


def _load_gene_symbol_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    if not path.exists():
        raise FileNotFoundError(f"Gene map not found: {path}")
    table = pd.read_csv(path, sep="\t", usecols=["gene_id", "gene_name"], dtype=str)
    table = table.dropna(subset=["gene_id", "gene_name"])
    return dict(zip(table["gene_id"].astype(str), table["gene_name"].astype(str)))


def _read_gene_ids_only(path: Path) -> list[str]:
    genes: list[str] = []
    with gzip.open(path, "rt") as fh:
        next(fh, None)  # header
        for line in fh:
            if line:
                genes.append(line.split("\t", 1)[0])
    return genes


def _var_from_gene_ids(gene_ids: list[str], gene_map: dict[str, str]) -> pd.DataFrame:
    """Create AnnData var with marker-friendly unique symbols where possible.

    Public GSE162170 RNA matrices use ENSEMBL ids, while the factory marker
    annotator matches marker dictionaries against ``adata.var_names``.  To keep
    annotation viable without losing provenance, use unique gene symbols as
    var_names when a one-to-one symbol is available and fall back to the ENSEMBL
    id for unmapped or duplicate symbols.  The original id is always preserved
    in ``var['gene_id']``.
    """
    seen_symbols: set[str] = set()
    var_names: list[str] = []
    symbols: list[str] = []
    fallback_count = 0
    duplicate_symbol_count = 0
    unmapped_count = 0
    for gene_id in map(str, gene_ids):
        symbol = str(gene_map.get(gene_id, "")).strip()
        symbols.append(symbol)
        if not symbol:
            unmapped_count += 1
            fallback_count += 1
            var_names.append(gene_id)
        elif symbol in seen_symbols:
            duplicate_symbol_count += 1
            fallback_count += 1
            var_names.append(gene_id)
        else:
            seen_symbols.add(symbol)
            var_names.append(symbol)
    var = pd.DataFrame(
        {
            "gene_id": list(map(str, gene_ids)),
            "gene_symbol": symbols,
            "var_name_source": [
                "gene_symbol" if name == symbol and symbol else "gene_id"
                for name, symbol in zip(var_names, symbols)
            ],
        },
        index=pd.Index(var_names, name="gene"),
    )
    var.attrs["gene_symbol_mapping"] = {
        "mapped_to_symbol_var_names": int(sum(var["var_name_source"] == "gene_symbol")),
        "fallback_to_gene_id": int(fallback_count),
        "unmapped_gene_ids": int(unmapped_count),
        "duplicate_symbols_fell_back": int(duplicate_symbol_count),
    }
    return var


def _stream_gene_by_cell_matrix(
    path: Path,
    *,
    chunksize: int,
    wanted_cells: list[str] | None = None,
    expected_genes: list[str] | None = None,
    max_genes: int | None = None,
) -> tuple[sparse.csr_matrix, list[str], list[str], int]:
    """Return cells × genes CSR matrix streamed from a genes × cells TSV.

    When ``expected_genes`` is provided, the source file may have the same genes
    in a different order.  In that mode the stream keeps only the requested
    genes, then reorders columns to exactly match ``expected_genes``.  This is
    required for AnnData layers, whose shape/order must match ``.var``.
    """
    chunks: list[sparse.csr_matrix] = []
    genes: list[str] = []
    cells: list[str] | None = None
    nnz = 0
    target_genes = [str(g) for g in expected_genes] if expected_genes is not None else None
    target_set = set(target_genes) if target_genes is not None else None

    for chunk in pd.read_csv(path, sep="\t", index_col=0, compression="gzip", chunksize=chunksize):
        if cells is None:
            cells = chunk.columns.tolist()
            if wanted_cells is None:
                wanted_cells = cells
        assert wanted_cells is not None
        if max_genes is not None and target_genes is None:
            remaining = max_genes - len(genes)
            if remaining <= 0:
                break
            if len(chunk) > remaining:
                chunk = chunk.iloc[:remaining, :]

        chunk_genes = [str(x) for x in chunk.index.tolist()]
        if target_set is not None:
            keep_positions = [i for i, gene in enumerate(chunk_genes) if gene in target_set]
            if not keep_positions:
                del chunk
                continue
            chunk = chunk.iloc[keep_positions, :]
            chunk_genes = [chunk_genes[i] for i in keep_positions]

        missing = [c for c in wanted_cells if c not in chunk.columns]
        if missing:
            raise ValueError(f"{path.name} missing {len(missing)} requested cells; first={missing[:3]}")
        sub = chunk.loc[:, wanted_cells].to_numpy(dtype=np.int32, copy=False)
        mat = sparse.csr_matrix(sub.T, dtype=np.int32)
        nnz += int(mat.nnz)
        chunks.append(mat)
        genes.extend(chunk_genes)
        del chunk, sub, mat
        if target_genes is not None and len(genes) == len(target_genes):
            break

    if cells is None:
        raise ValueError(f"No matrix rows read from {path}")
    if not chunks:
        raise ValueError(f"No matrix chunks retained from {path}")
    matrix = sparse.hstack(chunks, format="csr", dtype=np.int32)
    if target_genes is not None:
        gene_to_col = {gene: idx for idx, gene in enumerate(genes)}
        missing_genes = [gene for gene in target_genes if gene not in gene_to_col]
        if missing_genes:
            raise ValueError(
                f"{path.name} missing {len(missing_genes)} expected genes; first={missing_genes[:5]}"
            )
        matrix = matrix[:, [gene_to_col[gene] for gene in target_genes]].tocsr()
        genes = target_genes
    return matrix, genes, wanted_cells or cells, nnz


def _obs_from_metadata(path: Path, cells: Iterable[str]) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", compression="gzip")
    if "Cell.ID" not in meta.columns:
        raise ValueError(f"{path} missing Cell.ID column")
    meta = meta.set_index("Cell.ID", drop=True)
    cells = list(cells)
    missing = [c for c in cells if c not in meta.index]
    if missing:
        raise ValueError(f"Metadata missing {len(missing)} cells; first={missing[:3]}")
    obs = meta.loc[cells].copy()
    # Normalize common names expected by factory modules while preserving the
    # original Trevino columns.
    if "Sample.ID" in obs.columns and "sample" not in obs.columns:
        obs["sample"] = obs["Sample.ID"].astype(str)
    if "Age" in obs.columns and "development_stage" not in obs.columns:
        obs["development_stage"] = obs["Age"].astype(str)
    if "seurat_clusters" in obs.columns:
        obs["trevino_seurat_clusters"] = obs["seurat_clusters"].astype(str)
    obs["public_input_source"] = "GSE162170_scRNA_GEO_processed_counts"
    return obs


def build_public_rna_h5ad(
    input_dir: Path,
    output: Path,
    *,
    chunksize: int,
    max_genes: int | None,
    layers: str,
    layer_gene_policy: str,
    gene_map: Path | None,
) -> dict:
    input_dir = input_dir.resolve()
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    target_gene_ids: list[str] | None = None
    gene_selection_summary: dict[str, object] = {
        "policy": "counts_order",
        "reason": "counts-only or strict layer mode uses the counts matrix gene universe",
    }
    if layers == "all" and layer_gene_policy == "intersection":
        counts_gene_ids = _read_gene_ids_only(input_dir / RNA_COUNTS)
        layer_gene_sets = {
            name: set(_read_gene_ids_only(input_dir / filename))
            for name, filename in [
                ("spliced", RNA_SPLICED),
                ("unspliced", RNA_UNSPLICED),
                ("ambiguous", RNA_AMBIGUOUS),
            ]
        }
        common = set(counts_gene_ids)
        for geneset in layer_gene_sets.values():
            common &= geneset
        target_gene_ids = [gene for gene in counts_gene_ids if gene in common]
        if max_genes is not None:
            target_gene_ids = target_gene_ids[:max_genes]
        gene_selection_summary = {
            "policy": "intersection",
            "counts_gene_count": len(counts_gene_ids),
            "spliced_gene_count": len(layer_gene_sets["spliced"]),
            "unspliced_gene_count": len(layer_gene_sets["unspliced"]),
            "ambiguous_gene_count": len(layer_gene_sets["ambiguous"]),
            "common_gene_count": len(common),
            "selected_gene_count": len(target_gene_ids),
            "counts_genes_dropped": len(set(counts_gene_ids) - common),
            "reason": (
                "velocity layers have a strict subset of counts genes and a "
                "different order; use counts-order intersection for layer-aligned AnnData"
            ),
        }
    elif layers == "all" and layer_gene_policy != "strict":
        raise ValueError(f"Unsupported layer gene policy: {layer_gene_policy!r}")

    counts, gene_ids, cells, counts_nnz = _stream_gene_by_cell_matrix(
        input_dir / RNA_COUNTS,
        chunksize=chunksize,
        expected_genes=target_gene_ids,
        max_genes=max_genes,
    )
    obs = _obs_from_metadata(input_dir / RNA_METADATA, cells)
    symbol_map = _load_gene_symbol_map(gene_map)
    var = _var_from_gene_ids(gene_ids, symbol_map)

    adata = ad.AnnData(X=counts, obs=obs, var=var)
    layer_summaries: dict[str, dict] = {}
    if layers == "all":
        for layer, filename in [
            ("spliced", RNA_SPLICED),
            ("unspliced", RNA_UNSPLICED),
            ("ambiguous", RNA_AMBIGUOUS),
        ]:
            mat, layer_genes, layer_cells, nnz = _stream_gene_by_cell_matrix(
                input_dir / filename,
                chunksize=chunksize,
                wanted_cells=cells,
                expected_genes=gene_ids,
                max_genes=max_genes,
            )
            if layer_cells != cells or layer_genes != gene_ids:
                raise ValueError(f"{filename} did not align with counts")
            adata.layers[layer] = mat
            layer_summaries[layer] = {"filename": filename, "nnz": nnz}
    elif layers == "none":
        layer_summaries["__skipped__"] = {
            "reason": "counts-only viability lane; spliced/unspliced/ambiguous alignment audited separately",
            "files_deferred": [RNA_SPLICED, RNA_UNSPLICED, RNA_AMBIGUOUS],
        }
    else:
        raise ValueError(f"Unsupported layers mode: {layers!r}")

    adata.uns["wave5_public_input_boundary"] = {
        "raw_fastq_public_status": "NOT_FOUND_IN_GEO_BIOPROJECT",
        "source": "GSE162170 GEO supplementary processed matrices",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "max_genes": max_genes,
        "layers_mode": layers,
        "layer_gene_policy": layer_gene_policy,
        "gene_selection": gene_selection_summary,
        "gene_map": str(gene_map.resolve()) if gene_map is not None else None,
        "gene_symbol_mapping": dict(var.attrs.get("gene_symbol_mapping", {})),
    }
    adata.write_h5ad(output, compression="gzip")

    summary = {
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "input_dir": str(input_dir),
        "output_h5ad": str(output),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "X_nnz": counts_nnz,
        "layers_mode": layers,
        "layer_gene_policy": layer_gene_policy,
        "gene_selection": gene_selection_summary,
        "layers": layer_summaries,
        "gene_map": str(gene_map.resolve()) if gene_map is not None else None,
        "gene_symbol_mapping": dict(var.attrs.get("gene_symbol_mapping", {})),
        "obs_columns": list(map(str, adata.obs.columns)),
        "var_index_head": list(map(str, adata.var_names[:10])),
        "var_gene_id_head": list(map(str, adata.var["gene_id"].iloc[:10])),
    }
    output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--chunksize", type=int, default=500)
    ap.add_argument("--max-genes", type=int, default=None)
    ap.add_argument(
        "--gene-map",
        type=Path,
        default=None,
        help=(
            "Optional TSV with gene_id and gene_name columns. When provided, "
            "unique gene symbols become var_names for marker annotation while "
            "gene_id is retained in adata.var."
        ),
    )
    ap.add_argument(
        "--layers",
        choices=["all", "none"],
        default="all",
        help=(
            "Whether to include spliced/unspliced/ambiguous layers. Use 'none' "
            "for the first counts-only pipeline viability lane; velocity layers "
            "remain a separate audited figure-loader task."
        ),
    )
    ap.add_argument(
        "--layer-gene-policy",
        choices=["strict", "intersection"],
        default="strict",
        help=(
            "When --layers all, 'strict' requires velocity-layer genes to match "
            "counts exactly; 'intersection' builds a separate layer-aligned object "
            "on the counts-order intersection of counts/spliced/unspliced/ambiguous."
        ),
    )
    args = ap.parse_args()

    summary = build_public_rna_h5ad(
        args.input_dir,
        args.output,
        chunksize=args.chunksize,
        max_genes=args.max_genes,
        layers=args.layers,
        layer_gene_policy=args.layer_gene_policy,
        gene_map=args.gene_map,
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
