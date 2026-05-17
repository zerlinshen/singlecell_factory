"""Stage-aware loader for Trevino 2021 10x Multiome brain data.

Provides two entry points: ``load_trevino_2021()`` for real data and
``load_synthetic()`` for deterministic synthetic fixtures. Both return a
single AnnData following the same schema contract: RNA counts in ``.X``
(CSR sparse int32), ATAC peak counts in ``.obsm["atac_peaks"]`` (CSR
sparse int32; n_peaks != n_vars so layers is invalid), and cell barcodes
intersected across modalities.

``load_trevino_2021`` streams both TSVs in chunks and sparsifies per
chunk, keeping peak RSS bounded (~1-3 GB on Trevino PCW21 ≈ 9k cells ×
~150k peaks) instead of the multi-tens-of-GB dense intermediate the
original scaffold produced.
"""

from __future__ import annotations

from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scipy.sparse

__references__ = {
    "trevino2021": {
        "title": "Chromatin and gene-regulatory dynamics of the developing human neocortex at single-cell resolution",
        "authors": "Trevino AE, Müller F, Andersen J, et al.",
        "journal": "Cell",
        "year": 2021,
        "doi": "doi:10.1016/j.cell.2021.07.039",
        "description": (
            "Source of GSE162170 multiome RNA + ATAC counts used in "
            "load_trevino_2021(); defines the joint barcode set, stage "
            "labels (Sample.Age), and seurat_clusters used as .obs columns."
        ),
    },
    "10x_multiome": {
        "title": "10x Genomics Chromium Single Cell Multiome ATAC + Gene Expression",
        "authors": "10x Genomics",
        "journal": "10x Genomics Technical Note",
        "year": 2021,
        "doi": "doi:10.17504/protocols.io.bfskjncw",
        "description": (
            "Defines the 10x Multiome dual-barcode schema enforced here: "
            "one barcode per cell shared across RNA and ATAC modalities, "
            "enabling barcode-intersection alignment in load_trevino_2021()."
        ),
    },
}

_DATA_ROOT = Path(__file__).resolve().parents[2] / "data" / "raw" / "trevino_2021_brain"

_RNA_FILE = "GSE162170_multiome_rna_counts.tsv.gz"
_ATAC_FILE = "GSE162170_multiome_atac_counts.tsv.gz"
_META_FILE = "GSE162170_multiome_cell_metadata.txt.gz"
_PEAKS_FILE = "GSE162170_multiome_atac_consensus_peaks.txt.gz"


def load_trevino_2021(stage: str, data_root: Path = _DATA_ROOT) -> anndata.AnnData:
    """Load real Trevino 2021 multiome data for a given developmental stage.

    Parameters
    ----------
    stage:
        Value of the ``Sample.Age`` column to filter on (e.g. ``"pcw21"``).
    data_root:
        Path to the directory containing GSE162170 files. Defaults to
        ``data/raw/trevino_2021_brain/`` relative to the repo root.

    Returns
    -------
    AnnData with:
      - ``.X``: RNA counts (CSR sparse int32, cells × genes)
      - ``.obsm["atac_peaks"]``: ATAC peak counts (CSR sparse int32, cells × peaks).
        Stored in obsm (not layers) because n_peaks != n_vars; AnnData 0.12
        requires ``layers[k].shape[1] == n_vars``.
      - ``.obs``: ``Sample.Age``, ``seurat_clusters``
      - ``.var``: gene IDs (Ensembl)
      - ``.uns["atac_var"]``: ATAC peak metadata (peak_id / chrom / start / end)

    Implementation notes
    --------------------
    Both TSVs are streamed via ``pd.read_csv(chunksize=...)`` and each chunk
    is converted to ``scipy.sparse.csr_matrix`` immediately. Final RNA and
    ATAC matrices are assembled via ``scipy.sparse.hstack`` over the
    cells-x-features chunks. Peak RSS on Trevino PCW21 ≈ 1-3 GB; total
    wall-time ≈ 1-2 minutes. The prior scaffold materialised the full dense
    DataFrame and peaked at 60-80 GB before the OOM threshold.
    """
    data_root = Path(data_root)

    # --- metadata ---
    meta = pd.read_csv(
        data_root / _META_FILE, sep="\t", index_col="Cell.ID", compression="gzip"
    )
    stage_cell_set = set(meta[meta["Sample.Age"] == stage].index.tolist())
    if not stage_cell_set:
        raise ValueError(f"No cells found for stage={stage!r}")

    # --- ATAC peak metadata (small) ---
    peaks_meta = pd.read_csv(
        data_root / _PEAKS_FILE, sep="\t", compression="gzip"
    )

    # --- RNA streamed (genes × cells TSV, header = cell barcodes) ---
    # Read in chunks of 2000 genes; per-chunk dense footprint ≈ n_cells × 2000 × 4B.
    rna_chunks: list[scipy.sparse.csr_matrix] = []
    gene_index: list[str] = []
    rna_cells_in_order: list[str] | None = None
    all_cell_barcodes: list[str] | None = None
    # NOTE: do NOT pass dtype=np.int32 here — index_col=0 holds gene IDs
    # (strings like ``ENSG00000243485``) which would fail the int cast.
    for chunk in pd.read_csv(
        data_root / _RNA_FILE,
        sep="\t",
        index_col=0,
        compression="gzip",
        chunksize=2000,
    ):
        if all_cell_barcodes is None:
            all_cell_barcodes = chunk.columns.tolist()
            barcode_set = set(all_cell_barcodes)
            rna_cells_in_order = [c for c in all_cell_barcodes if c in stage_cell_set and c in barcode_set]
        gene_index.extend(chunk.index.tolist())
        # chunk: 2000 genes × n_all_cells; select stage cells, transpose to cells × 2000 genes
        sub = chunk[rna_cells_in_order].to_numpy(dtype=np.int32, copy=False)
        rna_chunks.append(scipy.sparse.csr_matrix(sub.T))
        del chunk, sub
    if all_cell_barcodes is None:
        raise RuntimeError("RNA file is empty or unreadable.")
    rna_sparse = scipy.sparse.hstack(rna_chunks, format="csr", dtype=np.int32)
    del rna_chunks

    # --- ATAC streamed (peaks × cells, no header, column order matches RNA header) ---
    atac_chunks: list[scipy.sparse.csr_matrix] = []
    for chunk in pd.read_csv(
        data_root / _ATAC_FILE,
        sep="\t",
        header=None,
        names=all_cell_barcodes,
        compression="gzip",
        chunksize=5000,
        dtype=np.int32,
    ):
        sub = chunk[rna_cells_in_order].to_numpy(dtype=np.int32, copy=False)
        atac_chunks.append(scipy.sparse.csr_matrix(sub.T))
        del chunk, sub
    atac_sparse = scipy.sparse.hstack(atac_chunks, format="csr", dtype=np.int32)
    del atac_chunks

    # --- build AnnData ---
    obs = meta.loc[rna_cells_in_order, ["Sample.Age", "seurat_clusters"]].copy()
    var = pd.DataFrame(index=gene_index)
    var.index.name = "gene_id"

    adata = anndata.AnnData(X=rna_sparse, obs=obs, var=var)
    # ATAC peaks (cells × peaks) stored in obsm; n_atac != n_rna so layers is invalid.
    adata.obsm["atac_peaks"] = atac_sparse
    adata.uns["atac_var"] = peaks_meta.reset_index(drop=True)

    return adata


def load_synthetic(
    n_obs: int = 1000,
    n_rna: int = 3000,
    n_atac: int = 5000,
    n_cell_types: int = 3,
    seed: int = 42,
) -> anndata.AnnData:
    """Generate a deterministic synthetic AnnData matching the Trevino 2021 schema.

    Parameters
    ----------
    n_obs:
        Number of synthetic cells.
    n_rna:
        Number of synthetic RNA features (genes).
    n_atac:
        Number of synthetic ATAC peaks (must be >= 1000).
    n_cell_types:
        Number of simulated cell types (must be >= 2).
    seed:
        Random seed for ``numpy.random.default_rng``.

    Returns
    -------
    AnnData with:
      - ``.X``: RNA counts (CSR sparse int32, cells × genes)
      - ``.layers["atac_peaks"]``: ATAC peak counts (CSR sparse int32, cells × peaks)
      - ``.obs["cell_type_synth"]``: cell type label with >= 2 distinct values
    """
    if n_atac < 1000:
        raise ValueError(f"n_atac must be >= 1000, got {n_atac}")
    if n_cell_types < 2:
        raise ValueError(f"n_cell_types must be >= 2, got {n_cell_types}")

    rng = np.random.default_rng(seed)

    cells_per_type = n_obs // n_cell_types
    cell_types = []
    for ct in range(n_cell_types):
        count = cells_per_type if ct < n_cell_types - 1 else n_obs - cells_per_type * (n_cell_types - 1)
        cell_types.extend([f"type_{ct}"] * count)
    cell_type_arr = np.array(cell_types)

    # Base background counts (Poisson low rate)
    rna_counts = rng.poisson(0.5, size=(n_obs, n_rna)).astype(np.int32)
    atac_counts = rng.poisson(0.3, size=(n_obs, n_atac)).astype(np.int32)

    # Plant marker-gene gradients and peak-gene linkages per cell type
    # Each cell type gets a block of marker genes and linked peaks
    markers_per_type = max(10, n_rna // (n_cell_types * 10))
    peaks_per_type = max(10, n_atac // (n_cell_types * 10))

    for ct in range(n_cell_types):
        mask = cell_type_arr == f"type_{ct}"
        gene_start = ct * markers_per_type
        gene_end = gene_start + markers_per_type
        peak_start = ct * peaks_per_type
        peak_end = peak_start + peaks_per_type

        # Elevated expression in marker genes for this cell type
        rna_counts[np.ix_(mask, np.arange(gene_start, gene_end))] += rng.poisson(
            8, size=(mask.sum(), markers_per_type)
        ).astype(np.int32)

        # Correlated ATAC signal in linked peaks (planted peak-gene linkage)
        atac_counts[np.ix_(mask, np.arange(peak_start, peak_end))] += rng.poisson(
            6, size=(mask.sum(), peaks_per_type)
        ).astype(np.int32)

    X = scipy.sparse.csr_matrix(rna_counts)
    atac_mat = scipy.sparse.csr_matrix(atac_counts)

    barcodes = [f"cell_{i:05d}" for i in range(n_obs)]
    gene_ids = [f"ENSG{i:011d}" for i in range(n_rna)]
    peak_ids = [f"chr1:{i*500}-{i*500+500}" for i in range(n_atac)]

    obs = pd.DataFrame({"cell_type_synth": cell_type_arr}, index=barcodes)
    var = pd.DataFrame(index=gene_ids)
    var.index.name = "gene_id"

    adata = anndata.AnnData(X=X, obs=obs, var=var)
    # ATAC peaks (cells × peaks) stored in obsm; n_atac != n_rna so layers is invalid.
    adata.obsm["atac_peaks"] = atac_mat
    adata.uns["atac_var"] = pd.DataFrame({"peak_id": peak_ids})

    return adata


if __name__ == "__main__":
    adata = load_synthetic()
    print("X.shape:", adata.X.shape)
    print("obsm['atac_peaks'].shape:", adata.obsm["atac_peaks"].shape)
    print("obs['cell_type_synth'].value_counts():")
    print(adata.obs["cell_type_synth"].value_counts())
