"""Hi-C / scHi-C contact map ingest (Wave 2B / P2.S18).

Memory-safe: stores contact matrices as scipy.sparse CSR (no dense
materialization). Per-cell scHi-C is optional — if not provided, a single
bulk/pseudo-bulk contact matrix is ingested.

Inputs (via cfg):
  --hic-contacts-path PATH      .cool, .mcool, or .tsv.gz contact pairs
  --hic-resolution-bp INT       bin size in bp (default 25_000 / 25 kb)
  --hic-chromsizes-path PATH    optional chrom.sizes (chrom, length); inferred when omitted

Outputs:
  adata.uns["hic_contact_matrix"]   scipy.sparse.csr_matrix of inter-bin contacts (bins x bins)
  adata.uns["hic_bins"]             DataFrame: bin_id, chrom, start, end
  adata.uns["hic_ingest_metadata"]  resolution, n_bins, total_contacts, sparsity
  runs/<run-id>/hic_ingest/hic_summary.json
"""
from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "Lieberman_Aiden_HiC_2009": {
        "title": "Comprehensive Mapping of Long-Range Interactions Reveals Folding Principles of the Human Genome",
        "authors": "Lieberman-Aiden et al.",
        "journal": "Science",
        "year": "2009",
        "doi": "10.1126/science.1181369",
        "description": "Original Hi-C methodology — contact-map binning conventions adopted here.",
    },
    "Abdennur_cooler_2020": {
        "title": "Cooler: scalable storage for Hi-C data and other genomically labeled arrays",
        "authors": "Abdennur, Mirny",
        "journal": "Bioinformatics",
        "year": "2020",
        "doi": "10.1093/bioinformatics/btz540",
        "description": "Canonical .cool / .mcool sparse storage format. This module reads cooler files via the optional `cooler` Python package when installed; falls back to TSV for portability.",
    },
}


def _load_chromsizes(path: Path | None) -> pd.DataFrame | None:
    if path is None or not path.exists():
        return None
    return pd.read_csv(path, sep="\t", header=None, names=["chrom", "length"],
                       usecols=[0, 1], dtype={"chrom": str, "length": int})


def _build_bins(chromsizes: pd.DataFrame, resolution: int) -> pd.DataFrame:
    chroms: list[str] = []
    starts: list[np.ndarray] = []
    for _, r in chromsizes.iterrows():
        s = np.arange(0, int(r["length"]), resolution, dtype=int)
        chroms.append(np.repeat(np.asarray(r["chrom"], dtype=object), s.shape[0]))
        starts.append(s)
    if not starts:
        return pd.DataFrame(columns=["bin_id", "chrom", "start", "end"])
    chrom_col = np.concatenate(chroms)
    start_col = np.concatenate(starts)
    # ``length`` is per-chromosome; align it to each bin so the end is clamped
    # to chromosome length exactly as the row-wise version did.
    length_per_chrom = chromsizes.set_index("chrom")["length"]
    length_col = length_per_chrom.reindex(chrom_col).to_numpy().astype(int)
    end_col = np.minimum(start_col + resolution, length_col)
    return pd.DataFrame({
        "bin_id": np.arange(start_col.shape[0], dtype=int),
        "chrom": chrom_col,
        "start": start_col.astype(int),
        "end": end_col.astype(int),
    })


def _load_contacts_tsv(path: Path) -> pd.DataFrame:
    """Read tab-separated contact pairs: chrom1, pos1, chrom2, pos2, count."""
    opener = gzip.open if str(path).endswith(".gz") else open
    return pd.read_csv(
        opener(path, "rt"),
        sep="\t",
        header=None,
        names=["chrom1", "pos1", "chrom2", "pos2", "count"],
        usecols=[0, 1, 2, 3, 4],
        dtype={"chrom1": str, "pos1": int, "chrom2": str, "pos2": int, "count": float},
    )


def _binize_contacts(contacts: pd.DataFrame, bins: pd.DataFrame, resolution: int) -> sp.csr_matrix:
    """Sum contact counts into bin × bin sparse matrix.

    Vectorized equivalent of a per-row ``iterrows`` accumulation: contacts whose
    endpoints both map to a known bin contribute ``(b1, b2, count)`` and, when
    ``b1 != b2``, also the mirrored ``(b2, b1, count)``. Duplicate ``(row, col)``
    pairs are summed by ``coo_matrix`` exactly as the original ``csr_matrix``
    constructor did, and the diagonal is never double-counted.
    """
    n_bins = int(bins.shape[0])
    # Bin lookup keyed by (chrom, pos // resolution); identical key construction
    # to the row-wise version (bin starts are multiples of ``resolution``).
    bin_key = bins["start"].to_numpy() // resolution
    lookup = pd.Series(
        bins["bin_id"].to_numpy().astype(np.int64),
        index=pd.MultiIndex.from_arrays([bins["chrom"].to_numpy(), bin_key]),
    )

    if contacts.shape[0] == 0:
        return sp.csr_matrix((n_bins, n_bins), dtype=np.float32)

    key1 = pd.MultiIndex.from_arrays(
        [contacts["chrom1"].to_numpy(), contacts["pos1"].to_numpy() // resolution]
    )
    key2 = pd.MultiIndex.from_arrays(
        [contacts["chrom2"].to_numpy(), contacts["pos2"].to_numpy() // resolution]
    )
    b1 = lookup.reindex(key1).to_numpy()  # NaN where bin missing
    b2 = lookup.reindex(key2).to_numpy()
    count = contacts["count"].to_numpy().astype(np.float64)

    valid = ~(np.isnan(b1) | np.isnan(b2))
    b1 = b1[valid].astype(np.int64)
    b2 = b2[valid].astype(np.int64)
    count = count[valid]

    off_diag = b1 != b2
    rows = np.concatenate([b1, b2[off_diag]])
    cols = np.concatenate([b2, b1[off_diag]])
    data = np.concatenate([count, count[off_diag]])

    coo = sp.coo_matrix((data, (rows, cols)), shape=(n_bins, n_bins), dtype=np.float32)
    return coo.tocsr()


class HiCIngestModule:
    """Hi-C contact map ingest. Stores contacts as scipy.sparse CSR."""

    name = "hic_ingest"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {
        "uns": ["hic_contact_matrix", "hic_bins", "hic_ingest_metadata"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata
        contacts_path = getattr(ctx.cfg, "hic_contacts_path", None)
        if not contacts_path:
            ctx.status(self.name, "skipped", "hic_contacts_path not set")
            ctx.metadata["hic_ingest_status"] = "skipped_no_input"
            return

        contacts_path = Path(contacts_path)
        if not contacts_path.exists():
            ctx.status(self.name, "skipped", f"contacts missing: {contacts_path}")
            ctx.metadata["hic_ingest_status"] = "skipped_path_missing"
            return

        resolution = int(getattr(ctx.cfg, "hic_resolution_bp", 25_000))
        chromsizes_path = getattr(ctx.cfg, "hic_chromsizes_path", None)

        if str(contacts_path).endswith(".hic"):
            ctx.status(
                self.name,
                "skipped",
                ".hic ingest is not interchangeable with .cool/.mcool/TSV; convert or supply a validated processed table",
            )
            ctx.metadata["hic_ingest_status"] = "skipped_hic_unsupported"
            return

        # Cool/.mcool support requires `cooler`; fall back to TSV
        if str(contacts_path).endswith((".cool", ".mcool")):
            try:
                import cooler
            except ImportError:
                ctx.status(self.name, "skipped", "cooler not installed; pass --hic-contacts-path as TSV to use the fallback")
                ctx.metadata["hic_ingest_status"] = "skipped_cooler_missing"
                return
            c = cooler.Cooler(str(contacts_path))
            mat = c.matrix(balance=False, sparse=True)[:]
            bins_df = c.bins()[:].reset_index().rename(columns={"index": "bin_id"})
            n_bins = mat.shape[0]
            total_contacts = float(mat.sum())
            input_format = "cooler_unbalanced"
        else:
            # TSV path: load contacts + build bins
            chromsizes = _load_chromsizes(Path(chromsizes_path) if chromsizes_path else None)
            contacts = _load_contacts_tsv(contacts_path)
            if chromsizes is None:
                # Infer chromsizes from observed positions
                chrom_max = (
                    pd.concat([contacts[["chrom1", "pos1"]].rename(columns={"chrom1": "chrom", "pos1": "pos"}),
                               contacts[["chrom2", "pos2"]].rename(columns={"chrom2": "chrom", "pos2": "pos"})])
                    .groupby("chrom", as_index=False)["pos"].max()
                )
                chromsizes = chrom_max.rename(columns={"pos": "length"})
                chromsizes["length"] = chromsizes["length"] + resolution  # add a tail bin
            bins_df = _build_bins(chromsizes, resolution)
            mat = _binize_contacts(contacts, bins_df, resolution)
            n_bins = mat.shape[0]
            total_contacts = float(contacts["count"].sum())
            input_format = "tsv_contact_pairs"

        adata.uns["hic_contact_matrix"] = mat
        adata.uns["hic_bins"] = bins_df
        sparsity = 1.0 - float(mat.nnz) / max(1, n_bins * n_bins)
        adata.uns["hic_ingest_metadata"] = {
            "resolution_bp": resolution,
            "n_bins": int(n_bins),
            "total_contacts": float(total_contacts),
            "sparsity": float(sparsity),
            "matrix_format": "csr_sparse",
            "input_format": input_format,
            "normalization": "raw_counts_or_unbalanced",
        }
        out_dir = ctx.run_dir / "hic_ingest"
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "hic_summary.json").write_text(
            json.dumps(adata.uns["hic_ingest_metadata"], indent=2), encoding="utf-8"
        )

        ctx.metadata["hic_ingest_status"] = "ok"
        ctx.metadata["hic_n_bins"] = int(n_bins)
        logger.info("%s: %d bins (resolution=%d bp), %.0f total contacts, sparsity=%.4f",
                    self.name, n_bins, resolution, total_contacts, sparsity)
