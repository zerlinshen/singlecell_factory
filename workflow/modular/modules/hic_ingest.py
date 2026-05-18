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
    rows = []
    bin_id = 0
    for _, r in chromsizes.iterrows():
        starts = np.arange(0, int(r["length"]), resolution, dtype=int)
        for s in starts:
            rows.append({
                "bin_id": bin_id,
                "chrom": r["chrom"],
                "start": int(s),
                "end": int(min(s + resolution, r["length"])),
            })
            bin_id += 1
    return pd.DataFrame(rows)


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
    """Sum contact counts into bin × bin sparse matrix."""
    bin_lookup: dict[tuple[str, int], int] = {}
    for _, r in bins.iterrows():
        bin_lookup[(r["chrom"], int(r["start"]) // resolution)] = int(r["bin_id"])

    n_bins = bins.shape[0]
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for _, c in contacts.iterrows():
        b1 = bin_lookup.get((c["chrom1"], int(c["pos1"]) // resolution))
        b2 = bin_lookup.get((c["chrom2"], int(c["pos2"]) // resolution))
        if b1 is None or b2 is None:
            continue
        rows.append(b1)
        cols.append(b2)
        data.append(float(c["count"]))
        # Symmetric: also store (b2, b1) when b1 != b2
        if b1 != b2:
            rows.append(b2)
            cols.append(b1)
            data.append(float(c["count"]))
    return sp.csr_matrix((data, (rows, cols)), shape=(n_bins, n_bins), dtype=np.float32)


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

        adata.uns["hic_contact_matrix"] = mat
        adata.uns["hic_bins"] = bins_df
        sparsity = 1.0 - float(mat.nnz) / max(1, n_bins * n_bins)
        adata.uns["hic_ingest_metadata"] = {
            "resolution_bp": resolution,
            "n_bins": int(n_bins),
            "total_contacts": float(total_contacts),
            "sparsity": float(sparsity),
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
