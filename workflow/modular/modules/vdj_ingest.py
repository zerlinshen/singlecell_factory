"""VDJ ingest (Wave 2B / P2.S15): cellranger VDJ output → adata.obs annotation.

Memory-safe (pandas DataFrame I/O only; no AnnData duplication).
Fail-soft: skips gracefully when no VDJ input is provided or input is empty.

Inputs (via cfg):
  --vdj-contigs-path PATH   path to filtered_contig_annotations.csv (cellranger-vdj)
  --vdj-clonotypes-path PATH (optional) clonotypes.csv from cellranger;
                              when provided, clonotype IDs come from this file
                              directly instead of being derived from contigs.

Outputs:
  adata.obs["clonotype_id"]      string; "" when cell has no VDJ contig
  adata.obs["chain_pairing"]     string; e.g. "TRA+TRB", "IGH+IGK", "ambiguous", "no_chains"
  adata.uns["vdj_clonotypes"]    DataFrame: clonotype_id, n_cells, v_genes, j_genes, cdr3_nts
  runs/<run-id>/vdj_ingest/vdj_summary.json
"""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "Sturm_scirpy_2020": {
        "title": "Scirpy: a Scanpy extension for analyzing single-cell T-cell receptor-sequencing data",
        "authors": "Sturm et al.",
        "journal": "Bioinformatics",
        "year": "2020",
        "doi": "10.1093/bioinformatics/btaa611",
        "description": "Reference framework for AnnData-native VDJ analysis. This module follows scirpy's chain-pairing + clonotype semantics without taking a hard dependency on the package (so the pipeline runs in scirpy-less environments).",
    },
    "Bagaev_VDJtools_2015": {
        "title": "VDJtools: Unifying Post-analysis of T Cell Receptor Repertoires",
        "authors": "Bagaev et al.",
        "journal": "PLoS Computational Biology",
        "year": "2015",
        "doi": "10.1371/journal.pcbi.1004503",
        "description": "Canonical VDJ repertoire analysis conventions adopted here for clonotype definition (V-gene + J-gene + CDR3 nt).",
    },
    "10x_VDJ_software": {
        "title": "Cell Ranger V(D)J pipeline",
        "authors": "10x Genomics",
        "journal": "Software documentation",
        "year": "2024",
        "doi": "https://support.10xgenomics.com/single-cell-vdj/software",
        "description": "Authoritative spec for filtered_contig_annotations.csv / clonotypes.csv schemas consumed here.",
    },
}


# Canonical chain pair signatures (TCR + BCR).
_PAIR_TR = ("TRA", "TRB")
_PAIR_TG = ("TRG", "TRD")
_PAIR_IGK = ("IGH", "IGK")
_PAIR_IGL = ("IGH", "IGL")


def _label_chain_pairing(chains: set[str]) -> str:
    """Return a canonical chain-pairing label for the chains observed on a barcode."""
    if not chains:
        return "no_chains"
    if chains == set(_PAIR_TR):
        return "TRA+TRB"
    if chains == set(_PAIR_TG):
        return "TRG+TRD"
    if chains == set(_PAIR_IGK):
        return "IGH+IGK"
    if chains == set(_PAIR_IGL):
        return "IGH+IGL"
    if len(chains) == 1:
        return f"orphan_{next(iter(chains))}"
    return "ambiguous"


def _derive_clonotype_id(rows: pd.DataFrame) -> str:
    """Stable clonotype id from sorted (v_gene, j_gene, cdr3_nt) tuples."""
    keys = sorted(
        (str(r.get("v_gene", "")), str(r.get("j_gene", "")), str(r.get("cdr3_nt", "") or r.get("cdr3", "")))
        for _, r in rows.iterrows()
    )
    digest = hashlib.sha1(repr(keys).encode("utf-8")).hexdigest()[:12]
    return f"clonotype_{digest}"


class VDJIngestModule:
    """Parse cellranger VDJ contigs into per-cell clonotype + chain-pairing annotations."""

    name = "vdj_ingest"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {
        "obs": ["clonotype_id", "chain_pairing"],
        "uns": ["vdj_clonotypes"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata
        contigs_path = getattr(ctx.cfg, "vdj_contigs_path", None)
        clonotypes_path = getattr(ctx.cfg, "vdj_clonotypes_path", None)

        if not contigs_path:
            ctx.status(self.name, "skipped", "vdj_contigs_path not set")
            ctx.metadata["vdj_ingest_status"] = "skipped_no_input"
            return

        contigs_path = Path(contigs_path)
        if not contigs_path.exists():
            ctx.status(self.name, "skipped", f"vdj_contigs_path missing: {contigs_path}")
            ctx.metadata["vdj_ingest_status"] = "skipped_path_missing"
            return

        logger.info("%s: parsing VDJ contigs from %s", self.name, contigs_path)
        contigs = pd.read_csv(contigs_path)
        # cellranger-vdj filtered_contig_annotations.csv columns of interest:
        # barcode, is_cell, productive, chain, v_gene, j_gene, cdr3, cdr3_nt, raw_clonotype_id
        keep_cols = {"barcode", "chain", "v_gene", "j_gene", "cdr3_nt", "cdr3", "is_cell", "productive", "raw_clonotype_id"}
        missing_cols = keep_cols - set(contigs.columns)
        if "barcode" in missing_cols or "chain" in missing_cols:
            ctx.status(self.name, "skipped", f"contigs missing required cols: {missing_cols}")
            ctx.metadata["vdj_ingest_status"] = "skipped_bad_schema"
            return

        # Productive cells only (per cellranger convention).
        if "is_cell" in contigs.columns:
            contigs = contigs[contigs["is_cell"].astype(str).str.lower().isin({"true", "1", "yes"})]
        if "productive" in contigs.columns:
            contigs = contigs[contigs["productive"].astype(str).str.lower().isin({"true", "1", "yes"})]

        # Chain-pairing per barcode.
        per_barcode_chains = contigs.groupby("barcode")["chain"].apply(lambda s: set(s.dropna().astype(str)))
        chain_pairing = per_barcode_chains.map(_label_chain_pairing)

        # Clonotype id: prefer cellranger's raw_clonotype_id if present, else derive from V+J+CDR3.
        clonotype_id_per_barcode: dict[str, str] = {}
        if clonotypes_path and Path(clonotypes_path).exists():
            clonotypes_df = pd.read_csv(clonotypes_path)
            if "clonotype_id" in clonotypes_df.columns and "frequency" in clonotypes_df.columns:
                # Map each barcode to its raw_clonotype_id via the contigs table.
                if "raw_clonotype_id" in contigs.columns:
                    clonotype_id_per_barcode = (
                        contigs.dropna(subset=["raw_clonotype_id"])
                        .drop_duplicates(subset=["barcode"])
                        .set_index("barcode")["raw_clonotype_id"]
                        .astype(str)
                        .to_dict()
                    )
        if not clonotype_id_per_barcode:
            # Derive from V+J+CDR3 hash per barcode (canonical VDJtools convention).
            for bc, group in contigs.groupby("barcode"):
                clonotype_id_per_barcode[bc] = _derive_clonotype_id(group)

        # Build per-clonotype summary.
        clonotype_df_rows: list[dict] = []
        for bc, cid in clonotype_id_per_barcode.items():
            grp = contigs[contigs["barcode"] == bc]
            v_genes = sorted(set(grp.get("v_gene", pd.Series(dtype=str)).dropna()))
            j_genes = sorted(set(grp.get("j_gene", pd.Series(dtype=str)).dropna()))
            cdr3_nts = sorted(set(grp.get("cdr3_nt", grp.get("cdr3", pd.Series(dtype=str))).dropna()))
            clonotype_df_rows.append({
                "clonotype_id": cid,
                "barcode": bc,
                "v_genes": ",".join(map(str, v_genes)),
                "j_genes": ",".join(map(str, j_genes)),
                "cdr3_nts": ",".join(map(str, cdr3_nts)),
            })
        clonotypes_per_barcode_df = pd.DataFrame(clonotype_df_rows)
        if not clonotypes_per_barcode_df.empty:
            clonotype_summary = (
                clonotypes_per_barcode_df
                .groupby("clonotype_id", as_index=False)
                .agg(
                    n_cells=("barcode", "count"),
                    v_genes=("v_genes", lambda s: ",".join(sorted({g for sub in s for g in sub.split(",") if g}))),
                    j_genes=("j_genes", lambda s: ",".join(sorted({g for sub in s for g in sub.split(",") if g}))),
                    cdr3_nts=("cdr3_nts", lambda s: ",".join(sorted({g for sub in s for g in sub.split(",") if g}))),
                )
                .sort_values("n_cells", ascending=False)
                .reset_index(drop=True)
            )
        else:
            clonotype_summary = pd.DataFrame(columns=["clonotype_id", "n_cells", "v_genes", "j_genes", "cdr3_nts"])

        # Map onto adata.obs: cells absent from VDJ get empty string.
        adata.obs["clonotype_id"] = adata.obs_names.to_series().map(
            lambda bc: clonotype_id_per_barcode.get(bc, "")
        ).astype(str)
        adata.obs["chain_pairing"] = adata.obs_names.to_series().map(
            lambda bc: chain_pairing.get(bc, "no_chains")
        ).astype(str)
        adata.uns["vdj_clonotypes"] = clonotype_summary

        # Write per-run summary JSON
        out_dir = ctx.run_dir / "vdj_ingest"
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "vdj_summary.json").write_text(
            json.dumps({
                "n_cells_with_vdj": int((adata.obs["clonotype_id"] != "").sum()),
                "n_unique_clonotypes": int(clonotype_summary.shape[0]),
                "chain_pairing_counts": adata.obs["chain_pairing"].value_counts().to_dict(),
                "contigs_path": str(contigs_path),
            }, indent=2, default=str),
            encoding="utf-8",
        )
        ctx.metadata["vdj_ingest_status"] = "ok"
        ctx.metadata["vdj_n_clonotypes"] = int(clonotype_summary.shape[0])
        ctx.metadata["vdj_n_cells_with_vdj"] = int((adata.obs["clonotype_id"] != "").sum())
        logger.info(
            "%s: %d cells with VDJ → %d unique clonotypes",
            self.name, ctx.metadata["vdj_n_cells_with_vdj"], ctx.metadata["vdj_n_clonotypes"],
        )
