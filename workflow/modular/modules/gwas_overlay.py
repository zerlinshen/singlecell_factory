"""GWAS catalog SNP overlay onto ATAC peaks (Wave 5.1 / MV3).

# STATUS: promote-on-reuse

Overlays NHGRI-EBI GWAS Catalog SNP associations onto the ATAC peak BED.
For each cell, computes a GWAS SNP burden score = sum of overlapping SNPs
in accessible peaks per cell (weighted by accessibility). Output is used
to render F-7 (UMAP coloured by GWAS SNP burden).

Inputs (via cfg):
  --gwas-catalog-path PATH        Path to NHGRI-EBI gwas_catalog_v1.0-associations TSV
                                    (downloaded ahead-of-time; SHA-recorded in
                                    data/external/gwas_catalog/MANIFEST.json)
  --gwas-trait-filter STR         Optional substring filter on DISEASE/TRAIT column
                                    (e.g., 'schizophrenia' to subset; default: all)
  --gwas-p-value-threshold FLOAT  Filter SNPs by P-VALUE column (default: 5e-8)

Outputs:
  <run-dir>/gwas_overlay/snp_peak_overlap.csv   peak_idx, chrom, pos, snp, trait, p_value
  <run-dir>/gwas_overlay/per_cell_burden.csv    cell, gwas_snp_burden
  adata.obs['gwas_snp_burden']                  per-cell GWAS SNP burden score (float)
  adata.uns['gwas_snp_peak_overlap']            DataFrame of overlap table
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse

from ..context import PipelineContext
from .._contract_violation import ModuleContractError

logger = logging.getLogger(__name__)


__references__ = {
    "Sollis_GWAS_Catalog_2023": {
        "title": "The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource",
        "authors": "Sollis E, Mosaku A, Abid A, et al.",
        "journal": "Nucleic Acids Research",
        "year": "2023",
        "doi": "10.1093/nar/gkac1010",
        "description": "NHGRI-EBI GWAS Catalog: source of SNP-trait associations consumed by this module.",
    },
    "Trevino_2021": {
        "title": "Chromatin and gene-regulatory dynamics of the developing human neocortex at single-cell resolution",
        "authors": "Trevino AE, Müller F, Andersen J, et al.",
        "journal": "Cell",
        "year": "2021",
        "doi": "10.1016/j.cell.2021.07.039",
        "description": "Source paper for F-7 GWAS overlay reproduction (Wave-5.6 MV3).",
    },
}


def _load_gwas_catalog(
    path: Path,
    trait_filter: str | None,
    p_threshold: float,
) -> pd.DataFrame:
    """Load the GWAS catalog TSV; keep CHR_ID, CHR_POS, SNPS, DISEASE/TRAIT, P-VALUE."""
    df = pd.read_csv(
        path, sep="\t", dtype=str, low_memory=False, on_bad_lines="skip",
    )
    needed = {"CHR_ID", "CHR_POS", "SNPS", "DISEASE/TRAIT", "P-VALUE"}
    missing = needed - set(df.columns)
    if missing:
        raise ModuleContractError(
            f"gwas_overlay: catalog missing required columns {sorted(missing)}; "
            f"got {list(df.columns)[:20]}"
        )
    df = df.rename(columns={"DISEASE/TRAIT": "trait", "P-VALUE": "p_value"})
    df = df[["CHR_ID", "CHR_POS", "SNPS", "trait", "p_value"]].copy()
    # Coerce types; drop rows with no chromosome/position.
    df["CHR_ID"] = df["CHR_ID"].astype(str).str.strip()
    df = df[df["CHR_ID"].notna() & (df["CHR_ID"] != "") & (df["CHR_ID"] != "nan")]
    df["CHR_POS"] = pd.to_numeric(df["CHR_POS"], errors="coerce")
    df = df.dropna(subset=["CHR_POS"])
    df["CHR_POS"] = df["CHR_POS"].astype(np.int64)
    df["p_value_num"] = pd.to_numeric(df["p_value"], errors="coerce")
    df = df[df["p_value_num"].fillna(1.0) <= p_threshold]
    if trait_filter:
        df = df[df["trait"].astype(str).str.contains(trait_filter, case=False, na=False)]
    df["chrom"] = df["CHR_ID"].apply(lambda c: c if c.startswith("chr") else f"chr{c}")
    return df.reset_index(drop=True)


def _overlap_snps_with_peaks(
    snps_df: pd.DataFrame,
    peaks_df: pd.DataFrame,
) -> pd.DataFrame:
    """Return DataFrame with columns: peak_idx, chrom, pos, snp, trait, p_value.

    peaks_df must have chrom/start/end columns; peak_idx is its row index.
    Per-chrom merge using sorted intervals + searchsorted (no full join).
    """
    records: list[dict] = []
    snps_by_chrom = {c: g for c, g in snps_df.groupby("chrom")}
    for chrom, peak_group in peaks_df.groupby("chrom"):
        if chrom not in snps_by_chrom:
            continue
        snp_group = snps_by_chrom[chrom]
        peak_starts = peak_group["start"].to_numpy()
        peak_ends = peak_group["end"].to_numpy()
        peak_idx_arr = peak_group.index.to_numpy()
        # Sort peaks by start
        order = np.argsort(peak_starts)
        ps_sorted = peak_starts[order]
        pe_sorted = peak_ends[order]
        pidx_sorted = peak_idx_arr[order]
        for _, snp in snp_group.iterrows():
            pos = int(snp["CHR_POS"])
            # find candidate peak with start <= pos
            j = np.searchsorted(ps_sorted, pos, side="right") - 1
            if 0 <= j < len(ps_sorted) and pe_sorted[j] >= pos:
                records.append({
                    "peak_idx": int(pidx_sorted[j]),
                    "chrom": chrom,
                    "pos": pos,
                    "snp": snp["SNPS"],
                    "trait": snp["trait"],
                    "p_value": snp["p_value"],
                })
    cols = ["peak_idx", "chrom", "pos", "snp", "trait", "p_value"]
    return pd.DataFrame(records, columns=cols)


def _per_cell_snp_burden(
    A: scipy.sparse.spmatrix,
    overlap_peak_idx: np.ndarray,
) -> np.ndarray:
    """Per-cell sum of accessibility over overlapping peaks. Cell axis stays sparse."""
    if len(overlap_peak_idx) == 0:
        return np.zeros(A.shape[0], dtype=np.float32)
    A_csc = A.tocsc()
    sub = A_csc[:, overlap_peak_idx]
    return np.asarray(sub.sum(axis=1)).ravel().astype(np.float32)


class GWASOverlayModule:
    """Overlay GWAS catalog SNPs onto ATAC peaks; emit per-cell SNP burden."""

    name = "gwas_overlay"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {
        "obsm": ["atac_peaks"],
        "uns": ["atac_var"],
    }
    provides_keys: dict[str, list[str]] = {
        "obs": ["gwas_snp_burden"],
        "uns": ["gwas_snp_peak_overlap"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")
        adata = ctx.adata
        cfg = ctx.cfg

        if "atac_peaks" not in adata.obsm:
            ctx.status(self.name, "skipped", "adata.obsm['atac_peaks'] absent")
            ctx.metadata["gwas_overlay_status"] = "skipped_no_peaks"
            return

        gwas_path = getattr(cfg, "gwas_catalog_path", None)
        if not gwas_path:
            ctx.status(self.name, "skipped",
                       "gwas_catalog_path not set — see data/external/gwas_catalog/MANIFEST.json")
            ctx.metadata["gwas_overlay_status"] = "skipped_no_catalog_path"
            return
        gwas_path = Path(gwas_path)
        if not gwas_path.exists():
            ctx.status(self.name, "skipped", f"GWAS catalog missing: {gwas_path}")
            ctx.metadata["gwas_overlay_status"] = "skipped_catalog_missing"
            return

        trait_filter = getattr(cfg, "gwas_trait_filter", None)
        p_threshold = float(getattr(cfg, "gwas_p_value_threshold", 5e-8))
        logger.info("gwas_overlay: loading catalog %s (trait_filter=%s, p<=%g)",
                    gwas_path, trait_filter, p_threshold)
        snps = _load_gwas_catalog(gwas_path, trait_filter, p_threshold)
        logger.info("gwas_overlay: %d SNPs after filter", len(snps))

        # Build peak BED from atac_var
        atac_var = adata.uns["atac_var"]
        if {"seqnames", "start", "end"}.issubset(atac_var.columns):
            peaks_df = atac_var[["seqnames", "start", "end"]].rename(
                columns={"seqnames": "chrom"}
            ).copy()
        elif "peak_id" in atac_var.columns:
            parts = atac_var["peak_id"].str.extract(r"^([\w.]+):(\d+)-(\d+)$")
            peaks_df = pd.DataFrame({
                "chrom": parts[0],
                "start": parts[1].astype(int),
                "end": parts[2].astype(int),
            })
        else:
            raise ModuleContractError(
                "gwas_overlay: atac_var must have (seqnames,start,end) or peak_id"
            )
        peaks_df = peaks_df.reset_index(drop=True)

        overlap = _overlap_snps_with_peaks(snps, peaks_df)
        logger.info("gwas_overlay: %d SNP-peak overlaps", len(overlap))

        out_dir = ctx.run_dir / "gwas_overlay"
        out_dir.mkdir(parents=True, exist_ok=True)
        overlap.to_csv(out_dir / "snp_peak_overlap.csv", index=False)

        # Per-cell burden
        unique_peak_idx = np.array(sorted(set(overlap["peak_idx"].tolist())), dtype=np.int64)
        burden = _per_cell_snp_burden(adata.obsm["atac_peaks"], unique_peak_idx)
        adata.obs["gwas_snp_burden"] = burden
        adata.uns["gwas_snp_peak_overlap"] = overlap

        per_cell_df = pd.DataFrame({
            "cell": adata.obs_names,
            "gwas_snp_burden": burden,
        })
        per_cell_df.to_csv(out_dir / "per_cell_burden.csv", index=False)

        summary = {
            "n_snps_loaded": int(len(snps)),
            "n_snp_peak_overlaps": int(len(overlap)),
            "n_unique_peaks_with_snp": int(len(unique_peak_idx)),
            "p_value_threshold": p_threshold,
            "trait_filter": trait_filter,
            "burden_mean": float(burden.mean()),
            "burden_max": float(burden.max()),
        }
        (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

        ctx.metadata["gwas_overlay_status"] = "ok"
        ctx.metadata["gwas_overlay_n_overlaps"] = int(len(overlap))
        logger.info("gwas_overlay: complete; %d cells scored", len(burden))
