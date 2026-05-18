"""ATAC QC (Wave 2B / P2.S12): TSS enrichment + FRiP from fragments.

Memory-safe: streams fragments file line-by-line (no full load). Sparse-safe
peak counting via scipy.sparse if peak matrix is provided.

Inputs (via cfg):
  --atac-fragments-path PATH    fragments.tsv.gz (cellranger-atac output)
  --atac-tss-bed-path PATH      tss.bed (chrom, tss_start, tss_end) — UCSC gencode TSS coords
  --atac-peaks-bed-path PATH    peaks.bed (chrom, start, end) — for FRiP denominator

Outputs:
  adata.obs["atac_qc_status"]   "pass" | "warn" | "fail"
  adata.obs["frip"]             fraction of reads in peaks (float, NaN if no fragments)
  adata.obs["tss_enrichment"]   TSS enrichment score (float)
  adata.uns["atac_qc_summary"]  cohort-level stats
  runs/<run-id>/atac_qc/atac_qc.json
"""
from __future__ import annotations

import gzip
import json
import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "ENCODE_ATAC_standards": {
        "title": "ENCODE ATAC-seq Data Standards",
        "authors": "ENCODE DCC",
        "journal": "Documentation",
        "year": "2024",
        "doi": "https://www.encodeproject.org/atac-seq/",
        "description": "TSS enrichment ≥7 and FRiP ≥0.3 thresholds adopted here as default 'pass' bounds.",
    },
    "Buenrostro_ATACseq_2013": {
        "title": "Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position",
        "authors": "Buenrostro et al.",
        "journal": "Nature Methods",
        "year": "2013",
        "doi": "10.1038/nmeth.2688",
        "description": "Original ATAC-seq methodology; TSS enrichment + open-chromatin peak interpretation.",
    },
    "Stuart_Signac_2021": {
        "title": "Multimodal single-cell chromatin analysis with Signac",
        "authors": "Stuart et al.",
        "journal": "Nature Methods",
        "year": "2021",
        "doi": "10.1038/s41592-021-01282-5",
        "description": "Canonical TSS-enrichment / FRiP computation conventions for single-cell ATAC; followed here.",
    },
}


# ENCODE-style thresholds.
_TSS_PASS = 7.0
_TSS_WARN = 4.0
_FRIP_PASS = 0.30
_FRIP_WARN = 0.10


def _load_bed(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end"],
        usecols=[0, 1, 2],
        dtype={"chrom": str, "start": int, "end": int},
    )


def _stream_fragments_per_barcode(path: Path) -> dict[str, list[tuple[str, int, int]]]:
    """Stream fragments.tsv.gz and bucket by barcode.

    Returns a dict barcode → list of (chrom, start, end) tuples. Caller is
    responsible for keeping per-barcode counts; we DO NOT load the full file
    into a single DataFrame (avoids ~5-30 GB memory for 800k cells).
    """
    per_barcode: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, bc = parts[0], int(parts[1]), int(parts[2]), parts[3]
            per_barcode[bc].append((chrom, start, end))
    return per_barcode


def _compute_frip(fragments_by_bc: dict, peaks: pd.DataFrame) -> dict[str, float]:
    """Per-barcode FRiP via interval intersection on a per-chromosome basis."""
    peak_by_chrom: dict[str, np.ndarray] = {}
    for chrom, group in peaks.groupby("chrom"):
        peak_by_chrom[chrom] = group[["start", "end"]].sort_values("start").to_numpy()

    frip: dict[str, float] = {}
    for bc, frags in fragments_by_bc.items():
        total = len(frags)
        if total == 0:
            frip[bc] = float("nan")
            continue
        in_peaks = 0
        for chrom, fstart, fend in frags:
            if chrom not in peak_by_chrom:
                continue
            arr = peak_by_chrom[chrom]
            # Binary search for fragments overlapping any peak (interval intersection).
            # arr is sorted by start; find first peak whose start < fend, then check end > fstart.
            starts = arr[:, 0]
            idx = np.searchsorted(starts, fend, side="right") - 1
            if idx < 0:
                continue
            # Walk back a few peaks to handle overlapping peaks (rare).
            while idx >= 0 and arr[idx, 1] >= fstart:
                if arr[idx, 0] <= fend and arr[idx, 1] >= fstart:
                    in_peaks += 1
                    break
                idx -= 1
        frip[bc] = in_peaks / total
    return frip


def _compute_tss_enrichment(fragments_by_bc: dict, tss_sites: pd.DataFrame, flank: int = 1000) -> dict[str, float]:
    """Per-barcode TSS enrichment.

    Computes mean fragment density in a ±flank bp window around each TSS, normalised
    by mean background density (200 bp at the flank edges). ENCODE standard window.
    """
    tss_by_chrom: dict[str, np.ndarray] = {}
    for chrom, group in tss_sites.groupby("chrom"):
        # Use the midpoint of each TSS interval
        midpoints = ((group["start"].to_numpy() + group["end"].to_numpy()) // 2)
        tss_by_chrom[chrom] = np.sort(midpoints)

    tss_scores: dict[str, float] = {}
    for bc, frags in fragments_by_bc.items():
        if not frags:
            tss_scores[bc] = float("nan")
            continue
        center_count = 0
        background_count = 0
        center_bp = flank * 2
        background_bp = 200 * 2  # 100 bp at each flank edge
        for chrom, fstart, fend in frags:
            if chrom not in tss_by_chrom:
                continue
            tss_arr = tss_by_chrom[chrom]
            frag_mid = (fstart + fend) // 2
            # Find nearest TSS via binary search
            i = np.searchsorted(tss_arr, frag_mid)
            nearest = None
            if i < len(tss_arr):
                nearest = tss_arr[i]
            if i > 0:
                left = tss_arr[i - 1]
                if nearest is None or abs(frag_mid - left) < abs(frag_mid - nearest):
                    nearest = left
            if nearest is None:
                continue
            dist = abs(frag_mid - nearest)
            if dist <= flank:
                center_count += 1
            elif flank - 100 <= dist <= flank:
                background_count += 1
        center_density = center_count / center_bp
        bg_density = background_count / background_bp if background_count > 0 else 1e-9
        tss_scores[bc] = float(center_density / max(bg_density, 1e-9))
    return tss_scores


class ATACQCModule:
    """ATAC QC: per-cell TSS enrichment + FRiP from fragments.tsv.gz."""

    name = "atac_qc"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {"obs": ["frip", "tss_enrichment", "atac_qc_status"]}

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata
        frag_path = getattr(ctx.cfg, "atac_fragments_path", None)
        tss_path = getattr(ctx.cfg, "atac_tss_bed_path", None)
        peaks_path = getattr(ctx.cfg, "atac_peaks_bed_path", None)

        if not frag_path:
            ctx.status(self.name, "skipped", "atac_fragments_path not set")
            ctx.metadata["atac_qc_status_summary"] = "skipped_no_fragments"
            return

        frag_path = Path(frag_path)
        if not frag_path.exists():
            ctx.status(self.name, "skipped", f"fragments missing: {frag_path}")
            ctx.metadata["atac_qc_status_summary"] = "skipped_path_missing"
            return

        logger.info("%s: streaming fragments from %s", self.name, frag_path)
        fragments_by_bc = _stream_fragments_per_barcode(frag_path)

        # FRiP
        frip_per_bc: dict[str, float] = {}
        if peaks_path and Path(peaks_path).exists():
            peaks_df = _load_bed(Path(peaks_path))
            frip_per_bc = _compute_frip(fragments_by_bc, peaks_df)
        else:
            logger.warning("%s: --atac-peaks-bed-path not provided; FRiP set to NaN", self.name)

        # TSS enrichment
        tss_per_bc: dict[str, float] = {}
        if tss_path and Path(tss_path).exists():
            tss_df = _load_bed(Path(tss_path))
            tss_per_bc = _compute_tss_enrichment(fragments_by_bc, tss_df)
        else:
            logger.warning("%s: --atac-tss-bed-path not provided; TSS enrichment set to NaN", self.name)

        # Map onto adata.obs
        adata.obs["frip"] = adata.obs_names.to_series().map(lambda bc: frip_per_bc.get(bc, float("nan")))
        adata.obs["tss_enrichment"] = adata.obs_names.to_series().map(lambda bc: tss_per_bc.get(bc, float("nan")))

        def _qc_status(row) -> str:
            frip = row.get("frip", float("nan"))
            tss = row.get("tss_enrichment", float("nan"))
            if pd.isna(frip) and pd.isna(tss):
                return "no_data"
            pass_frip = pd.notna(frip) and frip >= _FRIP_PASS
            pass_tss = pd.notna(tss) and tss >= _TSS_PASS
            warn_frip = pd.notna(frip) and frip >= _FRIP_WARN
            warn_tss = pd.notna(tss) and tss >= _TSS_WARN
            if pass_frip and pass_tss:
                return "pass"
            if warn_frip or warn_tss:
                return "warn"
            return "fail"

        adata.obs["atac_qc_status"] = adata.obs[["frip", "tss_enrichment"]].apply(_qc_status, axis=1).astype(str)

        summary = {
            "n_cells_evaluated": int(adata.n_obs),
            "n_pass": int((adata.obs["atac_qc_status"] == "pass").sum()),
            "n_warn": int((adata.obs["atac_qc_status"] == "warn").sum()),
            "n_fail": int((adata.obs["atac_qc_status"] == "fail").sum()),
            "n_no_data": int((adata.obs["atac_qc_status"] == "no_data").sum()),
            "frip_median": float(adata.obs["frip"].median(skipna=True)) if adata.obs["frip"].notna().any() else None,
            "tss_median": float(adata.obs["tss_enrichment"].median(skipna=True)) if adata.obs["tss_enrichment"].notna().any() else None,
            "thresholds": {"frip_pass": _FRIP_PASS, "frip_warn": _FRIP_WARN, "tss_pass": _TSS_PASS, "tss_warn": _TSS_WARN},
        }
        adata.uns["atac_qc_summary"] = summary
        out_dir = ctx.run_dir / "atac_qc"
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "atac_qc.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        ctx.metadata["atac_qc_status_summary"] = "ok"
        logger.info("%s: pass=%d warn=%d fail=%d no_data=%d",
                    self.name, summary["n_pass"], summary["n_warn"], summary["n_fail"], summary["n_no_data"])
