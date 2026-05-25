"""chromVAR TF motif activity scoring (Wave 5.1 / MV3 §3.2 Hybrid Option C).

# STATUS: promote-on-reuse (lightweight Python deviation; pychromvar upgrade path)

Computes per-cell TF motif activity ("deviation Z-scores") on the ATAC peak
matrix. Two execution modes:

- ``mode="pychromvar"`` — full chromVAR via pychromvar 0.0.4+. Activated only
  when both ``ctx.cfg.chromvar_genome_fasta`` and
  ``ctx.cfg.chromvar_motif_pwm`` are set AND exist. Calls
  ``pychromvar.add_peak_seq -> add_gc_bias -> get_bg_peaks -> match_motif ->
  compute_deviations``. Mirrors R chromVAR (Schep et al. 2017) ±GC-bias matching.

- ``mode="lightweight"`` (default, FASTA-free) — Hybrid Option C lightweight
  Z-score deviation. Reads a curated TF→peak annotation table
  (``ctx.cfg.tf_peak_annotation_path``); per TF, computes per-cell deviation
  as ``(observed - expected) / sqrt(expected_var)`` where observed is the
  sparse dot product of per-cell peak accessibility with the TF indicator
  vector and expected is derived from the global per-peak open-probability
  background. Cell axis is NEVER densified (sparse-axis discipline mirrors
  peak_to_gene.py).

The chosen-mode is recorded in ``adata.uns['chromvar_metadata']`` along with
all configuration knobs and Scenario 3 (MV3 rpy2 brittleness) divergence
flagging per plan v5.1 §3.6.

Inputs (via cfg):
  --chromvar-genome-fasta PATH          GRCh38 primary assembly FASTA
                                          (pychromvar mode only; absent → lightweight)
  --chromvar-motif-pwm PATH             JASPAR2024 PWM bundle (.meme or .jaspar)
                                          (pychromvar mode only)
  --tf-peak-annotation-path PATH        TSV with columns: tf_name, peak_id
                                          (lightweight mode)
  --chromvar-min-peaks-per-tf INT       Skip TFs with fewer than N annotated peaks
                                          (default 5)

Outputs:
  adata.obsm['chromvar_scores']         (n_cells, n_tfs) np.float32 — deviation Z-scores
  adata.uns['chromvar_tfs']             list[str] of TF symbols (matches obsm cols)
  adata.uns['chromvar_metadata']        dict — mode, version, scenario3 flags
  <run-dir>/chromvar/chromvar_scores.csv  per-cell × per-TF table
  <run-dir>/chromvar/tf_top_per_cluster.json  top TF per cluster (when leiden/seurat_clusters present)
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
    "Schep_chromVAR_2017": {
        "title": "chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data",
        "authors": "Schep AN, Wu B, Buenrostro JD, Greenleaf WJ",
        "journal": "Nature Methods",
        "year": "2017",
        "doi": "10.1038/nmeth.4401",
        "description": "Canonical chromVAR deviation Z-score formulation; lightweight mode here implements the deviation computation without per-cell GC bias matching.",
    },
    "Castro_JASPAR2024_2024": {
        "title": "JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles",
        "authors": "Rauluseviciute I, Riudavets-Puig R, et al.",
        "journal": "Nucleic Acids Research",
        "year": "2024",
        "doi": "10.1093/nar/gkad1059",
        "description": "TF motif PWM source consumed by pychromvar mode and by curated tf_peak_annotation derivation.",
    },
    "Bravo_pychromvar_2023": {
        "title": "pychromvar: a Python implementation of chromVAR",
        "authors": "scverse community",
        "journal": "scverse software (GitHub)",
        "year": "2023",
        "doi": "https://github.com/scverse/pychromvar",
        "description": "Python port of chromVAR consumed in mode=pychromvar.",
    },
    "Trevino_2021": {
        "title": "Chromatin and gene-regulatory dynamics of the developing human neocortex at single-cell resolution",
        "authors": "Trevino AE, Müller F, Andersen J, et al.",
        "journal": "Cell",
        "year": "2021",
        "doi": "10.1016/j.cell.2021.07.039",
        "description": "Source paper for F-3 chromVAR motif heatmap reproduction (Wave-5.6 MV3).",
    },
}


def _sparse_observed_per_tf(
    A_csc: scipy.sparse.csc_matrix,
    tf_peak_indices: np.ndarray,
) -> np.ndarray:
    """Sum per-cell accessibility across peaks annotated for one TF.

    A_csc: (n_cells, n_peaks) CSC sparse. tf_peak_indices: 1-D int array of
    peak column indices for the TF. Cell axis stays sparse: column-subset →
    sum along axis=1 yields (n_cells,) dense vector (cell-side dense is the
    intended output dimensionality, not a densification of the matrix).
    """
    if len(tf_peak_indices) == 0:
        return np.zeros(A_csc.shape[0], dtype=np.float64)
    sub = A_csc[:, tf_peak_indices]
    # axis=1 sum on (n_cells, k) sparse → (n_cells, 1) → ravel
    return np.asarray(sub.sum(axis=1)).ravel().astype(np.float64)


def _lightweight_deviation_scores(
    A_csr: scipy.sparse.csr_matrix,
    tf_to_peak_idx: dict[str, np.ndarray],
    min_peaks_per_tf: int,
) -> tuple[np.ndarray, list[str], dict]:
    """Compute per-cell deviation Z-score per TF on sparse peak matrix.

    For each TF t with peak set P_t:
      observed_c = sum_{p in P_t} A[c, p]
      expected_c = total_accessibility_c * (sum_{p in P_t} background_p) / sum_p background_p
                 ≈ (n_TF_peaks_normalized) * per-cell total accessibility
      var_c     = total_accessibility_c * (1 - effective_TF_fraction) * effective_TF_fraction
      z_c       = (observed_c - expected_c) / sqrt(var_c + 1e-9)

    Returns (Z (n_cells, n_tfs), tf_names, summary_dict).
    """
    n_cells, n_peaks = A_csr.shape
    # global per-peak open probability proxy: column sums normalized
    A_csc = A_csr.tocsc()
    col_sums = np.asarray(A_csc.sum(axis=0)).ravel().astype(np.float64)
    total_peak_mass = float(col_sums.sum())
    if total_peak_mass <= 0:
        raise ModuleContractError(
            "chromvar: ATAC peak matrix is empty (sum=0); cannot compute background."
        )
    bg_p = col_sums / total_peak_mass  # (n_peaks,), sums to 1

    # per-cell total accessibility (cell-axis dense, but n_cells-sized only)
    row_sums = np.asarray(A_csr.sum(axis=1)).ravel().astype(np.float64)

    kept_tfs: list[str] = []
    z_columns: list[np.ndarray] = []
    skipped_low_peak: list[str] = []
    for tf, peak_idx in tf_to_peak_idx.items():
        peak_idx = np.asarray(peak_idx, dtype=np.int64)
        # Bound-check
        peak_idx = peak_idx[(peak_idx >= 0) & (peak_idx < n_peaks)]
        if len(peak_idx) < min_peaks_per_tf:
            skipped_low_peak.append(tf)
            continue
        observed = _sparse_observed_per_tf(A_csc, peak_idx)  # (n_cells,)
        f_tf = float(bg_p[peak_idx].sum())                  # effective TF fraction
        expected = row_sums * f_tf
        var = row_sums * f_tf * (1.0 - f_tf) + 1e-9
        z = (observed - expected) / np.sqrt(var)
        z_columns.append(z.astype(np.float32))
        kept_tfs.append(tf)

    if not kept_tfs:
        raise ModuleContractError(
            f"chromvar: no TF passed min_peaks_per_tf={min_peaks_per_tf}; "
            f"check tf_peak_annotation_path content. {len(skipped_low_peak)} TFs skipped."
        )

    Z = np.column_stack(z_columns).astype(np.float32)
    summary = {
        "n_tfs_kept": len(kept_tfs),
        "n_tfs_skipped_low_peak": len(skipped_low_peak),
        "min_peaks_per_tf": int(min_peaks_per_tf),
        "n_cells": int(n_cells),
        "n_peaks": int(n_peaks),
        "z_mean": float(Z.mean()),
        "z_std": float(Z.std()),
        "z_abs_max": float(np.abs(Z).max()),
    }
    return Z, kept_tfs, summary


def _load_tf_peak_annotation(
    path: Path,
    peak_id_to_idx: dict[str, int],
) -> dict[str, np.ndarray]:
    """Load a curated TF→peak annotation TSV.

    Expected columns: tf_name (str), peak_id (str matching atac_var.peak_id).
    Returns dict[tf_name] -> np.ndarray of peak column indices.
    """
    df = pd.read_csv(path, sep="\t", dtype={"tf_name": str, "peak_id": str})
    if "tf_name" not in df.columns or "peak_id" not in df.columns:
        raise ModuleContractError(
            f"chromvar: tf_peak_annotation must have columns 'tf_name','peak_id' "
            f"(got {list(df.columns)})"
        )
    df["peak_idx"] = df["peak_id"].map(peak_id_to_idx)
    df = df.dropna(subset=["peak_idx"])
    df["peak_idx"] = df["peak_idx"].astype(int)
    tf_map: dict[str, np.ndarray] = {}
    for tf, grp in df.groupby("tf_name"):
        tf_map[str(tf)] = grp["peak_idx"].to_numpy()
    return tf_map


def _peak_id_index(adata) -> dict[str, int]:
    """Build peak_id → column index from adata.uns['atac_var']."""
    atac_var = adata.uns.get("atac_var")
    if atac_var is None:
        raise ModuleContractError(
            "chromvar: adata.uns['atac_var'] missing — run atac_ingest first."
        )
    if "peak_id" in atac_var.columns:
        ids = atac_var["peak_id"].astype(str).tolist()
    else:
        ids = [f"{c}:{s}-{e}" for c, s, e in zip(
            atac_var["seqnames"].astype(str),
            atac_var["start"].astype(int),
            atac_var["end"].astype(int),
        )]
    return {pid: i for i, pid in enumerate(ids)}


def _maybe_pychromvar(adata, cfg) -> tuple[np.ndarray, list[str], dict] | None:
    """Attempt pychromvar full chromVAR run. Returns None if prerequisites missing."""
    genome_fa = getattr(cfg, "chromvar_genome_fasta", None)
    motif_pwm = getattr(cfg, "chromvar_motif_pwm", None)
    if not genome_fa or not motif_pwm:
        return None
    genome_fa = Path(genome_fa)
    motif_pwm = Path(motif_pwm)
    if not genome_fa.exists() or not motif_pwm.exists():
        logger.info(
            "chromvar: pychromvar prerequisites missing (genome=%s exists=%s, "
            "motifs=%s exists=%s) — falling back to lightweight mode.",
            genome_fa, genome_fa.exists(), motif_pwm, motif_pwm.exists(),
        )
        return None
    try:
        import pychromvar as pc
        from anndata import AnnData
    except ImportError as exc:
        logger.info("chromvar: pychromvar not importable (%s) — lightweight fallback.", exc)
        return None
    # pychromvar expects an AnnData with X = peaks; we shim from obsm['atac_peaks'].
    # 2026-05-22 audit fix: pychromvar 0.0.4 splits adata.var_names by a single
    # delimiter (default "-") and reads positions [0],[1],[2] as chrom/start/end.
    # The 10x convention "chr1:1000-1500" splits to only two tokens and breaks
    # the indexing. Normalize the index to "{chrom}-{start}-{end}" so the
    # default delimiter works, regardless of how upstream produced peak_id.
    atac_var = adata.uns["atac_var"].copy()
    if {"chrom", "start", "end"}.issubset(atac_var.columns):
        canonical_names = (
            atac_var["chrom"].astype(str)
            + "-"
            + atac_var["start"].astype(int).astype(str)
            + "-"
            + atac_var["end"].astype(int).astype(str)
        )
        var_for_pychromvar = atac_var.copy()
        var_for_pychromvar.index = canonical_names.values
    elif "peak_id" in atac_var.columns:
        var_for_pychromvar = atac_var.set_index("peak_id")
    else:
        var_for_pychromvar = atac_var
    atac_ad = AnnData(
        X=adata.obsm["atac_peaks"],
        obs=adata.obs[[]].copy(),
        var=var_for_pychromvar,
    )
    pc.add_peak_seq(atac_ad, genome_file=str(genome_fa))
    pc.add_gc_bias(atac_ad)
    pc.get_bg_peaks(atac_ad)
    # Parse motifs (delegating to user-supplied loader; pychromvar accepts a list).
    motifs = _load_meme_motifs(motif_pwm)
    pc.match_motif(atac_ad, motifs=motifs)
    dev_ad = pc.compute_deviations(atac_ad)
    Z = np.asarray(dev_ad.X, dtype=np.float32)
    tfs = list(dev_ad.var_names)
    summary = {
        "n_tfs_kept": len(tfs),
        "pychromvar_version": pc.__version__,
        "genome_fasta": str(genome_fa),
        "motif_pwm": str(motif_pwm),
        "n_cells": int(Z.shape[0]),
        "z_mean": float(Z.mean()),
        "z_std": float(Z.std()),
    }
    return Z, tfs, summary


def _load_meme_motifs(path: Path):
    """Parse motifs from MEME-format file (JASPAR2024 text or MEME XML).

    Biopython's ``Bio.motifs.parse(handle, "MEME")`` expects MEME-XML output
    (from the ``meme`` discovery tool). JASPAR distributes motifs in
    MEME-text format which Biopython exposes as ``"minimal"``. We probe the
    file header to choose the right parser; downstream pychromvar tolerates
    a list of ``Bio.motifs.Motif`` from either source.

    pychromvar 0.0.4's ``match_motif`` requires each motif object to expose
    ``matrix_id`` and ``name`` attributes in JASPAR convention. Biopython's
    minimal-MEME parser only sets ``.name`` to the raw ``MOTIF`` line text,
    so we attach a JASPAR-style ``matrix_id`` (first whitespace token of the
    MOTIF line) and split-out short ``name`` (remaining text) before
    returning. This adapter is invisible to callers and only normalizes
    metadata, never PWM weights.
    """
    from Bio import motifs
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    if text.lstrip().lower().startswith("<?xml") or "<MEME " in text[:512]:
        fmt = "MEME"
    else:
        # 2026-05-22 audit fix: JASPAR2024_CORE_non-redundant_pfms_meme.txt is
        # MEME minimal text, not MEME XML; biopython names that parser
        # "minimal".
        fmt = "minimal"
    with open(path) as fh:
        parsed = list(motifs.parse(fh, fmt))
    # Normalize matrix_id / name for pychromvar's JASPAR-style expectation.
    for idx, m in enumerate(parsed):
        raw_name = getattr(m, "name", None) or f"motif_{idx}"
        tokens = str(raw_name).split(None, 1)
        matrix_id = tokens[0]
        short_name = tokens[1] if len(tokens) > 1 else matrix_id
        try:
            m.matrix_id = matrix_id  # type: ignore[attr-defined]
            m.name = short_name  # type: ignore[attr-defined]
        except Exception:
            # Some Motif implementations use __slots__; fall back to setattr-via-dict.
            try:
                m.__dict__["matrix_id"] = matrix_id
                m.__dict__["name"] = short_name
            except Exception:
                pass
    return parsed


class ChromVARModule:
    """Compute TF motif accessibility deviations (chromVAR-style).

    Default mode is the lightweight FASTA-free Z-score deviation (Hybrid C);
    upgrades transparently to pychromvar when genome + motifs are configured.
    """

    name = "chromvar"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {
        "obsm": ["atac_peaks"],
        "uns": ["atac_peaks", "atac_var"],
    }
    provides_keys: dict[str, list[str]] = {
        "obsm": ["chromvar_scores"],
        "uns": ["chromvar_tfs", "chromvar_metadata"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")
        adata = ctx.adata
        cfg = ctx.cfg

        if "atac_peaks" not in adata.obsm:
            ctx.status(self.name, "skipped", "adata.obsm['atac_peaks'] absent — run atac_ingest first")
            ctx.metadata["chromvar_status"] = "skipped_no_peaks"
            return

        if not scipy.sparse.issparse(adata.obsm["atac_peaks"]):
            raise ModuleContractError(
                "chromvar sparse-axis violation: adata.obsm['atac_peaks'] must be sparse."
            )

        A_csr = adata.obsm["atac_peaks"].tocsr()

        # Try full pychromvar path first
        full = _maybe_pychromvar(adata, cfg)
        mode = "pychromvar" if full is not None else "lightweight"

        if mode == "pychromvar":
            Z, tfs, summary = full
        else:
            ann_path = getattr(cfg, "tf_peak_annotation_path", None)
            if not ann_path:
                ctx.status(
                    self.name, "skipped",
                    "no tf_peak_annotation_path AND no pychromvar prerequisites — set one",
                )
                ctx.metadata["chromvar_status"] = "skipped_no_annotation"
                return
            ann_path = Path(ann_path)
            if not ann_path.exists():
                ctx.status(self.name, "skipped", f"tf_peak_annotation missing: {ann_path}")
                ctx.metadata["chromvar_status"] = "skipped_annotation_missing"
                return
            peak_id_to_idx = _peak_id_index(adata)
            tf_map = _load_tf_peak_annotation(ann_path, peak_id_to_idx)
            min_peaks = int(getattr(cfg, "chromvar_min_peaks_per_tf", 5))
            Z, tfs, summary = _lightweight_deviation_scores(A_csr, tf_map, min_peaks)

        adata.obsm["chromvar_scores"] = Z
        adata.uns["chromvar_tfs"] = list(tfs)
        adata.uns["chromvar_metadata"] = {
            "mode": mode,
            "summary": summary,
            "scenario3_divergence_documented": (mode == "lightweight"),
            "scenario3_rollback_path": (
                "Wave-6: set chromvar_genome_fasta + chromvar_motif_pwm to "
                "auto-switch to mode=pychromvar; re-run module."
                if mode == "lightweight" else None
            ),
        }

        # Persist tables
        out_dir = ctx.run_dir / "chromvar"
        out_dir.mkdir(parents=True, exist_ok=True)
        scores_df = pd.DataFrame(Z, index=adata.obs_names, columns=tfs)
        scores_df.to_csv(out_dir / "chromvar_scores.csv")

        # per-cluster aggregation when a clustering exists
        cluster_col = None
        for cand in ("seurat_clusters", "leiden"):
            if cand in adata.obs.columns:
                cluster_col = cand
                break
        if cluster_col is not None:
            scores_df["_cluster"] = adata.obs[cluster_col].astype(str).values
            cluster_means = (
                scores_df.groupby("_cluster", observed=True)
                .mean(numeric_only=True)
            )
            cluster_means.to_csv(out_dir / "chromvar_per_cluster.csv")
            top_per_cluster = {
                str(c): cluster_means.loc[c].nlargest(5).to_dict()
                for c in cluster_means.index
            }
            (out_dir / "tf_top_per_cluster.json").write_text(
                json.dumps(top_per_cluster, indent=2, default=float), encoding="utf-8",
            )

        (out_dir / "chromvar_summary.json").write_text(
            json.dumps({"mode": mode, **summary}, indent=2, default=float),
            encoding="utf-8",
        )

        # H-1 audit fix (2026-05-22): give lightweight vs strict modes
        # distinct top-level status strings so downstream consumers reading
        # the manifest cannot mistake the FASTA-free deviation lane (no
        # GC-matched background, no permutation null) for the canonical
        # Schep 2017 / pychromvar workflow.
        ctx.metadata["chromvar_status"] = (
            "ok_pychromvar" if mode == "pychromvar" else "ok_lightweight"
        )
        ctx.metadata["chromvar_mode"] = mode
        ctx.metadata["chromvar_method_actually_used"] = (
            "pychromvar_full_gc_matched_background"
            if mode == "pychromvar"
            else "lightweight_fasta_free_deviation_proxy"
        )
        ctx.metadata["chromvar_n_tfs"] = len(tfs)
        logger.info("chromvar: mode=%s, n_tfs=%d, n_cells=%d",
                    mode, len(tfs), Z.shape[0])
