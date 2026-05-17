"""Wave-5 peak-to-gene linkage v2 — vectorized permutation Pearson.

Replaces the sparse per-permutation row-reindex loop in
``workflow/modular/modules/peak_to_gene.py::compute_peak_gene_linkages``
with a dense BLAS matmul on z-scored matrices. Adds a protein-coding
biotype filter on candidate genes.

Why a separate script vs patching the factory module:
  - Factory module's permutation path is correct but slow (50+ min on
    467k peaks); patching it risks breaking the CI smoke that exists.
  - This script is the "Wave-5 v2 round" deliverable — once validated
    against Trevino S2F at overlap ≥ 0.70, the algorithmic core can be
    upstreamed into the factory module with a separate CI gate.

Wave-5 spec: .omc/specs/deep-interview-wave5-completion.md
Wave-5 plan: .omc/plans/wave5-completion-consensus-2026-05-16.md
Round v4 (this script): targets AC-VAL-3 overlap ≥ 0.70 via:
  1) protein-coding-only candidate filter (~halves pool)
  2) permutation FDR < 0.05 statistical gate (matches Trevino Methods)
  3) dense vectorized permutation Pearson (~10x faster than sparse path)
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scipy.sparse


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR within a candidate set."""
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranks = np.empty(n)
    ranks[order] = np.arange(1, n + 1)
    q = np.minimum.accumulate((p * n / ranks)[order][::-1])[::-1]
    out = np.empty(n)
    out[order] = q
    return np.minimum(out, 1.0)


def _log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {msg}", flush=True)


def compute_linkages_v2(
    adata: anndata.AnnData,
    biotype_tsv: Path,
    n_perms: int = 100,
    window_bp: int = 500_000,
    fdr_alpha: float = 0.05,
    protein_coding_only: bool = True,
    random_state: int = 42,
    progress_every: int = 25_000,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Run the optimized permutation pass.

    Returns
    -------
    top1000 : DataFrame  (peak_idx, gene, pearson, p_perm, fdr) — top 1000 by |pearson| among FDR-pass pairs
    significant : DataFrame  — all (peak, gene) pairs passing FDR < fdr_alpha
    summary : dict
    """
    n_cells = adata.n_obs
    n_peaks = adata.obsm["atac_peaks"].shape[1]

    # ---- 1. Load biotype + filter candidate gene indices ----
    biotypes = pd.read_csv(biotype_tsv, sep="\t")
    biotypes = biotypes.set_index("gene_id")
    var_idx = pd.Index(adata.var.index.astype(str))
    if protein_coding_only:
        coding_set = set(biotypes[biotypes["gene_type"] == "protein_coding"].index)
        gene_mask = var_idx.isin(coding_set)
    else:
        gene_mask = np.ones(adata.n_vars, dtype=bool)
    coding_gene_idx = np.where(gene_mask)[0]
    _log(f"  candidate gene pool: {len(coding_gene_idx)}/{adata.n_vars} "
         f"({'protein_coding only' if protein_coding_only else 'all genes'})")

    # ---- 2. Build windowed (peak, gene) candidate pairs ----
    gene_chrom = adata.var["chrom"].astype(str).to_numpy()
    gene_tss = adata.var["tss"].astype(np.int64).to_numpy()
    # Restrict gene lookup to coding genes
    coding_set_arr = np.zeros(adata.n_vars, dtype=bool)
    coding_set_arr[coding_gene_idx] = True

    gene_by_chrom: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom in np.unique(gene_chrom):
        mask = (gene_chrom == chrom) & coding_set_arr
        idxs = np.where(mask)[0]
        if len(idxs) == 0:
            continue
        order = np.argsort(gene_tss[idxs])
        gene_by_chrom[chrom] = (idxs[order], gene_tss[idxs][order])

    pmeta = adata.uns["atac_peaks"]
    chrom_col = "chrom" if "chrom" in pmeta.columns else "seqnames"
    pchrom = pmeta[chrom_col].astype(str).to_numpy()
    pmid = ((pmeta["start"].to_numpy(np.int64) + pmeta["end"].to_numpy(np.int64)) // 2)

    candidates: list[np.ndarray] = []  # per-peak candidate gene indices
    t0 = time.time()
    n_peaks_with_cand = 0
    for i in range(n_peaks):
        c = pchrom[i]
        if c not in gene_by_chrom:
            candidates.append(np.empty(0, dtype=np.int64))
            continue
        gidxs, gtss = gene_by_chrom[c]
        m = pmid[i]
        lo = np.searchsorted(gtss, m - window_bp)
        hi = np.searchsorted(gtss, m + window_bp)
        cand_i = gidxs[lo:hi]
        candidates.append(cand_i)
        if len(cand_i) > 0:
            n_peaks_with_cand += 1
    total_pairs = sum(len(c) for c in candidates)
    _log(f"  windowed pairing in {time.time()-t0:.1f}s; {total_pairs} candidate pairs "
         f"across {n_peaks_with_cand} peaks with ≥1 coding-gene candidate")

    # ---- 3. Pre-compute z-scored gene matrix (DENSE, cells × genes) ----
    t1 = time.time()
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.tocsc().astype(np.float32)
        X_dense = X.toarray()
    else:
        X_dense = np.asarray(X, dtype=np.float32)
    g_mean = X_dense.mean(axis=0)
    g_std = X_dense.std(axis=0)
    g_std[g_std == 0] = 1.0
    Z_gene = (X_dense - g_mean) / g_std  # (n_cells, n_genes)
    del X_dense
    _log(f"  z-scored gene matrix shape={Z_gene.shape} dtype={Z_gene.dtype} "
         f"({Z_gene.nbytes/1e9:.2f} GB) built in {time.time()-t1:.1f}s")

    # ---- 4. ATAC peak matrix as CSC for fast column access ----
    A_csc = adata.obsm["atac_peaks"].tocsc()
    n_cells_int = int(n_cells)

    # ---- 5. Permutation Pearson loop ----
    rng = np.random.default_rng(random_state)
    # Pre-generate one permutation index matrix; will reuse across peaks (each peak still gets
    # a fresh statistic because the OBSERVED Pearson differs; permuted r values are identically
    # distributed under H0 across peaks, so sharing the perm matrix is statistically valid).
    perm_idx = np.stack([rng.permutation(n_cells_int) for _ in range(n_perms)]).astype(np.int32)
    _log(f"  permutation index matrix shape={perm_idx.shape}")

    records: list[dict] = []
    t2 = time.time()
    n_skipped = 0
    for peak_idx in range(n_peaks):
        cand = candidates[peak_idx]
        if len(cand) == 0:
            continue
        # Densify + center the peak column
        a_col = A_csc[:, peak_idx].toarray().ravel().astype(np.float32)
        a_mean = a_col.mean()
        a_std = a_col.std()
        if a_std == 0:
            n_skipped += 1
            continue
        z_a = (a_col - a_mean) / a_std  # (n_cells,)

        Z_g = Z_gene[:, cand]  # (n_cells, n_cand)

        # Observed Pearson
        r_obs = z_a @ Z_g / n_cells_int  # (n_cand,)

        # Permuted Pearsons via fancy index then matmul
        z_a_perm = z_a[perm_idx]  # (n_perms, n_cells)
        r_perm = z_a_perm @ Z_g / n_cells_int  # (n_perms, n_cand)

        exceed = (np.abs(r_perm) >= np.abs(r_obs)).sum(axis=0)
        p_perm = (exceed + 1.0) / (n_perms + 1.0)
        fdr = _bh_fdr(p_perm)

        for local_i, gene_i in enumerate(cand):
            records.append({
                "peak_idx": int(peak_idx),
                "gene_idx": int(gene_i),
                "pearson": float(r_obs[local_i]),
                "p_perm": float(p_perm[local_i]),
                "fdr": float(fdr[local_i]),
            })

        if (peak_idx + 1) % progress_every == 0:
            elapsed = time.time() - t2
            est_total = elapsed * n_peaks / (peak_idx + 1)
            eta = est_total - elapsed
            _log(f"    peak {peak_idx+1}/{n_peaks}  elapsed={elapsed/60:.1f}m  ETA≈{eta/60:.1f}m  records={len(records)}")

    _log(f"  permutation pass done in {(time.time()-t2)/60:.1f} min  "
         f"(skipped {n_skipped} zero-std peaks, {len(records)} pairs)")

    df = pd.DataFrame.from_records(records)
    var_index = list(adata.var.index)
    df["gene"] = [var_index[g] for g in df["gene_idx"]]

    significant = df[df["fdr"] < fdr_alpha].copy()
    # Top-1000 by |pearson| ONLY among significant pairs (Trevino-style: stat gate first, then rank)
    abs_r = significant["pearson"].abs()
    top1000 = significant.iloc[abs_r.argsort()[::-1].values[:1000]].copy()

    summary = {
        "n_peaks_total":         int(n_peaks),
        "n_peaks_with_cand":     int(n_peaks_with_cand),
        "n_cand_pairs":          int(total_pairs),
        "n_perms":               int(n_perms),
        "window_bp":             int(window_bp),
        "fdr_alpha":             float(fdr_alpha),
        "protein_coding_only":   bool(protein_coding_only),
        "n_significant_pairs":   int(len(significant)),
        "n_top1000":             int(len(top1000)),
        "candidate_gene_pool":   int(len(coding_gene_idx)),
    }
    return top1000, significant, summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-h5ad", required=True, type=Path)
    ap.add_argument("--biotype-tsv", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--n-perms", type=int, default=100)
    ap.add_argument("--fdr-alpha", type=float, default=0.05)
    ap.add_argument("--no-coding-filter", action="store_true",
                    help="Disable the protein_coding-only candidate filter (debug).")
    ap.add_argument("--random-state", type=int, default=42)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    _log(f"loading {args.input_h5ad}")
    adata = anndata.read_h5ad(args.input_h5ad)
    _log(f"  n_obs={adata.n_obs} n_vars={adata.n_vars} atac_peaks={adata.obsm['atac_peaks'].shape}")

    top1000, sig, summary = compute_linkages_v2(
        adata,
        biotype_tsv=args.biotype_tsv,
        n_perms=args.n_perms,
        fdr_alpha=args.fdr_alpha,
        protein_coding_only=not args.no_coding_filter,
        random_state=args.random_state,
    )
    _log(f"summary: {summary}")

    top1000.to_csv(args.out_dir / "peak_to_gene_top1000_v2.csv", index=False)
    sig.to_csv(args.out_dir / "peak_to_gene_significant_v2.csv", index=False)
    (args.out_dir / "peak_to_gene_v2_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    _log(f"wrote {args.out_dir}/peak_to_gene_top1000_v2.csv ({len(top1000)} rows)")
    _log(f"wrote {args.out_dir}/peak_to_gene_significant_v2.csv ({len(sig)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
