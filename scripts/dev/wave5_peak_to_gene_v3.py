"""Wave-5 peak-to-gene v3 — adds kNN metacell aggregation before Pearson.

Diagnostic from v2: our top-1000 |Pearson| capped at 0.51, while Trevino
S2F's top-1000 reaches |corr|=1.0 (median 0.84). Single-cell Pearson on
~9k cells cannot produce |r|≈1.0 — the only mechanism is kNN smoothing
("metacell aggregation"), which is ArchR's `addPeak2GeneLinks` default
methodology and what Trevino 2021 used per their Methods.

v3 algorithm:
  1) Build kNN connectivity matrix K (n_cells, n_cells), k=100,
     normalized so each row sums to 1 (mean aggregation).
  2) Smooth gene matrix:  Z_gene_smooth = K @ Z_gene   (sparse @ dense)
  3) For each peak: smooth peak column: a_smooth = K @ A[:, p]
  4) Compute Pearson and permutation null on smoothed z-vectors.

Protein-coding filter + permutation FDR carried over from v2.

Wave-5 plan AC-VAL-3: top-1000 overlap ≥ 0.70 vs Trevino S2F.
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
from scipy.sparse import csr_matrix


def _bh_fdr(p: np.ndarray) -> np.ndarray:
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


def _log(m: str) -> None:
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {m}", flush=True)


def _build_knn_aggregator(
    adata: anndata.AnnData,
    use_rep: str,
    k: int,
) -> csr_matrix:
    """Build a row-normalized kNN connectivity matrix.

    Each row i has k nonzeros at column j iff j is among the k nearest
    neighbors of i in ``adata.obsm[use_rep]``. Row sums = 1 → smoothing
    by mean aggregation.
    """
    from sklearn.neighbors import NearestNeighbors
    X = np.asarray(adata.obsm[use_rep])
    n = X.shape[0]
    nn = NearestNeighbors(n_neighbors=k, algorithm="auto", n_jobs=-1).fit(X)
    _, idx = nn.kneighbors(X, return_distance=True)
    row = np.repeat(np.arange(n), k)
    col = idx.ravel()
    data = np.full(row.shape, 1.0 / k, dtype=np.float32)
    K = csr_matrix((data, (row, col)), shape=(n, n), dtype=np.float32)
    return K


def compute_linkages_v3(
    adata: anndata.AnnData,
    biotype_tsv: Path,
    n_perms: int = 100,
    window_bp: int = 500_000,
    fdr_alpha: float = 0.05,
    k_neighbors: int = 100,
    smoothing_rep: str = "X_wnn",
    random_state: int = 42,
    progress_every: int = 25_000,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    n_cells = adata.n_obs
    n_peaks = adata.obsm["atac_peaks"].shape[1]

    # ---- 1. Protein-coding filter + windowed candidate pairs ----
    biotypes = pd.read_csv(biotype_tsv, sep="\t").set_index("gene_id")
    coding_set = set(biotypes[biotypes["gene_type"] == "protein_coding"].index)
    var_idx = pd.Index(adata.var.index.astype(str))
    gene_mask = var_idx.isin(coding_set)
    coding_gene_idx = np.where(gene_mask)[0]
    _log(f"  candidate gene pool: {len(coding_gene_idx)}/{adata.n_vars} (protein_coding)")

    gene_chrom = adata.var["chrom"].astype(str).to_numpy()
    gene_tss = adata.var["tss"].astype(np.int64).to_numpy()
    coding_mask_arr = np.zeros(adata.n_vars, dtype=bool); coding_mask_arr[coding_gene_idx] = True
    gene_by_chrom = {}
    for chrom in np.unique(gene_chrom):
        mask = (gene_chrom == chrom) & coding_mask_arr
        idxs = np.where(mask)[0]
        if len(idxs) == 0:
            continue
        order = np.argsort(gene_tss[idxs])
        gene_by_chrom[chrom] = (idxs[order], gene_tss[idxs][order])

    pmeta = adata.uns["atac_peaks"]
    chrom_col = "chrom" if "chrom" in pmeta.columns else "seqnames"
    pchrom = pmeta[chrom_col].astype(str).to_numpy()
    pmid = ((pmeta["start"].to_numpy(np.int64) + pmeta["end"].to_numpy(np.int64)) // 2)

    candidates = []
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
    _log(f"  windowed pairing in {time.time()-t0:.1f}s; {total_pairs} pairs across {n_peaks_with_cand} peaks")

    # ---- 2. Build kNN aggregator on the multi-modal embedding ----
    t1 = time.time()
    if smoothing_rep not in adata.obsm:
        raise RuntimeError(f"adata.obsm['{smoothing_rep}'] missing; cannot build kNN aggregator")
    K = _build_knn_aggregator(adata, smoothing_rep, k_neighbors)
    _log(f"  kNN aggregator K shape={K.shape} nnz={K.nnz} (k={k_neighbors}, rep={smoothing_rep!r}) built in {time.time()-t1:.1f}s")

    # ---- 3. Smooth + z-score gene matrix ----
    t2 = time.time()
    X = adata.X.tocsr() if scipy.sparse.issparse(adata.X) else scipy.sparse.csr_matrix(adata.X)
    X_smooth = (K @ X).toarray().astype(np.float32)  # (n_cells, n_genes)
    g_mean = X_smooth.mean(axis=0)
    g_std = X_smooth.std(axis=0)
    g_std[g_std == 0] = 1.0
    Z_gene = (X_smooth - g_mean) / g_std
    del X_smooth
    _log(f"  smoothed+z-scored gene matrix shape={Z_gene.shape} dtype={Z_gene.dtype} "
         f"({Z_gene.nbytes/1e9:.2f} GB) built in {time.time()-t2:.1f}s")

    # ---- 4. Pre-smooth ATAC peak columns lazily via K @ A on the fly ----
    A_csc = adata.obsm["atac_peaks"].tocsc().astype(np.float32)
    # K @ A_csc is dense — costs n_cells * n_peaks * 4 bytes = 8981 * 467315 * 4 ~ 16.8 GB; too big.
    # Instead, smooth peak-by-peak: a_smooth = K @ A[:, p].toarray()  →  (n_cells,) float32.

    rng = np.random.default_rng(random_state)
    n_cells_int = int(n_cells)
    perm_idx = np.stack([rng.permutation(n_cells_int) for _ in range(n_perms)]).astype(np.int32)
    _log(f"  permutation index matrix shape={perm_idx.shape}")

    records: list[dict] = []
    t3 = time.time()
    n_skipped = 0
    K_dense_row_sum_inv = None  # not needed; K is row-normalized
    for peak_idx in range(n_peaks):
        cand = candidates[peak_idx]
        if len(cand) == 0:
            continue
        a_col = A_csc[:, peak_idx].toarray().ravel()
        # Smooth: K @ a_col (treat a_col as dense vector, output dense (n_cells,))
        a_smooth = (K @ a_col).astype(np.float32)
        a_mean = a_smooth.mean()
        a_std = a_smooth.std()
        if a_std == 0:
            n_skipped += 1
            continue
        z_a = (a_smooth - a_mean) / a_std

        Z_g = Z_gene[:, cand]

        r_obs = z_a @ Z_g / n_cells_int
        z_a_perm = z_a[perm_idx]
        r_perm = z_a_perm @ Z_g / n_cells_int

        exceed = (np.abs(r_perm) >= np.abs(r_obs)).sum(axis=0)
        p_perm = (exceed + 1.0) / (n_perms + 1.0)
        fdr = _bh_fdr(p_perm)

        for li, gi in enumerate(cand):
            records.append({
                "peak_idx": int(peak_idx),
                "gene_idx": int(gi),
                "pearson": float(r_obs[li]),
                "p_perm": float(p_perm[li]),
                "fdr": float(fdr[li]),
            })

        if (peak_idx + 1) % progress_every == 0:
            elapsed = time.time() - t3
            est_total = elapsed * n_peaks / (peak_idx + 1)
            _log(f"    peak {peak_idx+1}/{n_peaks}  elapsed={elapsed/60:.1f}m  ETA≈{(est_total-elapsed)/60:.1f}m  records={len(records)}")

    _log(f"  permutation pass done in {(time.time()-t3)/60:.1f} min (skipped {n_skipped}, {len(records)} pairs)")

    df = pd.DataFrame.from_records(records)
    var_index = list(adata.var.index)
    df["gene"] = [var_index[g] for g in df["gene_idx"]]

    significant = df[df["fdr"] < fdr_alpha].copy()
    abs_r = significant["pearson"].abs()
    top1000 = significant.iloc[abs_r.argsort()[::-1].values[:1000]].copy()

    summary = {
        "n_peaks_total": int(n_peaks),
        "n_peaks_with_cand": int(n_peaks_with_cand),
        "n_cand_pairs": int(total_pairs),
        "n_perms": int(n_perms),
        "window_bp": int(window_bp),
        "fdr_alpha": float(fdr_alpha),
        "k_neighbors": int(k_neighbors),
        "smoothing_rep": str(smoothing_rep),
        "n_significant_pairs": int(len(significant)),
        "n_top1000": int(len(top1000)),
        "candidate_gene_pool": int(len(coding_gene_idx)),
    }
    return top1000, significant, summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-h5ad", required=True, type=Path)
    ap.add_argument("--biotype-tsv", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--n-perms", type=int, default=100)
    ap.add_argument("--k-neighbors", type=int, default=100)
    ap.add_argument("--smoothing-rep", default="X_wnn")
    ap.add_argument("--fdr-alpha", type=float, default=0.05)
    ap.add_argument("--random-state", type=int, default=42)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    _log(f"loading {args.input_h5ad}")
    adata = anndata.read_h5ad(args.input_h5ad)
    _log(f"  n_obs={adata.n_obs} n_vars={adata.n_vars} atac_peaks={adata.obsm['atac_peaks'].shape}")

    top1000, sig, summary = compute_linkages_v3(
        adata,
        biotype_tsv=args.biotype_tsv,
        n_perms=args.n_perms,
        fdr_alpha=args.fdr_alpha,
        k_neighbors=args.k_neighbors,
        smoothing_rep=args.smoothing_rep,
        random_state=args.random_state,
    )
    _log(f"summary: {summary}")
    top1000.to_csv(args.out_dir / "peak_to_gene_top1000_v3.csv", index=False)
    sig.to_csv(args.out_dir / "peak_to_gene_significant_v3.csv", index=False)
    (args.out_dir / "peak_to_gene_v3_summary.json").write_text(json.dumps(summary, indent=2))
    _log(f"wrote top1000 ({len(top1000)}) and significant ({len(sig)}) to {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
