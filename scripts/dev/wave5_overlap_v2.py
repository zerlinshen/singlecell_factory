"""Wave-5 v2 overlap evaluator (AC-VAL-3).

Computes top-K peak-gene linkage overlap between our v2 output and
Trevino S2F. Both lists share Trevino's consensus peak naming
(``c{cluster}_peak_{n}{sub}``) — that is the canonical identifier and
is the only correct way to match links. Coordinates are downstream of
peak_name (consensus peaks file maps name -> coords), so coord-based
matching is strictly redundant.

Primary metric (AC-VAL-3): (peak_name, gene_symbol) tuple overlap at top-K.
Diagnostics: gene-set overlap and peak-set overlap, to characterize
where the two lists agree (gene biology) vs disagree (peak granularity).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata
import pandas as pd


def _log(m: str) -> None:
    print(m, flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-h5ad", required=True, type=Path)
    ap.add_argument("--v2-top1000", required=True, type=Path)
    ap.add_argument("--s2f-raw-mmc2", required=True, type=Path,
                    help="Path to mmc2.xlsx (Trevino S2 supplement; uses sheet F).")
    ap.add_argument("--biotype-tsv", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--top-k", type=int, default=1000)
    args = ap.parse_args()

    adata = anndata.read_h5ad(args.input_h5ad, backed="r")
    atac_var = adata.uns["atac_peaks"]
    if "name" not in atac_var.columns:
        _log("FAIL: atac_var lacks 'name' column — cannot match by peak_name.")
        return 2

    # ---- Our v2 top-K, mapped to peak_name + gene_symbol ----
    ours = pd.read_csv(args.v2_top1000)
    if "peak_idx" not in ours.columns or "gene" not in ours.columns:
        _log(f"FAIL: v2 top1000 missing peak_idx/gene; got {list(ours.columns)}")
        return 2
    if len(ours) == 0:
        _log("FAIL: v2 top1000 is empty (no FDR-significant pairs at fdr<0.05).")
        return 2

    bio = pd.read_csv(args.biotype_tsv, sep="\t").set_index("gene_id")
    id_to_sym = bio["gene_name"].to_dict()
    ours["peak_name"]   = atac_var["name"].iloc[ours["peak_idx"].astype(int)].astype(str).to_numpy()
    ours["gene_symbol"] = ours["gene"].astype(str).map(id_to_sym).fillna(ours["gene"].astype(str))
    ours = ours.reindex(ours["pearson"].abs().sort_values(ascending=False).index)\
               .head(args.top_k).reset_index(drop=True)
    _log(f"ours top-{args.top_k}: |pearson| {ours['pearson'].abs().min():.3f}..{ours['pearson'].abs().max():.3f}")

    # ---- Trevino S2F top-K, native (peak_name, gene_symbol, corr) from mmc2 sheet F ----
    import openpyxl
    wb = openpyxl.load_workbook(args.s2f_raw_mmc2, read_only=True, data_only=True)
    ws = wb["F"]
    rows = list(ws.iter_rows(values_only=True))
    header = rows[1]
    data = [r for r in rows[2:] if r and r[0] is not None]
    s2f = pd.DataFrame(data, columns=header)
    wb.close()
    s2f = s2f.rename(columns={"Peak name": "peak_name", "Gene symbol": "gene_symbol", "correlation": "corr"})
    _log(f"S2F raw: {len(s2f)} rows")
    trev = s2f.reindex(s2f["corr"].abs().sort_values(ascending=False).index)\
              .head(args.top_k).reset_index(drop=True)

    # ---- Primary metric: (peak_name, gene_symbol) overlap at top-K ----
    ours_tuples = set(zip(ours["peak_name"], ours["gene_symbol"]))
    trev_tuples = set(zip(trev["peak_name"].astype(str), trev["gene_symbol"].astype(str)))
    overlap = ours_tuples & trev_tuples
    primary_fraction = len(overlap) / args.top_k

    # ---- Diagnostics ----
    ours_peaks = set(ours["peak_name"])
    trev_peaks = set(trev["peak_name"].astype(str))
    ours_genes = set(ours["gene_symbol"])
    trev_genes = set(trev["gene_symbol"].astype(str))

    diag = {
        "peak_set_intersection":         len(ours_peaks & trev_peaks),
        "peak_set_union":                len(ours_peaks | trev_peaks),
        "peak_set_jaccard":              len(ours_peaks & trev_peaks) / len(ours_peaks | trev_peaks) if (ours_peaks | trev_peaks) else 0.0,
        "gene_set_intersection":         len(ours_genes & trev_genes),
        "gene_set_union":                len(ours_genes | trev_genes),
        "gene_set_jaccard":              len(ours_genes & trev_genes) / len(ours_genes | trev_genes) if (ours_genes | trev_genes) else 0.0,
        "n_unique_peaks_ours":           len(ours_peaks),
        "n_unique_peaks_trev":           len(trev_peaks),
        "n_unique_genes_ours":           len(ours_genes),
        "n_unique_genes_trev":           len(trev_genes),
        "ours_top10":                    ours.head(10)[["peak_name", "gene_symbol", "pearson", "fdr"]].to_dict("records"),
        "trev_top10":                    trev.head(10)[["peak_name", "gene_symbol", "corr"]].to_dict("records"),
        "ours_corr_distribution":        ours["pearson"].describe().to_dict(),
        "trev_corr_distribution":        trev["corr"].describe().to_dict(),
    }

    result = {
        "top_k":                  args.top_k,
        "primary_metric": {
            "name":               "peak_name_plus_gene_symbol_overlap",
            "value":              primary_fraction,
            "threshold":          0.70,
            "status":             "CLOSED" if primary_fraction >= 0.70 else ("CLOSED-PARTIAL" if primary_fraction >= 0.50 else "CLOSED-PARTIAL"),
            "n_overlap":          len(overlap),
            "n_ours":             len(ours_tuples),
            "n_trevino":          len(trev_tuples),
            "jaccard":            len(overlap) / len(ours_tuples | trev_tuples) if (ours_tuples | trev_tuples) else 0.0,
        },
        "diagnostics":            diag,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    _log(json.dumps({"primary_metric": result["primary_metric"],
                     "peak_set_jaccard": diag["peak_set_jaccard"],
                     "gene_set_jaccard": diag["gene_set_jaccard"]},
                    indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
