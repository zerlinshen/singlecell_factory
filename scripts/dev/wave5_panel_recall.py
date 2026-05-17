# STATUS: one-off, promote-on-reuse
#!/usr/bin/env python3
"""Wave-5 panel recall + top-50 hit rate computation.

Computes panel_recall and top50_hit_rate of a curated cortical panel against
the canonical peak-to-gene top-1000 CSV (v3). Pure tabular; no AnnData.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--top1000-csv", required=True, type=Path,
                   help="path to peak_to_gene_top1000_v3.csv")
    p.add_argument("--panel-json", required=True, type=Path,
                   help="path to wave5_cortical_panel_v1.json")
    p.add_argument("--ensg-symbol-tsv", required=True, type=Path,
                   help="path to gencode TSV with gene_id (ENSG) + gene_name (symbol)")
    p.add_argument("--out-json", required=True, type=Path,
                   help="output JSON path")
    p.add_argument("--top-n", type=int, default=1000,
                   help="N for top-N panel recall (default 1000)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.top1000_csv)
    required_cols = {"peak_idx", "gene_idx", "pearson", "p_perm", "fdr", "gene"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"top1000 CSV missing columns: {sorted(missing)}")

    df = df.assign(abs_pearson=df["pearson"].abs())
    df_sorted = df.sort_values("abs_pearson", ascending=False).reset_index(drop=True)
    df_topn = df_sorted.head(args.top_n).copy()
    df_top50 = df_sorted.head(50).copy()

    gencode = pd.read_csv(args.ensg_symbol_tsv, sep="\t")
    if "gene_id" not in gencode.columns or "gene_name" not in gencode.columns:
        raise ValueError(
            f"gencode TSV missing gene_id/gene_name columns: {list(gencode.columns)}"
        )
    ensg_to_symbol = dict(zip(gencode["gene_id"].astype(str),
                              gencode["gene_name"].astype(str)))

    def map_ensg(series: pd.Series) -> pd.Series:
        return series.astype(str).map(ensg_to_symbol)

    df_topn["gene_symbol"] = map_ensg(df_topn["gene"])
    df_top50["gene_symbol"] = map_ensg(df_top50["gene"])

    n_ensg_unmapped = int(df_topn["gene_symbol"].isna().sum())
    df_topn_mapped = df_topn.dropna(subset=["gene_symbol"]).copy()
    df_top50_mapped = df_top50.dropna(subset=["gene_symbol"]).copy()

    unique_topn_symbols = set(df_topn_mapped["gene_symbol"].unique().tolist())
    unique_top50_symbols = set(df_top50_mapped["gene_symbol"].unique().tolist())

    panel = json.loads(Path(args.panel_json).read_text())
    if not (isinstance(panel, dict) and "markers" in panel):
        raise ValueError(
            f"panel JSON shape unrecognized: expected {{markers: [...]}}, "
            f"got keys={list(panel.keys()) if isinstance(panel, dict) else type(panel)}"
        )
    panel_entries = panel["markers"]

    panel_pairs = [(str(e["gene_symbol"]), str(e.get("lineage", "")))
                   for e in panel_entries]
    panel_symbol_set = set(s for s, _ in panel_pairs)
    panel_total = len(panel_symbol_set)

    hits_symbols = panel_symbol_set & unique_topn_symbols
    top50_hits_symbols = panel_symbol_set & unique_top50_symbols
    panel_recall = len(hits_symbols) / panel_total if panel_total else 0.0
    top50_hit_rate = len(top50_hits_symbols) / 50.0

    panel_hits = sorted([(s, lin) for (s, lin) in panel_pairs if s in hits_symbols])
    panel_misses = sorted([(s, lin) for (s, lin) in panel_pairs if s not in hits_symbols])

    # top-10 unique gene symbols by max-|pearson| (gene-symbol-deduped).
    dedup = (df_topn_mapped
             .groupby("gene_symbol", as_index=False)["abs_pearson"].max()
             .sort_values("abs_pearson", ascending=False))
    top10_genes = dedup.head(10)["gene_symbol"].tolist()

    panel_sha = sha256_file(args.panel_json)

    out = {
        "panel_recall": panel_recall,
        "top50_hit_rate": top50_hit_rate,
        "panel_hits": panel_hits,
        "panel_misses": panel_misses,
        "top10_genes": top10_genes,
        "n_ensg_unmapped": n_ensg_unmapped,
        "n_top_n": int(len(df_topn)),
        "n_unique_symbols_top_n": int(len(unique_topn_symbols)),
        "panel_sha": panel_sha,
        "panel_total": panel_total,
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(out, indent=2))
    print(f"[wave5_panel_recall] wrote {args.out_json}")
    print(f"  panel_recall={panel_recall:.4f}  top50_hit_rate={top50_hit_rate:.4f}")
    print(f"  panel_hits={len(panel_hits)}/{panel_total}  n_ensg_unmapped={n_ensg_unmapped}")


if __name__ == "__main__":
    main()
