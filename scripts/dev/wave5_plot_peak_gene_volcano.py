# STATUS: one-off, promote-on-reuse
#!/usr/bin/env python3
"""Wave-5 peak-gene volcano-style plot.

Scatter of |pearson| vs -log10(fdr) over all significant linkages with
top-1000 highlighted orange and panel-gene hits highlighted red (star).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--all-significant-csv", required=True, type=Path,
                   help="path to peak_to_gene_significant_v3.csv (or fallback)")
    p.add_argument("--top1000-csv", required=True, type=Path,
                   help="path to peak_to_gene_top1000_v3.csv for highlight overlay")
    p.add_argument("--panel-json", required=True, type=Path,
                   help="path to wave5_cortical_panel_v1.json")
    p.add_argument("--ensg-symbol-tsv", required=True, type=Path,
                   help="path to gencode TSV (gene_id + gene_name)")
    p.add_argument("--out-png", required=True, type=Path,
                   help="output PNG path")
    return p.parse_args()


def load_panel_symbols(panel_json: Path) -> set[str]:
    panel = json.loads(panel_json.read_text())
    if not (isinstance(panel, dict) and "markers" in panel):
        raise ValueError(
            f"panel JSON shape unrecognized: expected {{markers: [...]}}, "
            f"got keys={list(panel.keys()) if isinstance(panel, dict) else type(panel)}"
        )
    return {str(e["gene_symbol"]) for e in panel["markers"]}


def main() -> None:
    args = parse_args()

    if args.all_significant_csv.exists():
        print(f"[wave5_plot_peak_gene_volcano] using all-significant CSV: {args.all_significant_csv}")
        df_all = pd.read_csv(args.all_significant_csv)
        source = "all_significant"
    else:
        print(f"[wave5_plot_peak_gene_volcano] all-significant CSV absent ({args.all_significant_csv}); "
              f"using top-1000 CSV exclusively")
        df_all = pd.read_csv(args.top1000_csv)
        source = "top1000_only"

    df_top = pd.read_csv(args.top1000_csv)

    gencode = pd.read_csv(args.ensg_symbol_tsv, sep="\t")
    ensg_to_symbol = dict(zip(gencode["gene_id"].astype(str),
                              gencode["gene_name"].astype(str)))
    df_top["gene_symbol"] = df_top["gene"].astype(str).map(ensg_to_symbol)

    panel_symbols = load_panel_symbols(args.panel_json)
    df_top_panel = df_top[df_top["gene_symbol"].isin(panel_symbols)].copy()

    fig, ax = plt.subplots(figsize=(9, 7))

    def neg_log10_fdr(s: pd.Series) -> np.ndarray:
        return -np.log10(s.astype(float).to_numpy() + 1e-300)

    x_all = df_all["pearson"].abs().to_numpy()
    y_all = neg_log10_fdr(df_all["fdr"])
    ax.scatter(x_all, y_all, s=8, c="lightgray", alpha=0.3,
               label=f"all ({source}, n={len(df_all)})")

    x_top = df_top["pearson"].abs().to_numpy()
    y_top = neg_log10_fdr(df_top["fdr"])
    ax.scatter(x_top, y_top, s=14, c="orange", alpha=0.7,
               label=f"top-1000 (n={len(df_top)})")

    if len(df_top_panel) > 0:
        x_panel = df_top_panel["pearson"].abs().to_numpy()
        y_panel = neg_log10_fdr(df_top_panel["fdr"])
        ax.scatter(x_panel, y_panel, s=120, c="red", marker="*",
                   edgecolors="black", linewidths=0.5,
                   label=f"panel hits (n={len(df_top_panel)})")
        for _, row in df_top_panel.iterrows():
            ax.annotate(
                str(row["gene_symbol"]),
                xy=(abs(float(row["pearson"])), float(-np.log10(float(row["fdr"]) + 1e-300))),
                xytext=(4, 4), textcoords="offset points",
                fontsize=8, color="darkred",
            )

    ax.set_xlabel("|pearson|")
    ax.set_ylabel("-log10(fdr + 1e-300)")
    ax.set_title("Wave-5 peak-gene linkages (volcano-style)")
    ax.legend(loc="best", fontsize=9)
    ax.grid(alpha=0.2)

    args.out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"[wave5_plot_peak_gene_volcano] wrote {args.out_png}")
    print(f"  source={source}  n_all={len(df_all)}  n_top1000={len(df_top)}  n_panel_hits={len(df_top_panel)}")


if __name__ == "__main__":
    main()
