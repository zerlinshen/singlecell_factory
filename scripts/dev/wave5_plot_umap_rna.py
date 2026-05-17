# STATUS: one-off, promote-on-reuse
#!/usr/bin/env python3
"""Wave-5 UMAP plot for RNA modality (Leiden-colored).

Loads canonical final_adata.h5ad, plots UMAP colored by Leiden cluster.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", required=True, type=Path,
                   help="path to final_adata.h5ad")
    p.add_argument("--obsm-key", default="X_umap",
                   help="obsm embedding key (default X_umap)")
    p.add_argument("--color-key", default="leiden",
                   help="obs column for coloring (default leiden)")
    p.add_argument("--out-png", required=True, type=Path,
                   help="output PNG path")
    p.add_argument("--title", default=None, type=str)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    adata = sc.read_h5ad(args.h5ad)
    if args.obsm_key not in adata.obsm:
        raise KeyError(
            f"obsm key {args.obsm_key!r} not found; available: {list(adata.obsm.keys())}"
        )
    if args.color_key not in adata.obs.columns:
        raise KeyError(
            f"obs column {args.color_key!r} not found; available: {list(adata.obs.columns)}"
        )

    args.out_png.parent.mkdir(parents=True, exist_ok=True)

    # Use sc.pl.embedding to honor obsm-key explicitly (basis = obsm_key without X_ prefix).
    basis = args.obsm_key
    if basis.startswith("X_"):
        basis = basis[2:]
    fig = sc.pl.embedding(
        adata,
        basis=basis,
        color=args.color_key,
        show=False,
        return_fig=True,
        title=args.title,
    )
    fig.savefig(args.out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    width, height = fig.get_size_inches() * fig.dpi
    print(f"[wave5_plot_umap_rna] wrote {args.out_png}")
    print(f"  approx pixel size: {int(width)} x {int(height)}")
    print(f"  obsm={args.obsm_key}  color={args.color_key}  n_obs={adata.n_obs}")


if __name__ == "__main__":
    main()
