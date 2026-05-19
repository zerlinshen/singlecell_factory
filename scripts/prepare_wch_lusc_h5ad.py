from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare WCH lung cancer atlas subset as h5ad")
    p.add_argument("--dataset-dir", required=True)
    p.add_argument("--disease", default="LUSC")
    p.add_argument("--stages", default="I,II,III")
    p.add_argument("--output", default="prepared_input.h5ad")
    p.add_argument("--chunksize", type=int, default=500)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    meta_path = dataset_dir / "meta.csv"
    expr_path = dataset_dir / "raw_data.csv.gz"
    out_path = dataset_dir / args.output
    wanted_stages = {s.strip() for s in args.stages.split(",") if s.strip()}

    meta = pd.read_csv(meta_path, low_memory=False)
    meta = meta[(meta["Disease"] == args.disease) & (meta["Stage"].isin(wanted_stages))].copy()
    if meta.empty:
        raise SystemExit("No cells matched requested disease/stages")
    meta = meta.drop_duplicates(subset=["Cells"]).set_index("Cells")

    with gzip.open(expr_path, "rt", newline="") as f:
        header = next(csv.reader(f))
    matrix_cells = header[1:]
    selected_cells = [c for c in matrix_cells if c in meta.index]
    if not selected_cells:
        raise SystemExit("No selected cells found in expression matrix header")

    first_col = pd.read_csv(expr_path, compression="gzip", nrows=0).columns[0]
    usecols = [first_col, *selected_cells]
    obs = meta.loc[selected_cells].copy()
    obs.index.name = None
    obs["sample"] = obs["SampleID"].astype(str)
    obs["patient"] = obs["PatientID"].astype(str)
    obs["stage"] = obs["Stage"].astype(str)
    obs["disease"] = obs["Disease"].astype(str)
    obs["tissue"] = obs["Tissue"].astype(str)
    obs["cell_type_hint"] = obs["CellName"].astype(str)

    blocks = []
    genes: list[str] = []
    for chunk in pd.read_csv(
        expr_path,
        compression="gzip",
        usecols=usecols,
        index_col=0,
        chunksize=args.chunksize,
    ):
        genes.extend(chunk.index.astype(str).tolist())
        arr = chunk.to_numpy(dtype=np.float32, copy=False)
        blocks.append(sparse.csr_matrix(arr.T))

    X = sparse.hstack(blocks, format="csr")
    var = pd.DataFrame(index=pd.Index(genes, name=None))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    adata.uns["source_dataset"] = "wch_lung_cancer_atlas"
    adata.uns["subset_disease"] = args.disease
    adata.uns["subset_stages"] = sorted(wanted_stages)
    adata.write_h5ad(out_path, compression="gzip")

    print(f"Saved {out_path}")
    print(f"cells={adata.n_obs} genes={adata.n_vars}")
    print(f"stages={sorted(obs['stage'].value_counts().to_dict().items())}")


if __name__ == "__main__":
    main()
