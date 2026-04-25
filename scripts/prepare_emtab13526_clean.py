from __future__ import annotations

from pathlib import Path
import gzip
import shutil

import anndata as ad
from anndata.experimental import concat_on_disk
import numpy as np
import pandas as pd
from scipy import io

ROOT = Path('/home/zerlinshen/singlecell_factory/data/raw/nc2024_nsclc_emtab13526')
SDRF = ROOT / 'E-MTAB-13526.sdrf.txt'
PARTS_ROOT = ROOT / '_prepared_parts_clean'


def read_features(path: Path) -> pd.DataFrame:
    with gzip.open(path, 'rt') as fh:
        df = pd.read_csv(fh, sep='\t', header=None)
    if df.shape[1] == 3:
        df.columns = ['gene_id', 'gene_name', 'feature_type']
    elif df.shape[1] == 2:
        df.columns = ['gene_id', 'gene_name']
        df['feature_type'] = 'Gene Expression'
    else:
        raise ValueError(f'Unexpected feature file format: {path}')
    return df


def read_barcodes(path: Path) -> pd.Series:
    with gzip.open(path, 'rt') as fh:
        return pd.read_csv(fh, sep='\t', header=None)[0]


def read_matrix(path: Path):
    with gzip.open(path, 'rb') as fh:
        return io.mmread(fh).tocsr().T.tocsr()


def build_sample_meta() -> pd.DataFrame:
    s = pd.read_csv(SDRF, sep='\t')
    keep = [
        'Source Name',
        'Characteristics[individual]',
        'Characteristics[disease]',
        'Characteristics[FACS]',
        'Characteristics[sampling site]',
        'Characteristics[sex]',
        'Characteristics[original source name]',
    ]
    meta = s[keep].drop_duplicates().rename(columns={
        'Source Name': 'sample',
        'Characteristics[individual]': 'patient',
        'Characteristics[disease]': 'disease',
        'Characteristics[FACS]': 'sorting',
        'Characteristics[sampling site]': 'sampling_site',
        'Characteristics[sex]': 'sex',
        'Characteristics[original source name]': 'original_source_name',
    })
    meta = meta.astype(str)
    meta['condition'] = np.where(meta['disease'].str.lower().eq('normal'), 'healthy_background', 'tumor')
    meta['tumor_type'] = np.where(meta['condition'].eq('tumor'), 'NSCLC', 'non_involved')
    meta['batch'] = meta['sample']
    return meta.sort_values('sample').reset_index(drop=True)


def make_var_out(var_ref: pd.DataFrame) -> pd.DataFrame:
    var_out = var_ref.copy()
    var_out.index = pd.Index(var_out['gene_name'].astype(str))
    var_out.index.name = None
    var_out['gene_symbol'] = var_out['gene_name'].astype(str)
    var_out = var_out.drop(columns=['gene_name'])
    return var_out


def write_part(sample: str, subset: str, mat, obs: pd.DataFrame, var_ref: pd.DataFrame) -> Path:
    var_out = make_var_out(var_ref)
    adata = ad.AnnData(X=mat, obs=obs, var=var_out)
    adata.var.index.name = None
    adata.var_names_make_unique()
    part_path = PARTS_ROOT / subset / f'{sample}.h5ad'
    print(f'WRITE_PART {sample} -> {part_path}', flush=True)
    adata.write_h5ad(part_path)
    return part_path


def main() -> None:
    sample_meta = build_sample_meta().set_index('sample')
    matrix_files = sorted(ROOT.glob('*-matrix.mtx.gz'))
    if not matrix_files:
        raise SystemExit('No matrix files found')

    if PARTS_ROOT.exists():
        shutil.rmtree(PARTS_ROOT)
    (PARTS_ROOT / 'tumor').mkdir(parents=True, exist_ok=True)
    (PARTS_ROOT / 'healthy_background').mkdir(parents=True, exist_ok=True)

    var_ref = None
    part_paths = {'tumor': [], 'healthy_background': []}

    for matrix_path in matrix_files:
        sample = matrix_path.name.replace('-matrix.mtx.gz', '')
        print(f'READING {sample}', flush=True)
        feature_path = ROOT / f'{sample}-features.tsv.gz'
        barcode_path = ROOT / f'{sample}-barcodes.tsv.gz'
        if not feature_path.exists() or not barcode_path.exists():
            raise FileNotFoundError(f'Missing feature/barcode files for {sample}')
        if sample not in sample_meta.index:
            raise KeyError(f'Sample {sample} missing from SDRF metadata')

        var = read_features(feature_path)
        if var_ref is None:
            var_ref = var.copy()
        else:
            if not var[['gene_id', 'gene_name']].equals(var_ref[['gene_id', 'gene_name']]):
                raise ValueError(f'Feature mismatch for sample {sample}')

        barcodes = read_barcodes(barcode_path).astype(str)
        mat = read_matrix(matrix_path)
        if mat.shape[0] != len(barcodes):
            raise ValueError(f'Barcode mismatch for {sample}: {mat.shape[0]} vs {len(barcodes)}')

        meta_row = sample_meta.loc[sample]
        cell_ids = pd.Index([f'{sample}:{bc}' for bc in barcodes], name='cell_id')
        obs = pd.DataFrame(index=cell_ids)
        obs['sample'] = sample
        obs['patient'] = meta_row['patient']
        obs['batch'] = meta_row['batch']
        obs['disease'] = meta_row['disease']
        obs['condition'] = meta_row['condition']
        obs['sorting'] = meta_row['sorting']
        obs['sampling_site'] = meta_row['sampling_site']
        obs['sex'] = meta_row['sex']
        obs['original_source_name'] = meta_row['original_source_name']
        obs['tumor_type'] = meta_row['tumor_type']

        subset = str(meta_row['condition'])
        part_path = write_part(sample, subset, mat, obs, var_ref)
        part_paths[subset].append(str(part_path))

    summary = sample_meta.reset_index()[['sample','patient','condition','sorting','sampling_site','sex','original_source_name']].drop_duplicates().sort_values(['condition','patient','sample'])
    summary.to_csv(ROOT / 'sample_summary.csv', index=False)

    for subset in ['tumor', 'healthy_background']:
        out_dir = ROOT / subset
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / 'prepared_input.zarr'
        if out_path.exists():
            shutil.rmtree(out_path)
        print(f'CONCAT_{subset.upper()} -> {out_path}', flush=True)
        concat_on_disk(part_paths[subset], out_path, axis=0, join='outer', max_loaded_elems=50_000_000)
        (out_dir / 'prepared_input.ready').write_text('ok\n')
        print(f'CONCAT_DONE_{subset.upper()}', flush=True)

if __name__ == '__main__':
    main()
