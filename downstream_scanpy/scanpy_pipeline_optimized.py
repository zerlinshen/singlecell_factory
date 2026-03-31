"""
Optimized Scanpy Pipeline for scRNA-seq Analysis
Enhanced with better QC, visualization, and modular structure.

Steps:
1. Read input (10x MTX, .h5, or .h5ad)
2. Comprehensive QC metrics (mito, ribo, doublets)
3. QC filtering with visualizations
4. Normalization, log-transform, HVG selection
5. Scaling, PCA (elbow-guided n_pcs), neighbors (adaptive k), UMAP
6. Leiden clustering
7. Marker gene detection
8. Cell type annotation (marker-based + optional CellTypist)
9. Comprehensive visualization output
"""

import argparse
import os
import logging
from datetime import datetime
import re
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Dict, List, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTPUT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "output"))
TIMESTAMP_PATTERN = re.compile(r'_\d{8}_\d{4,6}$')

# Default configuration
DEFAULT_CONFIG = {
    # QC thresholds
    'min_genes': 200,
    'max_genes': 6000,
    'max_mito_pct': 10.0,
    'max_ribo_pct': 50.0,
    'min_cells': 3,

    # Analysis parameters
    'n_top_genes': 2000,
    'n_pcs': 40,             # hard cap; actual PCs chosen by elbow (85% variance)
    'auto_n_pcs': True,      # Fix 2: use elbow-guided PCs instead of fixed 40
    'n_neighbors': None,     # Fix 1: None → adaptive k = round(log2(N)*2), capped 10-30
    'leiden_resolution': 0.8,
    'target_sum': 1e4,
    'random_state': 0,

    # Doublet detection
    'run_doublet_detection': True,
    'doublet_threshold': 0.25,
    'fail_on_doublet_error': True,

    # CellTypist reference annotation (Fix 4)
    'run_celltypist': True,
    'celltypist_model': 'Immune_All_Low.pkl',  # fine-grained immune model

    # Marker genes — Fix 3: expanded subtype resolution
    'pbmc_markers': {
        # T cells
        "Naive CD4 T":   ["IL7R", "CCR7", "TCF7", "SELL"],
        "Memory CD4 T":  ["IL7R", "S100A4", "CD44"],
        "Treg":          ["FOXP3", "IL2RA", "CTLA4", "TIGIT"],
        "CD8 T":         ["CD8A", "CD8B", "GZMB", "PRF1"],
        "NK":            ["GNLY", "NKG7", "FCGR3A", "KLRD1"],
        # B cells — Fix 3: subtype resolution matching 10x
        "Naive B":       ["MS4A1", "CD79A", "IGHD", "TCL1A", "IL4R"],
        "Memory B":      ["MS4A1", "CD79A", "CD27", "AIM2", "TNFRSF13B"],
        "Plasma B":      ["IGHG1", "IGKC", "MZB1", "SDC1", "PRDM1"],
        "B (Ig light)":  ["IGKV1-33", "IGLV4-69", "IGLC2", "IGLC3"],
        # Myeloid
        "CD14+ Mono":    ["LYZ", "CD14", "S100A8", "S100A9", "VCAN"],
        "FCGR3A+ Mono":  ["FCGR3A", "MS4A7", "RHOC"],
        "Macrophage":    ["CD68", "SPP1", "MRC1", "MARCO", "C1QA"],
        "DC":            ["FCER1A", "CST3", "CLEC9A", "CD1C"],
        # Other
        "Platelet":      ["PPBP", "PF4", "GP1BA"],
        # Tumour / stromal (used when running on non-PBMC datasets)
        "Tumour (SCC)":  ["KRT5", "KRT6A", "KRT14", "TP63", "SOX2"],
        "Fibroblast":    ["COL1A1", "FAP", "ACTA2", "VIM"],
        "Endothelial":   ["PECAM1", "VWF", "CLDN5"],
        "Mast cell":     ["TPSAB1", "CPA3", "KIT"],
    }
}


def build_timestamped_project_name(project_name: str) -> str:
    """Append a timestamp unless the name already ends with one."""
    normalized = re.sub(r'[^A-Za-z0-9._-]+', '_', project_name.strip()).strip('_')
    if not normalized:
        normalized = 'unnamed_project'
    if TIMESTAMP_PATTERN.search(normalized):
        return normalized
    return f"{normalized}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"


def setup_output_directories(
    project_name: str,
    base_output_dir: str = DEFAULT_OUTPUT_ROOT
) -> Dict[str, str]:
    """Create standardized project-based output directory structure.
    
    Structure:
        output/
        └── <project_name>/
            ├── data/           # Processed AnnData files
            ├── figures/        # All visualizations
            │   ├── qc/         # QC plots
            │   ├── analysis/   # Dimensionality reduction plots
            │   └── markers/    # Marker gene plots
            ├── tables/         # CSV results
            ├── logs/           # Pipeline logs and configs
            └── cache/          # Temporary/cache files (isolated per project)
    """
    # Main project directory
    project_dir = os.path.abspath(os.path.join(base_output_dir, project_name))
    
    dirs = {
        'root': project_dir,
        'data': os.path.join(project_dir, 'data'),
        'figures': os.path.join(project_dir, 'figures'),
        'qc': os.path.join(project_dir, 'figures', 'qc'),
        'analysis': os.path.join(project_dir, 'figures', 'analysis'),
        'markers': os.path.join(project_dir, 'figures', 'markers'),
        'tables': os.path.join(project_dir, 'tables'),
        'logs': os.path.join(project_dir, 'logs'),
        'cache': os.path.join(project_dir, 'cache')  # Project-specific cache
    }
    
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    
    logger.info(f"Project output directory: {project_dir}")
    return dirs


def load_data(
    tenx_dir: Optional[str] = None,
    input_h5ad: Optional[str] = None,
    tenx_var_names: str = 'gene_symbols',
    cache_dir: Optional[str] = None
) -> sc.AnnData:
    """Load data from 10x directory or h5ad file.
    
    Args:
        tenx_dir: Path to 10x filtered_feature_bc_matrix directory
        input_h5ad: Path to existing h5ad file
        tenx_var_names: Gene identifier type
        cache_dir: Project-specific cache directory (isolated per project)
    """
    if tenx_dir and os.path.isdir(tenx_dir):
        logger.info(f"Reading 10x MTX from: {tenx_dir}")
        
        # Use project-specific cache directory to avoid cross-project contamination
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
            # Change to cache directory for scanpy's cache
            original_dir = os.getcwd()
            os.chdir(cache_dir)
            try:
                adata = sc.read_10x_mtx(
                    tenx_dir,
                    var_names=tenx_var_names,
                    cache=True  # Cache in project-specific directory
                )
            finally:
                os.chdir(original_dir)
        else:
            adata = sc.read_10x_mtx(
                tenx_dir,
                var_names=tenx_var_names,
                cache=False  # No caching if no cache_dir specified
            )
            
    elif input_h5ad and os.path.isfile(input_h5ad):
        # Cell Ranger .h5 (filtered_feature_bc_matrix.h5) vs AnnData .h5ad
        if input_h5ad.endswith('.h5') and not input_h5ad.endswith('.h5ad'):
            logger.info(f"Reading Cell Ranger .h5 from: {input_h5ad}")
            adata = sc.read_10x_h5(input_h5ad)
            adata.var_names_make_unique()
        else:
            logger.info(f"Reading h5ad from: {input_h5ad}")
            adata = sc.read_h5ad(input_h5ad)
    else:
        raise ValueError("Must provide either --tenx_dir or --input")

    logger.info(f"Loaded dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata


def _gene_symbol_series(adata: sc.AnnData) -> pd.Series:
    """Return uppercase gene symbols aligned to adata.var_names."""
    if 'gene_symbols' in adata.var.columns:
        symbols = adata.var['gene_symbols'].astype(str)
    else:
        symbols = adata.var_names.astype(str).to_series(index=adata.var_names)
    return symbols.str.upper()


def _resolve_marker_var_names_from_raw(adata: sc.AnnData, markers: List[str]) -> List[str]:
    """Resolve marker symbols to raw.var_names for both gene_symbols/gene_ids modes."""
    if adata.raw is None:
        source = adata
    else:
        source = adata.raw.to_adata()
    if 'gene_symbols' in source.var.columns:
        symbol_series = pd.Series(source.var_names, index=source.var['gene_symbols'].astype(str))
        symbol_series = symbol_series[~symbol_series.index.duplicated(keep='first')]
        return [symbol_series[g] for g in markers if g in symbol_series.index]
    return [g for g in markers if g in source.var_names]


def calculate_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """Calculate comprehensive QC metrics."""
    symbols_upper = _gene_symbol_series(adata)

    # Mitochondrial genes (MT- prefix)
    adata.var['mt'] = symbols_upper.str.startswith('MT-').values

    # Ribosomal genes (RPS, RPL prefixes)
    adata.var['ribo'] = (
        symbols_upper.str.startswith('RPS') |
        symbols_upper.str.startswith('RPL')
    ).values

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo'],
        inplace=True,
        percent_top=(50, 100, 200, 500)
    )

    # Additional derived metrics
    adata.obs['log1p_n_genes_by_counts'] = np.log1p(adata.obs['n_genes_by_counts'])
    adata.obs['log1p_total_counts'] = np.log1p(adata.obs['total_counts'])
    adata.obs['n_counts_per_gene'] = adata.obs['total_counts'] / adata.obs['n_genes_by_counts']

    logger.info("QC metrics calculated")
    return adata


def suggest_qc_thresholds(adata: sc.AnnData, config: Dict) -> Dict:
    """Suggest QC thresholds based on data distribution (MAD method).
    
    Uses Median Absolute Deviation (MAD) for robust outlier detection.
    Returns suggested thresholds that may be less stringent than defaults.
    """
    suggestions = {}
    
    # Gene count thresholds using MAD
    n_genes = adata.obs['n_genes_by_counts']
    gene_median = n_genes.median()
    gene_mad = np.median(np.abs(n_genes - gene_median))
    
    # Suggest min_genes: median - 3*MAD or 5th percentile, whichever is higher
    suggested_min_genes = max(int(gene_median - 3 * gene_mad), int(n_genes.quantile(0.05)))
    suggestions['min_genes'] = max(suggested_min_genes, 100)  # At least 100 genes
    
    # Suggest max_genes: median + 5*MAD or 99th percentile for doublet detection
    suggested_max_genes = min(int(gene_median + 5 * gene_mad), int(n_genes.quantile(0.99)))
    suggestions['max_genes'] = suggested_max_genes
    
    # Mitochondrial threshold using percentile
    mito_pct = adata.obs['pct_counts_mt']
    # Suggest threshold at 95th percentile if most cells are healthy
    suggested_mito = min(float(mito_pct.quantile(0.95)), 20.0)  # Cap at 20%
    suggestions['max_mito_pct'] = max(suggested_mito, 5.0)  # At least 5%
    
    # Calculate expected cell retention with suggested thresholds
    mask = (
        (adata.obs['n_genes_by_counts'] > suggestions['min_genes']) &
        (adata.obs['n_genes_by_counts'] < suggestions['max_genes']) &
        (adata.obs['pct_counts_mt'] < suggestions['max_mito_pct'])
    )
    expected_retention = mask.sum() / len(adata) * 100
    suggestions['expected_retention_pct'] = expected_retention
    
    logger.info("=" * 60)
    logger.info("ADAPTIVE QC THRESHOLD SUGGESTIONS (MAD-based)")
    logger.info("=" * 60)
    logger.info(f"Current config -> min_genes: {config['min_genes']}, "
                f"max_mito: {config['max_mito_pct']}%")
    logger.info(f"Suggested      -> min_genes: {suggestions['min_genes']}, "
                f"max_genes: {suggestions['max_genes']}, "
                f"max_mito: {suggestions['max_mito_pct']:.1f}%")
    logger.info(f"Expected retention with suggested: {expected_retention:.1f}% "
                f"({mask.sum()}/{len(adata)} cells)")
    
    # Warning if current config is too stringent
    current_mask = (
        (adata.obs['n_genes_by_counts'] > config['min_genes']) &
        (adata.obs['n_genes_by_counts'] < config['max_genes']) &
        (adata.obs['pct_counts_mt'] < config['max_mito_pct'])
    )
    current_retention = current_mask.sum() / len(adata) * 100
    
    if current_retention < 20:
        logger.warning(f"⚠️  Current config retains only {current_retention:.1f}% of cells!")
        logger.warning(f"   Consider using: --min_genes={suggestions['min_genes']} "
                      f"--max_mito_pct={suggestions['max_mito_pct']:.1f}")
    
    logger.info("=" * 60)
    
    return suggestions


def run_doublet_detection(adata: sc.AnnData, config: Dict) -> sc.AnnData:
    """Run Scrublet doublet detection.
    
    Scrublet simulates doublets from the data and uses a nearest-neighbor
    classifier to identify potential doublets in the dataset.
    """
    if not config.get('run_doublet_detection', True):
        logger.info("Skipping doublet detection (disabled in config)")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        return adata
    
    try:
        # Run Scrublet
        logger.info("Running Scrublet doublet detection...")
        sc.pp.scrublet(
            adata,
            batch_key=None,
            threshold=config['doublet_threshold'],
            random_state=config.get('random_state', 0)
        )
        
        # Calculate doublet statistics
        n_doublets = adata.obs['predicted_doublet'].sum()
        doublet_pct = n_doublets / len(adata) * 100
        
        logger.info(f"Doublet detection complete:")
        logger.info(f"  Predicted doublets: {n_doublets} ({doublet_pct:.2f}%)")
        logger.info(f"  Mean doublet score: {adata.obs['doublet_score'].mean():.3f}")
        
        # Warning if doublet rate is unusually high
        if doublet_pct > 15:
            logger.warning(f"⚠️  High doublet rate detected ({doublet_pct:.1f}%)!")
            logger.warning(f"   Expected: 2-8% depending on cell loading")
        
        return adata
        
    except Exception as e:
        if config.get('fail_on_doublet_error', True):
            raise RuntimeError(f"Doublet detection failed and fail_on_doublet_error=True: {e}") from e
        logger.warning(f"Doublet detection failed: {e}")
        logger.warning("Continuing without doublet detection (--allow_doublet_fail)")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        return adata


def plot_qc_metrics(adata: sc.AnnData, qc_dir: str, config: Dict) -> None:
    """Plot comprehensive QC visualizations."""
    sc.settings.figdir = qc_dir
    sc.settings.autoshow = False

    # Violin plots of key metrics
    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']
    for metric in metrics:
        sc.pl.violin(adata, metric, jitter=0.4, size=1, save=f'_{metric}.png')
        plt.close()

    # Scatter plots
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_counts_vs_mito.png')
    plt.close()

    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_counts_vs_genes.png')
    plt.close()

    sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt', save='_genes_vs_mito.png')
    plt.close()

    # Histograms
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    axes[0].hist(adata.obs['n_genes_by_counts'], bins=50, alpha=0.7)
    axes[0].axvline(config['min_genes'], color='red', linestyle='--')
    axes[0].axvline(config['max_genes'], color='red', linestyle='--')
    axes[0].set_title('Number of genes per cell')

    axes[1].hist(adata.obs['total_counts'], bins=50, alpha=0.7, log=True)
    axes[1].set_title('Total counts per cell (log scale)')

    axes[2].hist(adata.obs['pct_counts_mt'], bins=50, alpha=0.7)
    axes[2].axvline(config['max_mito_pct'], color='red', linestyle='--')
    axes[2].set_title('Mitochondrial percentage')

    axes[3].hist(adata.obs['pct_counts_ribo'], bins=50, alpha=0.7)
    axes[3].set_title('Ribosomal percentage')

    plt.tight_layout()
    plt.savefig(os.path.join(qc_dir, 'qc_histograms.png'), dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"QC plots saved to: {qc_dir}")


def apply_qc_filtering(adata: sc.AnnData, config: Dict) -> Tuple[sc.AnnData, Dict]:
    """Apply QC filtering and track stats."""
    # Store initial stats
    initial_cells = adata.n_obs
    initial_genes = adata.n_vars
    initial_doublets = adata.obs['predicted_doublet'].sum() if 'predicted_doublet' in adata.obs.columns else 0

    # Filter cells by gene count, mito percentage, and doublets
    mask = (
        (adata.obs['n_genes_by_counts'] > config['min_genes']) &
        (adata.obs['n_genes_by_counts'] < config['max_genes']) &
        (adata.obs['pct_counts_mt'] < config['max_mito_pct']) &
        (adata.obs['pct_counts_ribo'] < config['max_ribo_pct'])
    )
    
    # Add doublet filtering if available
    if 'predicted_doublet' in adata.obs.columns:
        mask = mask & (~adata.obs['predicted_doublet'])
    
    adata = adata[mask, :].copy()

    if adata.n_obs == 0:
        logger.error("No cells remain after QC filtering. Check your thresholds.")
        raise ValueError("Empty dataset after QC filtering")

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=config['min_cells'])

    # Track filtering stats
    final_doublets_remaining = adata.obs['predicted_doublet'].sum() if 'predicted_doublet' in adata.obs.columns else 0
    
    filter_stats = {
        'initial_cells': initial_cells,
        'initial_genes': initial_genes,
        'initial_doublets': int(initial_doublets),
        'final_cells': adata.n_obs,
        'final_genes': adata.n_vars,
        'final_doublets': int(final_doublets_remaining),
        'cells_removed': initial_cells - adata.n_obs,
        'genes_removed': initial_genes - adata.n_vars,
        'doublets_removed': int(initial_doublets - final_doublets_remaining),
        'cells_passed_pct': (adata.n_obs / initial_cells) * 100,
        'genes_passed_pct': (adata.n_vars / initial_genes) * 100,
        'doublet_removal_rate_pct': (initial_doublets - final_doublets_remaining) / initial_doublets * 100 if initial_doublets > 0 else 0
    }

    logger.info(f"QC Filtering Results:")
    logger.info(f"  Initial: {initial_cells} cells, {initial_genes} genes, {initial_doublets} doublets")
    logger.info(f"  Removed: {filter_stats['cells_removed']} cells ({100-filter_stats['cells_passed_pct']:.1f}%), "
                f"{filter_stats['genes_removed']} genes, {filter_stats['doublets_removed']} doublets")
    logger.info(f"  Final: {adata.n_obs} cells, {adata.n_vars} genes")

    return adata, filter_stats


def normalize_and_scale(adata: sc.AnnData, config: Dict) -> sc.AnnData:
    """Normalize, log-transform, select HVGs, and scale data."""
    # Preserve raw integer counts before any normalization
    adata.layers['counts'] = adata.X.copy()

    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=config['target_sum'])
    sc.pp.log1p(adata)

    # Store all genes at log-normalized scale (used by rank_genes_groups)
    adata.raw = adata

    # Select highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        flavor='seurat',
        n_top_genes=config['n_top_genes'],
        batch_key=None
    )
    n_hvg = adata.var['highly_variable'].sum()

    # Subset to HVGs, then scale in-place (avoids redundant intermediate copy)
    adata = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata, max_value=10)

    logger.info(f"Normalization complete. {n_hvg} HVGs selected")
    return adata


def compute_elbow_pcs(adata: sc.AnnData, target_variance: float = 0.85, max_pcs: int = 40) -> int:
    """Fix 2: Choose n_pcs by finding the elbow where cumulative variance ≥ target_variance.

    Args:
        adata: AnnData after sc.tl.pca has been called.
        target_variance: Cumulative variance ratio threshold (default 0.85 = 85%).
        max_pcs: Hard upper cap on PCs to use.

    Returns:
        Number of PCs to use for neighbour graph.
    """
    var_ratio = adata.uns['pca']['variance_ratio']
    cumvar = np.cumsum(var_ratio)
    elbow = int(np.searchsorted(cumvar, target_variance)) + 1
    n_pcs = max(10, min(elbow, max_pcs, len(var_ratio)))
    logger.info(f"Elbow PCs: {n_pcs} ({cumvar[n_pcs-1]*100:.1f}% variance @ target {target_variance*100:.0f}%)")
    return n_pcs


def compute_adaptive_k(n_cells: int) -> int:
    """Fix 1: Dataset-adaptive k for neighbour graph, mirroring 10x log-scaling.

    k = round(log2(n_cells) * 2), capped between 10 and 30.
    Examples:
        1,000 cells  → k = 20
        2,500 cells  → k = 22
        10,000 cells → k = 27
        50,000 cells → k = 31 → capped at 30
    """
    import math
    if n_cells > 0:
        k = max(10, min(30, round(math.log2(n_cells) * 2)))
        logger.info(f"Adaptive k: {k} for {n_cells} cells (log2({n_cells})*2 = {math.log2(n_cells)*2:.1f})")
    else:
        k = 15
        logger.info(f"Adaptive k: {k} (fallback for {n_cells} cells)")
    return k


def run_dimred_and_clustering(adata: sc.AnnData, config: Dict) -> sc.AnnData:
    """Run PCA, build neighbor graph (adaptive k), UMAP, and Leiden clustering."""
    n_comps = config['n_pcs']

    # PCA — compute up to the configured cap
    sc.tl.pca(
        adata,
        n_comps=n_comps,
        svd_solver='arpack',
        random_state=config.get('random_state', 0)
    )

    # Fix 2: elbow-guided n_pcs (overrides cap if auto_n_pcs=True)
    if config.get('auto_n_pcs', True):
        n_pcs_use = compute_elbow_pcs(adata, target_variance=0.85, max_pcs=n_comps)
    else:
        n_pcs_use = min(n_comps, adata.obsm['X_pca'].shape[1])

    # Fix 1: adaptive k (use config value if explicitly set, else auto)
    if config.get('n_neighbors') is None:
        k = compute_adaptive_k(adata.n_obs)
    else:
        k = config['n_neighbors']
        logger.info(f"Using configured k={k}")

    # Store chosen values for logging/reproducibility
    adata.uns['pipeline_params'] = {
        'n_pcs_used': n_pcs_use,
        'n_neighbors_used': k,
        'leiden_resolution': config['leiden_resolution'],
    }

    # Neighbour graph
    sc.pp.neighbors(adata, n_neighbors=k, n_pcs=n_pcs_use)

    # UMAP
    sc.tl.umap(adata, random_state=config.get('random_state', 0))

    # Leiden clustering
    sc.tl.leiden(
        adata,
        resolution=config['leiden_resolution'],
        flavor='igraph',
        directed=False,
        random_state=config.get('random_state', 0)
    )

    n_clusters = adata.obs['leiden'].nunique()
    logger.info(f"Clustering: {n_clusters} Leiden clusters "
                f"(res={config['leiden_resolution']}, k={k}, n_pcs={n_pcs_use})")
    return adata


def plot_dimred_results(adata: sc.AnnData, analysis_dir: str) -> None:
    """Plot PCA and UMAP results."""
    sc.settings.figdir = analysis_dir
    sc.settings.autoshow = False

    # PCA variance — only plot as many PCs as were actually computed
    n_pcs_available = adata.obsm['X_pca'].shape[1]
    sc.pl.pca_variance_ratio(adata, log=True, save='.png', n_pcs=n_pcs_available)
    plt.close()

    # PCA plots
    sc.pl.pca(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], save='_qc.png')
    plt.close()

    sc.pl.pca(adata, color=['leiden'], save='_leiden.png')
    plt.close()

    # UMAP plots
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data', save='_leiden.png')
    plt.close()

    sc.pl.umap(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], save='_qc.png')
    plt.close()

    # UMAP coloured by doublet score (for QC)
    if 'doublet_score' in adata.obs.columns:
        sc.pl.umap(adata, color=['doublet_score'], save='_doublet_score.png')
        plt.close()

    logger.info(f"Dimred plots saved to: {analysis_dir}")


def find_marker_genes(
    adata: sc.AnnData,
    groupby: str = 'leiden',
    method: str = 'wilcoxon'
) -> sc.AnnData:
    """Find marker genes per cluster."""
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        use_raw=True,
        pts=True,
        tie_correct=True
    )
    logger.info("Marker genes identified")
    return adata


def plot_marker_genes(
    adata: sc.AnnData,
    markers_dir: str,
    marker_dict: Optional[Dict[str, List[str]]] = None
) -> None:
    """Plot marker gene visualizations."""
    sc.settings.figdir = markers_dir
    sc.settings.autoshow = False

    # Top marker genes per cluster
    sc.pl.rank_genes_groups(
        adata,
        n_genes=20,
        sharey=False,
        save='_top20.png'
    )
    plt.close()

    # Dotplot for top markers
    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=5,
        save='_top5.png'
    )
    plt.close()

    # Heatmap
    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=5,
        swap_axes=True,
        show_gene_labels=True,
        save='_top5.png'
    )
    plt.close()

    # Custom markers if provided
    if marker_dict:
        # Collect valid marker genes
        var_names = []
        seen = set()
        for _, genes in marker_dict.items():
            resolved = _resolve_marker_var_names_from_raw(adata, genes)
            for g in resolved:
                if g not in seen:
                    var_names.append(g)
                    seen.add(g)

        if var_names:
            # Dotplot of canonical markers
            fig = sc.pl.dotplot(
                adata,
                var_names=var_names,
                groupby='leiden',
                standard_scale='var',
                return_fig=True
            )
            fig.savefig(
                os.path.join(markers_dir, 'dotplot_canonical_markers.png'),
                bbox_inches='tight',
                dpi=150
            )
            plt.close()

            # Violin plots for key markers
            if len(var_names) > 0:
                # Plot in batches of 3
                for i in range(0, min(len(var_names), 12), 3):
                    sc.pl.violin(
                        adata,
                        keys=var_names[i:i+3],
                        groupby='leiden',
                        save=f'_canonical_markers_{i//3+1}.png'
                    )
                    plt.close()

    logger.info(f"Marker plots saved to: {markers_dir}")


def annotate_cell_types(
    adata: sc.AnnData,
    marker_dict: Dict[str, List[str]],
    cluster_key: str = 'leiden'
) -> Tuple[sc.AnnData, pd.DataFrame]:
    """Assign cluster and cell-level annotations by canonical marker-set scores."""
    if adata.raw is None:
        raise ValueError("adata.raw is required for marker-based annotation")
    raw = adata.raw.to_adata()
    expr = raw.to_df()
    expr[cluster_key] = adata.obs[cluster_key].values

    score_by_type = {}
    n_markers_used = {}
    for cell_type, markers in marker_dict.items():
        resolved = _resolve_marker_var_names_from_raw(adata, markers)
        if len(resolved) == 0:
            continue
        n_markers_used[cell_type] = len(resolved)
        score_by_type[cell_type] = expr.groupby(cluster_key, observed=False)[resolved].mean().mean(axis=1)

    if not score_by_type:
        adata.obs['cell_type'] = 'Unknown'
        adata.obs['cell_type_confidence'] = 0.0
        empty = pd.DataFrame(columns=['cluster', 'cell_type', 'score_best', 'score_second', 'confidence', 'n_markers'])
        return adata, empty

    score_df = pd.DataFrame(score_by_type).fillna(0.0)
    best_type = score_df.idxmax(axis=1)
    best_score = score_df.max(axis=1)
    second_score = score_df.apply(
        lambda row: row.sort_values(ascending=False).iloc[1] if len(row) > 1 else 0.0,
        axis=1
    )
    confidence = best_score - second_score

    cluster_annotation = pd.DataFrame({
        'cluster': score_df.index.astype(str),
        'cell_type': best_type.astype(str).values,
        'score_best': best_score.values,
        'score_second': second_score.values,
        'confidence': confidence.values,
    })
    cluster_annotation['n_markers'] = cluster_annotation['cell_type'].map(n_markers_used).fillna(0).astype(int)

    cluster_to_type = dict(zip(cluster_annotation['cluster'], cluster_annotation['cell_type']))
    cluster_to_conf = dict(zip(cluster_annotation['cluster'], cluster_annotation['confidence']))
    adata.obs['cell_type'] = adata.obs[cluster_key].astype(str).map(cluster_to_type).fillna('Unknown')
    adata.obs['cell_type_confidence'] = adata.obs[cluster_key].astype(str).map(cluster_to_conf).fillna(0.0).astype(float)

    adata.uns['cluster_annotation'] = cluster_annotation.set_index('cluster').to_dict(orient='index')
    logger.info("Cell-type annotation assigned from marker-set scoring")
    return adata, cluster_annotation


def run_celltypist(adata: sc.AnnData, config: Dict) -> sc.AnnData:
    """Fix 4: Reference-based cell type annotation using CellTypist.

    Downloads a pre-trained model and annotates each cell independently of
    the marker scoring approach. Results stored in:
        adata.obs['cell_type_celltypist']   — majority-vote per cell
        adata.obs['cell_type_celltypist_conf'] — confidence score
    """
    if not config.get('run_celltypist', True):
        logger.info("Skipping CellTypist (disabled in config)")
        return adata

    try:
        import celltypist
        from celltypist import models

        model_name = config.get('celltypist_model', 'Immune_All_Low.pkl')
        logger.info(f"Running CellTypist with model: {model_name}")

        # CellTypist expects log-normalised counts — use adata.raw
        if adata.raw is not None:
            ct_adata = adata.raw.to_adata()
            # Normalize if raw counts (not yet log-normalised)
            import scanpy as sc_inner
            sc_inner.pp.normalize_total(ct_adata, target_sum=1e4)
            sc_inner.pp.log1p(ct_adata)
        else:
            ct_adata = adata

        # Download model if not cached
        model = models.Model.load(model=model_name)

        predictions = celltypist.annotate(
            ct_adata,
            model=model,
            majority_voting=True,
            over_clustering=adata.obs['leiden']
        )

        adata.obs['cell_type_celltypist'] = predictions.predicted_labels['majority_voting'].values
        adata.obs['cell_type_celltypist_conf'] = predictions.probability_matrix.max(axis=1).values

        n_types = adata.obs['cell_type_celltypist'].nunique()
        logger.info(f"CellTypist annotation complete: {n_types} cell types identified")

        # Log comparison with marker-based annotation
        if 'cell_type' in adata.obs.columns:
            match = (adata.obs['cell_type'].astype(str) ==
                     adata.obs['cell_type_celltypist'].astype(str)).mean()
            logger.info(f"Agreement with marker scoring: {match*100:.1f}%")

    except ImportError:
        logger.warning("CellTypist not installed. Skipping. Install with: pip install celltypist")
        adata.obs['cell_type_celltypist'] = 'Not run'
        adata.obs['cell_type_celltypist_conf'] = 0.0
    except Exception as e:
        logger.warning(f"CellTypist failed: {e}. Continuing without it.")
        adata.obs['cell_type_celltypist'] = 'Failed'
        adata.obs['cell_type_celltypist_conf'] = 0.0

    return adata


def plot_cell_type_umap(adata: sc.AnnData, analysis_dir: str) -> None:
    """Plot UMAPs for marker-based and (if available) CellTypist annotations."""
    sc.settings.figdir = analysis_dir
    sc.settings.autoshow = False

    if 'cell_type' in adata.obs.columns:
        sc.pl.umap(adata, color=['cell_type'], legend_loc='on data', save='_cell_type.png')
        plt.close()

    if ('cell_type_celltypist' in adata.obs.columns and
            adata.obs['cell_type_celltypist'].iloc[0] not in ('Not run', 'Failed')):
        sc.pl.umap(adata, color=['cell_type_celltypist'],
                   legend_loc='on data', save='_cell_type_celltypist.png')
        plt.close()
        sc.pl.umap(adata, color=['cell_type_celltypist_conf'],
                   save='_celltypist_confidence.png')
        plt.close()

    logger.info(f"Cell-type UMAPs saved to: {analysis_dir}")


def generate_cellranger_comparison(
    cellranger_metrics_csv: Optional[str],
    filter_stats: Dict
) -> pd.DataFrame:
    """Parse Cell Ranger metrics_summary.csv and compare against 10x official expectations.

    Returns a DataFrame with metric / our_result / expected_range / status columns
    that is saved to tables/cellranger_vs_10x.csv.
    """
    # Official 10x expected ranges for PBMC 1k v3
    EXPECTED = {
        'Estimated Number of Cells':         ('900–1,300',   900,  1300),
        'Median Genes per Cell':              ('3,000–3,500', 3000, 3500),
        'Mean Reads per Cell':                ('50,000–60,000', 50000, 60000),
        'Valid Barcodes':                     ('≥95%',        95,   100),
        'Sequencing Saturation':              ('65–75%',      65,   75),
        'Reads Mapped Confidently to Transcriptome': ('75–85%', 75, 85),
        'Total Genes Detected':               ('25,000–27,000', 25000, 27000),
        'Median UMI Counts per Cell':         ('9,000–11,000', 9000, 11000),
    }

    rows = []

    if cellranger_metrics_csv and os.path.isfile(cellranger_metrics_csv):
        try:
            cr = pd.read_csv(cellranger_metrics_csv)
            # metrics_summary.csv has one row; each column is a metric
            for metric, (rng, lo, hi) in EXPECTED.items():
                if metric in cr.columns:
                    raw = str(cr[metric].iloc[0]).replace(',', '').replace('%', '').strip()
                    try:
                        val = float(raw)
                    except ValueError:
                        val = float('nan')
                    ok = lo <= val <= hi if not np.isnan(val) else False
                    rows.append({
                        'metric': metric,
                        'our_result': cr[metric].iloc[0],
                        'expected_range': rng,
                        'status': '✅ Pass' if ok else '⚠️  Check'
                    })
        except Exception as e:
            logger.warning(f"Could not parse Cell Ranger metrics: {e}")

    # Downstream filtering stats
    rows.append({
        'metric': 'Cells After QC Filtering',
        'our_result': filter_stats.get('final_cells', 'N/A'),
        'expected_range': '≥900',
        'status': '✅ Pass' if filter_stats.get('final_cells', 0) >= 900 else '⚠️  Check'
    })
    rows.append({
        'metric': 'Cell Retention Rate',
        'our_result': f"{filter_stats.get('cells_passed_pct', 0):.1f}%",
        'expected_range': '≥70%',
        'status': '✅ Pass' if filter_stats.get('cells_passed_pct', 0) >= 70 else '⚠️  Check'
    })

    df = pd.DataFrame(rows)
    logger.info("Cell Ranger vs 10x comparison table generated")
    return df


def save_results(
    adata: sc.AnnData,
    dirs: Dict[str, str],
    filter_stats: Optional[Dict] = None,
    config: Optional[Dict] = None,
    annotation_df: Optional[pd.DataFrame] = None,
    cellranger_metrics_csv: Optional[str] = None
) -> None:
    """Save all results to disk in standardized project structure."""
    # Save AnnData files to data/ directory
    adata.write(os.path.join(dirs['data'], 'processed.h5ad'))

    logger.info(f"Saved processed data: {os.path.join(dirs['data'], 'processed.h5ad')}")

    # Save marker genes
    de_df = sc.get.rank_genes_groups_df(adata, None)
    de_df.to_csv(os.path.join(dirs['tables'], 'marker_genes.csv'), index=False)

    # Save cluster counts
    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    cluster_counts.to_csv(os.path.join(dirs['tables'], 'cluster_counts.csv'))

    # Save cluster composition with cell counts and percentages
    cluster_summary = pd.DataFrame({
        'cluster': cluster_counts.index,
        'cell_count': cluster_counts.values,
        'percentage': (cluster_counts.values / cluster_counts.sum() * 100).round(2)
    })
    cluster_summary.to_csv(os.path.join(dirs['tables'], 'cluster_summary.csv'), index=False)

    if annotation_df is not None and not annotation_df.empty:
        annotation_df.to_csv(os.path.join(dirs['tables'], 'cluster_annotation.csv'), index=False)

    # Save QC stats
    if filter_stats:
        qc_stats_df = pd.DataFrame([filter_stats])
        qc_stats_df.to_csv(os.path.join(dirs['tables'], 'filtering_stats.csv'), index=False)

    # Save Cell Ranger vs 10x comparison table
    if filter_stats:
        comparison_df = generate_cellranger_comparison(cellranger_metrics_csv, filter_stats)
        comparison_df.to_csv(os.path.join(dirs['tables'], 'cellranger_vs_10x.csv'), index=False)
        logger.info(f"Saved Cell Ranger comparison: {os.path.join(dirs['tables'], 'cellranger_vs_10x.csv')}")

    # Save config
    if config:
        config_df = pd.DataFrame([config]).T
        config_df.columns = ['value']
        config_df.to_csv(os.path.join(dirs['tables'], 'pipeline_config.csv'))

        # Also save as JSON in logs for readability
        import json
        with open(os.path.join(dirs['logs'], 'pipeline_config.json'), 'w') as f:
            json.dump(config, f, indent=2, default=str)

    logger.info(f"All results saved to project directory: {dirs['root']}")


def run_pipeline(
    tenx_dir: Optional[str] = None,
    input_h5ad: Optional[str] = None,
    project_name: str = 'unnamed_project',
    base_output_dir: str = DEFAULT_OUTPUT_ROOT,
    config: Optional[Dict] = None,
    tenx_var_names: str = 'gene_symbols',
    cellranger_metrics_csv: Optional[str] = None
) -> None:
    """Run the complete pipeline with project-based output organization.

    Args:
        tenx_dir: Path to Cell Ranger output (filtered_feature_bc_matrix)
        input_h5ad: Path to existing h5ad file (alternative to tenx_dir)
        project_name: Name for this analysis project (creates output/<project_name>/)
        base_output_dir: Base directory for all outputs (default: ../output)
        config: Pipeline configuration dictionary
        tenx_var_names: Gene identifier type ('gene_symbols' or 'gene_ids')
        cellranger_metrics_csv: Optional path to Cell Ranger metrics_summary.csv
                                 for automated 10x comparison table
    """
    start_time = datetime.now()
    logger.info("=" * 70)
    logger.info(f"Starting scRNA-seq Analysis Pipeline")
    project_name = build_timestamped_project_name(project_name)
    logger.info(f"Project: {project_name}")
    logger.info(f"Started at: {start_time}")
    logger.info("=" * 70)

    # Use default config if not provided
    if config is None:
        config = DEFAULT_CONFIG.copy()

    # Setup output directories
    dirs = setup_output_directories(project_name, base_output_dir)
    
    # Setup file logging
    log_file = os.path.join(dirs['logs'], 'pipeline.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    try:
        logger.info(f"Logging to: {log_file}")

        # Load data (with project-specific cache to ensure isolation)
        adata = load_data(tenx_dir, input_h5ad, tenx_var_names, cache_dir=dirs['cache'])

        # QC
        adata = calculate_qc_metrics(adata)
        
        # Suggest adaptive QC thresholds based on data distribution
        suggest_qc_thresholds(adata, config)
        
        # Run doublet detection before filtering
        adata = run_doublet_detection(adata, config)
        
        # Plot QC metrics with current config
        plot_qc_metrics(adata, dirs['qc'], config)
        
        # Apply QC filtering
        adata, filter_stats = apply_qc_filtering(adata, config)

        # Save intermediate after QC
        adata.write(os.path.join(dirs['data'], 'after_qc.h5ad'))
        logger.info(f"Saved post-QC data: {os.path.join(dirs['data'], 'after_qc.h5ad')}")

        # Normalization and dimred
        adata = normalize_and_scale(adata, config)
        adata = run_dimred_and_clustering(adata, config)

        # Plot dimred
        plot_dimred_results(adata, dirs['analysis'])

        # Marker genes
        adata = find_marker_genes(adata)
        plot_marker_genes(adata, dirs['markers'], config['pbmc_markers'])
        adata, annotation_df = annotate_cell_types(adata, config['pbmc_markers'])

        # Fix 4: CellTypist reference-based annotation
        adata = run_celltypist(adata, config)

        # UMAPs coloured by marker-based + CellTypist annotations
        plot_cell_type_umap(adata, dirs['analysis'])

        # Auto-detect Cell Ranger metrics_summary.csv when tenx_dir is provided
        if cellranger_metrics_csv is None and tenx_dir:
            candidate = os.path.join(os.path.dirname(tenx_dir), 'metrics_summary.csv')
            if os.path.isfile(candidate):
                cellranger_metrics_csv = candidate
                logger.info(f"Auto-detected Cell Ranger metrics: {candidate}")

        # Save final results
        save_results(adata, dirs, filter_stats, config,
                     annotation_df=annotation_df,
                     cellranger_metrics_csv=cellranger_metrics_csv)

        # Completion summary
        end_time = datetime.now()
        duration = end_time - start_time
        
        logger.info("=" * 70)
        logger.info("PIPELINE COMPLETED SUCCESSFULLY")
        logger.info(f"Duration: {duration}")
        logger.info(f"Project directory: {dirs['root']}")
        logger.info(f"  - Data files: {dirs['data']}")
        logger.info(f"  - Figures: {dirs['figures']}")
        logger.info(f"  - Tables: {dirs['tables']}")
        logger.info(f"  - Logs: {dirs['logs']}")
        logger.info("=" * 70)
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}", exc_info=True)
        raise
    finally:
        # Remove file handler
        logger.removeHandler(file_handler)
        file_handler.close()


def main():
    parser = argparse.ArgumentParser(
        description='Optimized Scanpy scRNA-seq Pipeline - Project-based Output Organization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic run with project name (RECOMMENDED)
  python scanpy_pipeline_optimized.py --tenx_dir=../pbmc_1k_v3_count/outs/filtered_feature_bc_matrix --project=pbmc1k_default

  # From saved h5ad
  python scanpy_pipeline_optimized.py --input=previous_project/data/after_qc.h5ad --project=pbmc1k_reanalysis

  # Custom QC thresholds (permissive)
  python scanpy_pipeline_optimized.py --tenx_dir=... --project=pbmc1k_permissive --min_genes=100 --max_mito_pct=10
  
  # Skip doublet detection for faster execution
  python scanpy_pipeline_optimized.py --tenx_dir=... --project=pbmc1k_fast --skip_doublets
  
  # Stringent filtering with higher resolution clustering
  python scanpy_pipeline_optimized.py --tenx_dir=... --project=pbmc1k_stringent --min_genes=500 --max_mito_pct=3 --leiden_resolution=1.0

Output Structure:
  output/
  └── <project_name>/
      ├── data/           # Processed .h5ad files
      ├── figures/        # All visualizations
      │   ├── qc/         # QC plots
      │   ├── analysis/   # PCA/UMAP plots  
      │   └── markers/    # Marker gene plots
      ├── tables/         # CSV results
      └── logs/           # Pipeline logs and configs
        """
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--tenx_dir', help='Path to 10x filtered_feature_bc_matrix directory')
    input_group.add_argument('--input', help='Path to input .h5ad file')

    parser.add_argument('--tenx_var_names',
                       default='gene_symbols',
                       choices=['gene_symbols', 'gene_ids'],
                       help='Which column to use as gene names when reading 10x MTX')

    parser.add_argument('--project',
                       required=True,
                       help='Base project name for output organization; timestamp is appended automatically')
    parser.add_argument('--base_output_dir',
                       default=DEFAULT_OUTPUT_ROOT,
                       help=f'Base output directory (default: {DEFAULT_OUTPUT_ROOT})')

    # QC parameters
    qc_group = parser.add_argument_group('QC Thresholds')
    qc_group.add_argument('--min_genes', type=int,
                         help='Minimum genes per cell (default: 200)')
    qc_group.add_argument('--max_genes', type=int,
                         help='Maximum genes per cell (default: 6000)')
    qc_group.add_argument('--max_mito_pct', type=float,
                         help='Maximum mitochondrial percentage (default: 10.0)')
    qc_group.add_argument('--max_ribo_pct', type=float,
                         help='Maximum ribosomal percentage (default: 50.0)')
    qc_group.add_argument('--min_cells', type=int,
                         help='Minimum cells for gene filtering (default: 3)')

    # Doublet detection parameters
    doublet_group = parser.add_argument_group('Doublet Detection (Scrublet)')
    doublet_group.add_argument('--skip_doublets', action='store_true',
                              help='Skip doublet detection')
    doublet_group.add_argument('--doublet_threshold', type=float,
                              help='Doublet score threshold (default: 0.25)')
    doublet_group.add_argument('--allow_doublet_fail', action='store_true',
                              help='If set, continue with doublet_score=0 when scrublet fails')

    # Analysis parameters
    analysis_group = parser.add_argument_group('Analysis Parameters')
    analysis_group.add_argument('--n_top_genes', type=int,
                               help='Number of highly variable genes (default: 2000)')
    analysis_group.add_argument('--n_pcs', type=int, default=DEFAULT_CONFIG['n_pcs'],
                               help=f'Number of PCs to use (default: {DEFAULT_CONFIG["n_pcs"]})')
    analysis_group.add_argument('--n_neighbors', type=int, default=None,
                               help='Number of neighbors (default: None for adaptive)')
    analysis_group.add_argument('--leiden_resolution', type=float, default=DEFAULT_CONFIG['leiden_resolution'],
                               help=f'Leiden clustering resolution (default: {DEFAULT_CONFIG["leiden_resolution"]})')
    analysis_group.add_argument('--random_state', type=int,
                               help='Random seed for PCA/UMAP/Leiden (default: 0)')
    analysis_group.add_argument('--no_auto_n_pcs', action='store_true',
                               help='Disable elbow-guided PC selection; use --n_pcs directly')
    analysis_group.add_argument('--no_celltypist', action='store_true',
                               help='Skip CellTypist reference-based annotation')
    analysis_group.add_argument('--celltypist_model', type=str,
                               help='CellTypist model name (default: Immune_All_Low.pkl)')
    parser.add_argument('--cellranger_metrics_csv',
                       default=None,
                       help='Path to Cell Ranger metrics_summary.csv for automated 10x comparison '
                            '(auto-detected from --tenx_dir/../metrics_summary.csv if not set)')

    args = parser.parse_args()

    # Build config from defaults and CLI
    config = DEFAULT_CONFIG.copy()
    for key in config:
        if hasattr(args, key) and getattr(args, key) is not None:
            config[key] = getattr(args, key)
    
    # Handle doublet detection flag
    if args.skip_doublets:
        config['run_doublet_detection'] = False
    if args.doublet_threshold:
        config['doublet_threshold'] = args.doublet_threshold
    if args.allow_doublet_fail:
        config['fail_on_doublet_error'] = False

    # Handle Fix 2 / Fix 4 flags
    if args.no_auto_n_pcs:
        config['auto_n_pcs'] = False
    if args.no_celltypist:
        config['run_celltypist'] = False
    if args.celltypist_model:
        config['celltypist_model'] = args.celltypist_model

    # Run pipeline
    run_pipeline(
        tenx_dir=args.tenx_dir,
        input_h5ad=args.input,
        project_name=args.project,
        base_output_dir=args.base_output_dir,
        config=config,
        tenx_var_names=args.tenx_var_names,
        cellranger_metrics_csv=args.cellranger_metrics_csv
    )


if __name__ == "__main__":
    main()
