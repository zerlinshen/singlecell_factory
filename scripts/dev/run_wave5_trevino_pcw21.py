"""Wave-5 Trevino PCW21 driver — orchestrates the multiomics sub-DAG on the
cached AnnData produced by ``load_trevino_2021``.

The factory CLI (``workflow.modular.cli``) assumes a 10x cellranger output as
``--sample-root``, which is incompatible with Trevino 2021's published TSVs.
This driver loads the pre-built h5ad (RNA in ``.X``, ATAC in
``.obsm['atac_peaks']``), performs RNA preprocessing to produce
``obsm['X_pca']`` (required by ``multimodal_integration``), then runs the
Wave-5 sub-DAG modules directly:

    atac_lsi -> multimodal_integration(WNN, X_lsi) -> clustering
              -> annotation -> peak_to_gene -> trajectory

Outputs land under ``--project-root/runs/<run-id>/`` matching the factory
contract: ``final_adata.h5ad``, ``run_manifest.json``, ``module_status.csv``.

Wave-5 spec: ``.omc/specs/deep-interview-wave5-completion.md``
Wave-5 plan: ``.omc/plans/wave5-completion-consensus-2026-05-16.md``
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.modular.config import (
    CellRangerConfig,
    PipelineConfig,
)
from workflow.modular.context import PipelineContext
from workflow.modular.modules.atac_lsi import ATACLSIModule
from workflow.modular.modules.multimodal_integration import (
    MultimodalIntegrationConfig,
    MultimodalIntegrationModule,
)
from workflow.modular.modules.clustering import ClusteringModule
from workflow.modular.modules.peak_to_gene import PeakToGeneModule
from workflow.modular.modules.trajectory import TrajectoryModule


RUN_ID_RE = re.compile(r"^[0-9]{8}T[0-9]{4}Z-[0-9a-f]{7}$")


def _git_sha(repo: Path) -> str:
    head = repo / ".git" / "HEAD"
    if not head.exists():
        return "unknown"
    try:
        ref = head.read_text().strip()
        if ref.startswith("ref: "):
            ref_path = repo / ".git" / ref[5:]
            return ref_path.read_text().strip()[:12]
        return ref[:12]
    except OSError:
        return "unknown"


def _git_output(repo: Path, args: list[str]) -> str:
    try:
        return subprocess.run(
            ["git", *args],
            cwd=repo,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return ""


def _dirty_tree_metadata(repo: Path, *, allow_dirty: bool) -> dict[str, object]:
    """Return dirty-tree provenance and enforce the factory dirty gate."""
    status = _git_output(repo, ["status", "--short"])
    diff = _git_output(repo, ["diff", "--binary"])
    dirty = bool(status.strip())
    if dirty and not allow_dirty:
        raise SystemExit(
            "Factory git tree is dirty. Re-run with --allow-dirty to record "
            "diff/status hashes in run_manifest.json."
        )
    return {
        "factory_tree_dirty": dirty,
        "allow_dirty": allow_dirty,
        "diff_sha256": hashlib.sha256(diff.encode("utf-8")).hexdigest() if dirty else None,
        "git_status_short_sha256": (
            hashlib.sha256(status.encode("utf-8")).hexdigest() if dirty else None
        ),
        "git_status_short_count": len([line for line in status.splitlines() if line.strip()]),
    }


def _run_atac_lsi(ctx: PipelineContext) -> None:
    mod = ATACLSIModule()
    ctx.set_module_dir(mod.name)
    mod.run(ctx)


def _run_multimodal(ctx: PipelineContext, second_obsm_key: str = "X_lsi") -> None:
    cfg = MultimodalIntegrationConfig(
        engine="wnn",
        second_obsm_key=second_obsm_key,
        _explicitly_set=True,
    )
    mod = MultimodalIntegrationModule(config=cfg)
    ctx.set_module_dir(mod.name)
    mod.run(ctx)
    if cfg.wnn_obsm_output_key not in ctx.adata.obsm:
        status = ctx.adata.uns.get("multimodal_status", {})
        raise RuntimeError(
            "multimodal_integration completed without "
            f"{cfg.wnn_obsm_output_key!r}; status={status}. "
            "Set RSCRIPT_BIN to an R runtime with Seurat+arrow or fix the WNN driver."
        )


def _run_clustering_on_wnn(ctx: PipelineContext) -> None:
    """Run Leiden directly on the joint WNN embedding (X_wnn).

    ClusteringModule from the factory expects to operate on PCA + batch
    integration; here we want joint Leiden on X_wnn for US-W5-5 ARI. We
    compute the neighbor graph on X_wnn and run Leiden inline so the
    cluster labels reflect the multi-omic state. Also stores the same
    graph under the *default* key so downstream ``sc.tl.paga`` /
    ``sc.tl.dpt`` (which read the default neighbors slot) keep operating
    on the joint WNN manifold.
    """
    ctx.set_module_dir("clustering")
    if "X_wnn" not in ctx.adata.obsm:
        raise RuntimeError("X_wnn missing — multimodal_integration must run first")
    # Default key — used by trajectory's paga/diffmap/dpt
    sc.pp.neighbors(
        ctx.adata,
        use_rep="X_wnn",
        n_neighbors=ctx.cfg.clustering.n_neighbors,
        random_state=ctx.random_state,
    )
    sc.tl.leiden(
        ctx.adata,
        resolution=ctx.cfg.clustering.leiden_resolution,
        random_state=ctx.random_state,
        key_added="leiden",
    )
    sc.tl.umap(ctx.adata, random_state=ctx.random_state)
    ctx.adata.obsm["X_umap"] = ctx.adata.obsm["X_umap"].astype(np.float32)
    ctx.status("clustering", True, "Leiden on X_wnn (default neighbors key)")


def _run_peak_to_gene(ctx: PipelineContext) -> None:
    """Compute peak-to-gene links and top-1000 by |Pearson|.

    The factory's ``compute_peak_gene_linkages`` loops 100 sparse-row
    permutations per (peak, gene) pair, which is dominated by the
    ``A_col[perm_indices, :]`` reindex on a CSR matrix; on 467k peaks that
    is >1h of single-threaded Python. Wave-5 AC-VAL-3 ranks by
    ``|Pearson|`` (observation, not permutation p-value), so this driver
    computes observed Pearson once per (peak, gene) candidate pair and
    writes ``uns["peak_to_gene_top1000"]``. The ledger records the
    deviation: ``peak_to_gene_n_perms=0`` (no permutation null in this
    run) — overlap metric is unaffected; ``p_perm`` / ``fdr`` are not
    populated and any consumer that needs them must re-run with the
    upstream factory module.
    """
    import scipy.sparse
    import time

    ctx.set_module_dir("peak_to_gene")
    out_dir = ctx.run_dir / "peak_to_gene"

    adata = ctx.adata
    atac_peaks_meta = adata.uns["atac_peaks"]
    n_peaks = atac_peaks_meta.shape[0]
    if atac_peaks_meta.shape[0] != adata.obsm["atac_peaks"].shape[1]:
        raise RuntimeError(
            f"atac_peaks meta rows ({n_peaks}) != obsm atac_peaks cols ({adata.obsm['atac_peaks'].shape[1]})"
        )

    # 1. Nearest-TSS pairing within ±500 kb (same semantics as factory PeakToGeneModule).
    window = 500_000
    chrom_col = "chrom" if "chrom" in atac_peaks_meta.columns else "seqnames"
    gene_chrom = adata.var["chrom"].astype(str).values
    gene_tss = adata.var["tss"].astype(np.int64).values
    # group gene indices by chrom for fast windowed lookup
    from collections import defaultdict
    gene_by_chrom: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom in np.unique(gene_chrom):
        mask = gene_chrom == chrom
        idxs = np.where(mask)[0]
        order = np.argsort(gene_tss[idxs])
        gene_by_chrom[chrom] = (idxs[order], gene_tss[idxs][order])

    candidates: list[tuple[int, int]] = []  # (peak_idx, gene_idx)
    pmeta = atac_peaks_meta
    pmid = ((pmeta["start"].to_numpy(dtype=np.int64) + pmeta["end"].to_numpy(dtype=np.int64)) // 2)
    pchrom = pmeta[chrom_col].astype(str).to_numpy()
    t0 = time.time()
    for i in range(n_peaks):
        c = pchrom[i]
        if c not in gene_by_chrom:
            continue
        gidxs, gtss = gene_by_chrom[c]
        m = pmid[i]
        lo = np.searchsorted(gtss, m - window)
        hi = np.searchsorted(gtss, m + window)
        for j in gidxs[lo:hi]:
            candidates.append((i, int(j)))
    print(f"  nearest-TSS pairing done in {time.time()-t0:.1f}s; {len(candidates)} (peak,gene) candidates", flush=True)

    # 2. Compute observed Pearson per (peak, gene) candidate using sparse mean+stdev.
    # Use CSC for fast O(nnz_col) column slicing (CSR column slicing is O(n_cells)
    # and dominates wall-time for 462k peak accesses).
    A = adata.obsm["atac_peaks"].tocsc()
    X = adata.X.tocsc()
    n_cells = adata.n_obs

    # Precompute column means and stdevs (cell-axis) for both matrices.
    def _col_stats(M):
        s1 = np.asarray(M.sum(axis=0)).ravel()
        M2 = M.copy(); M2.data = M2.data ** 2
        s2 = np.asarray(M2.sum(axis=0)).ravel()
        mean = s1 / n_cells
        var = s2 / n_cells - mean ** 2
        var[var < 0] = 0
        return mean, np.sqrt(var)

    pA_mean, pA_std = _col_stats(A)
    pX_mean, pX_std = _col_stats(X)

    cand = np.asarray(candidates, dtype=np.int64)  # (n_candidates, 2)
    # Compute cross-moments E[XY] for each candidate pair via grouped dot products.
    # Strategy: per-peak block, take all candidate gene indices for that peak, do A_col.T @ X[:, gene_idxs].
    t1 = time.time()
    pearsons = np.zeros(len(cand), dtype=np.float32)
    write_pos = 0
    # group candidates by peak index
    cur_peak = -1
    block_start = 0
    n_peaks_with_any = int(np.unique(cand[:, 0]).size)
    print(f"  computing observed Pearson for {len(cand)} pairs across {n_peaks_with_any} peaks", flush=True)
    for k in range(len(cand) + 1):
        if k == len(cand) or cand[k, 0] != cur_peak:
            if cur_peak >= 0 and block_start < k:
                gene_idxs = cand[block_start:k, 1]
                A_col = A[:, cur_peak]                    # (n_cells, 1) sparse
                X_block = X[:, gene_idxs]                 # (n_cells, n_block) sparse
                # E[AX] = A_col.T @ X_block, shape (1, n_block)
                EAX = np.asarray((A_col.T @ X_block).todense()).ravel() / n_cells
                cov = EAX - pA_mean[cur_peak] * pX_mean[gene_idxs]
                denom = pA_std[cur_peak] * pX_std[gene_idxs]
                with np.errstate(invalid="ignore", divide="ignore"):
                    r = np.where(denom > 0, cov / denom, 0.0)
                pearsons[block_start:k] = r.astype(np.float32)
            if k < len(cand):
                cur_peak = int(cand[k, 0])
                block_start = k
        # else: same peak, continue
    print(f"  Pearson pass done in {time.time()-t1:.1f}s", flush=True)

    # 3. Build top-1000 by |Pearson|.
    abs_r = np.abs(pearsons)
    order = np.argsort(-abs_r)
    top_k = min(1000, len(order))
    top_idx = order[:top_k]
    var_index = list(adata.var.index)
    top_df = pd.DataFrame({
        "peak_idx": cand[top_idx, 0],
        "gene": [var_index[int(g)] for g in cand[top_idx, 1]],
        "pearson": pearsons[top_idx],
    })
    adata.uns["peak_to_gene_top1000"] = top_df
    adata.uns["peak_to_gene_linkages"] = top_df  # no FDR pass in this run; use top1000 as best-available linkages table

    # 4. Also produce the nearest-TSS link table that the factory module writes.
    out_dir.mkdir(parents=True, exist_ok=True)
    chrom_col_meta = chrom_col
    # Use a memory-efficient writer: one row per nearest-TSS-linked peak.
    nearest_records = []
    for i, (mids, _) in gene_by_chrom.items():
        pass  # unused placeholder
    # Re-derive nearest gene per peak for the CSV (cheap).
    for i in range(n_peaks):
        c = pchrom[i]
        if c not in gene_by_chrom:
            continue
        gidxs, gtss = gene_by_chrom[c]
        m = pmid[i]
        pos = np.searchsorted(gtss, m)
        best = None
        for p in (pos - 1, pos):
            if 0 <= p < len(gtss):
                d = abs(int(gtss[p]) - int(m))
                if d <= window and (best is None or d < best[0]):
                    best = (d, int(gidxs[p]))
        if best is not None:
            nearest_records.append({
                "peak_idx": i,
                "peak_chrom": c,
                "peak_start": int(pmeta["start"].iloc[i]),
                "peak_end":   int(pmeta["end"].iloc[i]),
                "gene": var_index[best[1]],
                "dist_bp": best[0],
            })
    link_df = pd.DataFrame(nearest_records)
    adata.uns["peak_to_gene"] = link_df
    link_df.to_csv(out_dir / "peak_to_gene.csv", index=False)
    (out_dir / "peak_to_gene_summary.json").write_text(
        json.dumps({
            "n_peaks_total":  int(n_peaks),
            "n_peaks_linked": int(len(link_df)),
            "n_pearson_pairs": int(len(cand)),
            "top1000_n": int(top_k),
            "window_bp":      window,
            "n_perms":        0,
            "deviation":      "compute_peak_gene_linkages permutation phase bypassed for wallclock; top1000 ranked by |Pearson| of observation only.",
        }, indent=2),
        encoding="utf-8",
    )
    top_df.to_csv(out_dir / "peak_to_gene_top1000.csv", index=False)
    ctx.metadata["peak_to_gene_n_perms"] = 0
    ctx.metadata["peak_to_gene_n_linked"] = int(len(link_df))
    ctx.status("peak_to_gene", True, f"top1000 by |Pearson| (no perms; ranked from {len(cand)} ±{window}bp pairs)")


def _run_trajectory(ctx: PipelineContext) -> None:
    mod = TrajectoryModule()
    ctx.set_module_dir(mod.name)
    mod.run(ctx)


def _preprocess_rna(adata: anndata.AnnData, n_pcs: int = 30, n_top_genes: int = 3000) -> None:
    """Minimal RNA pp to produce ``obsm['X_pca']`` (required by
    ``multimodal_integration``)."""
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat_v3", layer="counts")
    sc.tl.pca(adata, n_comps=n_pcs, random_state=42, mask_var="highly_variable")


def _build_ctx(project_root: Path, project_name: str, run_id: str | None = None) -> PipelineContext:
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%MZ")
    sha = _git_sha(Path(__file__).resolve().parents[2])
    run_id = run_id or f"{ts}-{sha[:7]}"
    run_dir = project_root / "runs" / run_id / "python"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Dummy cellranger config — required by PipelineConfig but unused in this driver.
    dummy_cr = CellRangerConfig(
        sample_root=Path("/dev/null"),
        outs_dir=Path("/dev/null"),
        sample_id="trevino_pcw21",
    )
    cfg = PipelineConfig(
        project=project_name,
        output_dir=project_root / "runs",
        cellranger=dummy_cr,
        optional_modules=[
            "atac_lsi",
            "multimodal_integration",
            "clustering",
            "peak_to_gene",
            "trajectory",
        ],
        random_state=42,
        gpu_mode="auto",
        checkpoint=True,
    )
    # PeakToGeneModule reads these as attributes via getattr(); set after
    # construction so we don't need to extend PipelineConfig's signature.
    cfg.gene_tss_bed_path = Path(
        "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/data/external/gencode_grch38_tss_by_ensg.bed"
    )
    cfg.peak_to_gene_window = 500_000  # Trevino 2021 Methods: ±500 kb window
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    ctx.metadata["run_id"] = run_id
    ctx.metadata["driver"] = "scripts/dev/run_wave5_trevino_pcw21.py"
    ctx.metadata["repo_sha_at_manifest_write"] = sha
    return ctx


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument(
        "--run-id",
        default=None,
        help=(
            "Optional project run id to populate. Format: "
            "YYYYMMDDTHHMMZ-<7 hex chars>. Existing source-truth outputs are not overwritten."
        ),
    )
    ap.add_argument(
        "--allow-dirty",
        action="store_true",
        help="Allow execution from a dirty factory tree and record diff/status hashes.",
    )
    ap.add_argument("--project", default="trevino_pcw21_v6_3")
    ap.add_argument(
        "--input-h5ad",
        required=True,
        type=Path,
        help="Cached AnnData from load_trevino_2021 (RNA in .X, ATAC in obsm['atac_peaks']).",
    )
    args = ap.parse_args()

    if args.run_id and not RUN_ID_RE.fullmatch(args.run_id):
        ap.error(f"--run-id must match {RUN_ID_RE.pattern}: {args.run_id}")

    repo_root = Path(__file__).resolve().parents[2]
    dirty_meta = _dirty_tree_metadata(repo_root, allow_dirty=args.allow_dirty)
    ctx = _build_ctx(args.project_root, args.project, run_id=args.run_id)
    ctx.metadata.update(dirty_meta)
    for source_truth in ("final_adata.h5ad", "run_manifest.json", "module_status.csv"):
        if (ctx.run_dir / source_truth).exists():
            raise SystemExit(
                f"Refusing to overwrite existing source-truth artifact: {ctx.run_dir / source_truth}"
            )
    log = lambda m: print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {m}", flush=True)

    log(f"run_id={ctx.metadata['run_id']}  run_dir={ctx.run_dir}")
    log(f"loading {args.input_h5ad}")
    ctx.adata = anndata.read_h5ad(args.input_h5ad)
    log(f"  n_obs={ctx.adata.n_obs} n_vars={ctx.adata.n_vars} atac_peaks={ctx.adata.obsm['atac_peaks'].shape}")

    # PeakToGeneModule reads adata.uns["atac_peaks"] expecting columns
    # {chrom, start, end}. Loader stores peak metadata as adata.uns["atac_var"]
    # with column names {seqnames, start, end, ...} matching the Trevino bed.
    # Reproject into the contract shape here so the module's contract is met.
    if "atac_peaks" not in ctx.adata.uns and "atac_var" in ctx.adata.uns:
        av = ctx.adata.uns["atac_var"]
        if "seqnames" in av.columns and "chrom" not in av.columns:
            av = av.rename(columns={"seqnames": "chrom"})
        ctx.adata.uns["atac_peaks"] = av
        log(f"  built adata.uns['atac_peaks'] (n_peaks={len(av)}) from loader's atac_var")

    # Inject `peak_id` (format "chr1:start-end") into atac_var so
    # compute_peak_gene_linkages can parse peak coordinates and use the
    # ±500 kb windowed candidate set instead of the high-variance fallback
    # (3,410 genes × 467k peaks × 100 perms is ~hours of compute).
    av = ctx.adata.uns["atac_var"]
    if "peak_id" not in av.columns:
        chrom_col = "chrom" if "chrom" in av.columns else "seqnames"
        av["peak_id"] = (
            av[chrom_col].astype(str)
            + ":"
            + av["start"].astype(int).astype(str)
            + "-"
            + av["end"].astype(int).astype(str)
        )
        ctx.adata.uns["atac_var"] = av
        log(f"  injected adata.uns['atac_var']['peak_id'] (e.g. {av['peak_id'].iloc[0]})")

    # Inject `chrom` + `tss` into adata.var so compute_peak_gene_linkages'
    # windowed mode activates (has_coords=True). var.index is bare ENSG;
    # match against the GENCODE TSS table built from the 10x GRCh38 GTF.
    tss_tsv = "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/data/external/gencode_grch38_tss_by_ensg.tsv"
    tss_df = pd.read_csv(tss_tsv, sep="\t").set_index("gene_id")
    chrom_map = tss_df["chrom"]
    tss_map = tss_df["tss"]
    # chrom is stored as int32 by the module after a chrom_map = {c: i for i, c in enumerate(...)}; use string here
    # and let the module hash chroms via the same enumerate path.
    ctx.adata.var["chrom"] = ctx.adata.var.index.map(chrom_map)
    ctx.adata.var["tss"] = ctx.adata.var.index.map(tss_map).astype("Int64")
    n_with_coords = int(ctx.adata.var["chrom"].notna().sum())
    log(f"  injected adata.var['chrom'/'tss']: {n_with_coords}/{ctx.adata.n_vars} genes have coords")
    # Drop genes without coords from candidate pool: set tss to NaN and chrom to None — the module's
    # `gene_chrom == peak_chrom` comparison will exclude them naturally only if we put them at a sentinel
    # chrom that never matches; cast to int32 chrom_hash later in module via enumerate(). Leaving NaN
    # would break the int32 cast. Replace NaN chrom with a sentinel that's chr-shaped but unused: "chrNA".
    ctx.adata.var.loc[ctx.adata.var["chrom"].isna(), "chrom"] = "chrNA"
    ctx.adata.var.loc[ctx.adata.var["tss"].isna(), "tss"] = -1
    ctx.adata.var["tss"] = ctx.adata.var["tss"].astype(np.int64)

    stage_starts: list[tuple[str, float]] = []

    def stage(name: str, fn):
        t0 = time.time()
        log(f">>> {name} START")
        try:
            fn(ctx)
            ok = True
            msg = "ok"
        except Exception as e:
            ok = False
            msg = f"{type(e).__name__}: {e}"
            log(f"!!! {name} FAILED: {msg}")
            traceback.print_exc()
        dt = time.time() - t0
        stage_starts.append((name, dt))
        log(f"<<< {name} END ({dt:.1f}s) status={'ok' if ok else 'FAILED'}")
        return ok

    if not stage("rna_pp + pca", lambda c: _preprocess_rna(c.adata)):
        return 2
    if not stage("atac_lsi", _run_atac_lsi):
        return 2
    if not stage("multimodal_integration", lambda c: _run_multimodal(c)):
        return 2
    if not stage("clustering(X_wnn)", _run_clustering_on_wnn):
        return 2
    if not stage("peak_to_gene", _run_peak_to_gene):
        return 2
    if not stage("trajectory", _run_trajectory):
        return 2

    # finalize
    final_path = ctx.run_dir / "final_adata.h5ad"
    log(f"writing {final_path}")
    ctx.adata.write_h5ad(final_path)
    manifest = {
        "run_id": ctx.metadata["run_id"],
        "project": args.project,
        "project_root": str(args.project_root),
        "run_dir": str(ctx.run_dir),
        "driver": ctx.metadata["driver"],
        "repo_sha_at_manifest_write": ctx.metadata.get("repo_sha_at_manifest_write"),
        "factory_tree_dirty": ctx.metadata.get("factory_tree_dirty"),
        "allow_dirty": ctx.metadata.get("allow_dirty"),
        "diff_sha256": ctx.metadata.get("diff_sha256"),
        "git_status_short_sha256": ctx.metadata.get("git_status_short_sha256"),
        "git_status_short_count": ctx.metadata.get("git_status_short_count"),
        "input_h5ad": str(args.input_h5ad),
        "modules_run": [name for name, _ in stage_starts],
        "stage_seconds": {name: round(dt, 2) for name, dt in stage_starts},
        "final_adata": str(final_path),
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "env_vars": {
            k: os.environ.get(k, "")
            for k in (
                "SC_REQUIRE_PROJECT_ROOT",
                "SC_MULTIMODAL_SECOND_OBSM",
                "SC_MEM_GUARD",
                "SC_CHECKPOINT_SCHEMA_VERSION",
            )
        },
    }
    (ctx.run_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    pd.DataFrame(ctx.module_status).to_csv(ctx.run_dir / "module_status.csv", index=False)
    ctx.flush_figures()
    log("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
