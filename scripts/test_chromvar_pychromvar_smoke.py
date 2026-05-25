"""Smoke test for the chromVAR pychromvar (canonical) path.

Verifies that the full-mode `ChromVARModule` (Schep 2017 / pychromvar with
GC-matched background) is invokable given:
  - GRCh38 genome FASTA (already present in repo)
  - JASPAR2024 CORE non-redundant motif PFMs (downloaded 2026-05-22)
  - sc10x_methods conda env (with pychromvar 0.0.4, biopython 1.87)

Acceptance:
  - chromvar_status == "ok_pychromvar"
  - chromvar_metadata["mode"] == "pychromvar"
  - chromvar_method_actually_used == "pychromvar_full_gc_matched_background"
  - n_tfs >= 1 (full validation criterion in ralplan is n_tfs > 50 on real data;
    synthetic 50-peak input will produce fewer motif matches, so we relax to >=1
    for this smoke; the real-data threshold lives in the ralplan acceptance gate)

Run with the sc10x_methods env:
  /home/zerlinshen/conda/envs/sc10x_methods/bin/python \\
    scripts/test_chromvar_pychromvar_smoke.py

Exits non-zero on failure.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


SUITE = Path("/home/zerlinshen/Bioinformatics Research Pipeline")
GENOME = SUITE / "singlecell_factory/ref/reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
PWM = SUITE / "singlecell_factory/data/references/JASPAR2024/JASPAR2024_CORE_non-redundant_pfms_meme.txt"


def _build_synthetic_atac(n_cells: int = 80, n_peaks: int = 60) -> ad.AnnData:
    """Build a minimal AnnData with a sparse peak matrix and atac_var coordinates.

    Peaks anchored on real chr1 coordinates with 500bp width so pychromvar's
    add_peak_seq + add_gc_bias have something coherent to read from the FASTA.
    """
    rng = np.random.default_rng(7)
    X_rna = sp.csr_matrix(rng.poisson(0.5, size=(n_cells, n_peaks)).astype("float32"))
    peak_mat = sp.csr_matrix((rng.binomial(1, 0.15, size=(n_cells, n_peaks))).astype("float32"))
    # Anchor peaks at chr1:1_000_000 onward, spaced 10 kb apart, 500 bp wide.
    coords = []
    start = 1_000_000
    for i in range(n_peaks):
        s = start + i * 10_000
        coords.append(("chr1", s, s + 500, f"chr1:{s}-{s + 500}"))
    var_atac = pd.DataFrame(
        coords, columns=["chrom", "start", "end", "peak_id"]
    )
    adata = ad.AnnData(
        X=X_rna,
        obs=pd.DataFrame(index=[f"cell{i:03d}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"gene{i:03d}" for i in range(n_peaks)]),
    )
    adata.obsm["atac_peaks"] = peak_mat
    adata.uns["atac_peaks"] = {"present": True}
    adata.uns["atac_var"] = var_atac
    return adata


class _DummyCfg:
    chromvar_genome_fasta = str(GENOME)
    chromvar_motif_pwm = str(PWM)
    tf_peak_annotation_path = None
    chromvar_min_peaks_per_tf = 3


class _DummyCtx:
    """Minimal PipelineContext shim for the chromvar module."""

    def __init__(self, run_dir: Path) -> None:
        self.adata = None
        self.cfg = _DummyCfg()
        self.run_dir = run_dir
        self.metadata: dict = {}
        self._status: list = []

    def status(self, mod: str, state: str, msg: str) -> None:
        self._status.append((mod, state, msg))
        print(f"  [status] {mod}={state}: {msg}")


def main() -> int:
    if not GENOME.exists():
        print(f"FAIL: GRCh38 FASTA missing at {GENOME}", file=sys.stderr)
        return 2
    if not PWM.exists():
        print(f"FAIL: JASPAR2024 PWM missing at {PWM}", file=sys.stderr)
        return 2

    sys.path.insert(0, str(SUITE / "singlecell_factory"))
    from workflow.modular.modules.chromvar import ChromVARModule

    with tempfile.TemporaryDirectory() as td:
        run_dir = Path(td)
        ctx = _DummyCtx(run_dir=run_dir)
        ctx.adata = _build_synthetic_atac()
        print("synthetic ATAC adata:", ctx.adata.shape, "peaks=", ctx.adata.obsm["atac_peaks"].shape)
        try:
            ChromVARModule().run(ctx)
        except Exception as exc:
            print(f"FAIL: ChromVARModule raised: {exc!r}", file=sys.stderr)
            import traceback; traceback.print_exc()
            return 3

        status = ctx.metadata.get("chromvar_status")
        mode = ctx.metadata.get("chromvar_mode")
        method = ctx.metadata.get("chromvar_method_actually_used")
        n_tfs = ctx.metadata.get("chromvar_n_tfs")
        print(json.dumps({
            "chromvar_status": status,
            "chromvar_mode": mode,
            "chromvar_method_actually_used": method,
            "chromvar_n_tfs": n_tfs,
        }, indent=2))

        if mode != "pychromvar":
            print(f"FAIL: expected mode=='pychromvar', got {mode}", file=sys.stderr)
            return 4
        if status != "ok_pychromvar":
            print(f"FAIL: expected status=='ok_pychromvar', got {status}", file=sys.stderr)
            return 5
        if method != "pychromvar_full_gc_matched_background":
            print(f"FAIL: expected chromvar_method_actually_used=='pychromvar_full_gc_matched_background', got {method}", file=sys.stderr)
            return 6
        if not isinstance(n_tfs, int) or n_tfs < 1:
            print(f"FAIL: expected n_tfs >= 1, got {n_tfs}", file=sys.stderr)
            return 7

    print("ALL CHROMVAR PYCHROMVAR SMOKE CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
