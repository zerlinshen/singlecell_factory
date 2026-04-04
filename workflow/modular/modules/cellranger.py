from __future__ import annotations

from pathlib import Path
import subprocess

import scanpy as sc

from ..context import PipelineContext


class CellRangerModule:
    """Mandatory module: validate/run Cell Ranger and load matrix."""

    name = "cellranger"
    required = True

    def run(self, ctx: PipelineContext) -> None:
        cfg = ctx.cfg.cellranger
        outs = cfg.outs_dir

        # If user requests force run or no previous Cell Ranger output exists, run `cellranger count`.
        if cfg.force_run or not outs.exists():
            if not cfg.run_if_missing:
                raise FileNotFoundError(f"Cell Ranger output not found: {outs}")
            self._run_cellranger(
                cfg.sample_root,
                cfg.fastq_dir,
                cfg.transcriptome_dir,
                cfg.sample_id,
                cfg.localcores,
                cfg.localmem,
            )

        if not outs.exists():
            raise FileNotFoundError(f"Cell Ranger output still missing after attempt: {outs}")

        adata = sc.read_10x_mtx(str(outs), var_names="gene_symbols", cache=False)
        adata.var_names_make_unique()
        ctx.adata = adata
        ctx.metadata["cellranger_outs"] = str(outs)
        ctx.metadata["raw_cells"] = int(adata.n_obs)
        ctx.metadata["raw_genes"] = int(adata.n_vars)

    @staticmethod
    def _run_cellranger(
        sample_root: Path,
        fastq_dir: Path | None,
        transcriptome_dir: Path | None,
        sample_id: str,
        localcores: int,
        localmem: int,
    ) -> None:
        if fastq_dir is None or transcriptome_dir is None:
            raise ValueError("Running Cell Ranger requires fastq_dir and transcriptome_dir.")
        cmd = [
            "cellranger",
            "count",
            f"--id={sample_id}",
            "--create-bam=true",
            f"--fastqs={fastq_dir}",
            f"--transcriptome={transcriptome_dir}",
            f"--localcores={localcores}",
            f"--localmem={localmem}",
        ]
        subprocess.run(cmd, check=True, cwd=sample_root)
