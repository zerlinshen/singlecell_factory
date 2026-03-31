from __future__ import annotations

from pathlib import Path
import shutil


def migrate_layout(root: Path) -> None:
    """Migrate raw inputs into `data/raw` and references into `ref`.

    Parameters
    ----------
    root
        Repository root path.
    """

    raw = root / "data" / "raw"
    ref = root / "ref"
    raw.mkdir(parents=True, exist_ok=True)
    ref.mkdir(parents=True, exist_ok=True)
    moves = [
        ("fastq", "data/raw/fastq"),
        ("pbmc_1k_v3_count", "data/raw/pbmc_1k_v3_count"),
        ("pbmc_1k_v3_count_backup", "data/raw/pbmc_1k_v3_count_backup"),
        ("lung_carcinoma_3k_count", "data/raw/lung_carcinoma_3k_count"),
        ("cellranger_runs", "data/raw/cellranger_runs"),
        ("reference", "ref/reference"),
    ]
    for src_rel, dst_rel in moves:
        src = root / src_rel
        dst = root / dst_rel
        if src.exists() and not src.is_symlink() and not dst.exists():
            shutil.move(str(src), str(dst))


if __name__ == "__main__":
    migrate_layout(Path("/home/zerlinshen/singlecell_factory"))
