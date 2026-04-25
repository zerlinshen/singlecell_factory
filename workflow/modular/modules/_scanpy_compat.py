from __future__ import annotations

from types import SimpleNamespace
from typing import Any


def _make_stub(exc: Exception) -> Any:
    """Return a lightweight scanpy-like stub for test monkeypatching."""
    return SimpleNamespace(
        read_10x_mtx=None,
        pp=SimpleNamespace(
            calculate_qc_metrics=None,
            filter_genes=None,
            normalize_total=None,
            log1p=None,
            scale=None,
            highly_variable_genes=None,
            neighbors=None,
            regress_out=None,
            combat=None,
        ),
        tl=SimpleNamespace(
            pca=None,
            umap=None,
            leiden=None,
            rank_genes_groups=None,
            score_genes=None,
            score_genes_cell_cycle=None,
            paga=None,
            diffmap=None,
            dpt=None,
        ),
        pl=SimpleNamespace(
            umap=None,
            violin=None,
            scatter=None,
            embedding=None,
            paga=None,
            rank_genes_groups_dotplot=None,
            rank_genes_groups_heatmap=None,
        ),
        get=SimpleNamespace(rank_genes_groups_df=None),
        external=SimpleNamespace(pp=SimpleNamespace(harmony_integrate=None)),
        settings=SimpleNamespace(n_jobs=1),
        _scanpy_stub=True,
        _scanpy_import_error=exc,
    )


def import_scanpy_or_stub() -> Any:
    """Import scanpy, or return a stub if import fails in this environment."""
    try:
        import scanpy as sc

        return sc
    except Exception as exc:  # pragma: no cover - exercised in environment-dependent tests
        return _make_stub(exc)


def has_api(sc_obj: Any, dotted_path: str) -> bool:
    """Check whether a scanpy object exposes a dotted API path."""
    cur = sc_obj
    for part in dotted_path.split("."):
        if not hasattr(cur, part):
            return False
        cur = getattr(cur, part)
    return cur is not None


def scanpy_import_error(sc_obj: Any) -> Exception | None:
    return getattr(sc_obj, "_scanpy_import_error", None)
