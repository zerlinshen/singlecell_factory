"""Apply the suite's publication theme to the pipeline's diagnostic figures.

Why this exists
---------------
``plotting_factory/theme/nature_high_impact.yaml`` is a carefully built
journal template (7 pt type, 89/183 mm columns, colourblind-safe palettes,
editable vector text, 300/600 dpi export profiles). Until now **no pipeline
figure ever saw it**: a grep for ``plotting_factory``/``theme_loader`` across
``workflow/modular`` returned nothing, and the ~16 figures every run emits were
drawn with raw scanpy/matplotlib defaults and saved at ``dpi=160`` — below the
300 dpi review floor, with a full four-sided axes box, oversized bold labels and
scanpy's default (non-colourblind-safe, silently wrapping) categorical palette.

That is why the run figures look unlike the governed figure packages: they were
never styled at all. This module closes that gap by applying the same tokens the
presentation layer uses, so a diagnostic figure and a manuscript figure share one
visual language and one dpi contract.

Scope and boundary
------------------
This is *styling only*. It sets matplotlib rcParams; it does not compute,
transform, or render anything, so the render-only ownership boundary is
unaffected — the pipeline still owns its diagnostic figures, it just stops
drawing them in matplotlib's defaults.

``singlecell_factory`` must remain runnable when the sibling presentation repo is
absent, so a missing theme is reported and recorded in the manifest rather than
raised: the run continues with matplotlib defaults and says so, instead of
silently claiming a styling it did not apply.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# <suite-root>/singlecell_factory/workflow/modular/_figure_theme.py -> <suite-root>
_SUITE_ROOT = Path(__file__).resolve().parents[3]
_THEME_PATH = _SUITE_ROOT / "plotting_factory" / "theme" / "nature_high_impact.yaml"

# Export profile whose dpi governs the run's raster diagnostics.
_RASTER_PROFILE = "review_raster"


def _rcparams(theme: dict[str, Any]) -> dict[str, Any]:
    """Mirror plotting_factory's ``rcparams_from_theme`` for the raster lane.

    Kept token-driven (no magic numbers) so the two surfaces cannot drift on a
    value; the presentation layer stays the single source of the tokens.
    """
    typography = theme["typography"]
    geometry = theme["geometry"]
    dpi = int(theme["export_profiles"][_RASTER_PROFILE]["dpi"])
    return {
        "font.family": "sans-serif",
        "font.sans-serif": [
            typography["font_family"], "Arimo", "Helvetica", "DejaVu Sans",
        ],
        "font.size": typography["base_size_pt"],
        "axes.titlesize": typography["title_size_pt"],
        "axes.labelsize": typography["base_size_pt"],
        "xtick.labelsize": typography["small_size_pt"],
        "ytick.labelsize": typography["small_size_pt"],
        "legend.fontsize": typography["small_size_pt"],
        "axes.linewidth": geometry["axis_line_width_pt"],
        "xtick.major.width": geometry["axis_line_width_pt"],
        "ytick.major.width": geometry["axis_line_width_pt"],
        "grid.linewidth": geometry["grid_line_width_pt"],
        "lines.linewidth": geometry["axis_line_width_pt"],
        # Publication conventions: open axes, no legend box, white canvas.
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        # Governs every savefig that does not pass an explicit dpi. The previous
        # hardcoded dpi=160 was below the 300 dpi review floor.
        "savefig.dpi": dpi,
        "figure.dpi": dpi,
        # Editable text if a diagnostic is ever exported to vector.
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }


def _load_theme() -> dict[str, Any] | None:
    try:
        import yaml
        return yaml.safe_load(_THEME_PATH.read_text(encoding="utf-8"))
    except Exception:
        return None


def journal_figure_size(column: str = "double", height_mm: float = 70.0) -> tuple[float, float]:
    """Figure size in inches at a real journal column width.

    A panel drawn at an arbitrary size is scaled by the typesetter to fit the
    column, and that scaling shrinks every font with it: a 12-inch canvas placed
    in a 183 mm column reduces 7 pt type to roughly 3.5 pt, well under the ~5 pt
    floor. Sizing to the column up front is what makes the typography tokens mean
    anything on the printed page.
    """
    theme = _load_theme()
    if theme is None:
        return (7.2, 5.0) if column == "double" else (3.5, 2.8)
    geometry = theme["geometry"]
    key = "double_column_width_mm" if column == "double" else "single_column_width_mm"
    return (float(geometry[key]) / 25.4, float(height_mm) / 25.4)


def qualitative_colors(token: str, n: int) -> list[str]:
    """``n`` distinct hex colours from a theme palette, interpolating past its size.

    Mirrors ``plotting_factory``'s ``theme_loader.qualitative_colors``: indexing a
    fixed palette modulo its length gives two categories the identical colour
    while the legend still lists both. Real runs reach 40+ Leiden clusters, so
    that regime is the norm here, not an edge case.
    """
    theme = _load_theme()
    palettes = (theme or {}).get("palettes", {})
    pal = palettes.get(token)
    if not isinstance(pal, list) or not pal or n <= 0:
        return []
    if n <= len(pal):
        return list(pal[:n])

    # Beyond the primary palette, widen the anchor set with the theme's OTHER
    # qualitative palettes before falling back to interpolation. Interpolating a
    # single 8-colour ramp to 40+ points yields near-identical neighbours — the
    # colour collision returns in a subtler form — whereas the theme collectively
    # supplies ~30 curated, colourblind-conscious hues to draw on first.
    seen: list[str] = []
    for name in (token, "cell_type_qualitative", "sample_qualitative",
                 "leiden_qualitative", "group_qualitative"):
        candidate = palettes.get(name)
        if not isinstance(candidate, list):
            continue
        for colour in candidate:
            if isinstance(colour, str) and colour.upper() not in {c.upper() for c in seen}:
                seen.append(colour)
    if n <= len(seen):
        return seen[:n]

    from matplotlib.colors import LinearSegmentedColormap, to_hex

    cmap = LinearSegmentedColormap.from_list(f"{token}_interp", seen or list(pal))
    return [to_hex(cmap(i / (n - 1))) for i in range(n)]


def apply_pipeline_figure_theme() -> dict[str, Any]:
    """Apply the suite publication theme to matplotlib; return provenance.

    Returns a dict recorded in the run manifest so a reader can tell whether a
    run's figures were themed, and by which token file.
    """
    provenance: dict[str, Any] = {"theme_path": str(_THEME_PATH)}
    try:
        import matplotlib as mpl
        import yaml
    except Exception as exc:  # pragma: no cover - matplotlib is a hard dep
        provenance.update(status="unavailable", reason=f"import failed: {exc}")
        return provenance

    if not _THEME_PATH.is_file():
        provenance.update(status="theme_not_found")
        logger.warning(
            "FIGURE_THEME: %s not found; run figures keep matplotlib defaults "
            "(no journal typography, open axes, or dpi floor). Recorded as "
            "figure_theme.status=theme_not_found.",
            _THEME_PATH,
        )
        return provenance

    try:
        theme = yaml.safe_load(_THEME_PATH.read_text(encoding="utf-8"))
        params = _rcparams(theme)
        mpl.rcParams.update(params)
    except Exception as exc:
        provenance.update(status="failed", reason=str(exc))
        logger.warning("FIGURE_THEME: could not apply %s: %s", _THEME_PATH, exc)
        return provenance

    provenance.update(
        status="applied",
        journal_style=theme.get("journal_style"),
        schema_version=theme.get("schema_version"),
        font_family=theme["typography"]["font_family"],
        base_size_pt=theme["typography"]["base_size_pt"],
        raster_dpi=params["savefig.dpi"],
    )
    logger.info(
        "FIGURE_THEME: applied %s (%s) — %s %gpt, raster dpi %d.",
        theme.get("journal_style"), _THEME_PATH.name,
        theme["typography"]["font_family"], theme["typography"]["base_size_pt"],
        params["savefig.dpi"],
    )
    return provenance
