"""Resolved software-version recorder for the Python producer.

Why this exists
---------------
``environment.yml`` / ``environment_gpu.yml`` deliberately leave most of the
scientific stack unpinned so conda can solve a working CUDA/BLAS combination on
whatever host the run lands on.  The cost of that freedom is that the manifest
alone could not identify *which* stack produced a run: two runs of the same git
SHA on the same input could disagree numerically and nothing in the run record
would explain it.  Recording the resolved versions closes that gap without
freezing the solver.

This is not hypothetical for this factory.  ``batch_correction._run_harmony_direct``
routes around a wrapper bug that exists in the *specific* pair
(scanpy 1.12, harmonypy 0.2.0), and ``rna_velocity._patch_scvelo_numpy2`` only
engages for numpy >= 2 with scVelo < 0.5.  Whether those version-conditional
code paths were active during a given run is unrecoverable unless the run
writes down the versions it actually imported.

Resolution strategy
-------------------
Versions come from :mod:`importlib.metadata`, i.e. the *installed distribution*
metadata, not from a module ``__version__`` attribute.  That choice matters:
metadata lookup neither imports nor initialises the package, so calling this at
manifest-write time cannot pull CUDA contexts or 100 MB of scanpy into a
process that did not otherwise need them.

Distribution names are not import names (``scikit-learn`` -> ``sklearn``), and
python-igraph has shipped its metadata under both ``igraph`` and
``python-igraph`` depending on how it was installed, so every entry carries a
tuple of candidate distribution names and the first hit wins.  As a last resort
a package that is *already imported* is read from ``sys.modules``; nothing new
is ever imported here.

Absent packages are recorded explicitly as ``None`` rather than omitted.  A
missing optional backend is itself provenance — "scrublet is absent" proves the
CPU Scrublet lane in ``doublet_detection`` could not have run.
"""

from __future__ import annotations

import platform
import sys
import warnings
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _distribution_version

# Packages whose version can move a number in a published figure. Each entry is
# (canonical_name, candidate_distribution_names, import_name).
#
# CORE: present in every environment that can run the mandatory pipeline
# (cellranger -> qc -> ambient_correction -> doublet_detection) plus clustering.
# A missing core package means the manifest is describing a broken environment
# and is worth surfacing loudly via ``unresolved``.
_CORE_PACKAGES: tuple[tuple[str, tuple[str, ...], str], ...] = (
    ("scanpy", ("scanpy",), "scanpy"),
    ("anndata", ("anndata",), "anndata"),
    ("numpy", ("numpy",), "numpy"),
    ("scipy", ("scipy",), "scipy"),
    ("scikit-learn", ("scikit-learn", "sklearn"), "sklearn"),
    ("pandas", ("pandas",), "pandas"),
    ("matplotlib", ("matplotlib",), "matplotlib"),
    # Leiden is the production clustering call
    # (``sc.tl.leiden(flavor="igraph")``); its two backends set the partition
    # quality and therefore the cluster count.
    ("leidenalg", ("leidenalg",), "leidenalg"),
    ("igraph", ("igraph", "python-igraph"), "igraph"),
)

# OPTIONAL: legitimately absent in a minimal or CPU-only install. Recorded as
# None when missing so a reader can tell "backend unavailable" apart from
# "recorder did not look".
_OPTIONAL_PACKAGES: tuple[tuple[str, tuple[str, ...], str], ...] = (
    ("rapids-singlecell", ("rapids-singlecell", "rapids_singlecell"), "rapids_singlecell"),
    ("harmonypy", ("harmonypy",), "harmonypy"),
    ("scrublet", ("scrublet",), "scrublet"),
)

CORE_PACKAGE_NAMES: tuple[str, ...] = tuple(name for name, _dists, _mod in _CORE_PACKAGES)
OPTIONAL_PACKAGE_NAMES: tuple[str, ...] = tuple(name for name, _dists, _mod in _OPTIONAL_PACKAGES)
RECORDED_PACKAGE_NAMES: tuple[str, ...] = CORE_PACKAGE_NAMES + OPTIONAL_PACKAGE_NAMES


def resolve_package_version(
    distribution_names: tuple[str, ...] | str,
    import_name: str | None = None,
) -> str | None:
    """Return the installed version of a package, or None when unresolvable.

    Tries each candidate distribution name in order, then falls back to the
    ``__version__`` of an already-imported module. Never imports anything and
    never raises: a provenance recorder that can abort a completed run would
    trade real results for bookkeeping.
    """
    if isinstance(distribution_names, str):
        distribution_names = (distribution_names,)

    for dist in distribution_names:
        try:
            resolved = _distribution_version(dist)
        except PackageNotFoundError:
            continue
        except Exception:  # corrupt or unreadable dist-info on disk
            continue
        if resolved:
            return str(resolved)

    # Fallback for packages installed without metadata (e.g. a source tree on
    # PYTHONPATH). Restricted to sys.modules so this stays import-free.
    # anndata >= 0.11 raises FutureWarning on ``__version__`` access; recording
    # provenance must not pollute a run's warning stream, so it is suppressed
    # here rather than at the caller.
    if import_name:
        module = sys.modules.get(import_name)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            attr = getattr(module, "__version__", None) if module is not None else None
        if attr:
            return str(attr)

    return None


def collect_software_versions() -> dict:
    """Return the ``software_versions`` manifest block for this interpreter.

    Shape::

        {
          "python": "3.13.5",
          "python_implementation": "CPython",
          "platform": "Linux-...-x86_64-with-glibc2.39",
          "executable": "/home/.../envs/sc_gpu/bin/python",
          "packages": {"scanpy": "1.12", ..., "scrublet": None},
          "core_packages": [...],
          "unresolved": ["scrublet"],
          "core_complete": true
        }

    ``executable`` is recorded because it identifies the conda env (sc_gpu vs
    sc10x) that the manifest's versions were read from, which the version
    strings alone do not.
    """
    packages: dict[str, str | None] = {}
    for name, dists, import_name in _CORE_PACKAGES + _OPTIONAL_PACKAGES:
        packages[name] = resolve_package_version(dists, import_name)

    unresolved = sorted(name for name, value in packages.items() if not value)

    return {
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "executable": sys.executable,
        "packages": packages,
        "core_packages": list(CORE_PACKAGE_NAMES),
        "unresolved": unresolved,
        "core_complete": all(packages.get(name) for name in CORE_PACKAGE_NAMES),
    }


if __name__ == "__main__":
    import json

    print(json.dumps(collect_software_versions(), indent=2))
