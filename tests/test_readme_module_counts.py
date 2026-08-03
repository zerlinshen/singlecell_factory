"""Guard against README module-count drift vs. the module catalog.

README.md states the mandatory/optional/total module counts by hand in two
places (the suite-intro sentence and the Features list). Both must match
``workflow/modular/module_catalog.py`` exactly, or a future module addition/
removal will silently rot the documented counts (this happened once already:
the catalog grew from 22/23 optional modules to 40 without either README
mention being updated).

This test is intentionally a drift *detector*, not a doc generator: it parses
the stated numbers out of README.md with a narrow regex and asserts them
against the catalog. If the README wording changes enough that the regex no
longer matches, the test fails loudly (via ``pytest.fail``) instead of
silently passing, so the numeric claim can never go unchecked.
"""

from __future__ import annotations

import re
from pathlib import Path

from workflow.modular.module_catalog import MANDATORY_MODULES, MODULE_SPECS, optional_module_names


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"


def _catalog_counts() -> tuple[int, int, int]:
    mandatory = len(MANDATORY_MODULES)
    optional = len(optional_module_names())
    total = len(MODULE_SPECS)
    assert mandatory + optional == total
    return mandatory, optional, total


def test_readme_intro_sentence_module_counts_match_catalog() -> None:
    mandatory, optional, total = _catalog_counts()
    text = README.read_text(encoding="utf-8")

    match = re.search(
        r"\*\*(\d+) mandatory QC modules \+ (\d+) optional analysis modules "
        r"\((\d+) total\)",
        text,
    )
    if match is None:
        raise AssertionError(
            "README.md suite-intro sentence no longer matches the expected "
            "'<N> mandatory QC modules + <N> optional analysis modules (<N> "
            "total)' phrasing. Update the regex in this test AND verify the "
            "restated counts against workflow/modular/module_catalog.py."
        )
    readme_mandatory, readme_optional, readme_total = (int(g) for g in match.groups())
    assert (readme_mandatory, readme_optional, readme_total) == (
        mandatory,
        optional,
        total,
    ), (
        f"README intro states {readme_mandatory} mandatory / {readme_optional} "
        f"optional / {readme_total} total modules, but "
        f"workflow/modular/module_catalog.py currently has {mandatory} "
        f"mandatory / {optional} optional / {total} total. Update README.md."
    )


def test_readme_features_section_module_counts_match_catalog() -> None:
    mandatory, optional, total = _catalog_counts()
    text = README.read_text(encoding="utf-8")

    mandatory_match = re.search(r"\*\*(\d+) mandatory modules\*\*", text)
    optional_match = re.search(
        r"\*\*(\d+) optional analysis modules\*\*[^\n]*\((\d+) total",
        text,
    )
    if mandatory_match is None or optional_match is None:
        raise AssertionError(
            "README.md Features section no longer matches the expected "
            "'**<N> mandatory modules**' / '**<N> optional analysis modules** "
            "... (<N> total' phrasing. Update the regex in this test AND "
            "verify the restated counts against "
            "workflow/modular/module_catalog.py."
        )
    readme_mandatory = int(mandatory_match.group(1))
    readme_optional, readme_total = (int(g) for g in optional_match.groups())
    assert (readme_mandatory, readme_optional, readme_total) == (
        mandatory,
        optional,
        total,
    ), (
        f"README Features section states {readme_mandatory} mandatory / "
        f"{readme_optional} optional / {readme_total} total modules, but "
        f"workflow/modular/module_catalog.py currently has {mandatory} "
        f"mandatory / {optional} optional / {total} total. Update README.md."
    )
