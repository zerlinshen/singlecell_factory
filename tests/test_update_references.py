from pathlib import Path

from scripts import update_references


def test_extract_references_accepts_null_doi(tmp_path):
    modules_dir = tmp_path / "modules"
    modules_dir.mkdir()
    (modules_dir / "pseudo_velocity.py").write_text(
        """
__references__ = {
    "project_local_proxy": {
        "title": "Project-local KNN-based pseudo-velocity proxy",
        "authors": "singlecell_factory contributors",
        "journal": "Internal documentation",
        "year": "2025",
        "doi": None,
        "description": "proxy method",
    },
}
""",
        encoding="utf-8",
    )

    refs = update_references.extract_references_from_modules(modules_dir)

    assert len(refs) == 1
    assert refs[0]["doi"] is None
    assert refs[0]["module"] == "pseudo_velocity"


def test_update_readme_renders_null_doi_as_na(tmp_path):
    readme = tmp_path / "README.md"
    readme.write_text(
        "### Complete Citation List\n\n"
        "| # | Reference | DOI | Used by |\n"
        "|---|---|---|---|\n"
        "\n---\n",
        encoding="utf-8",
    )
    refs = [
        {
            "authors": "singlecell_factory contributors",
            "journal": "Internal documentation",
            "year": "2025",
            "title": "Project-local KNN-based pseudo-velocity proxy",
            "doi": None,
            "module": "pseudo_velocity",
            "description": "proxy method",
        }
    ]

    changed = update_references.update_readme_references(Path(readme), refs)

    assert changed is True
    content = readme.read_text(encoding="utf-8")
    assert (
        "| 1 | singlecell_factory contributors et al., "
        "*Internal documentation*, 2025 | N/A | "
        "`pseudo_velocity` (proxy method) |"
    ) in content
