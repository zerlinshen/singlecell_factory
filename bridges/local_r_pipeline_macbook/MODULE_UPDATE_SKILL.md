# SCRA Pipeline Module Update Skill (Reusable checklist)

Use this checklist every time a new module is added or an existing module is changed.

- Update module list section in [README.md](README.md).
- Update or append rows in the **Module provenance and references** table.
- If new dependencies are required, add them to `main.R` `required_pkgs`.
- Add at least one usage example if behavior changes.
- Keep any output file naming changes backward-compatible when possible.
- If external package methods are introduced, include:
  - package/repository name
  - article or documentation reference
  - whether behavior is deterministic for reproducibility (set seed where needed)
- Run a short smoke review:
  - `Rscript main.R --help`
  - `Rscript main.R --help | grep -n "module name"` (when relevant)
- Run module maintenance check before finishing:
  - `Rscript scripts/check_module_updates.R`

This project rule is now part of the pipeline maintenance workflow.
