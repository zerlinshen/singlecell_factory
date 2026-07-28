"""Tests verifying NC2024 paper-aligned parameter defaults and cohort subset behavior.

Rewritten 2026-07-28. The paper's clustering geometry (15-PC Harmony space,
Leiden resolution 1.0) used to be asserted on the *dataclass defaults*, which is
how those values ended up shadowing every programmatic `PipelineConfig()` caller
while the CLI resolved the canonical 40 PCs / resolution 0.8. The paper values
are now asserted through the named `paper-15pc` scientific profile, and
`test_dataclass_defaults_match_canonical_profile` guards the divergence itself.
"""
from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData
from pathlib import Path


# ---------------------------------------------------------------------------
# 1. DoubletConfig: n_prin_comps=30, expected_doublet_rate=0.06
# ---------------------------------------------------------------------------
def test_doublet_config_paper_defaults():
    from workflow.modular.config import DoubletConfig
    cfg = DoubletConfig()
    assert cfg.n_prin_comps == 30, "Paper uses n_prin_comps=30 for Scrublet"
    assert cfg.expected_doublet_rate == 0.06, "Paper uses expected_doublet_rate=0.06"


# ---------------------------------------------------------------------------
# 2. ClusteringConfig carries the CANONICAL profile, not the paper's numbers.
#
# 40 PCs / resolution 0.8 is what the CLI has always resolved. The dataclass
# used to carry the paper's 15 / 1.0, so `PipelineConfig()` in a notebook or a
# test ran a different analysis than the identical CLI invocation, with nothing
# in the manifest saying so. Rationale for 40 over 15: rare populations in
# heterogeneous multi-batch tumour tissue carry signal past the first ~15
# components, and over-inclusion costs far less than truncation (Luecken & Theis
# 2019, Mol Syst Biol 15:e8746; Heumos 2023, Nat Rev Genet 24:550-572).
# ---------------------------------------------------------------------------
def test_clustering_config_canonical_defaults():
    from workflow.modular.config import ClusteringConfig
    cfg = ClusteringConfig()
    assert cfg.n_pcs == 40, "Programmatic default must be the canonical 40-PC space"
    assert cfg.leiden_resolution == 0.8, "Programmatic default must be canonical resolution=0.8"


# ---------------------------------------------------------------------------
# 2b. THE drift guard.
#
# Every dataclass default that has a canonical-profile counterpart must equal
# it. This is the test that would have caught the 15/1.0 vs 40/0.8 split at the
# moment it was introduced, and it catches the next one regardless of which
# side moves. The explicit mapping is also a coverage gate: adding a key to
# _CANONICAL_SCIENTIFIC_PARAMETERS without wiring it to a config field fails
# here rather than silently going unguarded.
# ---------------------------------------------------------------------------
def _canonical_to_dataclass_values():
    """Return {canonical_key: value read off the config dataclass defaults}."""
    from workflow.modular.config import ClusteringConfig, PipelineConfig

    clustering = ClusteringConfig()
    pipeline = PipelineConfig(project="x", output_dir=Path("/tmp"), cellranger=None)
    return {
        # CLI carries optional_modules as a comma-joined string; the dataclass
        # carries the parsed list. Compare in the CLI's shape.
        "optional_modules": ",".join(pipeline.optional_modules),
        "n_top_genes": clustering.n_top_genes,
        "n_pcs": clustering.n_pcs,
        "n_neighbors": clustering.n_neighbors,
        "leiden_resolution": clustering.leiden_resolution,
        "de_n_genes": pipeline.de_n_genes,
        "doublet_strategy": pipeline.doublet_strategy,
        "clustering_engine": pipeline.clustering_engine,
    }


def test_dataclass_defaults_match_canonical_profile():
    from workflow.modular.cli import _CANONICAL_SCIENTIFIC_PARAMETERS

    observed = _canonical_to_dataclass_values()

    assert set(observed) == set(_CANONICAL_SCIENTIFIC_PARAMETERS), (
        "A canonical scientific parameter has no mapped config dataclass field "
        "(or vice versa). Wire it into _canonical_to_dataclass_values so the "
        "programmatic and CLI entrypoints stay pinned to the same science."
    )
    assert observed == dict(_CANONICAL_SCIENTIFIC_PARAMETERS), (
        "config.py dataclass defaults diverged from the canonical scientific "
        "profile. A programmatic PipelineConfig() caller would silently run a "
        "different analysis than the same run launched through the CLI."
    )


# ---------------------------------------------------------------------------
# 2c. The CLI's argparse defaults must also resolve to canonical (third copy).
# ---------------------------------------------------------------------------
def test_cli_defaults_match_canonical_profile(monkeypatch):
    import workflow.modular.cli as cli_mod

    monkeypatch.setattr("sys.argv", ["prog", "--project", "x", "--sample-root", "/tmp"])
    args = cli_mod._apply_scientific_profile(cli_mod.parse_args())

    for field_name, canonical_value in cli_mod._CANONICAL_SCIENTIFIC_PARAMETERS.items():
        assert getattr(args, field_name) == canonical_value, (
            f"CLI default for {field_name} drifted from the canonical profile"
        )
    assert args.resolved_scientific_parameter_diff == {}


# ---------------------------------------------------------------------------
# 2d. The paper's clustering geometry stays reachable — through a NAMED profile.
#
# Source: docs/PUBLICATION_READY.md — Harmony on 15 PCs, Leiden resolution 1.0.
# ---------------------------------------------------------------------------
def test_paper_profile_reaches_paper_clustering_geometry(monkeypatch):
    import workflow.modular.cli as cli_mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project", "x",
            "--sample-root", "/tmp",
            "--scientific-profile", "paper-15pc",
            "--acknowledge-scientific-non-equivalence",
        ],
    )
    args = cli_mod._apply_scientific_profile(cli_mod.parse_args())

    assert args.n_pcs == 15, "paper-15pc must cluster on the 15-PC Harmony space"
    assert args.leiden_resolution == 1.0, "paper-15pc must use Leiden resolution=1.0"
    # The departure has to be auditable in the manifest, not just applied.
    assert args.resolved_scientific_parameter_diff == {
        "n_pcs": {"canonical": 40, "resolved": 15},
        "leiden_resolution": {"canonical": 0.8, "resolved": 1.0},
    }
    assert args.scientific_non_equivalence_acknowledged is True
    # Clustering geometry only: DE thresholds stay canonical/launcher-level.
    assert args.de_n_genes == 300
    assert args.n_top_genes == 3000


def test_paper_profile_requires_non_equivalence_acknowledgement(monkeypatch):
    import workflow.modular.cli as cli_mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project", "x",
            "--sample-root", "/tmp",
            "--scientific-profile", "paper-15pc",
        ],
    )
    with pytest.raises(SystemExit, match="not scientifically equivalent"):
        cli_mod._apply_scientific_profile(cli_mod.parse_args())


def _profile_names():
    from workflow.modular.cli import _SCIENTIFIC_PROFILE_OVERRIDES
    return list(_SCIENTIFIC_PROFILE_OVERRIDES)


@pytest.mark.parametrize("profile", _profile_names())
def test_every_implemented_profile_is_selectable(profile, monkeypatch):
    """Advertised choices and the override table must be the same list.

    A profile that exists in the table but is not accepted by --scientific-profile
    is unreachable; one accepted but missing from the table KeyErrors inside
    _apply_scientific_profile. Deriving `choices` from the table makes both
    impossible, and this test pins that wiring.
    """
    import workflow.modular.cli as cli_mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog", "--project", "x", "--sample-root", "/tmp",
            "--scientific-profile", profile,
            "--acknowledge-scientific-non-equivalence",
        ],
    )
    args = cli_mod._apply_scientific_profile(cli_mod.parse_args())
    assert args.scientific_profile == profile


def test_unknown_scientific_profile_is_rejected(monkeypatch):
    import workflow.modular.cli as cli_mod

    monkeypatch.setattr(
        "sys.argv",
        ["prog", "--project", "x", "--sample-root", "/tmp",
         "--scientific-profile", "paper-99pc"],
    )
    with pytest.raises(SystemExit):
        cli_mod.parse_args()


# ---------------------------------------------------------------------------
# 3. DEConfig: marker_correction default is 'benjamini-hochberg' (backward compat)
# ---------------------------------------------------------------------------
def test_de_config_default_correction():
    from workflow.modular.config import DEConfig
    cfg = DEConfig()
    assert cfg.marker_correction == "benjamini-hochberg"


# ---------------------------------------------------------------------------
# 4. PseudobulkDEConfig: padj_threshold=0.05, abs_logfc_threshold=1.0
# ---------------------------------------------------------------------------
def test_pseudobulk_de_config_paper_defaults():
    from workflow.modular.config import PseudobulkDEConfig
    cfg = PseudobulkDEConfig()
    assert cfg.padj_threshold == 0.05
    assert cfg.abs_logfc_threshold == 1.0


# ---------------------------------------------------------------------------
# 5. PipelineConfig.de_correction default is 'benjamini-hochberg'
# ---------------------------------------------------------------------------
def test_pipeline_config_de_correction_default():
    from workflow.modular.config import PipelineConfig
    cfg = PipelineConfig(project="x", output_dir=Path("/tmp"), cellranger=None)
    assert cfg.de_correction == "benjamini-hochberg"


# ---------------------------------------------------------------------------
# 6. CLI --cohort-subset parses correctly into cfg.cohort_subset
# ---------------------------------------------------------------------------
def test_cli_cohort_subset_parsed(monkeypatch):
    import workflow.modular.cli as cli_mod
    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project", "test",
            "--sample-root", "/tmp",
            "--cohort-subset", "disease=lung_adenocarcinoma,lung_squamous_cell_carcinoma",
        ],
    )
    args = cli_mod.parse_args()
    assert args.cohort_subset == ["disease=lung_adenocarcinoma,lung_squamous_cell_carcinoma"]


# ---------------------------------------------------------------------------
# 7. cohort_subset filter applied correctly in CellRangerModule
# ---------------------------------------------------------------------------
def test_cohort_subset_filter_applied(tmp_path):
    from workflow.modular.modules.cellranger import CellRangerModule
    from workflow.modular.config import PipelineConfig, CellRangerConfig
    from workflow.modular.context import PipelineContext

    adata = AnnData(np.ones((6, 3), dtype=float))
    adata.obs["disease"] = ["lung_adenocarcinoma", "lung_squamous_cell_carcinoma",
                            "normal", "normal", "lung_adenocarcinoma", "other"]
    adata.obs_names = [f"C{i}" for i in range(6)]
    adata.var_names = ["G1", "G2", "G3"]

    cfg = PipelineConfig(
        project="t",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        cohort_subset=["disease=lung_adenocarcinoma,lung_squamous_cell_carcinoma"],
    )
    ctx = PipelineContext(cfg=cfg, run_dir=tmp_path, figure_dir=tmp_path, table_dir=tmp_path)
    ctx.adata = adata

    mod = CellRangerModule()
    result = mod._apply_cohort_subset(adata, ctx)

    assert result.n_obs == 3
    assert set(result.obs["disease"].unique()) == {"lung_adenocarcinoma", "lung_squamous_cell_carcinoma"}
    assert ctx.metadata["cohort_subset_n_obs_before"] == 6
    assert ctx.metadata["cohort_subset_n_obs_after"] == 3


# ---------------------------------------------------------------------------
# 8. cohort_subset=None is a no-op
# ---------------------------------------------------------------------------
def test_cohort_subset_none_is_noop(tmp_path):
    from workflow.modular.modules.cellranger import CellRangerModule
    from workflow.modular.config import PipelineConfig, CellRangerConfig
    from workflow.modular.context import PipelineContext

    adata = AnnData(np.ones((4, 2), dtype=float))
    adata.obs_names = [f"C{i}" for i in range(4)]
    adata.var_names = ["G1", "G2"]

    cfg = PipelineConfig(
        project="t",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        cohort_subset=None,
    )
    ctx = PipelineContext(cfg=cfg, run_dir=tmp_path, figure_dir=tmp_path, table_dir=tmp_path)

    mod = CellRangerModule()
    result = mod._apply_cohort_subset(adata, ctx)

    assert result.n_obs == 4
    assert "cohort_subset_applied" not in ctx.metadata
