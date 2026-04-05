"""
Standard template for a new pipeline module.

Usage:
1. Copy this file to `workflow/modular/modules/my_new_module.py`.
2. Update the class name, module name, and docstring.
3. Add `__references__` at the module level for the citation manager to pick up.
4. Implement the `run` method, ensuring it uses the `PipelineContext` appropriately.
5. Add the module to `workflow/modular/pipeline.py` dependencies and registry.
"""
from __future__ import annotations

import logging

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from ..context import PipelineContext

logger = logging.getLogger(__name__)

# The reference manager (scripts/update_references.py) parses this dictionary.
# Add any relevant publications or tools used by this module.
__references__ = {
    "module_tool_reference": {
        "title": "A Great Tool for Single Cell Analysis",
        "authors": "Smith et al.",
        "journal": "Nature Methods",
        "year": "2024",
        "doi": "10.1038/s41592-024-0000-0",
        "description": "Used for fantastic analysis" # Optional brief description
    }
}

class MyNewModule:
    """
    Module: MyNewModule
    Description: Performs a fantastic analysis on the AnnData object.
    """

    name = "my_new_module"
    required = False # Set to True if this module must run in every pipeline execution

    def run(self, ctx: PipelineContext) -> None:
        """
        Executes the module's logic.

        Args:
            ctx (PipelineContext): The pipeline context containing the AnnData object,
                                   configuration, and output directory paths.
        """
        adata = ctx.adata
        if adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        logger.info(f"Running {self.name}...")

        # 1. Configuration Access
        # cfg = ctx.cfg.my_new_module_config
        
        # 2. Perform Analysis
        # Use GPU acceleration (rapids-singlecell) if applicable and configured
        # e.g. if getattr(ctx.cfg, "use_gpu", False):
        #          import rapids_singlecell as rsc
        #          rsc.tl.fantastic_analysis(adata)
        #      else:
        #          sc.tl.fantastic_analysis(adata)

        # 3. Store Metadata
        # Record useful statistics for the run_manifest.json
        ctx.metadata[f"{self.name}_metric"] = "value"

        # 4. Save Tables (CSV)
        # Save tabular results to the dedicated module table directory
        # df.to_csv(ctx.table_dir / f"{self.name}_results.csv", index=False)

        # 5. Save Figures (PNG)
        # Save plots to the dedicated module figure directory
        # sc.pl.umap(adata, color="fantastic_score", show=False)
        # plt.savefig(ctx.figure_dir / f"{self.name}_umap.png", dpi=160, bbox_inches="tight")
        # plt.close()

        logger.info(f"Finished {self.name}.")
