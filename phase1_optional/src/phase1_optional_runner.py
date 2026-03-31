import os
import runpy


if __name__ == "__main__":
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    target = os.path.join(
        root,
        "rna_velocity_pseudotime_analysis",
        "src",
        "rna_velocity_pseudotime_runner.py",
    )
    runpy.run_path(target, run_name="__main__")
