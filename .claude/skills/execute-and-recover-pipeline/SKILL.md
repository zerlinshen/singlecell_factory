---
name: execute-and-recover-pipeline
description: Execute modular pipeline runs and handle failures end-to-end with checkpoint/resume, evidence-based diagnosis, and concrete recovery actions. Use for run setup, module selection, run-time triage, failed module recovery, and post-run sanity checks.
---

Run and recover modular workflows with reproducible commands.

1. Validate run prerequisites (dependency preflight).
- Confirm sample root and expected matrix paths.
- Choose optional modules intentionally.
- Enable checkpointing by default for non-trivial runs.
- Environment check: verify conda env is activated and has core packages (scanpy, anndata, scipy, sklearn).
- Optional package check: for each requested module, verify its optional dependencies are importable (e.g., scvelo for rna_velocity, liana for cell_communication, infercnvpy for cnv_inference, cellrank for cell_fate, gseapy for pathway_analysis).
- GPU check: if rapids-singlecell is expected, verify CUDA is available (`torch.cuda.is_available()` or `cupy.cuda.runtime.getDeviceCount()`); warn if GPU modules requested but no GPU detected.
- Disk space check: verify output directory has sufficient free space (>1GB for small runs, >10GB for full pipeline with velocity).

2. Execute with reproducibility.
- Build a single explicit CLI command.
- Start with conservative parallelism when debugging.
- Persist command and environment assumptions in the summary.

3. Diagnose failures with evidence-first triage.
- Inspect `module_status.csv`, `run_manifest.json`, and module artifacts.
- Identify first failing module and nearest causal error.
- Distinguish root cause from downstream cascade errors.

4. Recover instead of restarting blindly.
- Prefer `--resume-from <failed-module>` using same project/module set.
- Provide exact recovery command and expected next checks.
- If checkpoint state is invalid, justify full rerun.

5. Close run with minimal quality gate.
- Confirm all modules reach expected status.
- Verify key output artifacts exist for requested analysis scope.
- Escalate to `/review-and-validate-quality` for full scientific/reporting checks.
