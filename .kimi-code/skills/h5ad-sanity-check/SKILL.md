---
name: h5ad-sanity-check
description: Use after a pipeline module writes an .h5ad output and before accepting the run or reporting success. Runs a structural probe (shape, sparsity, NaN/negatives, duplicate barcodes, obs/var/layers/uns inventory) and issues a verdict against the run contract. Catches silent corruption that logs and plots miss.
---

# h5ad-sanity-check — structural evidence gate for AnnData outputs

Processing principle made executable: **every pipeline output is validated before anyone
claims success.** This skill checks the data structure; `plot-qc-review` checks figures;
together they replace "the log said OK".

## Run the probe

Use a python that has `anndata` (a project env under `envs/`, or the env the run used).
**Prefer the env that wrote the file**: newer anndata writes encodings older installs
cannot read (`IOSpec ... no read method registered`). If the probe fails that way,
retry with the writer env before concluding anything; a file only its writer env can
read is itself a portability NOTE for multi-env suites, not proof of corruption.

```sh
<env-python> .kimi-code/skills/h5ad-sanity-check/h5ad_sanity.py <output.h5ad>
```

It prints JSON facts: shape, X type/dtype/density, NaN and negative counts, obs/var
columns, layers/obsm/uns inventory, duplicate barcode count, MT-prefix presence.
Backed-mode read, safe for 900k-cell files.

## Verdict rules (compare facts against `contracts/project_run_contract.yaml`)

FAIL (reject the output, route to `execute-and-recover-pipeline`):
- `n_obs == 0` or `n_vars == 0`
- `obs_name_duplicates > 0`
- `x_nan > 0` (or `x_nan_sampled > 0` for backed/dense arrays; unless the module documents NaN masking)
- An obs column, layer, or embedding (`obsm`) the contract requires is missing

NOTES (accept with caveats recorded in the run brief):
- `x_negative > 0` / `x_negative_sampled > 0` in a matrix expected to hold raw counts
- `x_density` / `x_density_sampled` ≈ 1.0 where a sparse counts matrix was expected
- `mt_gene_prefix_present == false` while percent-mt QC is part of the workflow

PASS: none of the above, and shape is consistent with the run's input expectations.

## Reporting

State the verdict, the JSON facts it rests on (quote the numbers), and the contract
fields compared. Never report "validated" without the probe output attached.
