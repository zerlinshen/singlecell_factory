# Mac Pull Recipe — NC2024 v2 Runs (2026-04-26)

Generated: 2026-04-26. Do NOT run on Mac automatically — user performs pull manually.

---

## Remote host

| Field | Value |
|---|---|
| Tailscale FQDN | `zerlinmacbook-pro` |
| Remote user | `zerlinshen` |
| Remote repo root | `/home/zerlinshen/singlecell_factory` |

---

## Packed assets

### Tumor cohort

| File | Size | SHA256 |
|---|---|---|
| `bundle.tgz` | 94M | `0def56a7514797c01f4ef66fb9647f266ce7a2c826402c1a3a61f896a3a00257` |
| `figures.tgz` | 1.1M | `7ced4988dc8561cefe84298f6cbd7cc40830fcb0d4196c400662c174acaca4ce` |
| `report.tgz` | 4.0K | `ec56d2f2e37ae55c6d1b3aefa73542543e240620592bc5352ea5d3f819c21e6c` |

Remote assets path: `/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2/mac_assets/`

### Background + healthy cohort

| File | Size | SHA256 |
|---|---|---|
| `bundle.tgz` | 996K | `edfb9ae23dc7222a7e746c1abf3472d8b63a38c3772658c1860ac6e53cfa721a` |
| `figures.tgz` | 372K | `5be86e15371ff85972b4b34582258e3719733fb202a479ea5053e6b5b0b8c67b` |
| `report.tgz` | 4.0K | `c2775e0e670c5b09ade8dcada7538d6a34dad26532d9af863dfd37f370a821fb` |

Remote assets path: `/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/mac_assets/`

---

## Pull commands (run ON THE MAC)

The pack script generated ready-to-run pull scripts on the remote. You can either
use those directly or use the scp commands below.

### Option A — use generated pull scripts

```bash
# Tumor
ssh zerlinshen@zerlinmacbook-pro \
  "bash /home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2/mac_assets/mac_pull_command.sh \
  zerlinshen@zerlinmacbook-pro ~/scrna_remote_pulls/nc2024_tumor_20260426_v2"

# B/H
ssh zerlinshen@zerlinmacbook-pro \
  "bash /home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/mac_assets/mac_pull_command.sh \
  zerlinshen@zerlinmacbook-pro ~/scrna_remote_pulls/nc2024_bh_20260426_v2"
```

### Option B — manual scp

```bash
# Tumor
REMOTE=zerlinshen@zerlinmacbook-pro
LOCAL_TUMOR=~/scrna_remote_pulls/nc2024_tumor_20260426_v2
mkdir -p "$LOCAL_TUMOR"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2/mac_assets/bundle.tgz" "$LOCAL_TUMOR/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2/mac_assets/bundle.tgz.sha256" "$LOCAL_TUMOR/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2/mac_assets/figures.tgz" "$LOCAL_TUMOR/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2/mac_assets/figures.tgz.sha256" "$LOCAL_TUMOR/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2/mac_assets/report.tgz" "$LOCAL_TUMOR/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_tumor_20260426_v2/mac_assets/report.tgz.sha256" "$LOCAL_TUMOR/"

# B/H
LOCAL_BH=~/scrna_remote_pulls/nc2024_bh_20260426_v2
mkdir -p "$LOCAL_BH"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/mac_assets/bundle.tgz" "$LOCAL_BH/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/mac_assets/bundle.tgz.sha256" "$LOCAL_BH/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/mac_assets/figures.tgz" "$LOCAL_BH/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/mac_assets/figures.tgz.sha256" "$LOCAL_BH/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/mac_assets/report.tgz" "$LOCAL_BH/"
scp "$REMOTE:/home/zerlinshen/singlecell_factory/results/nc2024_bh_20260426_v2/mac_assets/report.tgz.sha256" "$LOCAL_BH/"
```

---

## SHA256 verification (on Mac after pull)

```bash
# Verify tumor bundle
sha256sum -c ~/scrna_remote_pulls/nc2024_tumor_20260426_v2/bundle.tgz.sha256

# Verify BH bundle
sha256sum -c ~/scrna_remote_pulls/nc2024_bh_20260426_v2/bundle.tgz.sha256
```

---

## R analysis (run ON THE MAC after pull)

### Required R packages

```r
install.packages(c("Seurat", "arrow", "Matrix"))
```

### Load and render

```bash
# Tumor
Rscript multiomics_r_factory/main.R ~/scrna_remote_pulls/nc2024_tumor_20260426_v2/bundle.tgz

# B/H
Rscript multiomics_r_factory/main.R ~/scrna_remote_pulls/nc2024_bh_20260426_v2/bundle.tgz
```

Or extract first then point at the directory:

```bash
mkdir -p ~/scrna_remote_pulls/nc2024_tumor_20260426_v2/r_bundle
tar -xzf ~/scrna_remote_pulls/nc2024_tumor_20260426_v2/bundle.tgz \
  -C ~/scrna_remote_pulls/nc2024_tumor_20260426_v2/r_bundle

Rscript multiomics_r_factory/main.R \
  ~/scrna_remote_pulls/nc2024_tumor_20260426_v2/r_bundle
```

---

## R bundle structure (verified on remote)

Both bundles contain the full v2 schema:

```
r_bundle/
├── obs.parquet            — cell metadata (cell_type, leiden, batch, etc.)
├── obsm/
│   ├── X_pca.parquet      — 15-dim Harmony embedding
│   └── X_umap.parquet     — 2-dim UMAP coordinates
├── marker_expr.parquet    — per-cluster marker gene expression
├── bundle_manifest.json   — schema version, export args, source H5AD hash
├── bundle_manifest.tsv    — human-readable manifest
├── uns_summary.json       — pipeline metadata from adata.uns
└── README.md              — bundle usage notes
```

Note: counts matrix (`.mtx`) is not included in the v2 bundle by default — the
bundle is sized for R plotting/reporting, not full re-analysis. Use the remote
`final_adata.h5ad` if raw counts are needed.

---

## Bridge symlink status

Verified green at time of packaging:
```
OK: bridges/local_r_pipeline_macbook/R -> ../../../multiomics_r_factory/R
OK: bridges/local_r_pipeline_macbook/R_bundle -> ../../../multiomics_r_factory/R_bundle
PASS: Bridge symlink check succeeded.
```
