#!/usr/bin/env python3
"""
Calculate QC metrics for AnnData using scanpy.

Adds the following annotations:
- var: 'mt', 'ribo', 'hb' (boolean flags for gene types)
- obs: 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb'
"""

import os
import sys

# Fix numba caching issue in read-only containers
os.environ['NUMBA_CACHE_DIR'] = '/tmp/numba_cache'
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'
os.environ['XDG_CACHE_HOME'] = '/tmp/cache'

import platform
import importlib.metadata
import yaml
import scanpy as sc
import anndata as ad
import scipy.sparse
import numpy as np

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}" + ".h5ad"

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Calculating QC metrics for: {input_adata}")
print(f"Shape: {adata.shape}")

# Validate input data
n_cells, n_genes = adata.shape

if n_cells == 0:
    print("ERROR: AnnData has 0 cells/spots. Cannot calculate QC metrics.")
    sys.exit(1)

if n_genes == 0:
    print("ERROR: AnnData has 0 genes. Cannot calculate QC metrics.")
    sys.exit(1)

if n_genes < 10:
    print(f"WARNING: AnnData has only {n_genes} genes. This may indicate a problem with the input data.")

# Store raw counts layer if not already present
if "raw" not in adata.layers:
    adata.layers["raw"] = adata.X.copy()

# Ensure X is in sparse format for compatibility
if not scipy.sparse.issparse(adata.X):
    adata.X = scipy.sparse.csr_matrix(adata.X)

# Make var names unique before processing
adata.var_names_make_unique()

# Annotate gene types in var
# Mitochondrial genes (MT- prefix, case-insensitive)
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

# Ribosomal genes (RPS/RPL prefix, case-insensitive)
adata.var["ribo"] = adata.var_names.str.upper().str.match(r"^RP[SL]")

# Haemoglobin genes (HBA/HBB prefix, case-insensitive)
adata.var["hb"] = adata.var_names.str.upper().str.match(r"^HB[^P]")

# Print gene type counts
n_mt = adata.var["mt"].sum()
n_ribo = adata.var["ribo"].sum()
n_hb = adata.var["hb"].sum()

print("Gene type counts:")
print(f"  MT genes: {n_mt}")
print(f"  Ribo genes: {n_ribo}")
print(f"  Hb genes: {n_hb}")

# Determine which qc_vars to use (only those with at least one gene)
qc_vars = []
for var_name, count in [("mt", n_mt), ("ribo", n_ribo), ("hb", n_hb)]:
    if count > 0:
        qc_vars.append(var_name)
    else:
        print(f"WARNING: No {var_name} genes found in the dataset.")

# Determine percent_top based on number of genes
# Default is [50, 100, 200, 500], but we need to adjust if fewer genes
if n_genes >= 500:
    percent_top = [50, 100, 200, 500]
elif n_genes >= 200:
    percent_top = [50, 100, 200]
elif n_genes >= 100:
    percent_top = [50, 100]
elif n_genes >= 50:
    percent_top = [50]
elif n_genes >= 10:
    percent_top = [n_genes]
else:
    percent_top = []  # Use empty list instead of None

print(f"Using percent_top: {percent_top if percent_top else 'disabled'}")
print(f"Using qc_vars: {qc_vars if qc_vars else 'none'}")


def calculate_qc_manually(adata):
    """
    Manually calculate basic QC metrics when scanpy fails.
    """
    print("Calculating QC metrics manually...")

    X = adata.X
    if scipy.sparse.issparse(X):
        X_dense = X.toarray()
    else:
        X_dense = np.array(X)

    # Basic obs metrics
    adata.obs["total_counts"] = np.array(X.sum(axis=1)).flatten().astype(np.float64)
    adata.obs["n_genes_by_counts"] = np.array((X > 0).sum(axis=1)).flatten().astype(np.int64)

    # Basic var metrics
    adata.var["total_counts"] = np.array(X.sum(axis=0)).flatten().astype(np.float64)
    adata.var["n_cells_by_counts"] = np.array((X > 0).sum(axis=0)).flatten().astype(np.int64)
    adata.var["mean_counts"] = np.array(X.mean(axis=0)).flatten().astype(np.float64)

    # Calculate percentage for each gene type
    total_counts = adata.obs["total_counts"].values
    # Avoid division by zero
    total_counts_safe = np.where(total_counts == 0, 1, total_counts)

    for var_name in ["mt", "ribo", "hb"]:
        if var_name in adata.var.columns:
            var_mask = adata.var[var_name].values
            if var_mask.sum() > 0:
                var_counts = np.array(X[:, var_mask].sum(axis=1)).flatten().astype(np.float64)
            else:
                var_counts = np.zeros(adata.n_obs, dtype=np.float64)

            adata.obs[f"total_counts_{var_name}"] = var_counts
            pct = (var_counts / total_counts_safe * 100)
            # Set to 0 where total_counts was 0
            pct = np.where(total_counts == 0, 0.0, pct)
            adata.obs[f"pct_counts_{var_name}"] = pct.astype(np.float64)
        else:
            adata.obs[f"total_counts_{var_name}"] = 0.0
            adata.obs[f"pct_counts_{var_name}"] = 0.0

    print("Manual QC calculation completed.")


# Try scanpy's calculate_qc_metrics first
try:
    # scanpy requires qc_vars to be a list (can be empty) but not None
    # percent_top can be None or empty list to disable
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars if qc_vars else [],
        percent_top=percent_top if percent_top else None,
        inplace=True,
        log1p=False
    )
    print("Scanpy QC calculation completed successfully.")

except Exception as e:
    print(f"WARNING: Scanpy QC calculation failed: {e}")
    print("Falling back to manual calculation...")

    try:
        calculate_qc_manually(adata)
    except Exception as e2:
        print(f"ERROR: Manual QC calculation also failed: {e2}")
        sys.exit(1)

# Ensure all expected columns exist (add zeros if missing)
for var_name in ["mt", "ribo", "hb"]:
    if f"pct_counts_{var_name}" not in adata.obs.columns:
        adata.obs[f"pct_counts_{var_name}"] = 0.0
    if f"total_counts_{var_name}" not in adata.obs.columns:
        adata.obs[f"total_counts_{var_name}"] = 0.0

# Ensure basic columns exist
if "total_counts" not in adata.obs.columns:
    adata.obs["total_counts"] = np.array(adata.X.sum(axis=1)).flatten()
if "n_genes_by_counts" not in adata.obs.columns:
    adata.obs["n_genes_by_counts"] = np.array((adata.X > 0).sum(axis=1)).flatten()

# Print summary
print("QC metrics summary:")
print(f"  Median total counts: {adata.obs['total_counts'].median():.1f}")
print(f"  Median genes per spot: {adata.obs['n_genes_by_counts'].median():.1f}")
print(f"  Median MT %: {adata.obs['pct_counts_mt'].median():.2f}")
print(f"  Median Ribo %: {adata.obs['pct_counts_ribo'].median():.2f}")
print(f"  Median Hb %: {adata.obs['pct_counts_hb'].median():.2f}")

print(f"Total counts range: [{adata.obs['total_counts'].min():.0f}, {adata.obs['total_counts'].max():.0f}]")
print(f"Genes per spot range: [{adata.obs['n_genes_by_counts'].min():.0f}, {adata.obs['n_genes_by_counts'].max():.0f}]")

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with QC metrics to: {output_adata}")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
