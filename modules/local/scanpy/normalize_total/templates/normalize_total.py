#!/usr/bin/env python3
"""
Normalize total counts per cell/spot using scanpy.

Normalizes each cell/spot to have the same total count after normalization.
By default, uses median total counts as the target sum.
"""

import platform
import importlib.metadata
import yaml
import scanpy as sc
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"
target_sum_str = "${target_sum}"

# Parse target_sum parameter
if target_sum_str.lower() in ["null", "none", ""]:
    target_sum = None
else:
    target_sum = float(target_sum_str)

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Normalizing AnnData from: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Target sum: {target_sum if target_sum else 'median of total counts'}")

# Store pre-normalization statistics
median_counts_before = adata.obs["total_counts"].median() if "total_counts" in adata.obs else None
if median_counts_before:
    print(f"Median total counts before normalization: {median_counts_before:.2f}")

# Normalize total counts
sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)

# Store post-normalization statistics
if "total_counts" in adata.obs:
    # Recalculate after normalization
    adata.obs["total_counts_normalized"] = adata.X.sum(axis=1).A1 if hasattr(adata.X, 'A1') else adata.X.sum(axis=1)
    median_counts_after = adata.obs["total_counts_normalized"].median()
    print(f"Median total counts after normalization: {median_counts_after:.2f}")

# Store normalization info in uns
adata.uns["normalization"] = {
    "method": "normalize_total",
    "target_sum": target_sum if target_sum else "median",
}

# Write output
adata.write_h5ad(output_adata)
print(f"Written normalized AnnData to: {output_adata}")

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