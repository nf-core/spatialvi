#!/usr/bin/env python3
"""
Perform Principal Component Analysis (PCA) for dimensionality reduction.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import importlib.metadata
import yaml
import scanpy as sc
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"
n_comps = int("${n_principal_components}")
use_highly_variable = "${pca_use_highly_variable}".lower() == "true"

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Performing PCA on: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Number of components: {n_comps}")
print(f"Use highly variable genes: {use_highly_variable}")

# Check if HVGs are available
has_hvg = "highly_variable" in adata.var.columns
if use_highly_variable and has_hvg:
    n_hvgs = adata.var["highly_variable"].sum()
    print(f"Using {n_hvgs} highly variable genes for PCA")
elif use_highly_variable and not has_hvg:
    print("Warning: highly_variable not found in var, using all genes")
    use_highly_variable = False

# Perform PCA
sc.pp.pca(
    adata,
    n_comps=n_comps,
    use_highly_variable=use_highly_variable if has_hvg else False
)

# Print variance explained
variance_ratio = adata.uns["pca"]["variance_ratio"]
cumulative_variance = variance_ratio.cumsum()
print(f"Variance explained by first 10 PCs: {cumulative_variance[9]:.2%}")
print(f"Variance explained by first 20 PCs: {cumulative_variance[19]:.2%}")
print(f"Variance explained by all {n_comps} PCs: {cumulative_variance[-1]:.2%}")

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with PCA to: {output_adata}")

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
