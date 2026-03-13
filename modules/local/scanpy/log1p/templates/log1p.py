#!/usr/bin/env python3
"""
Apply log(1+x) transformation to the data matrix.
"""

import platform
import importlib.metadata
import yaml
import scanpy as sc
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Log-transforming AnnData from: {input_adata}")
print(f"Shape: {adata.shape}")

# Apply log1p transformation
sc.pp.log1p(adata)

# Update uns to track transformation
if "normalization" not in adata.uns:
    adata.uns["normalization"] = {}
adata.uns["normalization"]["log1p"] = True

print("Applied log(1+x) transformation")

# Write output
adata.write_h5ad(output_adata)
print(f"Written log-transformed AnnData to: {output_adata}")

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