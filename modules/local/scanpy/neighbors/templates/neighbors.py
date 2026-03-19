#!/usr/bin/env python3
"""
Compute a neighborhood graph of observations (cells/spots).
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
n_neighbors = int("${n_neighbors}")
n_pcs_str = "${n_pcs}"
use_rep_str = "${use_rep}"

# Parse optional parameters
n_pcs = None if n_pcs_str.lower() in ["null", "none", ""] else int(n_pcs_str)
use_rep = None if use_rep_str.lower() in ["null", "none", ""] else use_rep_str

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Computing neighbors for: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Number of neighbors: {n_neighbors}")
print(f"Number of PCs: {n_pcs if n_pcs else 'default'}")
print(f"Use representation: {use_rep if use_rep else 'X_pca'}")

# Compute neighbors
sc.pp.neighbors(
    adata,
    n_neighbors=n_neighbors,
    n_pcs=n_pcs,
    use_rep=use_rep
)

# Print summary
print("Computed neighbor graph")
print(f"  Connectivities shape: {adata.obsp['connectivities'].shape}")
print(f"  Distances shape: {adata.obsp['distances'].shape}")

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with neighbors to: {output_adata}")

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
