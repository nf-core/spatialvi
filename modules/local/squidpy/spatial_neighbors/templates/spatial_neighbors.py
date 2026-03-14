#!/usr/bin/env python3
"""
Compute spatial neighbors graph based on spatial coordinates.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import importlib.metadata
import yaml
import squidpy as sq
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"
coord_type = "${coord_type}"
n_neighbors = int("${n_neighbors}")

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Computing spatial neighbors for: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Coord type: {coord_type}")
print(f"Number of neighbors: {n_neighbors}")

# Check if spatial coordinates exist
if "spatial" not in adata.obsm:
    raise ValueError("Spatial coordinates not found in adata.obsm['spatial']")

# Compute spatial neighbors
sq.gr.spatial_neighbors(
    adata,
    coord_type=coord_type,
    n_neighs=n_neighbors
)

# Print summary
print("Computed spatial neighbor graph")
print(f"  Connectivities shape: {adata.obsp['spatial_connectivities'].shape}")
print(f"  Distances shape: {adata.obsp['spatial_distances'].shape}")

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with spatial neighbors to: {output_adata}")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "squidpy": importlib.metadata.version("squidpy"),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
