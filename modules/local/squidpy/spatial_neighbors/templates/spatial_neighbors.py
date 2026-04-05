#!/usr/bin/env python3
"""
Compute spatial neighbors graph based on spatial coordinates.

Creates a spatial connectivity graph where observations are connected
to their nearest neighbors in physical space. Results are stored in
adata.obsp as sparse matrices.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform

import anndata as ad
import squidpy as sq
import yaml


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
            "squidpy": importlib.metadata.version("squidpy")
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Compute spatial neighbors graph based on spatial coordinates."""
    # Template variables
    h5ad = "${adata}"
    coord_type = "${coord_type}"
    n_neighs = int("${n_neighs}")
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    print(f"Reading: {h5ad}")
    adata = ad.read_h5ad(h5ad)
    print(f"AnnData shape: {adata.shape}")
    print(f"Coord type: {coord_type}")
    print(f"Number of neighbors: {n_neighs}")

    if "spatial" not in adata.obsm:
        raise ValueError(
            "Spatial coordinates not found in adata.obsm['spatial']"
        )

    sq.gr.spatial_neighbors(adata, coord_type=coord_type, n_neighs=n_neighs)

    print("Computed spatial neighbor graph")
    print(f"Connectivities shape: {adata.obsp['spatial_connectivities'].shape}")
    print(f"Distances shape: {adata.obsp['spatial_distances'].shape}")

    adata.write_h5ad(output_adata)
    print(f"Written AnnData with spatial neighbors to: {output_adata}")

    write_versions(process_name)


if __name__ == "__main__":
    main()
