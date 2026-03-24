#!/usr/bin/env python3
"""
Compute spatial neighbors graph based on spatial coordinates.
"""

import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import yaml
import squidpy as sq
import anndata as ad


def read_adata(file_path):
    """Read an AnnData object from h5ad file."""
    print(f"Reading: {file_path}")
    adata = ad.read_h5ad(file_path)
    return adata


def validate_adata(adata):
    """Check that required data exists in the AnnData object."""
    if "spatial" not in adata.obsm:
        raise ValueError("Spatial coordinates not found in adata.obsm['spatial']")


def write_versions(process_name):
    """Write software versions to YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "squidpy": sq.__version__,
            "anndata": ad.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Compute spatial neighbors graph based on spatial coordinates."""
    input_adata = "${adata}"
    output_adata = "${prefix}.h5ad"
    coord_type = "${coord_type}"
    n_neighs = int("${n_neighs}")

    # Read AnnData
    adata = read_adata(input_adata)
    print(f"Shape: {adata.shape}")
    print(f"Coord type: {coord_type}")
    print(f"Number of neighbors: {n_neighs}")

    # Validate input
    validate_adata(adata)

    # Compute spatial neighbors
    sq.gr.spatial_neighbors(
        adata,
        coord_type=coord_type,
        n_neighs=n_neighs,
    )

    # Print summary
    print("Computed spatial neighbor graph")
    print(f"  Connectivities shape: {adata.obsp['spatial_connectivities'].shape}")
    print(f"  Distances shape: {adata.obsp['spatial_distances'].shape}")

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with spatial neighbors to: {output_adata}")

    write_versions("${task.process}")


if __name__ == "__main__":
    main()