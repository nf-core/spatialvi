#!/usr/bin/env python3
"""
Compute a neighborhood graph of observations (cells/spots).
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform
import yaml

import anndata as ad
import scanpy as sc


def compute_neighbors(adata, n_neighbors, n_pcs, use_rep):
    """Compute neighborhood graph for AnnData object."""
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

    return adata


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "scanpy": importlib.metadata.version("scanpy"),
            "anndata": importlib.metadata.version("anndata"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Compute neighborhood graph for an AnnData object."""

    # Template variables
    input_adata = "${adata}"
    n_neighbors = int("${n_neighbors}")
    n_pcs = int("${n_pcs}")
    use_rep = "${use_rep}"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    # Parse optional parameters
    use_rep = None if use_rep.lower() in ["null", "none", ""] else use_rep

    # Read AnnData
    print(f"Computing neighbors for: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Compute neighbors
    adata = compute_neighbors(adata, n_neighbors, n_pcs, use_rep)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with neighbors to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
