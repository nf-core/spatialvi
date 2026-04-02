#!/usr/bin/env python3
"""
Compute UMAP (Uniform Manifold Approximation and Projection) embedding.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform
import yaml

import anndata as ad
import scanpy as sc


def validate_neighbors(adata):
    """Validate that neighbor graph has been computed."""
    if "neighbors" not in adata.uns:
        raise ValueError("Neighbor graph not found. Run scanpy.pp.neighbors first.")


def compute_umap(adata, min_dist, spread, key_added):
    """Compute UMAP embedding for AnnData object."""
    print(f"Shape: {adata.shape}")
    print(f"Parameters: min_dist={min_dist}, spread={spread}")

    # Validate input
    validate_neighbors(adata)

    # Compute UMAP
    sc.tl.umap(
        adata,
        min_dist=min_dist,
        spread=spread,
        key_added=key_added
    )

    # Print summary
    print(f"UMAP embedding shape: {adata.obsm['X_umap'].shape}")
    print("UMAP coordinate ranges:")
    print(f"  UMAP1: [{adata.obsm['X_umap'][:, 0].min():.2f}, {adata.obsm['X_umap'][:, 0].max():.2f}]")
    print(f"  UMAP2: [{adata.obsm['X_umap'][:, 1].min():.2f}, {adata.obsm['X_umap'][:, 1].max():.2f}]")

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
    """Compute UMAP embedding for an AnnData object."""

    # Template variables
    input_adata = "${adata}"
    min_dist = float("${min_dist}")
    spread = float("${spread}")
    key_added = "${key_added}"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    # Read AnnData
    print(f"Computing UMAP for: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Compute UMAP
    adata = compute_umap(adata, min_dist, spread, key_added)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with UMAP to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
