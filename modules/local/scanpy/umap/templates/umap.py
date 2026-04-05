#!/usr/bin/env python3
"""
Compute UMAP (Uniform Manifold Approximation and Projection) embedding.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform

import anndata as ad
import scanpy as sc
import yaml


def compute_umap(adata, min_dist, spread, key_added):
    """Compute UMAP embedding for AnnData object."""
    print(f"AnnData shape: {adata.shape}")
    print(f"Parameters: min_dist={min_dist}, spread={spread}")

    if "neighbors" not in adata.uns:
        raise ValueError(
            "Neighbor graph not found; run scanpy.pp.neighbors first."
        )

    # Compute UMAP
    sc.tl.umap(
        adata,
        min_dist=min_dist,
        spread=spread,
        key_added=key_added
    )

    # Print summary
    print(f"UMAP embedding shape: {adata.obsm[key_added].shape}")
    print("UMAP coordinate ranges:")
    embedding = adata.obsm[key_added]
    print(f"  UMAP1: [{embedding[:, 0].min():.2f}, {embedding[:, 0].max():.2f}]")
    print(f"  UMAP2: [{embedding[:, 1].min():.2f}, {embedding[:, 1].max():.2f}]")

    return adata


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
            "scanpy": importlib.metadata.version("scanpy"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Compute UMAP embedding for an AnnData object."""
    # Template variables
    h5ad = "${adata}"
    min_dist = float("${min_dist}")
    spread = float("${spread}")
    key_added = "${key_added}"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    # Read AnnData
    print(f"Computing UMAP for: {h5ad}")
    adata = ad.read_h5ad(h5ad)

    # Compute UMAP
    adata = compute_umap(adata, min_dist, spread, key_added)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with UMAP to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
