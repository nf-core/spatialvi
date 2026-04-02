#!/usr/bin/env python3
"""
Perform Leiden clustering on the neighbor graph.
"""

import importlib.metadata
import platform
import yaml

import anndata as ad
import scanpy as sc


def perform_leiden(adata, resolution, key_added):
    """Perform Leiden clustering on AnnData object."""
    print(f"Shape: {adata.shape}")
    print(f"Resolution: {resolution}")
    print(f"Key added: {key_added}")

    # Validate input
    if "neighbors" not in adata.uns:
        raise ValueError("Neighbor graph not found. Run scanpy.pp.neighbors first.")

    # Perform Leiden clustering
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=key_added
    )

    # Store clustering parameters in uns
    n_clusters = adata.obs[key_added].nunique()
    adata.uns["leiden"] = {
        "resolution": resolution,
        "n_clusters": n_clusters,
    }

    # Print summary
    n_clusters = adata.obs[key_added].nunique()
    cluster_sizes = adata.obs[key_added].value_counts().sort_index()
    print(f"Found {n_clusters} clusters")
    print("Cluster sizes:")
    for cluster, size in cluster_sizes.items():
        print(f"  Cluster {cluster}: {size} spots ({size/adata.shape[0]*100:.1f}%)")

    return adata


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "scanpy": importlib.metadata.version("scanpy"),
            "anndata": importlib.metadata.version("anndata"),
            "leidenalg": importlib.metadata.version("leidenalg"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Perform Leiden clustering on an AnnData object."""

    # Template variables
    input_adata = "${adata}"
    resolution = float("${resolution}")
    key_added = "${key_added}"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    # Read AnnData
    print(f"Performing Leiden clustering on: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Perform clustering
    adata = perform_leiden(adata, resolution, key_added)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with clusters to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
