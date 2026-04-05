#!/usr/bin/env python3
"""
Perform Leiden clustering on the neighbor graph.

Leiden is a community detection algorithm that identifies clusters of
observations based on a pre-computed neighbor graph. Results are stored
in adata.obs.
"""

import importlib.metadata
import platform

import anndata as ad
import scanpy as sc
import yaml


def perform_leiden(adata, resolution, key_added):
    """
    Perform Leiden clustering on AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with neighbor graph computed.
    resolution : float
        Resolution parameter for clustering (higher = more clusters).
    key_added : str
        Key in adata.obs to store cluster labels.

    Returns
    -------
    AnnData
        AnnData with cluster labels in obs.
    """
    if "neighbors" not in adata.uns:
        raise ValueError("Neighbor graph not found; run sc.pp.neighbors first.")

    print(f"AnnData shape: {adata.shape}")
    print(f"Resolution: {resolution}")
    print(f"Key added: {key_added}")

    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)

    n_clusters = adata.obs[key_added].nunique()
    adata.uns["leiden"] = {
        "resolution": resolution,
        "n_clusters": n_clusters,
    }

    cluster_sizes = adata.obs[key_added].value_counts().sort_index()
    print(f"Found {n_clusters} clusters:")
    for cluster, size in cluster_sizes.items():
        pct = size / adata.shape[0] * 100
        print(f"  Cluster {cluster}: {size} obs ({pct:.1f}%)")

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
    h5ad = "${h5ad}"
    resolution = float("${resolution}")
    key_added = "${key_added}"
    output_h5ad = "${prefix}.h5ad"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)
    print(f"Performing Leiden clustering on: {h5ad}")

    adata = perform_leiden(adata, resolution=resolution, key_added=key_added)

    adata.write_h5ad(output_h5ad)
    print(f"Written AnnData with clusters to: {output_h5ad}")

    write_versions(process_name)


if __name__ == "__main__":
    main()
