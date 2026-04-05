#!/usr/bin/env python3
"""
Compute interaction matrix between clusters based on spatial neighbors.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform

import anndata as ad
import squidpy as sq
import yaml


def validate_adata(adata, cluster_key):
    """Check that required data exists in the AnnData object."""
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Column '{cluster_key}' not found in adata.obs")
    if "spatial_connectivities" not in adata.obsp:
        raise ValueError("Spatial connectivities not found; run squidpy.gr.spatial_neighbors first.")


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
    """Compute interaction matrix between clusters from spatial neighbors."""
    # Template variables
    h5ad = "${adata}"
    cluster_key = "${cluster_key}"
    output_h5ad = "${prefix}.h5ad"
    process_name = "${task.process}"

    print(f"Reading: {h5ad}")
    adata = ad.read_h5ad(h5ad)
    print(f"AnnData shape: {adata.shape}")
    print(f"Cluster key: {cluster_key}")

    validate_adata(adata, cluster_key)

    sq.gr.interaction_matrix(
        adata,
        cluster_key=cluster_key,
    )

    n_clusters = adata.obs[cluster_key].nunique()
    print(f"Computed interaction matrix for {n_clusters} clusters")
    print(f"Results stored in adata.uns['{cluster_key}_interactions']")

    adata.write_h5ad(output_h5ad)
    print(f"Written AnnData with interaction matrix to: {output_h5ad}")

    write_versions(process_name)


if __name__ == "__main__":
    main()
