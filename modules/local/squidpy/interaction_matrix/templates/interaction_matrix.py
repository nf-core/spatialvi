#!/usr/bin/env python3
"""
Compute interaction matrix between clusters based on spatial neighbors.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import logging
import pickle
import platform
from pathlib import Path

import anndata as ad
import squidpy as sq
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def validate_adata(adata, cluster_key):
    """Check that required data exists in the AnnData object."""
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Column '{cluster_key}' not found in adata.obs")
    if "spatial_connectivities" not in adata.obsp:
        raise ValueError("Spatial connectivities not found; run squidpy.gr.spatial_neighbors first.")


def write_pickle(data, slot, name):
    """Write data to a `<slot>/<name>.pkl` pickle file."""
    Path(slot).mkdir(exist_ok=True)
    with open(f"{slot}/{name}.pkl", "wb") as f:
        pickle.dump(data, f, protocol=5)


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
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    logger.info(f"Reading: {h5ad}")
    adata = ad.read_h5ad(h5ad)
    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Cluster key: {cluster_key}")

    validate_adata(adata, cluster_key)

    sq.gr.interaction_matrix(
        adata,
        cluster_key=cluster_key,
    )

    n_clusters = adata.obs[cluster_key].nunique()
    logger.info(f"Computed interaction matrix for {n_clusters} clusters")
    logger.info(f"Results stored in adata.uns['{cluster_key}_interactions']")

    uns_key = f"{cluster_key}_interactions"
    write_pickle(adata.uns[uns_key], "uns", uns_key)

    if write_adata:
        adata.write_h5ad(output_h5ad)
        logger.info(f"Written AnnData with interaction matrix to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
