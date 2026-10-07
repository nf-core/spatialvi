#!/usr/bin/env python3
"""
Perform Leiden clustering on the neighbor graph.

Leiden is a community detection algorithm that identifies clusters of
observations based on a pre-computed neighbor graph. Results are stored
in adata.obs.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

# Keep caches in the task's work directory, which is always writable and
# private to the task
os.environ["NUMBA_CACHE_DIR"] = os.path.join(os.getcwd(), ".cache", "numba")
os.environ["MPLCONFIGDIR"] = os.path.join(os.getcwd(), ".cache", "matplotlib")
os.environ["XDG_CACHE_HOME"] = os.path.join(os.getcwd(), ".cache")

import importlib.metadata
import logging
import pickle
import platform
from pathlib import Path

import anndata as ad
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


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

    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Resolution: {resolution}")
    logger.info(f"Key added: {key_added}")

    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=key_added,
        flavor="leidenalg",
        random_state=0
    )

    n_clusters = adata.obs[key_added].nunique()
    cluster_sizes = adata.obs[key_added].value_counts().sort_index()
    logger.info(f"Found {n_clusters} clusters:")
    for cluster, size in cluster_sizes.items():
        pct = size / adata.shape[0] * 100
        logger.info(f"  Cluster {cluster}: {size} obs ({pct:.1f}%)")

    return adata


def write_pickle(data, slot, name):
    """Write data to a `<slot>/<name>.pkl` pickle file."""
    Path(slot).mkdir(exist_ok=True)
    with open(f"{slot}/{name}.pkl", "wb") as f:
        pickle.dump(data, f, protocol=5)
    logger.info(f"Written slot data to: {slot}/{name}.pkl")


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
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    # `X` isn't used, so it stays on disk until the output is written
    adata = ad.read_h5ad(h5ad, backed="r")
    logger.info(f"Performing Leiden clustering on: {h5ad}")

    adata = perform_leiden(adata, resolution=resolution, key_added=key_added)

    write_pickle(adata.obs[[key_added]], "obs", key_added)
    write_pickle(adata.uns[key_added], "uns", key_added)

    if write_adata:
        adata.write_h5ad(output_h5ad)
        logger.info(f"Written AnnData with clusters to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
