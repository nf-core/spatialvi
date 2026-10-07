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
import squidpy as sq
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


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
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    logger.info(f"Reading: {h5ad}")
    # `X` isn't used, so it stays on disk until the output is written
    adata = ad.read_h5ad(h5ad, backed="r")
    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Coord type: {coord_type}")
    logger.info(f"Number of neighbors: {n_neighs}")

    if "spatial" not in adata.obsm:
        raise ValueError(
            "Spatial coordinates not found in adata.obsm['spatial']"
        )

    sq.gr.spatial_neighbors(adata, coord_type=coord_type, n_neighs=n_neighs)

    logger.info("Computed spatial neighbor graph")
    logger.info(f"Connectivities shape: {adata.obsp['spatial_connectivities'].shape}")
    logger.info(f"Distances shape: {adata.obsp['spatial_distances'].shape}")

    # Store `obsp` sparse matrices alongside an index in a dictionary
    for name in ["spatial_connectivities", "spatial_distances"]:
        obsp_dict = {"matrix": adata.obsp[name], "index": adata.obs_names}
        write_pickle(obsp_dict, "obsp", name)
    write_pickle(adata.uns["spatial_neighbors"], "uns", "spatial_neighbors")

    if write_adata:
        adata.write_h5ad(output_adata)
        logger.info(f"Written AnnData with spatial neighbors to: {output_adata}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
