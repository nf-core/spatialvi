#!/usr/bin/env python3
"""
Compute a neighborhood graph of observations.

The neighborhood graph is the basis for clustering and UMAP visualization.
It connects each observation to its nearest neighbors in the specified
representation space (typically PCA).
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
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def validate_representation(adata, use_rep):
    """Require an explicit, existing representation for the neighbor search."""
    if use_rep.lower() in ["", "none"]:
        raise ValueError(
            "`use_rep` is required: use 'X' for the data matrix or a key in "
            "`adata.obsm` (e.g. 'X_pca')"
        )
    if use_rep != "X" and use_rep not in adata.obsm:
        available = ", ".join(adata.obsm.keys()) or "none"
        raise ValueError(
            f"Representation '{use_rep}' not found in `adata.obsm` "
            f"(available: {available})"
        )


def compute_neighbors(adata, n_neighbors, n_pcs, use_rep):
    """
    Compute neighborhood graph for AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with PCA or other representation computed.
    n_neighbors : int
        Number of neighbors to use.
    n_pcs : int
        Number of dimensions of the representation to use (the first `n_pcs`
        columns); ignored when `use_rep` is 'X'.
    use_rep : str
        Representation to use: a key in `adata.obsm` (e.g. 'X_pca' or
        'X_harmony'), or 'X' to use the data matrix directly.

    Returns
    -------
    AnnData
        AnnData with neighbor graph in obsp.
    """
    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Number of neighbors: {n_neighbors}")
    logger.info(f"Number of PCs: {n_pcs}")
    logger.info(f"Representation: {use_rep}")

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
        random_state=0
    )

    logger.info("Computed neighbor graph:")
    logger.info(f"  Connectivities shape: {adata.obsp['connectivities'].shape}")
    logger.info(f"  Distances shape: {adata.obsp['distances'].shape}")

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
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Compute neighborhood graph for an AnnData object."""

    # Template variables
    h5ad = "${adata}"
    n_neighbors = int("${n_neighbors}")
    n_pcs = int("${n_pcs}")
    use_rep = "${use_rep}"
    output_h5ad = "${prefix}.h5ad"
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    logger.info(f"Reading: {h5ad}")
    adata = ad.read_h5ad(h5ad)

    validate_representation(adata, use_rep)

    logger.info(f"Computing neighbors for: {h5ad}")
    adata = compute_neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep
    )

    # Store `obsp` sparse matrices alongside an index in a dictionary
    for name in ["connectivities", "distances"]:
        obsp_dict = {"matrix": adata.obsp[name], "index": adata.obs_names}
        write_pickle(obsp_dict, "obsp", name)
    write_pickle(adata.uns["neighbors"], "uns", "neighbors")

    if write_adata:
        adata.write_h5ad(output_h5ad)
        logger.info(f"Written AnnData with neighbors to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
