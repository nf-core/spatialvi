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
import platform

import anndata as ad
import scanpy as sc
import yaml

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

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
        Number of principal components to use.
    use_rep : str or None
        Representation to use. If None, uses either `.X` when `.n_vars < 50` or
        `X_pca` otherwise.

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
        use_rep=use_rep
    )

    logger.info("Computed neighbor graph:")
    logger.info(f"  Connectivities shape: {adata.obsp['connectivities'].shape}")
    logger.info(f"  Distances shape: {adata.obsp['distances'].shape}")

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
    h5ad = "${adata}"
    n_neighbors = int("${n_neighbors}")
    n_pcs = int("${n_pcs}")
    use_rep = None if "${use_rep}".lower() in ["none", ""] else "${use_rep}"
    output_h5ad = "${prefix}.h5ad"
    process_name = "${task.process}"

    logger.info(f"Computing neighbors for: {h5ad}")
    adata = ad.read_h5ad(h5ad)

    adata = compute_neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep
    )

    adata.write_h5ad(output_h5ad)
    logger.info(f"Written AnnData with neighbors to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
