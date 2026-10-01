#!/usr/bin/env python3
"""
Integrate AnnData objects using Scanorama.

Scanorama is an algorithm for integrating single-cell data from multiple
batches by identifying and merging shared cell types across datasets.
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
import pandas as pd
import scanpy.external as sce
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def integrate_scanorama(adata, key, adjusted_basis):
    """
    Integrate observations using Scanorama.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    key : str
        Column in adata.obs containing batch/sample labels.
    adjusted_basis : str
        Name of the obsm key to store the integrated embedding.

    Returns
    -------
    AnnData
        AnnData with integrated embedding in obsm.
    """
    if key not in adata.obs.columns:
        raise ValueError(f"Integration key '{key}' not found in adata.obs.")

    if "X_pca" not in adata.obsm:
        raise ValueError(
            "PCA not found in adata.obsm; run PCA before integration."
        )

    n_batches = adata.obs[key].nunique()
    logger.info(f"Integrating {n_batches} batches using key: {key}")

    sce.pp.scanorama_integrate(adata, key=key, adjusted_basis=adjusted_basis)

    return adata


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
            "scanorama": importlib.metadata.version("scanorama"),
            "scanpy": importlib.metadata.version("scanpy"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Integrate observations in an AnnData object using Scanorama."""

    # Template variables
    h5ad = "${h5ad}"
    key = "${key}"
    adjusted_basis = "${embedding_added}"
    output_h5ad = "${prefix}.h5ad"
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)
    logger.info(f"Read AnnData object: {h5ad}")
    logger.info(f"AnnData shape: {adata.shape}")

    adata = integrate_scanorama(adata, key=key, adjusted_basis=adjusted_basis)

    # `obsm` needs an added index before writing to pickle
    df_obsm = pd.DataFrame(adata.obsm[adjusted_basis])
    df_obsm.index = adata.obs_names
    write_pickle(df_obsm, "obsm", adjusted_basis)

    if write_adata:
        adata.write_h5ad(output_h5ad)
        logger.info(f"Written integrated AnnData to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
