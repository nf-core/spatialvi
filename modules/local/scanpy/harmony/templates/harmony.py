#!/usr/bin/env python3
"""
Integrate AnnData objects using Harmony.

Harmony is an algorithm for integrating single-cell data from multiple
batches or experiments by removing batch effects while preserving
biological variation.
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
import pandas as pd
import scanpy.external as sce
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def integrate_harmony(adata, key, basis, adjusted_basis):
    """
    Integrate observations using Harmony.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with PCA computed.
    key : str
        Column in adata.obs containing batch/sample labels.
    basis : str
        Key in adata.obsm of the embedding to integrate, _e.g._ a PCA.

    Returns
    -------
    AnnData
        AnnData with integrated embedding in obsm[adjusted_basis].
    """
    if key not in adata.obs.columns:
        raise ValueError(f"Integration key '{key}' not found in adata.obs.")

    if basis not in adata.obsm:
        raise ValueError(
            f"Embedding '{basis}' not found in adata.obsm; run PCA before "
            "integration."
        )

    n_batches = adata.obs[key].nunique()
    logger.info(f"Integrating {n_batches} batches using key: {key}")

    sce.pp.harmony_integrate(
        adata,
        key=key,
        basis=basis,
        adjusted_basis=adjusted_basis,
        random_state=0
    )

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
            "anndata": importlib.metadata.version("anndata"),
            "harmonypy": importlib.metadata.version("harmonypy"),
            "scanpy": importlib.metadata.version("scanpy"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Integrate observations in an AnnData object using Harmony."""

    # Template variables
    h5ad = "${h5ad}"
    key = "${key}"
    basis = "${basis}"
    adjusted_basis = "${embedding_added}"
    output_h5ad = "${prefix}.h5ad"
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    # `X` isn't used, so it stays on disk until the output is written
    adata = ad.read_h5ad(h5ad, backed="r")
    logger.info(f"AnnData shape: {adata.shape}")

    adata = integrate_harmony(
        adata,
        key=key,
        basis=basis,
        adjusted_basis=adjusted_basis
    )

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
