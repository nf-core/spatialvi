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

import importlib.metadata
import logging
import platform

import anndata as ad
import scanpy.external as sce
import yaml

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def integrate_harmony(adata, key, adjusted_basis):
    """
    Integrate observations using Harmony.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with PCA computed.
    key : str
        Column in adata.obs containing batch/sample labels.

    Returns
    -------
    AnnData
        AnnData with integrated embedding in obsm[adjusted_basis].
    """
    if key not in adata.obs.columns:
        raise ValueError(f"Integration key '{key}' not found in adata.obs.")

    if "X_pca" not in adata.obsm:
        raise ValueError(
            "PCA not found in adata.obsm; run PCA before integration."
        )

    n_batches = adata.obs[key].nunique()
    logger.info(f"Integrating {n_batches} batches using key: {key}")

    sce.pp.harmony_integrate(adata, key=key, adjusted_basis=adjusted_basis)

    return adata


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
    adjusted_basis = "${embedding_added}"
    output_h5ad = "${prefix}.h5ad"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)
    logger.info(f"AnnData shape: {adata.shape}")

    adata = integrate_harmony(adata, key=key, adjusted_basis=adjusted_basis)

    adata.write_h5ad(output_h5ad)
    logger.info(f"Written integrated AnnData to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
