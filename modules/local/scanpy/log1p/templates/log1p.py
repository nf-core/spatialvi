#!/usr/bin/env python3
"""
Apply log(1+x) transformation to the data matrix.
"""

import importlib.metadata
import logging
import platform

import anndata as ad
import scanpy as sc
import yaml

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def log_transform(adata):
    """
    Apply log1p transformation to AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (typically normalized counts).

    Returns
    -------
    AnnData
        Log-transformed AnnData.
    """
    logger.info(f"AnnData shape: {adata.shape}")

    sc.pp.log1p(adata)

    if "normalization" not in adata.uns:
        adata.uns["normalization"] = {}
    adata.uns["normalization"]["log1p"] = True

    logger.info("Applied log(1+x) transformation")

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
    """Apply log1p transformation to an AnnData object."""

    # Template variables
    h5ad = "${adata}"
    output_h5ad = "${prefix}.h5ad"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)
    logger.info(f"Log-transforming AnnData from: {h5ad}")

    adata = log_transform(adata)

    adata.write_h5ad(output_h5ad)
    logger.info(f"Written log-transformed AnnData to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
