#!/usr/bin/env python3
"""
Normalize total counts per observation using scanpy.

Normalizes each observation to have the same total count after normalization.
By default, uses median total counts as the target sum.
"""

import importlib.metadata
import logging
import platform

import anndata as ad
import scanpy as sc
import yaml

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

def normalize_adata(adata, target_sum):
    """
    Normalize total counts per observation.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    target_sum : int or None
        Target total counts per observation. If None, uses median.

    Returns
    -------
    AnnData
        Normalized AnnData.
    """
    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Target sum: {target_sum if target_sum else 'median of total counts'}")

    # Calculate pre-normalization statistics
    if "total_counts" in adata.obs:
        median_before = adata.obs["total_counts"].median()
        logger.info(f"Median total counts before normalization: {median_before:.2f}")

    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)

    # Calculate post-normalization statistics
    if "total_counts" in adata.obs:
        if hasattr(adata.X, 'A1'):
            adata.obs["total_counts_normalized"] = adata.X.sum(axis=1).A1
        else:
            adata.obs["total_counts_normalized"] = adata.X.sum(axis=1)
        median_counts_after = adata.obs["total_counts_normalized"].median()
        logger.info(f"Median total counts after normalization: {median_counts_after:.2f}")

    adata.uns["normalization"] = {
        "method": "normalize_total",
        "target_sum": target_sum if target_sum else "median",
    }

    return adata

def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
            "scanpy": importlib.metadata.version("scanpy"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)

def main():
    """Normalize total counts in an AnnData object."""

    # Template variables
    h5ad = "${adata}"
    target_sum = None if int("${target_sum}") == 0 else int("${target_sum}")
    output_h5ad = "${prefix}.h5ad"
    process_name = "${task.process}"

    logger.info(f"Normalizing AnnData from: {h5ad}")
    adata = ad.read_h5ad(h5ad)

    adata = normalize_adata(adata, target_sum=target_sum)

    adata.write_h5ad(output_h5ad)
    logger.info(f"Written normalized AnnData to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
