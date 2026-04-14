#!/usr/bin/env python3
"""
Identify highly variable genes (HVGs) in the dataset.

Highly variable genes are genes that show significant variation across
observations, indicating they may be biologically relevant. These genes
are typically used for downstream dimensionality reduction and clustering.
"""

# Required for numba caching in read-only containers
import os
os.environ["NUMBA_CACHE_DIR"] = "/tmp/numba_cache"
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib"
os.environ["XDG_CACHE_HOME"] = "/tmp/cache"

import importlib.metadata
import logging
import platform

import anndata as ad
import numpy as np
import scanpy as sc
import yaml

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def validate_adata(adata):
    """
    Validate that AnnData has sufficient data for HVG selection.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    Raises
    ------
    ValueError
        If AnnData has 0 observations or 0 genes.

    Returns
    -------
    tuple
        Number of observations and genes.
    """
    n_obs, n_genes = adata.shape

    if n_obs == 0:
        raise ValueError("AnnData has 0 observations.")

    if n_genes == 0:
        raise ValueError("AnnData has 0 genes.")

    return n_obs, n_genes


def mark_all_genes_hvg(adata, flavor):
    """
    Mark all genes as highly variable when too few genes exist.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    flavor : str
        HVG selection flavor used.

    Returns
    -------
    AnnData
        AnnData with all genes marked as highly variable.
    """
    n_genes = adata.shape[1]

    logger.warning("Too few genes for meaningful HVG selection.")
    logger.info("Marking all genes as highly variable.")

    adata.var["highly_variable"] = True
    adata.var["highly_variable_rank"] = np.arange(n_genes)
    adata.var["means"] = np.array(adata.X.mean(axis=0)).flatten()
    adata.var["dispersions"] = np.zeros(n_genes)
    adata.var["dispersions_norm"] = np.zeros(n_genes)

    adata.uns["hvg"] = {
        "flavor": flavor,
        "n_hvgs_requested": n_genes,
        "n_hvgs_found": n_genes,
        "warning": f"Only {n_genes} genes available, all marked as HVG",
    }

    return adata


def find_highly_variable_genes(adata, n_top_genes, flavor):
    """
    Identify highly variable genes in the dataset.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_top_genes : int
        Number of highly variable genes to select.
    flavor : str
        Method for HVG selection (e.g., "seurat", "cell_ranger").

    Returns
    -------
    AnnData
        AnnData with HVG annotations in var.
    """
    n_obs, n_genes = validate_adata(adata)

    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"HVGs requested: {n_top_genes}")
    logger.info(f"Flavor: {flavor}")

    # Adjust n_top_genes if necessary
    if n_top_genes >= n_genes:
        logger.warning(f"Requested {n_top_genes} HVGs "
               "but only {n_genes} genes available.")
        return mark_all_genes_hvg(adata, flavor)

    try:
        sc.pp.highly_variable_genes(
            adata,
            flavor=flavor,
            n_top_genes=n_top_genes,
            inplace=True,
        )
    except ValueError as e:
        if "Bin edges must be unique" in str(e):
            logger.warning("Binning failed due to low gene variance.")
            return mark_all_genes_hvg(adata, flavor)
        raise

    adata.var["highly_variable"] = adata.var["highly_variable"].astype(bool)
    n_hvgs_found = adata.var["highly_variable"].sum()

    adata.uns["hvg"] = {
        "flavor": flavor,
        "n_hvgs_requested": n_top_genes,
        "n_hvgs_found": int(n_hvgs_found),
    }

    logger.info(f"Identified {n_hvgs_found} highly variable genes")
    logger.info(f"Percentage of genes: {n_hvgs_found / n_genes * 100:.1f}%")

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
    """Identify highly variable genes in an AnnData object."""

    # Template variables
    h5ad = "${h5ad}"
    n_top_genes = int("${n_hvgs}")
    flavor = "${flavor}"
    output_h5ad = "${prefix}.h5ad"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)
    logger.info(f"Finding highly variable genes in: {h5ad}")

    adata = find_highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor)

    adata.write_h5ad(output_h5ad)
    logger.info(f"Written AnnData with HVG annotations to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
