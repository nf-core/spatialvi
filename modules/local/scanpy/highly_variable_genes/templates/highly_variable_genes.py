#!/usr/bin/env python3
"""
Identify highly variable genes (HVGs) in the dataset.

Highly variable genes are genes that show significant variation across
observations, indicating they may be biologically relevant. These genes
are typically used for downstream dimensionality reduction and clustering.
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
import numpy as np
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))

# The `var` columns that `sc.pp.highly_variable_genes` adds
HVG_COLUMNS = ["highly_variable", "means", "dispersions", "dispersions_norm"]
HVG_BATCH_COLUMNS = ["highly_variable_nbatches", "highly_variable_intersection"]


def hvg_columns(batch_key):
    """Return the `var` columns that HVG selection adds."""
    return HVG_COLUMNS + (HVG_BATCH_COLUMNS if batch_key else [])


def mark_all_genes_hvg(adata, flavor, batch_key):
    """
    Mark all genes as highly variable when too few genes exist.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    flavor : str
        HVG selection flavor used.
    batch_key : str
        Column in `adata.obs` with the batches, or an empty string.

    Returns
    -------
    AnnData
        AnnData with all genes marked as highly variable.
    """
    n_genes = adata.shape[1]

    logger.warning("Too few genes for meaningful HVG selection.")
    logger.info("Marking all genes as highly variable.")

    # Add the same columns and `uns` entry as scanpy does; every gene counts as
    # selected in every batch
    values = {
        "highly_variable": True,
        "means": np.array(adata.X.mean(axis=0)).flatten(),
        "dispersions": np.zeros(n_genes),
        "dispersions_norm": np.zeros(n_genes),
    }
    if batch_key:
        values["highly_variable_nbatches"] = adata.obs[batch_key].nunique()
        values["highly_variable_intersection"] = True
    for col in hvg_columns(batch_key):
        adata.var[col] = values[col]
    adata.uns["hvg"] = {"flavor": flavor}

    return adata


def log_batches(adata, batch_key):
    """Check that the batch column exists, and log the size of each batch."""
    if batch_key not in adata.obs.columns:
        raise ValueError(
            f"Batch key '{batch_key}' not found in `adata.obs`; available "
            f"columns: {', '.join(adata.obs.columns)}"
        )
    batch_sizes = adata.obs[batch_key].value_counts()
    logger.info(f'Batches in `adata.obs["{batch_key}"]`: {len(batch_sizes)}')
    for batch, n_obs in batch_sizes.items():
        logger.info(f"  {batch}: {n_obs} observations")


def find_highly_variable_genes(adata, n_top_genes, flavor, batch_key):
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
    batch_key : str
        Column in `adata.obs` to select HVGs per batch, or an empty string.

    Returns
    -------
    AnnData
        AnnData with HVG annotations in var.
    """

    # Validate AnnData
    n_obs, n_var = adata.shape
    if n_obs == 0:
        raise ValueError("AnnData has 0 observations.")
    if n_var == 0:
        raise ValueError("AnnData has 0 variables.")

    allowed_flavors = ["seurat", "cell_ranger"]
    if flavor not in allowed_flavors:
        raise ValueError(
            f"Unsupported flavor '{flavor}'; use one of: "
            f"{', '.join(allowed_flavors)}"
        )

    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"HVGs requested: {n_top_genes}")
    logger.info(f"Flavor: {flavor}")
    if batch_key:
        log_batches(adata, batch_key)

    # Adjust n_top_genes if necessary
    if n_top_genes >= n_var:
        logger.warning(
            f"Requested {n_top_genes} HVGs but only {n_var} genes available."
        )
        return mark_all_genes_hvg(adata, flavor, batch_key)

    try:
        sc.pp.highly_variable_genes(
            adata,
            flavor=flavor,
            n_top_genes=n_top_genes,
            batch_key=batch_key or None,
            inplace=True,
        )
    except ValueError as e:
        if "Bin edges must be unique" in str(e):
            logger.warning("Binning failed due to low gene variance.")
            return mark_all_genes_hvg(adata, flavor, batch_key)
        raise

    adata.var["highly_variable"] = adata.var["highly_variable"].astype(bool)
    n_hvgs_found = adata.var["highly_variable"].sum()

    logger.info(f"Identified {n_hvgs_found} highly variable genes")
    logger.info(f"Percentage of genes: {n_hvgs_found / n_var * 100:.1f}%")

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
    """Identify highly variable genes in an AnnData object."""

    # Template variables
    h5ad = "${h5ad}"
    n_top_genes = int("${n_hvgs}")
    flavor = "${flavor}"
    batch_key = "${batch_key}"
    output_h5ad = "${prefix}.h5ad"
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)
    logger.info(f"Finding highly variable genes in: {h5ad}")

    adata = find_highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        batch_key=batch_key
    )

    write_pickle(
        adata.var[hvg_columns(batch_key)],
        "var",
        "highly_variable_genes"
    )
    write_pickle(adata.uns["hvg"], "uns", "hvg")

    if write_adata:
        adata.write_h5ad(output_h5ad)
        logger.info(f"Written AnnData with HVG annotations to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
