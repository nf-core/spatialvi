#!/usr/bin/env python3
"""
Perform Principal Component Analysis (PCA) for dimensionality reduction.

PCA reduces the dimensionality of the data while preserving the most
important variation. The results are stored in obsm["X_pca"] and are
used for downstream neighbor computation and visualization.
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
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def pca_keys(key_added):
    """
    Return the `obsm`, `varm` and `uns` keys that `sc.pp.pca` writes to: the
    defaults for `X_pca`, otherwise `key_added` for all three.
    """
    if key_added == "X_pca":
        return "X_pca", "PCs", "pca"
    return key_added, key_added, key_added


def log_variance_summary(adata, n_comps, uns_key):
    """Print summary of variance explained by principal components."""
    variance_ratio = adata.uns[uns_key]["variance_ratio"]
    cumulative_variance = variance_ratio.cumsum()

    logger.info("Variance explained:")
    for n in [10, 20, 50]:
        if n <= n_comps:
            logger.info(f"  First {n} PCs: {cumulative_variance[n - 1]:.2%}")

    logger.info(f"  All {n_comps} PCs: {cumulative_variance[-1]:.2%}")


def perform_pca(adata, n_comps, use_highly_variable, key_added):
    """
    Perform PCA on AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_comps : int
        Number of principal components to compute.
    use_highly_variable : bool
        Whether to use only highly variable genes.
    key_added : str
        Key for the results; `X_pca` keeps scanpy's default keys.

    Returns
    -------
    AnnData
        AnnData with PCA results in obsm[key_added].
    """
    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Number of components: {n_comps}")
    logger.info(f"Use highly variable genes: {use_highly_variable}")

    has_hvg = "highly_variable" in adata.var.columns

    if use_highly_variable and not has_hvg:
        raise ValueError("Highly variable genes not found in `adata.var`.")

    if use_highly_variable and has_hvg:
        n_hvgs = adata.var["highly_variable"].sum()
        logger.info(f"Using {n_hvgs} highly variable genes for PCA")

    # Without `key_added`, scanpy uses its default keys (see `pca_keys`)
    key_args = {} if key_added == "X_pca" else {"key_added": key_added}
    sc.pp.pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=use_highly_variable and has_hvg,
        random_state=0,
        **key_args,
    )

    _, _, uns_key = pca_keys(key_added)
    log_variance_summary(adata, n_comps, uns_key)

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
    """Perform PCA on an AnnData object."""

    # Template variables
    h5ad = "${adata}"
    n_comps = int("${n_pcs}")
    use_highly_variable = "${use_highly_variable}".lower() == "true"
    key_added = "${key_added}"
    output_h5ad = "${prefix}.h5ad"
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)
    logger.info(f"Performing PCA on: {h5ad}")

    adata = perform_pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=use_highly_variable,
        key_added=key_added
    )

    # `obsm` and `varm` need an added index before writing to pickle
    obsm_key, varm_key, uns_key = pca_keys(key_added)
    df_obsm = pd.DataFrame(adata.obsm[obsm_key])
    df_obsm.index = adata.obs_names
    df_varm = pd.DataFrame(adata.varm[varm_key])
    df_varm.index = adata.var_names
    write_pickle(df_obsm, "obsm", obsm_key)
    write_pickle(df_varm, "varm", varm_key)
    write_pickle(adata.uns[uns_key], "uns", uns_key)

    if write_adata:
        adata.write_h5ad(output_h5ad)
        logger.info(f"Written AnnData with PCA to: {output_h5ad}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
