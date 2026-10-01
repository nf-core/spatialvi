#!/usr/bin/env python3
"""
Compute UMAP (Uniform Manifold Approximation and Projection) embedding.
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
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def compute_umap(adata, min_dist, spread, key_added):
    """Compute UMAP embedding for AnnData object."""
    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Parameters: min_dist={min_dist}, spread={spread}")

    if "neighbors" not in adata.uns:
        raise ValueError(
            "Neighbor graph not found; run scanpy.pp.neighbors first."
        )

    # Compute UMAP
    sc.tl.umap(
        adata,
        min_dist=min_dist,
        spread=spread,
        key_added=key_added
    )

    # Print summary
    logger.info(f"UMAP embedding shape: {adata.obsm[key_added].shape}")
    logger.info("UMAP coordinate ranges:")
    embedding = adata.obsm[key_added]
    logger.info(f"  UMAP1: [{embedding[:, 0].min():.2f}, {embedding[:, 0].max():.2f}]")
    logger.info(f"  UMAP2: [{embedding[:, 1].min():.2f}, {embedding[:, 1].max():.2f}]")

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
            "scanpy": importlib.metadata.version("scanpy"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Compute UMAP embedding for an AnnData object."""

    # Template variables
    h5ad = "${adata}"
    min_dist = float("${min_dist}")
    spread = float("${spread}")
    key_added = "${key_added}"
    output_adata = "${prefix}.h5ad"
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    # Read AnnData
    logger.info(f"Computing UMAP for: {h5ad}")
    adata = ad.read_h5ad(h5ad)

    # Compute UMAP
    adata = compute_umap(adata, min_dist, spread, key_added)

    # `obsm` needs an added index before writing to pickle
    df_obsm = pd.DataFrame(adata.obsm[key_added])
    df_obsm.index = adata.obs_names
    write_pickle(df_obsm, "obsm", key_added)
    write_pickle(adata.uns[key_added], "uns", key_added)

    # Write output
    if write_adata:
        adata.write_h5ad(output_adata)
        logger.info(f"Written AnnData with UMAP to: {output_adata}")

    # Write versions
    write_versions(process_name)

if __name__ == "__main__":
    main()
