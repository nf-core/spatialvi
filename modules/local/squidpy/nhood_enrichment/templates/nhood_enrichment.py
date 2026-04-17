#!/usr/bin/env python3
"""
Compute neighborhood enrichment by permutation test.

Neighborhood enrichment analysis tests whether clusters are spatially
co-localized more or less frequently than expected by chance.
Results are stored in adata.uns.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"
os.environ["KMP_INIT_AT_FORK"] = "FALSE"

import importlib.metadata
import logging
import platform

import anndata as ad
import squidpy as sq
import yaml

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def validate_adata(adata, cluster_key):
    """Check that required data exists in the AnnData object."""
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Column '{cluster_key}' not found in adata.obs")
    if "spatial_connectivities" not in adata.obsp:
        raise ValueError(
            "Spatial connectivities not found; "
            "run squidpy.gr.spatial_neighbors first."
        )


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
            "squidpy": importlib.metadata.version("squidpy")
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Compute neighborhood enrichment by permutation test."""

    # Template variables
    h5ad = "${adata}"
    cluster_key = "${cluster_key}"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    logger.info(f"Reading: {h5ad}")
    adata = ad.read_h5ad(h5ad)
    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Cluster key: {cluster_key}")

    validate_adata(adata, cluster_key)

    sq.gr.nhood_enrichment(
        adata,
        cluster_key=cluster_key,
    )

    n_clusters = adata.obs[cluster_key].nunique()
    logger.info(f"Computed neighborhood enrichment for {n_clusters} clusters")
    logger.info(f"Results stored in adata.uns['{cluster_key}_nhood_enrichment']")

    adata.write_h5ad(output_adata)
    logger.info(f"Written AnnData with neighborhood enrichment to: {output_adata}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
