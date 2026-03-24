#!/usr/bin/env python3
"""
Compute neighborhood enrichment by permutation test.
"""

import os
os.environ["KMP_AFFINITY"] = "disabled"
os.environ["KMP_INIT_AT_FORK"] = "FALSE"

import platform
import yaml
import squidpy as sq
import anndata as ad


def read_adata(file_path):
    """Read an AnnData object from h5ad file."""
    print(f"Reading: {file_path}")
    adata = ad.read_h5ad(file_path)
    return adata


def validate_adata(adata, cluster_key):
    """Check that required data exists in the AnnData object."""
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Column '{cluster_key}' not found in adata.obs")
    if "spatial_connectivities" not in adata.obsp:
        raise ValueError("Spatial connectivities not found. Run squidpy.gr.spatial_neighbors first.")


def write_versions(process_name):
    """Write software versions to YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "squidpy": sq.__version__,
            "anndata": ad.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Compute neighborhood enrichment by permutation test."""
    input_adata = "${adata}"
    output_adata = "${prefix}.h5ad"
    cluster_key = "${cluster_key}"

    # Read AnnData
    adata = read_adata(input_adata)
    print(f"Shape: {adata.shape}")
    print(f"Cluster key: {cluster_key}")

    # Validate input
    validate_adata(adata, cluster_key)

    # Compute neighborhood enrichment
    sq.gr.nhood_enrichment(
        adata,
        cluster_key=cluster_key,
    )

    # Print summary
    n_clusters = adata.obs[cluster_key].nunique()
    print(f"Computed neighborhood enrichment for {n_clusters} clusters")
    print(f"Results stored in adata.uns['{cluster_key}_nhood_enrichment']")

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with neighborhood enrichment to: {output_adata}")

    write_versions("${task.process}")


if __name__ == "__main__":
    main()