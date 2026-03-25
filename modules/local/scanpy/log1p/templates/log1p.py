#!/usr/bin/env python3
"""
Apply log(1+x) transformation to the data matrix.
"""

import importlib.metadata
import platform
import yaml

import anndata as ad
import scanpy as sc


def log_transform(adata):
    """Apply log1p transformation to AnnData object."""
    print(f"Shape: {adata.shape}")

    # Apply log1p transformation
    sc.pp.log1p(adata)

    # Update uns to track transformation
    if "normalization" not in adata.uns:
        adata.uns["normalization"] = {}
    adata.uns["normalization"]["log1p"] = True

    print("Applied log(1+x) transformation")

    return adata


def write_versions(process_name):
    """Write software versions to YAML file."""
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
    input_adata = "${adata}"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    # Read AnnData
    print(f"Log-transforming AnnData from: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Transform
    adata = log_transform(adata)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written log-transformed AnnData to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
