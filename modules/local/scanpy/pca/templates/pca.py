#!/usr/bin/env python3
"""
Perform Principal Component Analysis (PCA) for dimensionality reduction.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform
import yaml

import anndata as ad
import scanpy as sc


def perform_pca(adata, n_comps, use_highly_variable):
    """Perform PCA on AnnData object."""
    print(f"Shape: {adata.shape}")
    print(f"Number of components: {n_comps}")
    print(f"Use highly variable genes: {use_highly_variable}")

    # Check if HVGs are available
    has_hvg = "highly_variable" in adata.var.columns

    if use_highly_variable and has_hvg:
        n_hvgs = adata.var["highly_variable"].sum()
        print(f"Using {n_hvgs} highly variable genes for PCA")
    elif use_highly_variable and not has_hvg:
        print("Warning: highly_variable not found in var, using all genes")
        use_highly_variable = False

    # Perform PCA
    sc.pp.pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=use_highly_variable if has_hvg else False
    )

    # Print summary of variance explained
    variance_ratio = adata.uns["pca"]["variance_ratio"]
    cumulative_variance = variance_ratio.cumsum()
    print(f"Variance explained by first 10 PCs: {cumulative_variance[9]:.2%}")
    print(f"Variance explained by first 20 PCs: {cumulative_variance[19]:.2%}")
    print(f"Variance explained by all {n_comps} PCs: {cumulative_variance[-1]:.2%}")

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
    """Perform PCA on an AnnData object."""

    # Template variables
    input_adata = "${adata}"
    n_comps = int("${n_pcs}")
    use_highly_variable = "${use_highly_variable}".lower() == "true"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    # Read AnnData
    print(f"Performing PCA on: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Perform PCA
    adata = perform_pca(adata, n_comps, use_highly_variable)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with PCA to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
