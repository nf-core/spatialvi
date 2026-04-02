#!/usr/bin/env python3
"""
Identify highly variable genes (HVGs) in the dataset.
"""

# Fix numba caching issue in read-only containers
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp/numba_cache'
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'
os.environ['XDG_CACHE_HOME'] = '/tmp/cache'

import importlib.metadata
import platform
import sys
import yaml

import anndata as ad
import numpy as np
import scanpy as sc


def validate_adata(adata):
    """Validate that AnnData has sufficient data for HVG selection."""
    n_cells, n_genes = adata.shape

    if n_cells == 0:
        print("ERROR: AnnData has 0 cells/spots.")
        sys.exit(1)

    if n_genes == 0:
        print("ERROR: AnnData has 0 genes.")
        sys.exit(1)

    return n_cells, n_genes


def mark_all_genes_hvg(adata, n_genes, flavor, n_requested):
    """Mark all genes as highly variable when too few genes exist."""
    print(f"WARNING: AnnData has only {n_genes} genes. This is too few for meaningful HVG selection.")
    print("Marking all genes as highly variable.")

    adata.var["highly_variable"] = True
    adata.var["highly_variable_rank"] = np.arange(n_genes)
    adata.var["means"] = np.array(adata.X.mean(axis=0)).flatten()
    adata.var["dispersions"] = np.zeros(n_genes)
    adata.var["dispersions_norm"] = np.zeros(n_genes)

    adata.uns["hvg"] = {
        "flavor": flavor,
        "n_highly_variable_genes": n_requested,
        "n_hvgs_found": int(n_genes),
        "warning": f"Only {n_genes} genes available, all marked as HVG"
    }

    return adata


def find_highly_variable_genes(adata, n_highly_variable_genes, flavor):
    """Identify highly variable genes in the dataset."""
    print(f"Shape: {adata.shape}")
    print(f"Number of top genes requested: {n_highly_variable_genes}")
    print(f"Flavor: {flavor}")

    n_cells, n_genes = validate_adata(adata)

    # Adjust n_highly_variable_genes if necessary
    actual_n_hvgs = min(n_highly_variable_genes, n_genes)
    if actual_n_hvgs < n_highly_variable_genes:
        print(f"WARNING: Requested {n_highly_variable_genes} HVGs but only {n_genes} genes available.")
        print(f"Adjusting to select top {actual_n_hvgs} genes.")

    # Get HVGs
    try:
        sc.pp.highly_variable_genes(
            adata,
            flavor=flavor,
            n_top_genes=actual_n_hvgs,
            inplace=True
        )
    except ValueError as e:
        if "Bin edges must be unique" in str(e):
            print("WARNING: Binning failed - marking all genes as HVG.")
            return mark_all_genes_hvg(adata, n_genes, flavor, n_highly_variable_genes)
        raise
    n_found = adata.var["highly_variable"].sum()
    print(f"Identified {n_found} highly variable genes")

    # Store parameters in uns
    adata.uns["hvg"] = {
        "flavor": flavor,
        "n_highly_variable_genes": n_highly_variable_genes,
        "actual_n_hvgs": actual_n_hvgs,
        "n_hvgs_found": int(actual_n_hvgs),
    }

    # Ensure highly_variable column exists and is boolean
    if "highly_variable" not in adata.var.columns:
        adata.var["highly_variable"] = True
    else:
        adata.var["highly_variable"] = adata.var["highly_variable"].astype(bool)

    # Print summary
    n_hvgs_final = adata.var["highly_variable"].sum()
    print(f"Final HVG count: {n_hvgs_final}")
    print(f"Percentage of genes: {n_hvgs_final / n_genes * 100:.1f}%")

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
    input_adata = "${adata}"
    output_adata = "${prefix}.h5ad"
    n_highly_variable_genes = int("${n_hvgs}")
    flavor = "${flavor}"
    process_name = "${task.process}"

    # Read AnnData
    print(f"Finding highly variable genes in: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Find HVGs
    adata = find_highly_variable_genes(adata, n_highly_variable_genes, flavor)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with HVG annotations to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
