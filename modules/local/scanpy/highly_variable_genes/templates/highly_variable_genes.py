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


def select_hvgs_by_mean_expression(adata, n_hvgs):
    """Fallback method: select top genes by mean expression."""
    print("Falling back to marking top genes by mean expression as HVG...")

    mean_expr = np.array(adata.X.mean(axis=0)).flatten()
    top_indices = np.argsort(mean_expr)[::-1][:n_hvgs]

    adata.var["highly_variable"] = False
    adata.var.iloc[top_indices, adata.var.columns.get_loc("highly_variable")] = True
    adata.var["means"] = mean_expr
    adata.var["highly_variable_rank"] = np.nan
    adata.var.loc[adata.var["highly_variable"], "highly_variable_rank"] = np.arange(n_hvgs)

    print(f"Selected top {n_hvgs} genes by mean expression")

    return adata, n_hvgs


def find_hvgs_scanpy(adata, flavor, n_hvgs):
    """Find highly variable genes using scanpy methods."""
    try:
        sc.pp.highly_variable_genes(
            adata,
            flavor=flavor,
            n_top_genes=n_hvgs,
            inplace=True
        )
        n_found = adata.var["highly_variable"].sum()
        print(f"Identified {n_found} highly variable genes")
        return adata, flavor, n_found

    except Exception as e:
        print(f"WARNING: HVG selection with flavor '{flavor}' failed: {e}")
        print("Attempting with flavor 'cell_ranger'...")

        try:
            sc.pp.highly_variable_genes(
                adata,
                flavor="cell_ranger",
                n_top_genes=n_hvgs,
                inplace=True
            )
            n_found = adata.var["highly_variable"].sum()
            print(f"Identified {n_found} highly variable genes with cell_ranger flavor")
            return adata, "cell_ranger", n_found

        except Exception as e2:
            print(f"WARNING: HVG selection with 'cell_ranger' also failed: {e2}")
            return None, None, None


def find_highly_variable_genes(adata, n_highly_variable_genes, flavor):
    """Identify highly variable genes in the dataset."""
    print(f"Shape: {adata.shape}")
    print(f"Number of top genes requested: {n_highly_variable_genes}")
    print(f"Flavor: {flavor}")

    n_cells, n_genes = validate_adata(adata)

    # Handle edge case: too few genes
    if n_genes < 10:
        return mark_all_genes_hvg(adata, n_genes, flavor, n_highly_variable_genes)

    # Adjust n_highly_variable_genes if necessary
    actual_n_hvgs = min(n_highly_variable_genes, n_genes)

    if actual_n_hvgs < n_highly_variable_genes:
        print(f"WARNING: Requested {n_highly_variable_genes} HVGs but only {n_genes} genes available.")
        print(f"Adjusting to select top {actual_n_hvgs} genes.")

    # Try scanpy HVG methods
    result = find_hvgs_scanpy(adata, flavor, actual_n_hvgs)

    if result[0] is not None:
        adata, used_flavor, n_hvgs = result
    else:
        # Fallback to mean expression method
        adata, n_hvgs = select_hvgs_by_mean_expression(adata, actual_n_hvgs)
        used_flavor = "mean_expression_fallback"

    # Store parameters in uns
    adata.uns["hvg"] = {
        "flavor": used_flavor,
        "n_highly_variable_genes": n_highly_variable_genes,
        "actual_n_hvgs": actual_n_hvgs,
        "n_hvgs_found": int(n_hvgs),
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
