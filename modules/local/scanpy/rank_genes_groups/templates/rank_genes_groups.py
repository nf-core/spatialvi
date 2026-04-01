#!/usr/bin/env python3
"""
Rank genes for characterizing groups (differential expression analysis).
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform
import yaml

import anndata as ad
import scanpy as sc


def rank_genes(adata, groupby, method):
    """Rank genes by groups for differential expression analysis."""
    print(f"Shape: {adata.shape}")
    print(f"Groupby: {groupby}")
    print(f"Method: {method}")

    # Validate input
    if groupby not in adata.obs.columns:
        raise ValueError(f"Column '{groupby}' not found in adata.obs")

    # Rank genes by groups
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method
    )

    # Print summary
    n_groups = adata.obs[groupby].nunique()
    print(f"Computed DEGs for {n_groups} groups")

    # Print top genes per group
    for group in adata.obs[groupby].unique():
        genes = sc.get.rank_genes_groups_df(adata, group=str(group))
        top_genes = genes.head(5)["names"].tolist()
        print(f"  Group {group} top 5: {', '.join(top_genes)}")

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
    """Rank genes for characterizing groups in an AnnData object."""

    # Template variables
    input_adata = "${adata}"
    groupby = "${groupby}"
    method = "${method}"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    # Read AnnData
    print(f"Performing differential expression analysis on: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Rank genes
    adata = rank_genes(adata, groupby, method)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with DEGs to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
