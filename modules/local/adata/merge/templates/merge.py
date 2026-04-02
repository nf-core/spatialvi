#!/usr/bin/env python3
"""
Merge multiple AnnData objects into one.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform
import yaml
from pathlib import Path

import anndata as ad
import pandas as pd


def add_var(adata, adata_list, join='inner'):
    """
    Adds `.var` back into a merged AnnData object from the original list of
    multiple AnnData objects.
    """
    merged_var = pd.concat([adata.var for adata in adata_list], join=join)
    merged_var = merged_var[~merged_var.index.duplicated()]
    adata.var = merged_var.loc[adata.var_names]
    print('Preserved `.var` data')
    return adata


def add_spatial(adata, adata_list):
    """
    Adds `.uns['spatial']` back into a merged AnnData objects from the original
    list of multiple AnnData objects.
    """
    merged_spatial = {}
    for adata_orig in adata_list:
        if 'spatial' in adata_orig.uns:
            merged_spatial.update(adata_orig.uns['spatial'])
    if merged_spatial:
        adata.uns['spatial'] = merged_spatial
        print('Preserved `.uns["spatial"]` data')
    return adata


def merge_adata(adata_list,
                keys,
                join='inner',
                label='library_id',
                preserve_var=True,
                preserve_spatial=True):
    """
    Merge multiple AnnData objects into one. Can optionally preserve both `.var`
    and `.uns['spatial']` for the final merged object.
    """

    print(f"Merging {len(adata_list)} AnnData objects")
    adata = ad.concat(
        adata_list,
        join=join,
        label=label,
        keys=keys,
        index_unique="-"
    )

    if preserve_var:
        adata = add_var(adata, adata_list, join)

    if preserve_spatial:
        adata = add_spatial(adata, adata_list)

    adata.uns['merge'] = {
        'n_samples': len(keys),
        'join': join,
        'label': label,
        'preserve_var': preserve_var,
        'preserve_spatial': preserve_spatial,
    }
    print(f"Final merged AnnData {adata}")

    return adata

def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Merge multiple AnnData objects into one."""

    # Template variables
    h5ad_files = "${h5ad}".split()
    join = "${join}"
    label = "${label}"
    preserve_var = "${preserve_var}" == "true"
    preserve_spatial = "${preserve_spatial}" == "true"
    output_file = "${prefix}.h5ad"
    process_name = "${task.process}"

    adata_list = []
    for h5ad in h5ad_files:
        adata = ad.read_h5ad(h5ad)
        adata_list.append(adata)
        print(f"Read AnnData object {adata}")

    sample_names = [Path(h5ad).stem for h5ad in h5ad_files]
    adata_integrated = merge_adata(
        adata_list,
        sample_names,
        join=join,
        label=label,
        preserve_var=preserve_var,
        preserve_spatial=preserve_spatial
    )

    adata_integrated.write_h5ad(output_file)
    print(f"Written merged AnnData to: {output_file}")

    write_versions(process_name)


if __name__ == "__main__":
    main()
