#!/usr/bin/env python3
"""
Integrate AnnData objects using Scanorama.
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
import scanorama
import scanpy as sc


def integrate_scanorama(adata_list, labels):
    """Integrate multiple samples using Scanorama."""

    # Perform integration across all AnnData objects
    print(f"Integrating {len(adata_list)} AnnData objects")
    adatas_corrected = scanorama.correct_scanpy(
        adata_list,
        return_dimred=True
    )

    # Merge all `.obs` and `.X` into one AnnData object
    adata = sc.concat(
        adatas_corrected,
        label="library_id",
        uns_merge="unique",
        keys=labels,
        index_unique="-"
    )

    # `.var` is dropped during the previous merge, so we add them back manually
    merged_var = pd.concat([adata.var for adata in adata_list], join="outer")
    merged_var = merged_var[~merged_var.index.duplicated()]
    adata.var = merged_var.loc[adata.var_names]

    print(f"Final integrated AnnData shape: {adata.shape}")

    return adata


def write_versions(process_name):
    """Write software versions to YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
            "scanorama": importlib.metadata.version("scanorama"),
            "scanpy": importlib.metadata.version("scanpy"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Integrate multiple AnnData objects into one."""

    # Template variables
    h5ad_files = "${h5ad}".split()
    output_file = "${prefix}.h5ad"
    process_name = "${task.process}"

    adata_list = []
    for h5ad in h5ad_files:
        print(f"Reading {h5ad} AnnData object")
        adata = ad.read_h5ad(h5ad)
        adata_list.append(adata)

    labels = [Path(f).stem for f in h5ad_files]
    adata_integrated = integrate_scanorama(adata_list, labels)

    adata_integrated.write_h5ad(output_file)
    print(f"Written integrated AnnData to: {output_file}")

    write_versions(process_name)


if __name__ == "__main__":
    main()
