#!/usr/bin/env python3
"""
Integrate AnnData objects using Harmony.
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
import scanpy as sc
import scanpy.external as sce


def integrate_harmony(adata_list, sample_names, sample_col='library_id'):
    """Integrate multiple samples using Harmony."""

    print(f"Concatenating {len(adata_list)} AnnData objects")

    # Merge all `.obs` and `.X` into one AnnData object
    adata = sc.concat(
        adata_list,
        label=sample_col,
        uns_merge="unique",
        keys=sample_names,
        join='outer',
        index_unique="-"
    )

    # `.var` is dropped during the previous merge, so we add them back manually
    merged_var = pd.concat([adata.var for adata in adata_list], join="outer")
    merged_var = merged_var[~merged_var.index.duplicated()]
    adata.var = merged_var.loc[adata.var_names]

    # Integrate with Harmony
    print(f'Integrating {len(adata_list)} samples with Harmony')
    sce.pp.harmony_integrate(
        adata,
        key=sample_col,
        adjusted_basis="X_harmony"
    )
    print(f"Final integrated AnnData shape: {adata.shape}")

    return adata


def write_versions(process_name):
    """Write software versions to YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
            "harmony": importlib.metadata.version("harmonypy"),
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

    sample_names = [Path(h5ad).stem for h5ad in h5ad_files]
    adata_integrated = integrate_harmony(adata_list, sample_names)

    adata_integrated.write_h5ad(output_file)
    print(f"Written integrated AnnData to: {output_file}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
