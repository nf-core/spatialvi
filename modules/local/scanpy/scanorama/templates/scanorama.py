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

import anndata as ad
import scanpy.external as sce
import scipy.sparse as sp


def integrate_scanorama(adata, key='library_id'):
    """Integrate multiple samples using Scanorama."""

    # Convert to CSR format as required by Scanorama, if data is CSC
    if sp.issparse(adata.X) and not isinstance(adata.X, sp.csr_matrix):
        print('Converting to CSR data format')
        adata.X = adata.X.tocsr()

    sce.pp.scanorama_integrate(
        adata,
        key=key,
        adjusted_basis='X_scanorama'
    )

    print(f"Final integrated AnnData shape: {adata.shape}")

    return adata


def write_versions(process_name):
    """Write software versions to a YAML file."""
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
    h5ad = "${h5ad}"
    output_file = "${prefix}.h5ad"
    process_name = "${task.process}"

    print(f"Read AnnData object {h5ad}")
    adata = ad.read_h5ad(h5ad)

    adata_integrated = integrate_scanorama(adata)

    adata_integrated.write_h5ad(output_file)
    print(f"Written integrated AnnData to: {output_file}")

    write_versions(process_name)


if __name__ == "__main__":
    main()
