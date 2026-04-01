#!/usr/bin/env python3
"""
Extract a legacy AnnData object from a SpatialData object. The legacy format
includes spatial coordinates and images compatible with scanpy.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform
import yaml

import scipy.sparse
import spatialdata
from spatialdata_io.experimental import to_legacy_anndata


def find_table_name(sdata, sample_id_clean):
    """Find the appropriate table name in a SpatialData object."""
    expected_name = f"{sample_id_clean}_table"
    if expected_name in sdata.tables:
        return expected_name
    fallback_name = list(sdata.tables.keys())[0]
    print(f"Warning: Expected table '{expected_name}' not found, using '{fallback_name}'")
    return fallback_name


def find_coordinate_system(sdata, sample_id_clean):
    """Find the appropriate coordinate system in a SpatialData object."""
    expected_name = f"{sample_id_clean}_downscaled_hires"
    if expected_name in sdata.coordinate_systems:
        return expected_name
    fallback_name = list(sdata.coordinate_systems)[0]
    print(f"Warning: Expected coordinate system not found, using '{fallback_name}'")
    return fallback_name


def extract_to_legacy_anndata(sdata, table_name, coord_system):
    """Convert SpatialData to the legacy AnnData format."""
    print(f"Using table: {table_name}")
    print(f"Using coordinate system: {coord_system}")
    adata = to_legacy_anndata(
        sdata,
        coordinate_system=coord_system,
        table_name=table_name,
        include_images=True
    )
    return adata


def ensure_sparse_csc(adata):
    """Convert `.X` matrix to CSC sparse format for compatibility."""
    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csc_matrix(adata.X)
    elif not scipy.sparse.isspmatrix_csc(adata.X):
        adata.X = scipy.sparse.csc_matrix(adata.X)
    return adata


def add_metadata(adata, sample_id_clean, table_name, coord_system):
    """Add raw counts layer and metadata to an AnnData object."""
    if 'raw' not in adata.layers:
        adata.layers['raw'] = adata.X.copy()
    # Store additional metadata in `.uns`
    adata.uns['sample_id'] = sample_id_clean
    adata.uns['table_name'] = table_name
    adata.uns['coordinate_system'] = coord_system
    return adata


def write_versions(process_name):
    """Write software versions to YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
            "spatialdata": importlib.metadata.version("spatialdata"),
            "spatialdata-io": importlib.metadata.version("spatialdata-io"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Extract legacy AnnData from SpatialData."""
    # Template variables
    input_sdata = "${sdata}"
    sample_id = "${meta.id}"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    # Sample ID must only contain alphanumerics, underscores and dashes
    sample_id_clean = "".join(
        filter(lambda x: x.isalnum() or x in ["_", "-"], sample_id)
    )

    print(f"Processing SpatialData from: {input_sdata}")
    sdata = spatialdata.read_zarr(input_sdata)

    table_name = find_table_name(sdata, sample_id_clean)
    coord_system = find_coordinate_system(sdata, sample_id_clean)

    adata = extract_to_legacy_anndata(sdata, table_name, coord_system)

    adata = ensure_sparse_csc(adata)
    adata = add_metadata(adata, sample_id_clean, table_name, coord_system)

    print(f"Sample ID: {sample_id_clean}")
    print(f"Tables found: {list(sdata.tables.keys())}")
    print(f"Coordinate systems: {list(sdata.coordinate_systems)}")
    print(f"AnnData shape: {adata.shape}")
    print(f"Spatial keys: {list(adata.uns.get('spatial', {}).keys())}")

    adata.write_h5ad(output_adata)
    print(f"Written legacy AnnData to: {output_adata}")

    write_versions(process_name)


if __name__ == "__main__":
    main()
