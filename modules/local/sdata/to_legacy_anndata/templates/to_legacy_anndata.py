#!/usr/bin/env python3
"""
Extract a legacy AnnData object from a SpatialData object.
The legacy format includes spatial coordinates and images compatible with scanpy.
"""

import platform
import importlib.metadata
import yaml
import scipy.sparse
import spatialdata
from spatialdata_io.experimental import to_legacy_anndata

# Parameters from Nextflow
input_sdata = "${sdata}"
output_adata = "${prefix}.h5ad"
sample_id = "${meta.id}"

# Keep only alphanumeric characters, underscores, and hyphens in the sample ID
sample_id_clean = "".join(
    filter(lambda x: x.isalnum() or x in ["_", "-"], sample_id)
)

# Read SpatialData object
sdata = spatialdata.read_zarr(input_sdata)

print(f"Processing SpatialData from: {input_sdata}")
print(f"Sample ID: {sample_id_clean}")
print(f"Tables found: {list(sdata.tables.keys())}")
print(f"Coordinate systems: {list(sdata.coordinate_systems)}")

# Find the appropriate table and coordinate system
# Convention: table name is "{sample_id}_table", coordinate system is "{sample_id}_downscaled_hires"
table_name = f"{sample_id_clean}_table"
coord_system = f"{sample_id_clean}_downscaled_hires"

# Fallback: use first available table/coordinate system if convention doesn't match
if table_name not in sdata.tables:
    table_name = list(sdata.tables.keys())[0]
    print(f"Warning: Expected table '{sample_id_clean}_table' not found, using '{table_name}'")

if coord_system not in sdata.coordinate_systems:
    coord_system = list(sdata.coordinate_systems)[0]
    print(f"Warning: Expected coordinate system not found, using '{coord_system}'")

print(f"Using table: {table_name}")
print(f"Using coordinate system: {coord_system}")

# Convert to legacy AnnData
adata = to_legacy_anndata(
    sdata,
    coordinate_system=coord_system,
    table_name=table_name,
    include_images=True
)

# Convert X matrix to CSC sparse format for compatibility
if not scipy.sparse.issparse(adata.X):
    adata.X = scipy.sparse.csc_matrix(adata.X)
elif not scipy.sparse.isspmatrix_csc(adata.X):
    adata.X = scipy.sparse.csc_matrix(adata.X)

# Store raw counts layer
if 'raw' not in adata.layers:
    adata.layers['raw'] = adata.X.copy()

# Store metadata for later use
adata.uns['sample_id'] = sample_id_clean
adata.uns['table_name'] = table_name
adata.uns['coordinate_system'] = coord_system

print(f"AnnData shape: {adata.shape}")
print(f"Spatial keys: {list(adata.uns.get('spatial', {}).keys())}")

# Write AnnData
adata.write_h5ad(output_adata)
print(f"Written legacy AnnData to: {output_adata}")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "spatialdata": importlib.metadata.version("spatialdata"),
        "spatialdata-io": importlib.metadata.version("spatialdata-io"),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)