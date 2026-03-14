#!/usr/bin/env python3
"""
Update a SpatialData object's table with a processed AnnData object.
Also updates associated spatial elements to match filtered observations.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

# Fix numba caching issue in read-only containers
os.environ['NUMBA_CACHE_DIR'] = '/tmp/numba_cache'
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'
os.environ['XDG_CACHE_HOME'] = '/tmp/cache'

import sys
import shutil
import platform
import importlib.metadata
import yaml
import anndata as ad
import numpy as np
import spatialdata
from spatialdata.models import TableModel

# Parameters from Nextflow
input_sdata = "${sdata}"
input_adata = "${adata}"
output_sdata = "${prefix}" + "_updated.zarr"
sample_id = "${meta.id}"

# Keep only alphanumeric characters, underscores, and hyphens in the sample ID
sample_id_clean = "".join(
    filter(lambda x: x.isalnum() or x in ["_", "-"], sample_id)
)

print(f"Input SpatialData: {input_sdata}")
print(f"Input AnnData: {input_adata}")
print(f"Output SpatialData: {output_sdata}")
print(f"Sample ID: {sample_id_clean}")

# Ensure input and output are different
if os.path.abspath(input_sdata) == os.path.abspath(output_sdata):
    print("ERROR: Input and output paths are the same!")
    sys.exit(1)

# Remove output directory if it exists
if os.path.exists(output_sdata):
    print(f"Removing existing output directory: {output_sdata}")
    shutil.rmtree(output_sdata)

# Read inputs
sdata = spatialdata.read_zarr(input_sdata)
adata = ad.read_h5ad(input_adata)

print(f"Available tables: {list(sdata.tables.keys())}")
print(f"Available shapes: {list(sdata.shapes.keys())}")

# Get table name from adata.uns or use convention
if "table_name" in adata.uns:
    table_name = adata.uns["table_name"]
elif f"{sample_id_clean}_table" in sdata.tables:
    table_name = f"{sample_id_clean}_table"
elif len(sdata.tables) > 0:
    table_name = list(sdata.tables.keys())[0]
else:
    print("ERROR: No tables found in SpatialData object")
    sys.exit(1)

print(f"Updating table: {table_name}")
print(f"Original table shape: {sdata.tables[table_name].shape}")
print(f"New AnnData shape: {adata.shape}")

# Get the original table to extract metadata
original_table = sdata.tables[table_name]

# Extract region and instance_key from original table's uns
spatialdata_attrs = original_table.uns.get("spatialdata_attrs", {})
region = spatialdata_attrs.get("region", None)
region_key = spatialdata_attrs.get("region_key", None)
instance_key = spatialdata_attrs.get("instance_key", None)

print(f"Original region: {region}")
print(f"Original region_key: {region_key}")
print(f"Original instance_key: {instance_key}")

# Handle region being a numpy array (convert to list/string)
if isinstance(region, np.ndarray):
    region = region.tolist()
if isinstance(region, list) and len(region) == 1:
    region = region[0]

# Ensure instance_key column exists in the new AnnData
if instance_key and instance_key not in adata.obs.columns:
    if instance_key in original_table.obs.columns:
        # Match by index
        common_idx = adata.obs.index.intersection(original_table.obs.index)
        if len(common_idx) == len(adata.obs.index):
            adata.obs[instance_key] = original_table.obs.loc[adata.obs.index, instance_key]
        else:
            print(f"WARNING: Could not match all indices for {instance_key}")
    else:
        print(f"WARNING: instance_key '{instance_key}' not found in original table")

# Ensure region_key column exists
if region_key and region_key not in adata.obs.columns:
    if region_key in original_table.obs.columns:
        common_idx = adata.obs.index.intersection(original_table.obs.index)
        if len(common_idx) == len(adata.obs.index):
            adata.obs[region_key] = original_table.obs.loc[adata.obs.index, region_key]
    elif region:
        # Set region_key to the region value for all cells
        region_value = region if isinstance(region, str) else (region[0] if region else None)
        if region_value:
            adata.obs[region_key] = region_value

# Create a properly formatted table using TableModel
try:
    if region and instance_key:
        new_table = TableModel.parse(
            adata,
            region=region,
            region_key=region_key,
            instance_key=instance_key
        )
        print("Created new table using TableModel.parse()")
    else:
        # Fallback: copy uns from original
        print("WARNING: Missing region/instance_key metadata, copying from original")
        adata.uns["spatialdata_attrs"] = spatialdata_attrs.copy()
        new_table = adata

except Exception as e:
    print(f"WARNING: TableModel.parse failed: {e}")
    print("Copying uns from original table...")
    adata.uns["spatialdata_attrs"] = spatialdata_attrs.copy()
    new_table = adata

# Delete old table and add new one
del sdata.tables[table_name]

try:
    sdata.tables[table_name] = new_table
    print(f"Successfully updated table '{table_name}'")
except Exception as e:
    print(f"WARNING: Failed to set table via sdata.tables: {e}")
    print("Attempting to set table using internal method...")

    try:
        # Try setting directly on the internal dict
        if hasattr(sdata, '_tables'):
            sdata._tables[table_name] = new_table
            print("Set table using _tables dict")
        else:
            raise AttributeError("No _tables attribute found")
    except Exception as e2:
        print(f"ERROR: All methods failed: {e2}")
        sys.exit(1)

# Optionally update associated spatial elements to match filtered observations
if region and adata.shape[0] < original_table.shape[0]:
    print(f"Observations were filtered: {original_table.shape[0]} -> {adata.shape[0]}")

    region_name = region if isinstance(region, str) else (region[0] if isinstance(region, list) else None)

    if region_name and region_name in sdata.shapes:
        try:
            matched_element, _ = spatialdata.match_element_to_table(
                sdata,
                element_name=region_name,
                table_name=table_name
            )
            sdata.shapes[region_name] = matched_element[region_name]
            print(f"Updated shapes '{region_name}' to match filtered table")
        except Exception as e:
            print(f"WARNING: Could not update shapes: {e}")
else:
    print("No filtering detected or no region to update")

# Write updated SpatialData to new location
try:
    sdata.write(output_sdata)
    print(f"Written updated SpatialData to: {output_sdata}")
except Exception as e:
    print(f"ERROR: Failed to write SpatialData: {e}")
    sys.exit(1)

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "spatialdata": importlib.metadata.version("spatialdata"),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
