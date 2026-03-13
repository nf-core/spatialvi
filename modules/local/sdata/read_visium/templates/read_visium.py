#!/usr/bin/env python3
"""
Read Visium or Visium HD data from Space Ranger output into SpatialData format.
"""

import os
import sys
import shutil

import platform
import importlib.metadata
import yaml
import spatialdata_io

# Parameters from Nextflow
spaceranger_dir = "${spaceranger_dir}"
output_sdata = "${prefix}" + ".zarr"
sample_id = "${meta.id}"
hd_bin_size = int("${hd_bin_size}")

# Keep only alphanumeric characters, underscores, and hyphens in the sample ID
sample_id_clean = "".join(
    filter(lambda x: x.isalnum() or x in ["_", "-"], sample_id)
)

print(f"Reading data from: {spaceranger_dir}")
print(f"Sample ID: {sample_id_clean}")
print(f"Output: {output_sdata}")

# Check if this is Visium HD data
is_visium_hd = os.path.isdir(os.path.join(spaceranger_dir, "binned_outputs"))

print(f"Visium HD: {is_visium_hd}")

# Read Visium data
if is_visium_hd:
    print(f"Reading Visium HD data with bin size: {hd_bin_size}")
    
    # Copy feature_slice.h5 file with sample ID prefix (required by spatialdata_io)
    feature_slice_src = os.path.join(spaceranger_dir, "feature_slice.h5")
    feature_slice_dst = os.path.join(spaceranger_dir, f"{sample_id_clean}_feature_slice.h5")
    
    if os.path.exists(feature_slice_src) and not os.path.exists(feature_slice_dst):
        print(f"Copying {feature_slice_src} to {feature_slice_dst}")
        shutil.copyfile(feature_slice_src, feature_slice_dst)
    
    sdata = spatialdata_io.visium_hd(
        spaceranger_dir,
        bin_size=[hd_bin_size],
        dataset_id=sample_id_clean,
    )
    table_name = f"square_{hd_bin_size:03d}um"
else:
    print("Reading standard Visium data")
    sdata = spatialdata_io.visium(
        spaceranger_dir,
        counts_file="raw_feature_bc_matrix.h5",
        dataset_id=sample_id_clean,
    )
    table_name = "table"

print(f"SpatialData object created")
print(f"Tables: {list(sdata.tables.keys())}")
print(f"Shapes: {list(sdata.shapes.keys())}")
print(f"Images: {list(sdata.images.keys())}")

# Check if the table exists
if table_name not in sdata.tables:
    print(f"ERROR: Expected table '{table_name}' not found in SpatialData")
    print(f"Available tables: {list(sdata.tables.keys())}")
    sys.exit(1)

# Print table info before renaming
adata = sdata.tables[table_name]
print(f"Table '{table_name}' shape: {adata.shape}")
print(f"Number of genes: {adata.n_vars}")
print(f"Number of cells/spots: {adata.n_obs}")

# Remove sample_id metadata from table uns (if present)
if sample_id_clean in adata.uns.keys():
    print(f"Removing '{sample_id_clean}' from table uns")
    del adata.uns[sample_id_clean]

# Make var_names unique
adata.var_names_make_unique()

# Rename table to include sample ID
new_table_name = f"{sample_id_clean}_table"
if table_name != new_table_name:
    print(f"Renaming table: '{table_name}' -> '{new_table_name}'")
    sdata.tables[new_table_name] = adata
    del sdata.tables[table_name]

# Verify final table
final_table = sdata.tables[new_table_name]
print(f"Final table '{new_table_name}' shape: {final_table.shape}")

# Write SpatialData to zarr
print(f"Writing SpatialData to: {output_sdata}")
sdata.write(output_sdata, overwrite=True)
print("Done!")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "spatialdata": importlib.metadata.version("spatialdata"),
        "spatialdata-io": importlib.metadata.version("spatialdata-io"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)