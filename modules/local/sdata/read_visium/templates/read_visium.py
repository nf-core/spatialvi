#!/usr/bin/env python3
"""
Read Visium or Visium HD data from Space Ranger output into SpatialData format.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform
import shutil
import sys
import yaml

import spatialdata_io


def read_visium_hd(spaceranger_dir, sample_id_clean, hd_bin_size):
    """Read Visium HD data."""
    print(f"Reading Visium HD data with bin size: {hd_bin_size}")

    # Add sample ID to feature_slice, if not present
    feature_slice_src = os.path.join(
        spaceranger_dir, "feature_slice.h5"
    )
    feature_slice_dst = os.path.join(
        spaceranger_dir, f"{sample_id_clean}_feature_slice.h5"
    )
    if os.path.exists(feature_slice_src) \
        and not os.path.exists(feature_slice_dst):
            print(f"Copying {feature_slice_src} to {feature_slice_dst}")
            shutil.copyfile(feature_slice_src, feature_slice_dst)

    sdata = spatialdata_io.visium_hd(
        spaceranger_dir,
        bin_size=[hd_bin_size],
        dataset_id=sample_id_clean,
    )
    table_name = f"square_{hd_bin_size:03d}um"

    return sdata, table_name


def read_visium_standard(spaceranger_dir, sample_id_clean):
    """Read standard Visium data."""
    print("Reading standard Visium data")
    sdata = spatialdata_io.visium(
        spaceranger_dir,
        counts_file="raw_feature_bc_matrix.h5",
        dataset_id=sample_id_clean,
    )
    table_name = "table"
    return sdata, table_name


def read_visium_data(spaceranger_dir, sample_id_clean, hd_bin_size):
    """Read Visium or Visium HD data from Space Ranger output."""
    print(f"Reading data from: {spaceranger_dir}")
    print(f"Sample ID: {sample_id_clean}")

    is_hd_data = os.path.isdir(os.path.join(spaceranger_dir, "binned_outputs"))
    print(f"Visium HD: {is_hd_data}")

    if is_hd_data:
        sdata, table_name = read_visium_hd(
            spaceranger_dir,
            sample_id_clean,
            hd_bin_size
        )
    else:
        sdata, table_name = read_visium_standard(
            spaceranger_dir,
            sample_id_clean
        )

    print("SpatialData object created")
    print(f"Tables: {list(sdata.tables.keys())}")
    print(f"Shapes: {list(sdata.shapes.keys())}")
    print(f"Images: {list(sdata.images.keys())}")

    return sdata, table_name


def validate_table(sdata, table_name, sample_id_clean):
    """Validate table exists and process it."""
    # Check if the table exists
    if table_name not in sdata.tables:
        print(f"ERROR: Expected table '{table_name}' not found in SpatialData")
        print(f"Available tables: {list(sdata.tables.keys())}")
        sys.exit(1)

    # Get table and print info
    adata = sdata.tables[table_name]
    print(f"Table '{table_name}' shape: {adata.shape}")
    print(f"Number of genes: {adata.n_vars}")
    print(f"Number of cells/spots: {adata.n_obs}")

    # Remove `sample_id` metadata from table uns, if present
    if sample_id_clean in adata.uns.keys():
        print(f"Removing '{sample_id_clean}' from table uns")
        del adata.uns[sample_id_clean]

    adata.var_names_make_unique()

    # Rename table to include sample ID
    new_table_name = f"{sample_id_clean}_table"
    if table_name != new_table_name:
        print(f"Renaming table: '{table_name}' -> '{new_table_name}'")
        sdata.tables[new_table_name] = adata
        del sdata.tables[table_name]

    # Validate final table
    final_table = sdata.tables[new_table_name]
    print(f"Final table '{new_table_name}' shape: {final_table.shape}")

    return sdata


def write_versions(process_name):
    """Write software versions to YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "spatialdata": importlib.metadata.version("spatialdata"),
            "spatialdata-io": importlib.metadata.version("spatialdata-io"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Read Visium data into SpatialData format."""
    # Template variables
    spaceranger_dir = "${spaceranger_dir}"
    sample_id = "${meta.id}"
    hd_bin_size = int("${hd_bin_size}")
    output_sdata = "${prefix}.zarr"
    process_name = "${task.process}"

    # Sample ID must only contain alphanumerics, underscores and dashes
    sample_id_clean = "".join(
        filter(lambda x: x.isalnum() or x in ["_", "-"], sample_id)
    )

    sdata, table_name = read_visium_data(
        spaceranger_dir,
        sample_id_clean,
        hd_bin_size
    )

    sdata = validate_table(sdata, table_name, sample_id_clean)

    print(f"Writing SpatialData to: {output_sdata}")
    sdata.write(output_sdata, overwrite=True)

    write_versions(process_name)


if __name__ == "__main__":
    main()
