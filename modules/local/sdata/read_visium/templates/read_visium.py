#!/usr/bin/env python3
"""
Read Visium or Visium HD data from Space Ranger output into SpatialData format.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import logging
import platform
import re
import shutil

import spatialdata_io
import yaml

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def read_visium_hd(spaceranger_dir, sample_id_clean, hd_bin_size):
    """Read Visium HD data."""
    logger.info(f"Reading Visium HD data with bin size: {hd_bin_size}")

    # Add sample ID to feature_slice, if not present
    feature_slice_src = os.path.join(
        spaceranger_dir, "feature_slice.h5"
    )
    feature_slice_dst = os.path.join(
        spaceranger_dir, f"{sample_id_clean}_feature_slice.h5"
    )
    if (os.path.exists(feature_slice_src)
    and not os.path.exists(feature_slice_dst)):
        logger.info(f"Copying {feature_slice_src} to {feature_slice_dst}")
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
    logger.info("Reading standard Visium data")
    sdata = spatialdata_io.visium(
        spaceranger_dir,
        counts_file="raw_feature_bc_matrix.h5",
        dataset_id=sample_id_clean,
    )
    table_name = "table"
    return sdata, table_name


def read_visium_data(spaceranger_dir, sample_id_clean, hd_bin_size):
    """Read Visium or Visium HD data from Space Ranger output."""
    logger.info(f"Reading data from: {spaceranger_dir}")
    logger.info(f"Sample ID: {sample_id_clean}")

    is_hd_data = os.path.isdir(os.path.join(spaceranger_dir, "binned_outputs"))
    logger.info(f"Visium HD: {is_hd_data}")

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

    logger.info("SpatialData object created:")
    logger.info(f"  Tables: {list(sdata.tables.keys())}")
    logger.info(f"  Shapes: {list(sdata.shapes.keys())}")
    logger.info(f"  Images: {list(sdata.images.keys())}")

    return sdata, table_name


def validate_table(sdata, table_name, sample_id_clean):
    """Validate table exists and process it."""
    if table_name not in sdata.tables:
        raise ValueError(
            f"Expected table '{table_name}' not found in SpatialData. "
            f"Available tables: {list(sdata.tables.keys())}"
        )

    # Get table and print info
    adata = sdata.tables[table_name]
    logger.info(f"Table '{table_name}' shape: {adata.shape}")
    logger.info(f"Number of genes: {adata.n_vars}")
    logger.info(f"Number of observations: {adata.n_obs}")

    # Remove `sample_id` metadata from table uns, if present
    if sample_id_clean in adata.uns.keys():
        logger.info(f"Removing '{sample_id_clean}' from table uns")
        del adata.uns[sample_id_clean]

    adata.var_names_make_unique()

    # Rename table to include sample ID
    new_table_name = f"{sample_id_clean}_table"
    if table_name != new_table_name:
        logger.info(f"Renaming table: '{table_name}' -> '{new_table_name}'")
        sdata.tables[new_table_name] = adata
        del sdata.tables[table_name]

    # Validate final table
    final_table = sdata.tables[new_table_name]
    logger.info(f"Final table '{new_table_name}' shape: {final_table.shape}")

    return sdata


def write_versions(process_name):
    """Write software versions to a YAML file."""
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
    sample_id_clean = re.sub(r"[^a-zA-Z0-9_-]", "", sample_id)

    sdata, table_name = read_visium_data(
        spaceranger_dir,
        sample_id_clean,
        hd_bin_size
    )

    sdata = validate_table(sdata, table_name, sample_id_clean)

    logger.info(f"Writing SpatialData to: {output_sdata}")
    sdata.write(output_sdata, overwrite=True)

    write_versions(process_name)

if __name__ == "__main__":
    main()
