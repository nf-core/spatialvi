#!/usr/bin/env python

# Load packages
import spatialdata
import argparse

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Merge SpatialData objects")
    parser.add_argument("files", nargs="+", help="List of SpatialData files to merge")
    parser.add_argument("output", help="Output file name")
    args = parser.parse_args()

    # Read all zarr SpatialData files/folders
    sdatas = []
    for file in args.files:
        sdata = spatialdata.read_zarr(file)
        sdatas.append(sdata)

    # Merge the data
    output_sdata = spatialdata.concatenate(
        sdatas,
        region_key=None,
        instance_key=None,
        concatenate_tables=False,
        obs_names_make_unique=True,
        modify_tables_inplace=False,
    )

    # Save the concatenated data
    output_sdata.write(
        args.output,
        overwrite=True
    )