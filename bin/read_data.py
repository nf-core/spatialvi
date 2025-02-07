#!/usr/bin/env python

# Load packages
import argparse
import os
import shutil

import spatialdata_io

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Load spatial transcriptomics data from MTX matrices and aligned images."
    )
    parser.add_argument(
        "--SRCountDir",
        metavar="SRCountDir",
        type=str,
        default=None,
        help="Input directory with Spaceranger data.",
    )
    parser.add_argument(
        "--output_sdata",
        metavar="output_sdata",
        type=str,
        default=None,
        help="Output spatialdata zarr path.",
    )
    parser.add_argument(
        "--sampleID",
        metavar="sampleID",
        type=str,
        default=None,
        help="Sample ID.",
    )
    parser.add_argument(
        "--visium_hd",
        action='store_true',
        help="Visium HD data.",
    )
    parser.add_argument(
        "--bin_size",
        metavar="bin_size",
        type=int,
        default=16,
        help="Bin size in micrometers.",
    )

    args = parser.parse_args()

    # Keep only alphanumeric characters, underscores, and hyphens in the sample ID
    args.sampleID = "".join(
        filter(lambda x: x.isalnum() or x in ["_", "-"], args.sampleID)
    )

    # Read Visium data
    if args.visium_hd:
        # Copy file f"{args.SRCountDir}/feature_slice.h5" to f"{args.sampleID}_feature_slice.h5"
        shutil.copyfile(
            os.path.join(args.SRCountDir, "feature_slice.h5"),
            os.path.join(args.SRCountDir, f"{args.sampleID}_feature_slice.h5")
        )
        spatialdata = spatialdata_io.visium_hd(
            args.SRCountDir,
            bin_size=[args.bin_size],
            dataset_id=args.sampleID,
        )
        table_name = f'square_{args.bin_size:03d}um'
    else:
        spatialdata = spatialdata_io.visium(
            args.SRCountDir,
            counts_file="raw_feature_bc_matrix.h5",
            dataset_id=args.sampleID,
        )
        table_name = 'table'

    # Remove sampleID metadata from table:
    if args.sampleID in spatialdata.tables[table_name].uns.keys():
        del spatialdata.tables[table_name].uns[args.sampleID]

    # Rename table into sample id
    spatialdata.tables[f'{args.sampleID}_table'] = spatialdata.tables[table_name]
    del spatialdata.tables[table_name]

    # Rename var_names to unique:
    spatialdata.tables[f'{args.sampleID}_table'].var_names_make_unique()

    # Write raw spatialdata to file
    spatialdata.write(args.output_sdata, overwrite=True)
