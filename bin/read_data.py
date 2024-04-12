#!/usr/bin/env python

# Load packages
import argparse

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
    # TODO Add argument with meta.id for dataset_id

    args = parser.parse_args()

    # Read Visium data
    spatialdata = spatialdata_io.visium(
        args.SRCountDir, counts_file="raw_feature_bc_matrix.h5", dataset_id="visium"
    )

    # Write raw spatialdata to file
    spatialdata.write(args.output_sdata, overwrite=True)
