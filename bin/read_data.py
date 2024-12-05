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
    parser.add_argument(
        "--sampleID",
        metavar="sampleID",
        type=str,
        default=None,
        help="Sample ID.",
    )

    args = parser.parse_args()

    # Keep only alphanumeric characters, underscores, and hyphens in the sample ID
    args.sampleID = "".join(
        filter(lambda x: x.isalnum() or x in ["_", "-"], args.sampleID)
    )

    # Read Visium data
    spatialdata = spatialdata_io.visium(
        args.SRCountDir,
        counts_file="raw_feature_bc_matrix.h5",
        dataset_id=args.sampleID,
    )

    # Remove sampleID metadata from table:
    if args.sampleID in spatialdata.tables['table'].uns.keys():
        del spatialdata.tables['table'].uns[args.sampleID]
    
    # Rename table into sample id
    spatialdata.tables[f'{args.sampleID}_table'] = spatialdata.tables['table']
    del spatialdata.tables['table']

    # Write raw spatialdata to file
    spatialdata.write(args.output_sdata, overwrite=True)
