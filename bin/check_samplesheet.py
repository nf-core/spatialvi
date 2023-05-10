#!/usr/bin/env python

import argparse
import errno
import os
import sys


def parse_args(argv=None):
    Description = "Check contents of nf-core/spatialtranscriptomics samplesheet."
    Epilog = "Example usage: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("file_in", help="Input samplesheet file.")
    parser.add_argument("file_out", help="Output validated samplesheet file.")
    parser.add_argument(
        "--is_raw_data", action="store_true", help="Whether input is raw data to be processed by SpaceRanger."
    )
    return parser.parse_args(argv)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out, is_raw_data):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.
    Validate the general shape of the table, expected columns and each row.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet.
        file_out (pathlib.Path): Where the validated samplesheet should be created.
        is_raw_data (boolean): Whether the input is raw spatial data to be processed by SpaceRanger.

    The following structure is expected:
        sample,tissue_positions_list,tissue_lowres_image,tissue_hires_image,scale_factors,barcodes,features,matrix
        SAMPLE,TISSUE_POSITIONS_LIST.csv,TISSUE_LOWRES_IMAGE.png,TISSUE_HIGHRES_IMAGE.png,SCALEFACTORS_JSON.json,BARCODES.tsv.gz,FEATURES.tsv.gz,MATRIX.mtx.gz

    For a complete example see:
    https://data.githubusercontent.com/nf-core/test-datasets/spatialtranscriptomics/testdata/test-dataset-subsampled/samplesheet.csv
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:
        # Get cols and header depending on samplesheet type
        if is_raw_data:
            MIN_COLS = 4
            HEADER = ["sample", "fastq_dir", "tissue_hires_image", "slide", "area"]
        else:
            MIN_COLS = 7
            HEADER = [
                "sample",
                "tissue_positions_list",
                "tissue_lowres_image",
                "tissue_hires_image",
                "scale_factors",
                "barcodes",
                "features",
                "matrix",
            ]

        # Check header
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        # Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            # Check sample name entries
            if is_raw_data:
                (
                    sample,
                    fastq_dir,
                    tissue_hires_image,
                    slide,
                    area,
                ) = lspl[: len(HEADER)]
                sample_info = [
                    fastq_dir,
                    tissue_hires_image,
                    slide,
                    area,
                ]
            else:
                (
                    sample,
                    tissue_positions_list,
                    tissue_lowres_image,
                    tissue_hires_image,
                    scale_factors,
                    barcodes,
                    features,
                    matrix,
                ) = lspl[: len(HEADER)]
                sample_info = [
                    tissue_positions_list,
                    tissue_lowres_image,
                    tissue_hires_image,
                    scale_factors,
                    barcodes,
                    features,
                    matrix,
                ]
            sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            # Create sample mapping dictionary
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    # Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(HEADER) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):
                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join(["{}".format(sample, idx + 1)] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(argv=None):
    args = parse_args(argv)
    check_samplesheet(args.file_in, args.file_out, args.is_raw_data)


if __name__ == "__main__":
    sys.exit(main())
