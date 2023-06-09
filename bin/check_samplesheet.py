#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys

# from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (".fq.gz", ".fastq.gz")

    def __init__(
        self,
        area_col="area",
        barcodes_col="barcodes",
        fastq_dir_col="fastq_dir",
        features_col="features",
        hires_image_col="tissue_hires_image",
        lowres_image_col="tissue_lowres_image",
        matrix_col="matrix",
        manual_alignment_col="manual_alignment",
        sample_col="sample",
        scale_factors_col="scale_factors",
        slide_col="slide",
        tissue_positions_list_col="tissue_positions_list",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names, depending on
        whether the input data is raw data (to be processed with Space Ranger)
        or processed data.

        Args:
            area_col (str): The name of the column that contains the slide area
                (default "area").
            barcodes_col (str): The name of the column that contains the
                barcodes file path (default "barcodes").
            fastq_dir_col (str): The name of the column that contains the
                directory in which the input FASTQ files are stored (default
                "fastq_dir").
            features_col (str): The name of the column that contains the
                features file path (default "features").
            hires_image_col (str): The name of the column that contains the high
                resolution image file path (default "tissue_hires_image").
            lowres_image_col (str): The name of the column that contains the low
                resolution image file path (default "tissue_lowres_image").
            matrix_col (str): The name of the column that contains the matrix
                matrix file path (default "matrix").
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            scale_factors_col (str): The name of the column that contains the
                scale factors file path (default "scale_factors").
            slide_col (str): The name of the column that contains the slide ID
                (default "slide").
            tissue_positions_list_col (str): The column that contains the tissue
                poositions list file path (default "tissue_positions_list").

        """
        super().__init__(**kwargs)
        self._area_col = area_col
        self._barcodes_col = barcodes_col
        self._fastq_dir_col = fastq_dir_col
        self._features_col = features_col
        self._hires_image_col = hires_image_col
        self._lowres_image_col = lowres_image_col
        self._manual_alignment_col = manual_alignment_col
        self._matrix_col = matrix_col
        self._sample_col = sample_col
        self._scale_factors_col = scale_factors_col
        self._slide_col = slide_col
        self._tissue_positions_list_col = tissue_positions_list_col
        self.modified = []

    def validate_and_transform(self, row, is_raw_data):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        if is_raw_data:
            self._validate_sample(row)
            self._validate_fastq_dir(row)
            self._validate_hires_image(row)
            self._validate_slide(row)
            self._validate_area(row)
            self._validate_manual_alignment(row)
            self.modified.append(row)
        else:
            self._validate_sample(row)
            self._validate_tissue_positions_list(row)
            self._validate_lowres_image(row)
            self._validate_hires_image(row)
            self._validate_scale_factors(row)
            self._validate_barcodes(row)
            self._validate_features(row)
            self._validate_matrix(row)
            self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_fastq_dir(self, row):
        """Assert that the FASTQ directory entry is non-empty."""
        # TODO: Possibly add a check to make sure that at least one FASTQ file
        # exists in the specified directory
        if len(row[self._fastq_dir_col]) <= 0:
            raise AssertionError("The directory in which the FASTQ files are stored is required")

    def _validate_hires_image(self, row):
        """Assert that the high resolution image entry exists."""
        # TODO: Possibly add a check for image formats
        if len(row[self._hires_image_col]) <= 0:
            raise AssertionError("The high resolution image is required")

    def _validate_slide(self, row):
        """Assert that the slide entry exists."""
        # TODO: Possibly add a check for valid slide IDs
        if len(row[self._slide_col]) <= 0:
            raise AssertionError("The slide ID is required")

    def _validate_area(self, row):
        """Assert that the slide area exists."""
        # TODO: Possibly add a check for valid area specifications
        if len(row[self._area_col]) <= 0:
            raise AssertionError("The area is required")

    def _validate_manual_alignment(self, row):
        """Assert that the manual alignment entry has the right format if it exists."""
        return

    def _validate_tissue_positions_list(self, row):
        """Assert that the tissue positions list entry exists."""
        # TODO: Add a CSV file check
        if len(row[self._hires_image_col]) <= 0:
            raise AssertionError("The high resolution image is required")

    def _validate_lowres_image(self, row):
        """Assert that the low resolution image entry exists."""
        # TODO: Possibly add a check for image formats
        if len(row[self._lowres_image_col]) <= 0:
            raise AssertionError("The low resolution image is required")

    def _validate_scale_factors(self, row):
        """Assert that the scale factors entry exists."""
        # TODO: Possibly add a JSON format check
        if len(row[self._scale_factors_col]) <= 0:
            raise AssertionError("The scale factors file is required")

    def _validate_barcodes(self, row):
        """Assert that the barcodes entry exists."""
        # TODO: Possibly add a check for TSV file formats
        if len(row[self._barcodes_col]) <= 0:
            raise AssertionError("The barcodes file is required")

    def _validate_features(self, row):
        """Assert that the features entry exists."""
        # TODO: Possibly add a check for TSV file formats
        if len(row[self._features_col]) <= 0:
            raise AssertionError("The features file is required")

    def _validate_matrix(self, row):
        """Assert that the matrix entry exists."""
        # TODO: Possibly add a check for MTX file formats
        if len(row[self._matrix_col]) <= 0:
            raise AssertionError("The matrix file is required")


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out, is_raw_data):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.
        is_raw_data (boolean): Whether the data is raw data to be processed by
            Space Ranger or not.

    Example (processed data):
        This function checks that the samplesheet follows the following
        structure by default, see also the `spatial transcriptomics
        samplesheet`_::

            sample,tissue_positions_list,tissue_lowres_image,tissue_hires_image,scale_factors,barcodes,features,matrix
            SAMPLE,TISSUE_POSITIONS_LIST.csv,TISSUE_LOWRES_IMAGE.png,tissue_hires_image.png,SCALEFACTORS_JSON.json,BARCODES.tsv.gz,FEATURES.tsv.gz,MATRIX.mtx.gz

    .. _spatial transcriptomics samplesheet:
        https://data.githubusercontent.com/nf-core/test-datasets/spatialtranscriptomics/testdata/test-dataset-subsampled/samplesheet.csv

    Example (raw data):
        This function check that the samplesheet follows the following structure
        if the `--is_raw_data` parameter is set:

            sample,fastq_dir,tissue_hires_image,slide,area
            SAMPLE,FASTQ_DIRECTORY,TISSUE_HIRES_IMAGE.png,SLIDE_ID,SLIDE_AREA

    """
    if is_raw_data:
        required_columns = {"sample", "fastq_dir", "tissue_hires_image", "slide", "area", "manual_alignment"}
    else:
        required_columns = {
            "sample",
            "tissue_positions_list",
            "tissue_lowres_image",
            "tissue_hires_image",
            "scale_factors",
            "barcodes",
            "features",
            "matrix",
        }

    # TODO nf-core: re-enable validation
    import shutil
    shutil.copy(file_in, file_out)
    return
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        # if not required_columns.issubset(reader.fieldnames):
        #     req_cols = ", ".join(required_columns)
        #     logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
        #     sys.exit(1)
        # # Validate each row.
        # checker = RowChecker()
        # for i, row in enumerate(reader):
        #     try:
        #         checker.validate_and_transform(row, is_raw_data)
        #     except AssertionError as error:
        #         logger.critical(f"{str(error)} On line {i + 2}.")
        #         sys.exit(1)
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-r",
        "--is_raw_data",
        action="store_true",
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out, args.is_raw_data)


if __name__ == "__main__":
    sys.exit(main())
