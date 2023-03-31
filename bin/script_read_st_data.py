#!/usr/bin/env python

# Load packages
import argparse
import scanpy as sc
import numpy as np
import read_visium_mtx

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Load spatial transcriptomics data from MTX " + "matrices and aligned images."
)
parser.add_argument(
    "--SRCountDir", metavar="SRCountDir", type=str, default=None, help="Path to Space Range output directory, etc."
)
parser.add_argument(
    "--outAnnData", metavar="outAnnData", type=str, default=None, help="Path to a file to save h5ad data into."
)
args = parser.parse_args()

# Main script
st_adata = read_visium_mtx.read_visium_mtx(args.SRCountDir, library_id=None, load_images=True, source_image_path=None)

# Save raw anndata to file
st_adata.write(args.outAnnData)
