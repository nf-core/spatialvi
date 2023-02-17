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
parser.add_argument("--outSTCounts", metavar="outSTCounts", type=str, default=None, help="Name of the NPZ file.")
parser.add_argument("--minCounts", metavar="minCounts", type=int, default=1, help="Min counts per spot.")
parser.add_argument("--minCells", metavar="minCells", type=int, default=1, help="Min cells per gene.")
args = parser.parse_args()

# Main script
st_adata = read_visium_mtx.read_visium_mtx(args.SRCountDir, library_id=None, load_images=True, source_image_path=None)

st_adata.var_names_make_unique()
sc.pp.filter_cells(st_adata, min_counts=args.minCounts)
sc.pp.filter_genes(st_adata, min_cells=args.minCells)

# Save raw anndata to file
st_adata.write(args.outAnnData)

# Save counts anndata to file
X = np.array(st_adata[st_adata.obs["in_tissue"] == 1].X.todense()).T
np.savez_compressed(args.outSTCounts, X)
