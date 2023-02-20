#!/usr/bin/env python

# Load packages
import os
import argparse
import scanpy as sc
import numpy as np
import pandas as pd

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Load scRNA-seq data from MTX.")
parser.add_argument(
    "--SRCountDir", metavar="SRCountDir", type=str, default=None, help="Path to Space Range output directory, etc."
)
parser.add_argument(
    "--outAnnData", metavar="outAnnData", type=str, default=None, help="Path to a file to save h5ad data into."
)
parser.add_argument("--outSCCounts", metavar="outSCCounts", type=str, default=None, help="Name of the NPZ file.")
parser.add_argument("--minCounts", metavar="minCounts", type=int, default=1, help="Min counts per spot.")
parser.add_argument("--minGenes", metavar="minGenes", type=int, default=1, help="Min genes per spot.")
parser.add_argument("--minCells", metavar="minCells", type=int, default=1, help="Min cells per gene.")
args = parser.parse_args()

matrixFilename = os.path.join(args.SRCountDir, "raw_feature_bc_matrix", "matrix.mtx.gz")
featuresFiles = os.path.join(args.SRCountDir, "raw_feature_bc_matrix", "features.tsv.gz")
barcodesFilename = os.path.join(args.SRCountDir, "raw_feature_bc_matrix", "barcodes.tsv.gz")

sc_adata = sc.read_mtx(matrixFilename).T
genes = pd.read_csv(featuresFiles, header=None, sep="\t")
print(genes)
if len(genes.columns) == 1:
    gs = genes[0]
else:
    gs = genes[1]
sc_adata.var_names = gs
sc_adata.var["gene_symbols"] = gs.values
sc_adata.obs_names = pd.read_csv(barcodesFilename, header=None)[0]
print(sc_adata.var)

sc_adata.var_names_make_unique()
sc.pp.filter_cells(sc_adata, min_counts=args.minCounts)
sc.pp.filter_genes(sc_adata, min_cells=args.minCells)

# Save raw anndata to file
sc_adata.write(args.outAnnData)

# Save counts anndata to file
X = np.array(sc_adata.X.todense()).T
np.savez_compressed(args.outSCCounts, X)
