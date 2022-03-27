#!/opt/conda/bin/python

# Load packages
import os
import sys
import argparse
import scanpy as sc
import numpy as np
import pandas as pd


# Parse command-line arguments
parser = argparse.ArgumentParser(description='Load single cell traqnscriptomics data from MTX.')

parser.add_argument('--outsPath', metavar='outspath', type=str, default=None, help='Path to Space Range outs directory, etc.')
parser.add_argument('--saveFile', metavar='savefile', type=str, default=None, help='Path to a file to save h5ad data into.')
parser.add_argument('--npCountsOutputName', metavar='npzoutput', type=str, default=None, help='Name of the NPZ file.')

parser.add_argument('--minCounts', metavar='cutoff', type=int, default=1, help='Min counts per spot.')
parser.add_argument('--minGenes', metavar='cutoff', type=int, default=1, help='Min genes per spot.')
parser.add_argument('--minCells', metavar='cutoff', type=int, default=1, help='Min cells per gene.')

args = parser.parse_args()


# Main script
countsFile = ''
files = os.listdir(args.outsPath)

for fname in files:
    if 'matrix.mtx' in fname:
        countsFile = fname
        break

if countsFile == '':
    for fname in files:
        if '.h5ad' in fname:
            countsFile = fname
            break

if '.h5ad' in countsFile:
    sc_adata = sc.read_h5ad(args.outsPath + '/' + countsFile)
else:
    sc_adata = sc.read_mtx(args.outsPath +'matrix.mtx.gz').T
    genes = pd.read_csv(args.outsPath + 'features.tsv.gz', header=None, sep='\t')
    print(genes)
    if len(genes.columns)==1:
        gs = genes[0]
    else:
        gs = genes[1]
    sc_adata.var_names = gs
    sc_adata.var['gene_symbols'] = gs.values
    sc_adata.obs_names = pd.read_csv(args.outsPath + 'barcodes.tsv.gz', header=None)[0]
    print(sc_adata.var)

sc_adata.var_names_make_unique()
sc.pp.filter_cells(sc_adata, min_counts=args.minCounts)
sc.pp.filter_genes(sc_adata, min_cells=args.minCells)

if not os.path.exists(os.path.dirname(args.saveFile)):
    os.makedirs(os.path.dirname(args.saveFile))

sc_adata.write(args.saveFile)

X = np.array(sc_adata.X.todense()).T
np.savez_compressed(os.path.dirname(args.saveFile) + '/' + args.npCountsOutputName, X)

exit(0)
