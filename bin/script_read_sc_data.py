#!/opt/conda/bin/python

import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd

outsPath = sys.argv[1]
saveFile = sys.argv[2]

print('outsPath:', '\t', outsPath)
print('saveFile:', '\t', saveFile)

countsFile = ''
files = os.listdir(outsPath)

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
    sc_adata = sc.read_h5ad(outsPath + '/' + countsFile)
else:
    sc_adata = sc.read_mtx(outsPath +'matrix.mtx.gz').T
    genes = pd.read_csv(outsPath + 'features.tsv.gz', header=None, sep='\t')
    print(genes)
    if len(genes.columns)==1:
        gs = genes[0]
    else:
        gs = genes[1]
    sc_adata.var_names = gs
    sc_adata.var['gene_symbols'] = gs.values
    sc_adata.obs_names = pd.read_csv(outsPath + 'barcodes.tsv.gz', header=None)[0]
    print(sc_adata.var)

sc_adata.var_names_make_unique()
sc.pp.filter_cells(sc_adata, min_counts=1)
sc.pp.filter_genes(sc_adata, min_cells=1)

if not os.path.exists(os.path.dirname(saveFile)):
    os.makedirs(os.path.dirname(saveFile))

sc_adata.write(saveFile)

X = np.array(sc_adata.X.todense()).T
np.savez_compressed(os.path.dirname(saveFile) + '/sc_adata_counts.npz', X)
