#!/opt/conda/bin/python

import os
import sys
import scanpy as sc
import numpy as np

outsPath = sys.argv[1]
saveFile = sys.argv[2]

print(outsPath)
print(saveFile)

sc_adata = sc.read_10x_mtx(outsPath)
sc_adata.var_names_make_unique()
sc.pp.filter_cells(sc_adata, min_counts=1)
sc.pp.filter_genes(sc_adata, min_cells=1)

if not os.path.exists(os.path.dirname(saveFile)):
    os.makedirs(os.path.dirname(saveFile))

sc_adata.write(saveFile)

X = np.array(sc_adata.X.todense()).T
np.savez_compressed(os.path.dirname(saveFile) + '/sc_adata_counts.npz', X)
