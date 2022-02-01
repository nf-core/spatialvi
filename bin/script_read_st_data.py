#!/opt/conda/bin/python

import os
import sys
import scanpy as sc
import numpy as np

outsPath = sys.argv[1]
saveFile = sys.argv[2]
countsFile = 'filtered_feature_bc_matrix.h5'

for fname in os.listdir(outsPath):
    if 'filtered_feature_bc_matrix.h5' in fname:
        countsFile = fname
        break

print(outsPath)
print(countsFile)
print(saveFile)

st_adata = sc.read_visium(outsPath,
                          count_file=countsFile,
                          library_id=None,
                          load_images=True,
                          source_image_path=None)

st_adata.var_names_make_unique()
sc.pp.filter_cells(st_adata, min_counts=1)
sc.pp.filter_genes(st_adata, min_cells=1)

if not os.path.exists(os.path.dirname(saveFile)):
    os.makedirs(os.path.dirname(saveFile))

st_adata.write(saveFile)

X = np.array(st_adata[st_adata.obs['in_tissue']==1].X.todense()).T
np.savez_compressed(os.path.dirname(saveFile) + '/st_adata_counts_in_tissue.npz', X)
