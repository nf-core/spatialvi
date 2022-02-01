#!/opt/conda/bin/python

import sys
import os
import scanpy as sc
import numpy as np
import pandas as pd
import SpatialDE

# See more settings at:
# https://scanpy.readthedocs.io/en/stable/generated/scanpy._settings.ScanpyConfig.html
sc.settings.figdir = sys.argv[1]

if not os.path.exists(sys.argv[1] + 'show/'):
    os.makedirs(sys.argv[1] + 'show/')

st_adata = sc.read(sys.argv[1] + '/' + sys.argv[2])
print(st_adata.shape)

counts = pd.DataFrame(st_adata.X.todense(), columns=st_adata.var_names, index=st_adata.obs_names)
coord = pd.DataFrame(st_adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=st_adata.obs_names)

df_results = SpatialDE.run(coord, counts)

df_results.index = df_results["g"]
df_results = df_results.sort_values("qval", ascending=True)

df_results.to_csv(sys.argv[1] + '/stSpatialDE.csv')

keys = df_results.index.values[:15]
sc.pl.spatial(st_adata, img_key="hires", color=keys, alpha=0.7, save='/st_SpatialDE.png', ncols=5)


