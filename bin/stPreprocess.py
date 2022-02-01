#!/opt/conda/bin/python

import sys
import os
import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scipy.stats
from scipy.sparse import csr_matrix

# See more settings at:
# https://scanpy.readthedocs.io/en/stable/generated/scanpy._settings.ScanpyConfig.html
sc.settings.figdir = sys.argv[1]

if not os.path.exists(sys.argv[1] + 'show/'):
    os.makedirs(sys.argv[1] + 'show/')

f_temp = np.load(sys.argv[1] + '/' + sys.argv[2] + '_factors.npz')
f_temp = f_temp[list(f_temp.keys())[0]]
print(f_temp.shape)

st_adata = sc.read(sys.argv[1] + '/' + sys.argv[3])
print(st_adata.shape)

#np.array(st_adata[st_adata.obs['in_tissue']==1].X.todense()).T

st_adata.obs['norm_factors'] = pd.Series(index=st_adata.obs[st_adata.obs['in_tissue']==1].index, data=f_temp).reindex(st_adata.obs.index)

mito = pd.read_csv(sys.argv[4], index_col=['Symbol', 'MCARTA2_LIST'], delimiter='\t')['EnsemblGeneID']
mito = mito.xs(1, level='MCARTA2_LIST').sort_index().reset_index()
print(mito)

st_adata.var["mt"] = st_adata.var_names.isin(mito['Symbol'])
sc.pp.calculate_qc_metrics(st_adata, qc_vars=["mt"], inplace=True)

plt.rcParams["figure.figsize"] = (6, 6)

keys = ["in_tissue", "pct_counts_mt", "total_counts", "n_genes_by_counts"]
st_adata_in = st_adata[st_adata.obs['in_tissue']==1]
sc.pl.spatial(st_adata_in, img_key="hires", color=keys, save='/st_QC_in.png')

keys = ["pct_counts_mt", "total_counts", "n_genes_by_counts"]
st_adata_out = st_adata[st_adata.obs['in_tissue']!=1]
sc.pp.filter_cells(st_adata_out, min_counts=500)
sc.pp.filter_cells(st_adata_out, min_genes=250)
sc.pp.filter_genes(st_adata_out, min_cells=1)
sc.pl.spatial(st_adata_out, img_key="hires", color=keys, save='/st_QC_out.png')

def histplotQC(se_data, bins, ax):
    ax.hist(se_data, density=True, bins=bins, color='navy', alpha=0.3)
    kde = scipy.stats.gaussian_kde(se_data)
    xx = np.linspace(min(se_data), max(se_data), 300)
    ax.set_xlabel(se_data.name)
    ax.set_ylabel('Density')
    ax.plot(xx, kde(xx), color='crimson')
    ax.set_xlim([0, ax.get_xlim()[1]])
    return
    
fig, axs = plt.subplots(1, 5, figsize=(20, 4))
histplotQC(st_adata.obs["total_counts"], bins=40, ax=axs[0])
histplotQC(st_adata.obs["total_counts"][st_adata.obs["total_counts"] < 10000], bins=40, ax=axs[1])
histplotQC(st_adata.obs["n_genes_by_counts"], bins=40, ax=axs[2])
histplotQC(st_adata.obs["n_genes_by_counts"][st_adata.obs["n_genes_by_counts"] < 4000], bins=40, ax=axs[3])
histplotQC(st_adata.obs["pct_counts_mt"], bins=40, ax=axs[4])
fig.tight_layout()
fig.savefig(sys.argv[1] + '/st_histogrtam_all.png', facecolor='white')

fig, axs = plt.subplots(1, 5, figsize=(20, 4))
histplotQC(st_adata[st_adata.obs['in_tissue']==1].obs["total_counts"], bins=40, ax=axs[0])
histplotQC(st_adata[st_adata.obs['in_tissue']==1].obs["total_counts"][st_adata[st_adata.obs['in_tissue']==1].obs["total_counts"] < 10000], bins=40, ax=axs[1])
histplotQC(st_adata[st_adata.obs['in_tissue']==1].obs["n_genes_by_counts"], bins=40, ax=axs[2])
histplotQC(st_adata[st_adata.obs['in_tissue']==1].obs["n_genes_by_counts"][st_adata[st_adata.obs['in_tissue']==1].obs["n_genes_by_counts"] < 4000], bins=40, ax=axs[3])
histplotQC(st_adata[st_adata.obs['in_tissue']==1].obs["pct_counts_mt"], bins=40, ax=axs[4])
fig.tight_layout()
fig.savefig(sys.argv[1] + '/st_histogrtam_in.png', facecolor='white')
plt.close(fig)

# TO DO: 'sc_histogrtam_out.png'

# Remove spots outside tissue
st_adata = st_adata[st_adata.obs['in_tissue']==1]
print(st_adata.shape)

# Effect of normalization by size factors
fig, ax = plt.subplots(figsize=(5, 5))
display_cutoff = 10**5
se = pd.Series(np.array(st_adata.X.sum(axis=1)).T[0])
se = se[se<display_cutoff]
print('Number of spots displayed:', se.shape)
se.hist(bins=100, alpha=0.75, ax=ax)
ax.set_xlim(0, display_cutoff);
st_adata_c = st_adata[st_adata.obs['in_tissue']==1].copy()
st_adata_c.X = csr_matrix(st_adata.X / st_adata.obs['norm_factors'].values[:, None])
se = pd.Series(np.array(st_adata_c.X.sum(axis=1)).T[0])
se = se[se<display_cutoff]
print('Number of spots displayed:', se.shape)
se.hist(bins=100, alpha=0.75, ax=ax)
ax.set_xlim(0, display_cutoff);
fig.savefig(sys.argv[1] + '/st_histogram_with_without_normalization.png', facecolor='white', dpi=300);
plt.close(fig)

# Save raw filtered data
st_adata.write(sys.argv[1] + '/st_adata_plain.h5ad')

# Save normalized data
st_adata.X = csr_matrix(st_adata.X / st_adata.obs['norm_factors'].values[:, None])
sc.pp.log1p(st_adata)
st_adata.write(sys.argv[1] + '/st_adata_norm.h5ad')

# Save to open in R
np.savez_compressed(sys.argv[1] + '/st_adata_X.npz', st_adata.X.T.todense())
st_adata.var.to_csv(sys.argv[1] + '/st_adata.var.csv')
st_adata.obs.to_csv(sys.argv[1] + '/st_adata.obs.csv')
