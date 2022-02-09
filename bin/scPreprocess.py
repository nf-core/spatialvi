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

sc_adata = sc.read(sys.argv[1] + '/' + sys.argv[3])
print(sc_adata.shape)

sc_adata.obs['norm_factors'] = pd.Series(index=sc_adata.obs.index, data=f_temp)

mito = pd.read_csv(sys.argv[4], index_col=['Symbol', 'MCARTA2_LIST'], delimiter='\t')['EnsemblGeneID']
mito = mito.xs(1, level='MCARTA2_LIST').sort_index().reset_index()
print(mito)

sc_adata.var["mt"] = sc_adata.var_names.isin(mito['Symbol'])
sc.pp.calculate_qc_metrics(sc_adata, qc_vars=["mt"], inplace=True)
print(sc_adata)

print(sc_adata.var)

plt.rcParams["figure.figsize"] = (6, 6)

def histplotQC(se_data, bins, ax):
    print(se_data.shape)
    ax.hist(se_data, density=True, bins=bins, color='navy', alpha=0.3)
    kde = scipy.stats.gaussian_kde(se_data)
    xx = np.linspace(min(se_data), max(se_data), 300)
    ax.set_xlabel(se_data.name)
    ax.set_ylabel('Density')
    ax.plot(xx, kde(xx), color='crimson')
    ax.set_xlim([0, ax.get_xlim()[1]])
    return

fig, axs = plt.subplots(1, 5, figsize=(20, 4))
histplotQC(sc_adata.obs["total_counts"], bins=40, ax=axs[0])
histplotQC(sc_adata.obs["total_counts"][sc_adata.obs["total_counts"] < 5000], bins=40, ax=axs[1])
histplotQC(sc_adata.obs["n_genes_by_counts"], bins=40, ax=axs[2])
histplotQC(sc_adata.obs["n_genes_by_counts"][sc_adata.obs["n_genes_by_counts"] < 2000], bins=40, ax=axs[3])
histplotQC(sc_adata.obs["pct_counts_mt"], bins=40, ax=axs[4])
fig.tight_layout()
fig.savefig(sys.argv[1] + '/sc_histogrtam_all.png', facecolor='white')

# Filter cells and genes
sc.pp.filter_cells(sc_adata, min_counts=500)
sc.pp.filter_cells(sc_adata, min_genes=250)
sc.pp.filter_genes(sc_adata, min_cells=1)
print(sc_adata.shape)

fig, axs = plt.subplots(1, 5, figsize=(20, 4))
histplotQC(sc_adata.obs["total_counts"], bins=40, ax=axs[0])
histplotQC(sc_adata.obs["total_counts"][sc_adata.obs["total_counts"] < 5000], bins=40, ax=axs[1])
histplotQC(sc_adata.obs["n_genes_by_counts"], bins=40, ax=axs[2])
histplotQC(sc_adata.obs["n_genes_by_counts"][sc_adata.obs["n_genes_by_counts"] < 2000], bins=40, ax=axs[3])
histplotQC(sc_adata.obs["pct_counts_mt"], bins=40, ax=axs[4])
fig.tight_layout()
fig.savefig(sys.argv[1] + '/sc_histogrtam_filtered.png', facecolor='white')

# Effect of normalization by size factors
fig, ax = plt.subplots(figsize=(5, 5))
display_cutoff = 5*10**3
se = pd.Series(np.array(sc_adata.X.sum(axis=1)).T[0])
se = se[se<display_cutoff]
print('Number of cells displayed:', se.shape)
se.hist(bins=100, alpha=0.75, ax=ax)
ax.set_xlim(0, display_cutoff);
sc_adata_c = sc_adata.copy()
sc_adata_c.X = csr_matrix(sc_adata.X / sc_adata.obs['norm_factors'].values[:, None])
se = pd.Series(np.array(sc_adata_c.X.sum(axis=1)).T[0])
se = se[se<display_cutoff]
print('Number of cells displayed:', se.shape)
se.hist(bins=100, alpha=0.75, ax=ax)
ax.set_xlim(0, display_cutoff);
fig.savefig(sys.argv[1] + '/sc_histogram_with_without_normalization.png', facecolor='white', dpi=300);
plt.close(fig)

# Save raw filtered data
sc_adata.write(sys.argv[1] + '/sc_adata_plain.h5ad')

# Save normalized data
sc_adata.X = csr_matrix(sc_adata.X / sc_adata.obs['norm_factors'].values[:, None])
sc.pp.log1p(sc_adata)
sc_adata.write(sys.argv[1] + '/sc_adata_norm.h5ad')

# Save to open in R
np.savez_compressed(sys.argv[1] + '/sc_adata_X.npz', sc_adata.X.T.todense())
sc_adata.var.to_csv(sys.argv[1] + '/sc_adata.var.csv')
sc_adata.obs.to_csv(sys.argv[1] + '/sc_adata.obs.csv')
