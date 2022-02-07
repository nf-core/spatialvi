#!/opt/conda/bin/python

import sys
import os
import scanpy as sc
import numpy as np
import pandas as pd
import scanorama
from sklearn.manifold import TSNE
from umap import UMAP
from matplotlib import pyplot as plt
from matplotlib import cm
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import cdist

def label_transfer(dist, labels):
    class_prob = pd.get_dummies(labels).to_numpy().T @ dist
    class_prob /= np.linalg.norm(class_prob, 2, axis=0)
    return (class_prob.T - class_prob.min(1)) / np.ptp(class_prob, axis=1)
   
# See more settings at:
# https://scanpy.readthedocs.io/en/stable/generated/scanpy._settings.ScanpyConfig.html
sc.settings.figdir = sys.argv[1]
sc.set_figure_params(dpi_save=300, facecolor='white')

for dirName in ['show/', 'umap/', 'umap_density_clusters_/', 'violin/']:
    if not os.path.exists(sys.argv[1] + dirName):
        os.makedirs(sys.argv[1] + dirName)

#### Label transfer
#1. To topics from single cells
#2. To spots from single cells
#3. Compare spots transfer to deconvolution results

#Task: Determine cluster composition of LDA topics <br>
#Challenge: topics are intertwined -> hard to match the clusters <br>

#Task: direct label transfer from SC clusters to ST spots <br>
#Prerequisite: Do integration via Scanorama or BBKNN <br>
#Post-task: compare deconvolution results to STdeconvolve proportions (correlate proportions) <br>

# Read pre-processed data
st_adata = sc.read(sys.argv[1] + 'st_adata_norm.h5ad')
sc_adata = sc.read(sys.argv[1] + 'sc_adata_norm.h5ad')
df_st = pd.DataFrame(index=st_adata.var.index, data=st_adata.X.T.todense(), columns=st_adata.obs.index)
df_sc = pd.DataFrame(index=sc_adata.var.index, data=sc_adata.X.T.todense(), columns=sc_adata.obs.index)
ind = df_sc.index.intersection(df_st.index)
df_st = df_st.loc[ind]
df_sc = df_sc.loc[ind]
print(df_st.shape, df_sc.shape)


# Read STdeconvolve results
se_clusters_SC = pd.read_csv(sys.argv[1] + 'STdeconvolve_sc_cluster_ids.csv', index_col=0).rename({'x': 'cluster_id'}, axis=1)['cluster_id']
df_theta = pd.read_csv(sys.argv[1] + 'STdeconvolve_prop_norm.csv', index_col=0)
df_beta = pd.read_csv(sys.argv[1] + 'STdeconvolve_beta_norm.csv', index_col=0)
hvg_SC = pd.read_csv(sys.argv[1] + 'STdeconvolve_sc_pca_feature_loadings.csv', index_col=0).index
df_markers_SC = pd.read_csv(sys.argv[1] + 'STdeconvolve_sc_cluster_markers.csv', index_col=0)


# Simple comparison of clusters to topics
df_beta_st = df_beta.loc[df_beta.index.intersection(df_st.index)]
df_beta_st.columns = 'X ' + df_beta_st.columns.astype(str)
df_sc_cl = df_sc.copy().loc[df_beta.index.intersection(df_st.index)]
df_sc_cl.columns = pd.MultiIndex.from_arrays(se_clusters_SC.reset_index().values.T, names=['cell', 'cluster'])
df_sc_cl = df_sc_cl.groupby(level='cluster', axis=1).mean()
df_sc_cl.columns = 'Cluster ' + df_sc_cl.columns.astype(str)
# Similar approach to label transfer
df = df_sc_cl.T @ df_beta_st / df_sc_cl.shape[0]
df /= np.linalg.norm(df, 2, axis=0)
df = (df.T - df.min(1)) / df.values.ptp(1)
df.T.round(3) #.style.background_gradient(axis=None)

# Integrating with BBKNN
# This is only to see the UMAP of SC and ST together
# Not using it with label transfer
if True: 
    st_adata = sc.read(sys.argv[1] + 'st_adata_norm.h5ad')
    sc_adata = sc.read(sys.argv[1] + 'sc_adata_norm.h5ad')
    adata_all = st_adata.concatenate(sc_adata, batch_categories=['st', 'sc'])
    sc.pp.highly_variable_genes(adata_all, flavor="seurat", n_top_genes=2000)
    sc.pp.pca(adata_all, use_highly_variable=True)
    sc.pp.neighbors(adata_all)
    sc.tl.umap(adata_all)
    sc.pl.umap(adata_all, color=['batch'], palette=sc.pl.palettes.vega_20_scanpy, save='/scanpy_UMAP_st_sc.png')
    df_before = pd.DataFrame(adata_all.obsp['connectivities'].todense(), index=adata_all.obs.batch, columns=adata_all.obs.batch).xs('st', axis=0).xs('sc', axis=1)
    sc.external.pp.bbknn(adata_all, batch_key='batch')
    sc.tl.umap(adata_all)
    sc.pl.umap(adata_all, color=['batch'], palette=sc.pl.palettes.vega_20_scanpy, save='/scanpy_UMAP_st_sc_BBKNN.png')
    df_after = pd.DataFrame(adata_all.obsp['connectivities'].todense(), index=adata_all.obs.batch, columns=adata_all.obs.batch).xs('st', axis=0).xs('sc', axis=1)

# Correlation similarity of gene expression in PC space
# Integrating with Scanorama
integrated = scanorama.integrate([df_st.values.T, df_sc.values.T], [df_st.index.values, df_sc.index.values])[0]
umap_st = UMAP(random_state=2).fit_transform(integrated[0]).T
umap_sc = UMAP(random_state=2).fit_transform(integrated[1]).T
fig, ax = plt.subplots(1, 1, figsize=(6, 6))
ax.scatter(umap_st[0], umap_st[1], alpha=0.95, s=2, label='st')
ax.scatter(umap_sc[0], umap_sc[1], alpha=0.95, s=2, label='sc')
plt.legend()
fig.savefig(sys.argv[1] + '/scanorama_UMAP_st_sc.png', facecolor='white', dpi=300);
plt.close(fig)

df_similarity = pd.DataFrame(data=1 - cdist(integrated[0], integrated[1], metric='correlation').T).fillna(0)
v_in_pc = label_transfer(df_similarity.values, se_clusters_SC.values)

# Plot individual histograms
if True:
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    pd.DataFrame(v_in_pc).hist(bins=50, ax=ax);
    fig.savefig(sys.argv[1] + '/LT_individual_histograms.png', facecolor='white', dpi=300);
    plt.close(fig)

# Plot combined histogram
if True:
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    for i in range(v_in_pc.shape[1]):
        pd.DataFrame(v_in_pc)[i].hist(bins=100, alpha=0.5, ax=ax)
    fig.savefig(sys.argv[1] + '/LT_individual_histograms_combined.png', facecolor='white', dpi=300);
    plt.close(fig)

# Add LT results to adata 
dfb = pd.DataFrame(v_in_pc)
dfb.columns = 'LT PC ' + dfb.columns.astype(str)
dfb.index = st_adata.obs.index.values
for col in dfb.columns:
    st_adata.obs[col] = dfb[col]
    
#dfc = pd.DataFrame(v_in_bbknn)
#dfc.columns = 'LT BBKNN ' + dfc.columns.astype(str)
#dfc.index = st_adata.obs.index.values
#for col in dfc.columns:
#    st_adata.obs[col] = dfc[col]

# Add LDA proportions to adata
df_theta = pd.read_csv(sys.argv[1] + 'STdeconvolve_prop_norm.csv', index_col=0)
df_theta.columns = 'Topic LDA ' + df_theta.columns.astype(str)
for col in df_theta.columns:
    st_adata.obs[col] = df_theta[col]   

# Add NMF proportions to adata
df_theta = pd.read_csv(sys.argv[1] + 'SPOTlight_prop_norm.csv', index_col=0).set_index('barcodes').drop('res_ss', axis=1)
df_theta.columns = 'Topic NMF ' + df_theta.columns.astype(str)
for col in df_theta.columns:
    st_adata.obs[col] = df_theta[col]

# Make plots with proportions
# These plots have a uniform color-mapping of scanpy
if True:
    plt.rcParams["figure.figsize"] = (5, 5)
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic LDA ')]
    sc.pl.spatial(st_adata, img_key="hires", color=keys, ncols=10, save='/Topics_LDA_spatial.png')
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic NMF ')]
    sc.pl.spatial(st_adata, img_key="hires", color=keys, ncols=10, save='/Topics_NMF_spatial.png')
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('LT PC ')]
    sc.pl.spatial(st_adata, img_key="hires", color=keys, ncols=10, save='/Topics_LT_PC_spatial.png')
    #keys = st_adata.var.index.intersection(all_markers).values
    #sc.pl.spatial(st_adata, img_key="hires", color=keys, ncols=9, save='/All_Markers_spatial.png')


# Clsutering and Selection
sc.pp.pca(st_adata)
sc.pp.neighbors(st_adata)
sc.tl.umap(st_adata)
sc.tl.leiden(st_adata, key_added="clusters", resolution=0.4)

# Make plots of UMAP of ST spots clusters
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(st_adata, color=["clusters", "total_counts", "n_genes_by_counts"], wspace=0.4, save='/UMAP_st_spots_clusters.png')
sc.tl.embedding_density(st_adata, basis='umap', groupby='clusters')
sc.pl.embedding_density(st_adata, groupby='clusters', ncols=10, save='/UMAP_clusters_embedding_density.png')

# Plot LDA cell topics proportions in UMAP
if True:
    plt.rcParams["figure.figsize"] = (4, 4)
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic LDA ')]
    sc.pl.umap(st_adata, color=keys, wspace=0.4, ncols=5, save='/UMAP_LDA_topics.png')
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic NMF ')]
    sc.pl.umap(st_adata, color=keys, wspace=0.4, ncols=5, save='/UMAP_NMF_topics.png')
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('LT PC ')]
    sc.pl.umap(st_adata, color=keys, wspace=0.4, ncols=5, save='/UMAP_LT_PC.png')

# Make plots of spatial ST spots clusters
if True:
    plt.rcParams["figure.figsize"] = (10, 10)
    sc.pl.spatial(st_adata, img_key="hires", color=["clusters"], save='/Clusters_scanpy_spatial.png') # groups=['1', '2', '3']

# Make violin plots of topic proportions of ST spots by clusters
if True:
    plt.rcParams["figure.figsize"] = (3.5, 3.5)
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic LDA ')]
    sc.pl.violin(st_adata, keys, jitter=0.4, groupby='clusters', rotation=0, save='/violin_topics_LDA.png')
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic NMF ')]
    sc.pl.violin(st_adata, keys, jitter=0.4, groupby='clusters', rotation=0, save='/violin_topics_NMF.png')
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('LT PC ')]
    sc.pl.violin(st_adata, keys, jitter=0.4, groupby='clusters', rotation=0, save='/violin_topics_LT_PC.png')

exit(0)