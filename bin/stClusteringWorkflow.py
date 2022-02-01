#!/opt/conda/bin/python

import sys
import os
import scanpy as sc
import numpy as np
import pandas as pd

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
st_adata = sc.read('st_adata_norm.h5ad')
sc_adata = sc.read('sc_adata_norm.h5ad')
df_st = pd.DataFrame(index=st_adata.var.index, data=st_adata.X.T.todense(), columns=st_adata.obs.index)
df_sc = pd.DataFrame(index=sc_adata.var.index, data=sc_adata.X.T.todense(), columns=sc_adata.obs.index)
ind = df_sc.index.intersection(df_st.index)
df_st = df_st.loc[ind]
df_sc = df_sc.loc[ind]
print(df_st.shape, df_sc.shape)


pddir = 'c:/Projects/A_ST/forPipelineDev/R_e_STdeconvolve_norm/'
se_clusters_SC = pd.read_csv(pddir + 'sc_cluster_ids.csv', index_col=0).rename({'x': 'cluster_id'}, axis=1)['cluster_id']
df_theta = pd.read_csv(pddir + 'theta_norm.csv', index_col=0)
df_beta = pd.read_csv(pddir + 'beta_norm.csv', index_col=0)
hvg_SC = pd.read_csv(pddir + 'sc_pca_feature_loadings.csv', index_col=0).index
df_markers_SC = pd.read_csv(pddir + 'sc_cluster_markers.csv', index_col=0)

df_beta_st = df_beta.loc[df_beta.index.intersection(df_st.index)]
df_beta_st.columns = 'X ' + df_beta_st.columns.astype(str)

df_sc_cl = df_sc.copy().loc[df_beta.index.intersection(df_st.index)]
df_sc_cl.columns = pd.MultiIndex.from_arrays(se_clusters_SC.reset_index().values.T, names=['cell', 'cluster'])
df_sc_cl = df_sc_cl.groupby(level='cluster', axis=1).mean()
df_sc_cl.columns = 'Cluster ' + df_sc_cl.columns.astype(str)

# df = pd.concat([df_sc_cl, df_beta_st], axis=1).corr()[df_sc_cl.columns].loc[df_beta_st.columns]
df = df_sc_cl.T @ df_beta_st / df_sc_cl.shape[0]
df /= np.linalg.norm(df, 2, axis=0)
df = (df.T - df.min(1)) / df.values.ptp(1)
df.T.round(3).style.background_gradient(axis=None)

def label_transfer(dist, labels):
    class_prob = pd.get_dummies(labels).to_numpy().T @ dist
    class_prob /= np.linalg.norm(class_prob, 2, axis=0)
    return (class_prob.T - class_prob.min(1)) / np.ptp(class_prob, axis=1)
    
# Integrating with BBKNN
st_adata = sc.read('st_adata_norm.h5ad')
sc_adata = sc.read('sc_adata_norm.h5ad')
adata_all = st_adata.concatenate(sc_adata, batch_categories=['st', 'sc'])
sc.pp.highly_variable_genes(adata_all, flavor="seurat", n_top_genes=2000)
sc.pp.pca(adata_all, use_highly_variable=True)
sc.pp.neighbors(adata_all)
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=['batch'], palette=sc.pl.palettes.vega_20_scanpy)
df_before = pd.DataFrame(adata_all.obsp['connectivities'].todense(), index=adata_all.obs.batch, columns=adata_all.obs.batch).xs('st', axis=0).xs('sc', axis=1)

sc.external.pp.bbknn(adata_all, batch_key='batch')
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=['batch'], palette=sc.pl.palettes.vega_20_scanpy)
df_after = pd.DataFrame(adata_all.obsp['connectivities'].todense(), index=adata_all.obs.batch, columns=adata_all.obs.batch).xs('st', axis=0).xs('sc', axis=1)

# Integrating with Scanorama
integrated = scanorama.integrate([df_st.values.T, df_sc.values.T], [df_st.index.values, df_sc.index.values])[0]
umap_st = UMAP(random_state=2).fit_transform(integrated[0]).T
umap_sc = UMAP(random_state=2).fit_transform(integrated[1]).T
plt.scatter(umap_st[0], umap_st[1], alpha=0.95, s=2, label='st')
plt.scatter(umap_sc[0], umap_sc[1], alpha=0.95, s=2, label='sc')
plt.legend()

# Correlation similarity of gene expression in PC space
df_similarity = pd.DataFrame(data=1 - cdist(integrated[0], integrated[1], metric='correlation').T).fillna(0)
v_in_pc = label_transfer(df_similarity.values, se_clusters_SC.values)

dfb = pd.DataFrame(v_in_pc)
dfb.columns = 'LT PC ' + dfb.columns.astype(str)
dfb.index = st_adata.obs.index.values
for col in dfb.columns:
    st_adata.obs[col] = dfb[col]
    
dfc = pd.DataFrame(v_in_bbknn)
dfc.columns = 'LT BBKNN ' + dfc.columns.astype(str)
dfc.index = st_adata.obs.index.values
for col in dfc.columns:
    st_adata.obs[col] = dfc[col]

pddir = 'c:/Projects/A_ST/forPipelineDev/R_e_STdeconvolve_norm/'
df_theta = pd.read_csv(pddir + 'theta_norm.csv', index_col=0)
df_theta.columns = 'Topic ' + df_theta.columns.astype(str)
for col in df_theta.columns:
    st_adata.obs[col] = df_theta[col]
    
pddir = 'c:/Projects/A_ST/forPipelineDev/R_e_NMF_norm/'
df_theta = pd.read_csv(pddir + 'prop_norm.csv', index_col=0).set_index('barcodes').drop('res_ss', axis=1)
df_theta.columns = 'Topic ' + df_theta.columns.astype(str)
for col in df_theta.columns:
    st_adata.obs[col] = df_theta[col]
    
plt.rcParams["figure.figsize"] = (5, 5)

sc.pl.spatial(st_adata, img_key="hires", color=["Topic %s" % (i + 1) for i in range(9)], ncols=10) # LDA

sc.pl.spatial(st_adata, img_key="hires", color=["Topic X%s" % (i) for i in range(10)], ncols=10) # NMF

sc.pl.spatial(st_adata, img_key="hires", color=["LT PC %s" % (i + 1) for i in range(9)], ncols=10) # LT

sc.pl.spatial(st_adata, img_key="hires", color=st_adata.var.index.intersection(all_markers).values, ncols=9)




# Clsutering and Selection Workflow
sc.pp.pca(st_adata)
sc.pp.neighbors(st_adata)
sc.tl.umap(st_adata)

sc.tl.leiden(st_adata, key_added="clusters", resolution=0.4)

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(st_adata, color=["clusters", "total_counts", "n_genes_by_counts"], wspace=0.4)

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.embedding_density(st_adata, groupby='clusters', ncols=10)

pddir = 'c:/Projects/A_ST/forPipelineDev/R_e_STdeconvolve_norm/'
df_theta = pd.read_csv(pddir + 'theta_norm.csv', index_col=0)
df_theta.columns = 'Topic ' + df_theta.columns.astype(str)
for col in df_theta.columns:
    st_adata.obs[col] = df_theta[col]

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(st_adata, color=["Topic 1", "Topic 2", "Topic 3", "Topic 4", "Topic 5", "Topic 6", "Topic 7", "Topic 8", "Topic 9"], wspace=0.4, ncols=5)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.spatial(st_adata, img_key="hires", color=["clusters"]) # groups=['1', '2', '3']

plt.rcParams["figure.figsize"] = (5, 5)
sc.pl.spatial(st_adata, img_key="hires", color=["Topic %s" % (i + 1) for i in range(9)], ncols=5)

plt.rcParams["figure.figsize"] = (3.5, 3.5)

sc.pl.violin(st_adata, ["Topic %s" % (i + 1) for i in range(0,5)], jitter=0.4, groupby='clusters', rotation=0)
sc.pl.violin(st_adata, ["Topic %s" % (i + 1) for i in range(5,9)], jitter=0.4, groupby='clusters', rotation=0)


st_adata_sub = st_adata[(st_adata.obsm["spatial"][:, 0] > 3000) & 
                        (st_adata.obsm["spatial"][:, 1] > 8000), :].copy()

plt.rcParams["figure.figsize"] = (5, 5)
sc.pl.spatial(st_adata_sub, img_key="hires", color=["clusters"]) # groups=['1', '2', '3']


