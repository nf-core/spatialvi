#!/usr/local/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

set.seed(123)

library(Matrix)
library(data.table)
library(Seurat)
#library(SeuratData)
library(dplyr)
library(gt) ###########
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)

normDataDir <- args[1]

library(reticulate)
np <- import("numpy")

filename <- list.files(path=args[2], pattern="filtered_feature_bc_matrix.h5")[1]
print(list.files(path=args[2], pattern="filtered_feature_bc_matrix.h5"))
print(list.files(path=args[2], pattern="filtered_feature_bc_matrix.h5")[1])
print(filename)
print(args[2])

se_st <- Seurat::Load10X_Spatial(data.dir = args[2], filename = filename)
matrix_st <- np$load(paste0(normDataDir, 'st_adata_X.npz'))[['arr_0']]
st_genes <- read.csv(paste0(normDataDir, 'st_adata.var.csv'))$X
st_obs <- read.csv(paste0(normDataDir, 'st_adata.obs.csv'))$X
rownames(matrix_st) <- st_genes
colnames(matrix_st) <- st_obs
se_st@assays$Spatial@counts <- as(100*matrix_st, "sparseMatrix")
se_st@assays$Spatial@data <- as(100*matrix_st, "sparseMatrix")

#se_sc <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(data.dir = 'c:/Projects/A_ST/sc mouse kidney/SRX3436301/outs/filtered_feature_bc_matrix/'))
matrix_sc <- np$load(paste0(normDataDir, 'sc_adata_X.npz'))[['arr_0']]
sc_genes <- read.csv(paste0(normDataDir, 'sc_adata.var.csv'))$X
sc_obs <- read.csv(paste0(normDataDir, 'sc_adata.obs.csv'))$X
rownames(matrix_sc) <- sc_genes
colnames(matrix_sc) <- sc_obs
se_sc <- Seurat::CreateSeuratObject(counts = as(100*matrix_sc, "sparseMatrix"))


# Cluster sc data
se_sc <- Seurat::FindVariableFeatures(se_sc, verbose = FALSE)
se_sc <- Seurat::ScaleData(se_sc, verbose = FALSE)
se_sc <- Seurat::RunPCA(se_sc, verbose = FALSE)
se_sc <- Seurat::RunUMAP(se_sc, dims = 1:30, verbose = FALSE)
se_sc <- Seurat::FindNeighbors(se_sc)
se_sc <- Seurat::FindClusters(se_sc, resolution=0.3)
Seurat::DimPlot(se_sc, group.by = "seurat_clusters", label = TRUE) + Seurat::NoLegend()
ggsave(paste0(args[1], "SPOTlight_sc_clusters.png"), dpi=600, scale=0.5, width=8, height=8, units='in')

cluster_markers_all <- Seurat::FindAllMarkers(object = se_sc, assay = NULL, slot = "data", verbose = TRUE, test.use = "wilcox", only.pos = TRUE)

spotlight_ls <- SPOTlight::spotlight_deconvolution(
  se_sc = se_sc,				# Single cell dataset
  counts_spatial = se_st@assays$Spatial@counts,	# Spatial dataset count
  clust_vr = "seurat_clusters", 		# Variable in sc_sc containing the cell-type annotation
  cluster_markers = cluster_markers_all, 	# Dataframe with the marker genes
  cl_n = 100, 					# Number of cells per cell type to use
  hvg = 3000, 					# Number of HVG to use
  ntop = NULL, 					# Number of top marker genes to use (by default all)
  transf = "uv", 				# Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", 				# Factorization method
  min_cont = 0)					# Remove those cells contributing to a spot below a certain threshold 

saveRDS(object = spotlight_ls, file = paste0(args[1], "SPOTlight_ls_mk_normed.rds"))

# spotlight_ls <- readRDS(file = paste0(args[1], "SPOTlight_ls_mk_normed.rds"))

# NMF topic profiles
nmf_mod <- spotlight_ls[[1]]
h <- NMF::coef(nmf_mod[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(h = h, train_cell_clust = nmf_mod[[2]])
topic_profile_plts[[2]] + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0), axis.text = ggplot2::element_text(size = 12))
ggsave(paste0(args[1], "SPOTlight_st_topic_profiles.png"), dpi=300, scale=0.5, width=10, height=8, units='in')

# Add deconvolution results to se_st meta.data
decon_mtrx <- spotlight_ls[[2]]
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx)!="res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rm(decon_mtrx_sub)
rownames(decon_mtrx) <- colnames(se_st) 
decon_df <- decon_mtrx %>% data.frame() %>% tibble::rownames_to_column("barcodes")
se_st@meta.data <- se_st@meta.data %>% tibble::rownames_to_column("barcodes") %>% dplyr::left_join(decon_df, by = "barcodes") %>% tibble::column_to_rownames("barcodes")

# Individual cell types on image
Seurat::SpatialFeaturePlot(object = se_st, features = colnames(decon_df)[-1][-length(colnames(decon_df))+1], alpha = c(0.1, 1), min.cutoff=0, max.cutoff=0.3, crop = FALSE, pt.size.factor=1.0)
ggsave(paste0(args[1], "SPOTlight_st_prop.png"), dpi=500, scale=1.25, width=8, height=8, units='in')


cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
cell_types_all <- paste("X", cell_types_all, sep = "")
SPOTlight::spatial_scatterpie(se_obj = se_st, cell_types_all = cell_types_all,
                              img_path = paste0(args[2], "spatial/tissue_lowres_image.png"),
                              cell_types_interest = NULL, slice = NULL, scatterpie_alpha = 1, pie_scale = 0.35)
# [c(1,7,3,4,5,6,2,8,9)]
# rainbow(length(cell_types_all))
ggsave(paste0(args[1], "SPOTlight_st_scatterpies.png"), dpi=600, scale=1.0, width=8, height=8, units='in')

# Remove cell types not predicted to be on the tissue
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]]
decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]
colnames(decon_mtrx_sub) <- paste("X", colnames(decon_mtrx_sub), sep = "")
decon_cor <- cor(decon_mtrx_sub)
p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)

# Visualize
ggcorrplot::ggcorrplot(corr = decon_cor, p.mat = p.mat[[1]], hc.order = TRUE, type = "full", insig = "blank",
  lab = TRUE, outline.col = "lightgrey", method = "square", colors = c("#6D9EC1", "white", "#E46726"),
  title = "Cell type proportions correlation", 
  legend.title = "Correlation\n(Pearson)") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(angle = 90), axis.text = ggplot2::element_text(size = 18, vjust = 0.5))
ggsave(paste0(args[1], "SPOTlight_st_prop_corr.png"), dpi=600, scale=0.75, width=8, height=8, units='in')

write.csv(decon_df, file=paste0(args[1], "SPOTlight_prop_norm.csv"))
#write.csv(cbind(cluster_id=nmf_mod[[2]], h_ds), file=paste0(args[1], "SPOTlight_beta_norm.csv"))

write.csv(se_sc@active.ident, file=paste0(args[1], "SPOTlight_sc_cluster_ids.csv"))
write.csv(se_sc@reductions[["pca"]]@cell.embeddings, file=paste0(args[1], "SPOTlight_sc_pca.csv"))
write.csv(se_sc@reductions[["pca"]]@feature.loadings, file=paste0(args[1], "SPOTlight_sc_pca_feature_loadings.csv"))
write.csv(cluster_markers_all, file=paste0(args[1], "SPOTlight_sc_cluster_markers.csv"))

quit(status=0)