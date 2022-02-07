#!/usr/local/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggplot2)
library(Seurat)
library(STdeconvolve)
library(corrplot)

library(reticulate)
np <- import("numpy")
normDataDir <- args[1]

filename <- list.files(path=args[2], pattern="filtered_feature_bc_matrix.h5")[1]
print(filename)
print(args[2])

#se_st <- Seurat::Load10X_Spatial(data.dir = args[2], filename = filename)
image <- Read10X_Image(image.dir=file.path(args[2], 'spatial'), filter.matrix=TRUE)
m <- Read10X(paste0(args[2], "raw_feature_bc_matrix/"), gene.column=2)
m <- m[,row.names(image@coordinates)]
m <- m[,colSums(m)>0]
se_st <- CreateSeuratObject(counts=m, assay="Spatial")
image <- image[Cells(x=se_st)]
DefaultAssay(object=image) <- "Spatial"
se_st[["slice1"]] <- image

print(dim(se_st))
print(sum(colSums(se_st@assays$Spatial@counts)==0))

matrix_st <- np$load(paste0(normDataDir, 'st_adata_X.npz'))[['arr_0']]
st_genes <- read.csv(paste0(normDataDir, 'st_adata.var.csv'))$X
st_obs <- read.csv(paste0(normDataDir, 'st_adata.obs.csv'))$X
rownames(matrix_st) <- st_genes
colnames(matrix_st) <- st_obs
se_st@assays$Spatial@counts <- as(matrix_st, "sparseMatrix")
se_st@assays$Spatial@data <- as(matrix_st, "sparseMatrix")
se_st@assays$Spatial@counts <- round(100*se_st@assays$Spatial@counts)

print(dim(se_st))
print(sum(colSums(se_st@assays$Spatial@counts)==0))

corpus <- restrictCorpus(se_st@assays$Spatial@counts, removeAbove=1.0, removeBelow = 0.05)
corpus <- corpus + 1
print(dim(as.matrix(corpus)))
print(sum(colSums(as.matrix(corpus))==0))

ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(8, 9, by = 1))
optLDA <- optimalModel(models = ldas, opt = "min")
saveRDS(object = optLDA, file = paste0(args[1], "STdeconvolve_optLDA.rds"))

optLDA <- readRDS(file = paste0(args[1], "STdeconvolve_optLDA.rds"))

results <- getBetaTheta(optLDA, t(as.matrix(corpus)))
deconProp <- results$theta
posk <- se_st@images[["slice1"]]@coordinates[c("imagerow", "imagecol")]
names(posk) <- c('y', 'x')
posk$x <- posk$x * se_st@images[["slice1"]]@scale.factors[["lowres"]]
posk$y <- 600 - posk$y * se_st@images[["slice1"]]@scale.factors[["lowres"]]
vizAllTopics(deconProp, posk, r=2.85, lwd=0, overlay=se_st@images[["slice1"]]@image)
ggsave(paste0(args[1], "STdeconvolve_st_scatterpies.png"), dpi=600, scale=1.0, width=8, height=8, units='in')


# Individual topics proportions spatial overlay
decon_df <- deconProp %>% data.frame() %>% tibble::rownames_to_column("barcodes")

print(colnames(decon_df))
print(dim(decon_df))

se_st@meta.data <- se_st@meta.data %>% tibble::rownames_to_column("barcodes") %>% dplyr::left_join(decon_df, by = "barcodes") %>% tibble::column_to_rownames("barcodes")
#saveRDS(object = se_st, file = paste0(args[1], "SpatialFeaturePlot_se_st.rds"))
#saveRDS(object = colnames(decon_df)[-1], file = paste0(args[1], "SpatialFeaturePlot_cols.rds"))
Seurat::SpatialFeaturePlot(object = se_st, features = colnames(decon_df)[-1], alpha = c(0.1, 1), min.cutoff=0, max.cutoff=0.3, crop = FALSE, pt.size.factor=1.0)
ggsave(paste0(args[1], "STdeconvolve_st_prop.png"), dpi=300, scale=2.0, width=8, height=8, units='in')

# Cell type proportions correlation
decon_mtrx_sub <- deconProp[, colnames(deconProp)[which(colnames(deconProp) != "res_ss")]]
decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]
colnames(decon_mtrx_sub) <- paste("X", colnames(decon_mtrx_sub), sep = "")
decon_cor <- cor(decon_mtrx_sub)
p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)
ggcorrplot::ggcorrplot(corr = decon_cor, p.mat = p.mat[[1]], hc.order = TRUE, type = "full", insig = "blank",
  lab = TRUE, outline.col = "lightgrey", method = "square", colors = c("#6D9EC1", "white", "#E46726"),
  title = "Cell type proportions correlation", legend.title = "Correlation\n(Pearson)") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(angle = 90), axis.text = ggplot2::element_text(size = 18, vjust = 0.5))
ggsave(paste0(args[1], "STdeconvolve_st_prop_corr.png"), dpi=600, scale=0.75, width=8, height=8, units='in')

write.csv(results$theta, file=paste0(args[1], "STdeconvolve_prop_norm.csv"))
write.csv(t(results$beta), file=paste0(args[1], "STdeconvolve_beta_norm.csv"))



##### For PCA and clustering only
#se_sc <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(data.dir = './SRX3436301/outs/filtered_feature_bc_matrix/'))
matrix_sc <- np$load(paste0(normDataDir, 'sc_adata_X.npz'))[['arr_0']]
sc_genes <- read.csv(paste0(normDataDir, 'sc_adata.var.csv'))
sc_obs <- read.csv(paste0(normDataDir, 'sc_adata.obs.csv'))
rownames(matrix_sc) <- get(colnames(sc_genes)[1], sc_genes)
colnames(matrix_sc) <- get(colnames(sc_obs)[1], sc_obs)
se_sc <- Seurat::CreateSeuratObject(counts = as(100*matrix_sc, "sparseMatrix"))

se_sc <- Seurat::FindVariableFeatures(se_sc, verbose = FALSE)
se_sc <- Seurat::ScaleData(se_sc, verbose = FALSE)
se_sc <- Seurat::RunPCA(se_sc, verbose = FALSE)
se_sc <- Seurat::RunUMAP(se_sc, dims = 1:30, verbose = FALSE)
se_sc <- Seurat::FindNeighbors(se_sc)
se_sc <- Seurat::FindClusters(se_sc, resolution=0.3)
Seurat::DimPlot(se_sc, group.by = "seurat_clusters", label = TRUE) + Seurat::NoLegend()
ggsave(paste0(args[1], "STdeconvolve_sc_clusters.png"), dpi=600, scale=0.5, width=8, height=8, units='in')
cluster_markers_all <- Seurat::FindAllMarkers(object = se_sc, assay = NULL, slot = "data", verbose = TRUE, test.use = "wilcox", only.pos = TRUE)


write.csv(se_sc@active.ident, file=paste0(args[1], "STdeconvolve_sc_cluster_ids.csv"))
write.csv(se_sc@reductions[["pca"]]@cell.embeddings, file=paste0(args[1], "STdeconvolve_sc_pca.csv"))
write.csv(se_sc@reductions[["pca"]]@feature.loadings, file=paste0(args[1], "STdeconvolve_sc_pca_feature_loadings.csv"))
write.csv(cluster_markers_all, file=paste0(args[1], "STdeconvolve_sc_cluster_markers.csv"))

quit(status=0)