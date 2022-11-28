#!/usr/local/bin/Rscript

# Load packages
library(argparse)
library(dplyr)
library(ggplot2)
library(Seurat)
library(STdeconvolve)
library(corrplot)
library(reticulate)

# Parse command-line arguments
parser <- ArgumentParser()
args <- parser$add_argument_group("Agruments", "required and optional arguments")
args$add_argument("--outsPath", help = "Path to data", metavar = "dir", required = TRUE)
args$add_argument("--nameX", default = "st_adata_X.npz", help = "Path to X", metavar = "file", required = FALSE)
args$add_argument("--nameVar", default = "st_adata.var.csv", help = "Path to features metadata", metavar = "file", required = FALSE)
args$add_argument("--nameObs", default = "st_adata.obs.csv", help = "Path to observation metadata", metavar = "file", required = FALSE)
args$add_argument("--SCnameX", default = "sc_adata_X.npz", help = "Path to X", metavar = "file", required = FALSE)
args$add_argument("--SCnameVar", default = "sc_adata.var.csv", help = "Path to features metadata", metavar = "file", required = FALSE)
args$add_argument("--SCnameObs", default = "sc_adata.obs.csv", help = "Path to observation metadata", metavar = "file", required = FALSE)
args$add_argument("--fileh5", default = "raw_feature_bc_matrix.h5", help = "File HDF5", metavar = "file", required = FALSE)
args$add_argument("--outsSubDir", default = "raw_feature_bc_matrix/", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--mtxGeneColumn", default = 2, type = "integer", help = "columns index", metavar = "col", required = FALSE)
args$add_argument("--countsFactor", default = 100, type = "integer", help = "factor", metavar = "factor", required = FALSE)
args$add_argument("--corpusRemoveAbove", default = 1.0, type = "double", help = "factor", metavar = "factor", required = FALSE)
args$add_argument("--corpusRemoveBelow", default = 0.05, type = "double", help = "factor", metavar = "factor", required = FALSE)
args$add_argument("--LDAminTopics", default = 8, type = "integer", help = "factor", metavar = "factor", required = FALSE)
args$add_argument("--LDAmaxTopics", default = 9, type = "integer", help = "factor", metavar = "factor", required = FALSE)
args$add_argument("--LDAsaveFile", default = "STdeconvolve_optLDA.rds", help = "File to save LDA RDS", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveScatterpiesName", default = "STdeconvolve_st_scatterpies.png", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveScatterpiesSize", default = 2.85, type = "double", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveFeaturesSizeFactor", default = 1.0, type = "double", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveFeaturesName", default = "STdeconvolve_st_prop.png", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveCorrName", default = "STdeconvolve_st_prop_corr.png", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolvePropNormName", default = "STdeconvolve_prop_norm.csv", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveBetaNormName", default = "STdeconvolve_beta_norm.csv", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveSCclustersName", default = "STdeconvolve_sc_clusters.png", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveSCclusterIds", default = "STdeconvolve_sc_cluster_ids.csv", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveSCpca", default = "STdeconvolve_sc_pca.csv", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveSCloadings", default = "STdeconvolve_sc_pca_feature_loadings.csv", help = "dir", metavar = "file", required = FALSE)
args$add_argument("--STdeconvolveSCclusterMarkers", default = "STdeconvolve_sc_cluster_markers.csv", help = "dir", metavar = "file", required = FALSE)
args <- parser$parse_args()

# Main script
set.seed(123)
np <- import("numpy")

filename <- list.files(path=args$outsPath, pattern = args$fileh5)[1]
print(args$outsPath)

if (!is.na(filename)) {
    print(filename)
    se_st <- Seurat::Load10X_Spatial(data.dir = args$outsPath, filename = filename)
} else {
    image <- Read10X_Image(image.dir = file.path(args$outsPath, 'spatial'), filter.matrix = TRUE)
    m <- Read10X(paste0(args$outsPath, args$outsSubDir), gene.column = args$mtxGeneColumn)
    m <- m[, row.names(image@coordinates)]
    m <- m[, colSums(m) > 0]
    se_st <- CreateSeuratObject(counts = m, assay = "Spatial")
    image <- image[Cells(x = se_st)]
    DefaultAssay(object = image) <- "Spatial"
    se_st[["slice1"]] <- image
}

matrix_st <- np$load(paste0(args$nameX))[['arr_0']]
st_genes <- read.csv(paste0(args$nameVar))$X
st_obs <- read.csv(paste0(args$nameObs))$X
rownames(matrix_st) <- st_genes
colnames(matrix_st) <- st_obs
se_st@assays$Spatial@counts <- as(matrix_st, "sparseMatrix")
se_st@assays$Spatial@data <- as(matrix_st, "sparseMatrix")
se_st@assays$Spatial@counts <- round((args$countsFactor) * se_st@assays$Spatial@counts)

corpus <- restrictCorpus(se_st@assays$Spatial@counts,
    removeAbove = args$corpusRemoveAbove,
    removeBelow = args$corpusRemoveBelow
)
corpus <- corpus + 1
print(dim(as.matrix(corpus)))
print(sum(colSums(as.matrix(corpus)) == 0))

ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(args$LDAminTopics, args$LDAmaxTopics, by = 1))
optLDA <- optimalModel(models = ldas, opt = "min")
saveRDS(object = optLDA, file = args$LDAsaveFile)

optLDA <- readRDS(file = args$LDAsaveFile)

results <- getBetaTheta(optLDA, t(as.matrix(corpus)))
deconProp <- results$theta
posk <- se_st@images[["slice1"]]@coordinates[c("imagerow", "imagecol")]
names(posk) <- c('y', 'x')
posk$x <- posk$x * se_st@images[["slice1"]]@scale.factors[["lowres"]]
posk$y <- dim(se_st@images[["slice1"]])[1] -
    posk$y * se_st@images[["slice1"]]@scale.factors[["lowres"]]
vizAllTopics(deconProp,
             posk,
             lwd     = 0,
             overlay = se_st@images[["slice1"]]@image,
             r       = args$STdeconvolveScatterpiesSize)
ggsave(args$STdeconvolveScatterpiesName,
    dpi    = 600,
    scale  = 1.0,
    width  = 8,
    height = 8,
    units  = 'in'
)

# Individual topics proportions spatial overlay
decon_df <- deconProp %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")

print(colnames(decon_df))
print(dim(decon_df))

se_st@meta.data <- se_st@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
Seurat::SpatialFeaturePlot(object         = se_st,
    features       = colnames(decon_df)[-1],
    alpha          = c(0.1, 1),
    min.cutoff     = 0,
    max.cutoff     = 0.3,
    crop           = FALSE,
    pt.size.factor = args$STdeconvolveFeaturesSizeFactor
)
ggsave(args$STdeconvolveFeaturesName,
    dpi    = 300,
    scale  = 2.0,
    width  = 8,
    height = 8,
    units  = 'in'
)

# Cell type proportions correlation
decon_mtrx_sub <- deconProp[, colnames(deconProp)[which(colnames(deconProp) != "res_ss")]]
decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]
colnames(decon_mtrx_sub) <- paste("X", colnames(decon_mtrx_sub), sep = "")
decon_cor <- cor(decon_mtrx_sub)
p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)
ggcorrplot::ggcorrplot(
        corr         = decon_cor,
        p.mat        = p.mat[[1]],
        hc.order     = TRUE,
        type         = "full",
        insig        = "blank",
        lab          = TRUE,
        outline.col  = "lightgrey",
        method       = "square",
        colors       = c("#6D9EC1", "white", "#E46726"),
        title        = "Cell type proportions correlation",
        legend.title = "Correlation\n(Pearson)") +
    ggplot2::theme(plot.title = ggplot2::element_text(size  = 22, hjust = 0.5, face  = "bold")) +
    ggplot2::theme(legend.text  = ggplot2::element_text(size  = 12)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size  = 15)) +
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90)) +
    ggplot2::theme(axis.text    = ggplot2::element_text(size  = 18, vjust = 0.5))
ggsave(args$STdeconvolveCorrName,
    dpi    = 600,
    scale  = 0.75,
    width  = 8,
    height = 8,
    units  = 'in'
)

write.csv(results$theta, file = args$STdeconvolvePropNormName)
write.csv(t(results$beta), file = args$STdeconvolveBetaNormName)

##### For PCA and clustering only
matrix_sc <- np$load(paste0(args$SCnameX))[['arr_0']]
sc_genes <- read.csv(paste0(args$SCnameVar))
sc_obs <- read.csv(paste0(args$SCnameObs))
rownames(matrix_sc) <- get(colnames(sc_genes)[1], sc_genes)
colnames(matrix_sc) <- get(colnames(sc_obs)[1], sc_obs)
se_sc <- Seurat::CreateSeuratObject(counts = as((args$countsFactor) * matrix_sc, "sparseMatrix"))

se_sc <- Seurat::FindVariableFeatures(se_sc, verbose = FALSE)
se_sc <- Seurat::ScaleData(se_sc, verbose = FALSE)
se_sc <- Seurat::RunPCA(se_sc, verbose = FALSE)
se_sc <- Seurat::RunUMAP(se_sc, dims = 1:30, verbose = FALSE)
se_sc <- Seurat::FindNeighbors(se_sc)
se_sc <- Seurat::FindClusters(se_sc, resolution = 0.3)
Seurat::DimPlot(se_sc, group.by = "seurat_clusters", label = TRUE) +
    Seurat::NoLegend()
ggsave(args$STdeconvolveSCclustersName,
    dpi    = 600,
    scale  = 0.5,
    width  = 8,
    height = 8,
    units  = 'in'
)
cluster_markers_all <- Seurat::FindAllMarkers(
    object   = se_sc,
    assay    = NULL,
    verbose  = TRUE,
    only.pos = TRUE,
    slot     = "data",
    test.use = "wilcox"
)

write.csv(se_sc@active.ident, file = args$STdeconvolveSCclusterIds)
write.csv(se_sc@reductions[["pca"]]@cell.embeddings, file = args$STdeconvolveSCpca)
write.csv(se_sc@reductions[["pca"]]@feature.loadings, file = args$STdeconvolveSCloadings)
write.csv(cluster_markers_all, file = args$STdeconvolveSCclusterMarkers)
