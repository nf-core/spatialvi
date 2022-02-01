#!/usr/local/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
#args = c("c:/Projects/A_ST/forPipelineDev/")

# #### Spatial clustering:
# The latent cluster is modeled to depend on three parameters: cluster mean, fixed precision matrix and scaling factor. <br>
# The parameters are set to follow the priors. <br>
# Clusters are initialized using a non-spatial clustering method (mclust). <br>
# Iteratively and sequentially, each of the three parameters is updated via MCMC (Gibbs sampling). <br>
# The cluster membership is updated via MCMC (Metropolis–Hastings algorithm) by taking into account both the likelihood and spatial prior information. <br>

# #### Enhancement of clustering resolution:
# The procedure reassigns the total expression (in the PC space) within a spot to its constituent subspots by leveraging spatial information. The latent expression of each subspot is initialized to be expression of the original spot, and updated via MCMC (the Metropolis–Hastings algorithm). The cluster membership is then determined similarly to the spatial clustering above. <br>

# #### Mapping gene expression:
# A regression model is trained for each gene of interest, and used to predict gene expression from the high-resolution PCs estimated using enhanced-resolution clustering. <br>

# #### Tuning number of clusters:
# Elbow of negative pseudo-log-likelihood curve, that estimates how well the model fit the data for each number of clusters. <br>

set.seed(123)

library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)

library(reticulate)
np <- import("numpy")

normDataDir <- args[1]

# Load gene names
rowData <- read.csv(paste0(normDataDir, "st_adata.var.csv"))$X
row_df <- as.data.frame(rowData)
row.names(row_df) <- rowData

# Load coordinates of spots
st_obs_all <- read.csv(paste0(normDataDir, "st_adata.obs.csv"))
col_df <- as.data.frame(st_obs_all$X)
row.names(col_df) <- st_obs_all$X
col_df$imagerow <- st_obs_all$array_row
col_df$imagecol <- st_obs_all$array_col
colnames(col_df) <- c('id', 'imagerow', 'imagecol')
col_df$row <- st_obs_all$array_row
col_df$col <- st_obs_all$array_col


# Load counts data
count.data <- np$load(paste0(normDataDir, "st_adata_X.npz"))[['arr_0']] * 100.
colnames(count.data) <- st_obs_all$X
rownames(count.data) <- rowData

dsp <- SingleCellExperiment(assays=list(counts=count.data), rowData=row_df, colData=col_df)
dsp <- spatialPreprocess(dsp, platform="Visium", n.PCs=7, n.HVGs=2000, log.normalize=TRUE)
    
dsp <- qTune(dsp, qs=seq(2, 10), platform="Visium", d=7)
qPlot(dsp)
ggsave(paste0(normDataDir, "st_bayes_qtune.png"), dpi=600, scale=0.55, width=8, height=8, units="in")
    
dsp <- spatialCluster(dsp, q=5, platform="Visium", d=7, init.method="mclust", model="t", gamma=2, nrep=1000, burn.in=100, save.chain=FALSE)
    
clusterPlot(dsp, palette=c("purple", "cyan", "blue", "yellow", "red"), color=NA) + theme_bw() + xlab("Column") + ylab("Row") + labs(fill="BayesSpace\ncluster", title="Spatial clustering")
ggsave(paste0(normDataDir, "st_bayes_clusters.png"), dpi=600, scale=0.75, width=8, height=8, units="in")

dsp.enhanced <- spatialEnhance(dsp, q=5, platform="Visium", d=7, model="t", gamma=2, jitter_prior=0.3, jitter_scale=3.5, nrep=1000, burn.in=100)
    
clusterPlot(dsp.enhanced, palette=c("purple", "cyan", "blue", "yellow", "red"), color=NA) + theme_bw() + xlab("Column") + ylab("Row") + labs(fill="BayesSpace\ncluster", title="Spatial clustering")
ggsave(paste0(normDataDir, "st_bayes_clusters_enhanced.png"), dpi=600, scale=0.75, width=8, height=8, units="in") 



ncol <- 6
markers <- dsp@rowRanges@partitioning@NAMES[0:ncol]

dsp.enhanced <- enhanceFeatures(dsp.enhanced, dsp, feature_names=markers, nrounds=0)
enhanced.plots <- purrr::map(markers, function(x) featurePlot(dsp.enhanced, x))
spot.plots <- purrr::map(markers, function(x) featurePlot(dsp, x, color=NA))
patchwork::wrap_plots(c(spot.plots, enhanced.plots), ncol=ncol)
ggsave(paste0(normDataDir, "st_bayes_original_and_enhanced.png"), dpi=600, scale=1.25, width=2.3*ncol, height=4.5, limitsize=FALSE, units="in")



df_subspot_cluster_and_coord <- subset(as.data.frame(dsp.enhanced@colData), select=-c(imagecol, imagerow))
write.csv(df_subspot_cluster_and_coord, file=paste0(args[1], "bayes_subspot_cluster_and_coord.csv"))

df_enhanced_markers <- as.data.frame(dsp.enhanced@assays@data@listData[["logcounts"]][markers,])
write.csv(df_enhanced_markers, file=paste0(args[1], "bayes_enhanced_markers.csv"))

df_spot_cluster <- subset(as.data.frame(dsp@colData), select=-c(imagecol, imagerow))
write.csv(df_spot_cluster, file=paste0(args[1], "bayes_spot_cluster.csv"))

quit(status=0)