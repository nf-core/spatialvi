#!/usr/bin/env Rscript

# Load packages
library(argparse)
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)
library(reticulate)

# Parse command-line arguments
parser <- ArgumentParser()
args <- parser$add_argument_group("Arguments",
                                  "required and optional arguments")
args$add_argument("--nameX",
                  default  = "st_adata_X.npz",
                  help     = "Path to X",
                  metavar  = "file",
                  required = FALSE)
args$add_argument("--nameVar",
                  default  = "st_adata.var.csv",
                  help     = "Path to features metadata",
                  metavar  = "file",
                  required = FALSE)
args$add_argument("--nameObs",
                  default  = "st_adata.obs.csv",
                  help     = "Path to observation metadata",
                  metavar  = "file",
                  required = FALSE)
args$add_argument("--countsFactor",
                  default  = 100,
                  help     = "factor",
                  metavar  = "factor",
                  required = FALSE)
args$add_argument("--numberHVG",
                  default  = 2000,
                  type     = "integer",
                  help     = "factor",
                  metavar  = "factor",
                  required = FALSE)
args$add_argument("--numberPCs",
                  default  = 7,
                  type     = "integer",
                  help     = "factor",
                  metavar  = "factor",
                  required = FALSE)
args$add_argument("--minClusters",
                  default  = 2,
                  type     = "integer",
                  help     = "factor",
                  metavar  = "factor",
                  required = FALSE)
args$add_argument("--maxClusters",
                  default  = 10,
                  type     = "integer",
                  help     = "factor",
                  metavar  = "factor",
                  required = FALSE)
args$add_argument("--optimalQ",
                  default  = 5,
                  type     = "integer",
                  help     = "factor",
                  metavar  = "factor",
                  required = FALSE)
args$add_argument("--STplatform",
                  default  = "Visium",
                  help     = "Technology grid",
                  metavar  = "factor",
                  required = FALSE)
args$add_argument("--qtuneSaveName",
                  default  = "st_bayes_qtune.png",
                  help     = "file name",
                  metavar  = "file",
                  required = FALSE)
args$add_argument("--bayesClustersName",
                  default  = "st_bayes_clusters.png",
                  help     = "file name",
                  metavar  = "file",
                  required = FALSE)
args$add_argument("--bayesClustersEnhancedName",
                  default  = "st_bayes_clusters_enhanced.png",
                  help     = "file name",
                  metavar  = "file",
                  required = FALSE)
args$add_argument("--bayesOriginalEnhancedFeatures",
                  default  = "st_bayes_original_and_enhanced.png",
                  help     = "file name",
                  metavar  = "file",
                  required = FALSE)
args$add_argument("--bayesEnhancedMarkers",
                  default  = "bayes_enhanced_markers.csv",
                  help     = "file name",
                  metavar  = "file",
                  required = FALSE)
args$add_argument("--bayesSubspotCoord",
                  default  = "bayes_subspot_cluster_and_coord.csv",
                  help     = "file name",
                  metavar  = "file",
                  required = FALSE)
args$add_argument("--bayesSpotCluster",
                  default  = "bayes_spot_cluster.csv",
                  help     = "file name",
                  metavar  = "file",
                  required = FALSE)
args <- parser$parse_args()

# #### Spatial clustering:
# The latent cluster is modeled to depend on three parameters: cluster mean,
# fixed precision matrix and scaling factor.  The parameters are set to follow
# the priors. Clusters are initialized using a non-spatial clustering method
# (mclust). Iteratively and sequentially, each of the three parameters is
# updated via MCMC (Gibbs sampling). The cluster membership is updated via MCMC
# (Metropolis–Hastings algorithm) by taking into account both the likelihood
# and spatial prior information.
#
# #### Enhancement of clustering resolution:
# The procedure reassigns the total expression (in the PC space) within a spot
# to its constituent subspots by leveraging spatial information. The latent
# expression of each subspot is initialized to be expression of the original
# spot, and updated via MCMC (the Metropolis–Hastings algorithm). The cluster
# membership is then determined similarly to the spatial clustering above.
#
# #### Mapping gene expression:
# A regression model is trained for each gene of interest, and used to predict
# gene expression from the high-resolution PCs estimated using
# enhanced-resolution clustering.
#
# #### Tuning number of clusters:
# Elbow of negative pseudo-log-likelihood curve, that estimates how well the
# model fit the data for each number of clusters.

# Main script
set.seed(123)
np <- import("numpy")

# Load gene names
rowData <- read.csv(paste0(args$nameVar))$X
row_df <- as.data.frame(rowData)
row.names(row_df) <- rowData

# Load coordinates of spots
st_obs_all <- read.csv(paste0(args$nameObs))
col_df <- as.data.frame(st_obs_all$X)
row.names(col_df) <- st_obs_all$X
col_df$imagerow <- st_obs_all$array_row
col_df$imagecol <- st_obs_all$array_col
colnames(col_df) <- c('id', 'imagerow', 'imagecol')
col_df$row <- st_obs_all$array_row
col_df$col <- st_obs_all$array_col

# Load counts data
count.data <- np$load(paste0(args$nameX))[['arr_0']] * args$countsFactor
colnames(count.data) <- st_obs_all$X
rownames(count.data) <- rowData

print(args$countsFactor)
print(args$STplatform)
print(args$numberPCs)
print(args$numberHVG)
print(args$minClusters)

dsp <- SingleCellExperiment(assays  = list(counts=count.data),
                            rowData = row_df,
                            colData = col_df)
dsp <- spatialPreprocess(dsp,
                         platform      = args$STplatform,
                         n.PCs         = args$numberPCs,
                         n.HVGs        = args$numberHVG,
                         log.normalize = TRUE)

dsp <- qTune(dsp,
             qs       = seq(args$minClusters, args$maxClusters),
             d        = args$numberPCs,
             platform = "Visium")
qPlot(dsp)
ggsave(paste0(args$qtuneSaveName),
       dpi    = 600,
       scale  = 0.55,
       width  = 8,
       height = 8,
       units  = "in")

# TODO: Need to determine optimal q!
dsp <- spatialCluster(dsp,
                      q           = args$optimalQ,
                      platform    = args$STplatform,
                      d           = args$numberPCs,
                      init.method = "mclust",
                      model       = "t",
                      gamma       = 2,
                      nrep        = 1000,
                      burn.in     = 100,
                      save.chain  = FALSE)

clusterPlot(dsp,
            palette = c("purple", "cyan", "blue", "yellow", "red"),
            color   = NA) +
    theme_bw() +
    xlab("Column") +
    ylab("Row") +
    labs(fill = "BayesSpace\ncluster", title = "Spatial clustering")
ggsave(paste0(args$bayesClustersName),
       dpi    = 600,
       scale  = 0.75,
       width  = 8,
       height = 8,
       units  = "in")

dsp.enhanced <- spatialEnhance(dsp,
                               q            = args$optimalQ,
                               platform     = "Visium",
                               d            = args$numberPCs,
                               model        = "t",
                               gamma        = 2,
                               jitter_prior = 0.3,
                               jitter_scale = 3.5,
                               nrep         = 1000,
                               burn.in      = 100)

clusterPlot(dsp.enhanced,
            palette = c("purple", "cyan", "blue", "yellow", "red"),
            color   = NA) +
    theme_bw() +
    xlab("Column") +
    ylab("Row") +
    labs(fill = "BayesSpace\ncluster", title = "Spatial clustering")
ggsave(paste0(args$bayesClustersEnhancedName),
       dpi    = 600,
       scale  = 0.75,
       width  = 8,
       height = 8,
       units  = "in")

# As of current implementation: take first 6 genes and enhance/plot them
ncol <- 6
markers <- dsp@rowRanges@partitioning@NAMES[0:ncol]

dsp.enhanced <- enhanceFeatures(dsp.enhanced,
                                dsp,
                                feature_names = markers,
                                nrounds       = 0)
enhanced.plots <- purrr::map(markers, function(x) featurePlot(dsp.enhanced, x))
spot.plots <- purrr::map(markers, function(x) featurePlot(dsp, x, color = NA))
patchwork::wrap_plots(c(spot.plots, enhanced.plots), ncol = ncol)
ggsave(paste0(args$bayesOriginalEnhancedFeatures),
       dpi       = 600
       scale     = 1.25
       width     = 2.3*ncol
       height    = 4.5
       limitsize = FALSE
       units     = "in")

df_subspot_cluster_and_coord <- subset(as.data.frame(dsp.enhanced@colData),
                                       select = -c(imagecol, imagerow))
write.csv(df_subspot_cluster_and_coord,
          file = paste0(args$bayesSubspotCoord))

df_enhanced_markers <- as.data.frame(dsp.enhanced@assays@data@listData[["logcounts"]][markers,])
write.csv(df_enhanced_markers,
          file = paste0(args$bayesEnhancedMarkers))

df_spot_cluster <- subset(as.data.frame(dsp@colData),
                          select = -c(imagecol, imagerow))
write.csv(df_spot_cluster,
          file = paste0(args$bayesSpotCluster))
