#!/usr/local/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

library(reticulate)
np <- import("numpy")
matrix_st <- np$load(paste0(args[1], args[2], '.npz'))[['arr_0']]
print(dim(matrix_st))

library(SpatialExperiment)
library(scran)
spe <- SpatialExperiment(list(counts=matrix_st))
sfs <- calculateSumFactors(spe, cluster=quickCluster(spe))
print(length(sfs))

np$savez_compressed(paste0(args[1], args[2], '_factors.npz'), sfs)
