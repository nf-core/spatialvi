#!/usr/local/bin/Rscript

# Load packages
library(argparse)
library(reticulate)
library(SpatialExperiment)
library(scran)


# Parse command-line arguments
parser <- ArgumentParser()

args <- parser$add_argument_group("Agruments", "required and optional arguments")

args$add_argument("--filePath", help = "Path to npz counts file", metavar="dir", required=TRUE)
args$add_argument("--npCountsOutputName", help = "Name of npz counts file", metavar="file", required=TRUE)
args$add_argument("--npFactorsOutputName", help = "Name of npz factors file", metavar="file", required=TRUE)

args <- parser$parse_args()


# Main script
np <- import("numpy")
matrix_st <- np$load(paste0(args$filePath, args$npCountsOutputName))[['arr_0']]
print(dim(matrix_st))

spe <- SpatialExperiment(list(counts=matrix_st))
sfs <- calculateSumFactors(spe, cluster=quickCluster(spe))
print(length(sfs))

np$savez_compressed(paste0(args$filePath, args$npFactorsOutputName), sfs)

quit(status=0)
