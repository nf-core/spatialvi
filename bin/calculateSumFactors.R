#!/usr/bin/env Rscript

# Load packages
library(argparse)
library(reticulate)
library(SpatialExperiment)
library(scran)

# Parse command-line arguments
parser <- ArgumentParser()

args <- parser$add_argument_group("Agruments", "required and optional arguments")
args$add_argument("--npCountsOutputName", required = TRUE, metavar = "file", help = "Name of npz counts file")
args$add_argument("--npFactorsOutputName", required = TRUE, metavar = "file", help = "Name of npz factors file")
args <- parser$parse_args()

# Main script
np <- import("numpy")
matrix_st <- np$load(args$npCountsOutputName)[['arr_0']]
print(dim(matrix_st))

spe <- SpatialExperiment(list(counts=matrix_st))
sfs <- calculateSumFactors(spe, cluster=quickCluster(spe))
print(length(sfs))

# Save to file
np$savez_compressed(args$npFactorsOutputName, sfs)
