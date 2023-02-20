#!/usr/bin/env python

# Load packages
import argparse
import scanpy as sc
import pandas as pd
import SpatialDE


# Parse command-line arguments
parser = argparse.ArgumentParser(description="Preprocess single cell transcriptomics data.")
parser.add_argument("--fileName", metavar="file", type=str, default="st_adata_norm.h5ad", help="File name.")
parser.add_argument("--saveFileName", metavar="file", type=str, default="stSpatialDE.csv", help="File name.")
parser.add_argument("--savePlotName", metavar="plot", type=str, default="stSpatialDE.png", help="File name.")
parser.add_argument("--plotTopHVG", metavar="number", type=int, default=15, help="File name.")
parser.add_argument("--numberOfColumns", metavar="number", type=int, default=5, help="File name.")
args = parser.parse_args()

# Main script
# See more settings at:
# https://scanpy.readthedocs.io/en/stable/generated/scanpy._settings.ScanpyConfig.html
st_adata = sc.read(args.fileName)
print(st_adata.shape)

counts = pd.DataFrame(st_adata.X.todense(), columns=st_adata.var_names, index=st_adata.obs_names)
coord = pd.DataFrame(st_adata.obsm["spatial"], columns=["x_coord", "y_coord"], index=st_adata.obs_names)

df_results = SpatialDE.run(coord, counts)

df_results.index = df_results["g"]
df_results = df_results.sort_values("qval", ascending=True)

df_results.to_csv(args.saveFileName)

# Plotting top most-HVG
keys = df_results.index.values[: args.plotTopHVG]
sc.pl.spatial(st_adata, img_key="hires", color=keys, alpha=0.7, save=args.savePlotName, ncols=args.numberOfColumns)
