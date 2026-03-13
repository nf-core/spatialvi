#!/usr/bin/env python3
"""
Rank genes for characterizing groups (differential expression analysis).
"""

import platform
import importlib.metadata
import yaml
import scanpy as sc
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"
groupby = "${groupby}"
method = "${method}"

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Performing differential expression analysis on: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Groupby: {groupby}")
print(f"Method: {method}")

# Check if groupby column exists
if groupby not in adata.obs.columns:
    raise ValueError(f"Column '{groupby}' not found in adata.obs")

# Suppress verbose output
sc.settings.verbosity = 0

# Rank genes by groups
sc.tl.rank_genes_groups(
    adata,
    groupby=groupby,
    method=method
)

# Print summary
n_groups = adata.obs[groupby].nunique()
print(f"Computed DEGs for {n_groups} groups")

# Get top genes per group
top_genes = {}
for group in adata.obs[groupby].unique():
    genes = sc.get.rank_genes_groups_df(adata, group=str(group))
    top_genes[group] = genes.head(5)["names"].tolist()
    print(f"  Group {group} top 5: {', '.join(top_genes[group])}")

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with DEGs to: {output_adata}")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)