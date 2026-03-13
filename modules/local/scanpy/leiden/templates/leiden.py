#!/usr/bin/env python3
"""
Perform Leiden clustering on the neighbor graph.
"""

import platform
import importlib.metadata
import yaml
import scanpy as sc
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"
resolution = float("${resolution}")
key_added = "${key_added}"

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Performing Leiden clustering on: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Resolution: {resolution}")
print(f"Key added: {key_added}")

# Check if neighbors have been computed
if "neighbors" not in adata.uns:
    raise ValueError("Neighbor graph not found. Run scanpy.pp.neighbors first.")

# Perform Leiden clustering
sc.tl.leiden(
    adata,
    resolution=resolution,
    key_added=key_added
)

# Print summary
n_clusters = adata.obs[key_added].nunique()
cluster_sizes = adata.obs[key_added].value_counts().sort_index()

print(f"Found {n_clusters} clusters")
print("Cluster sizes:")
for cluster, size in cluster_sizes.items():
    print(f"  Cluster {cluster}: {size} spots ({size/adata.shape[0]*100:.1f}%)")

# Store clustering parameters in uns
adata.uns["leiden"] = {
    "resolution": resolution,
    "n_clusters": n_clusters,
}

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with clusters to: {output_adata}")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "anndata": importlib.metadata.version("anndata"),
        "leidenalg": importlib.metadata.version("leidenalg"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)