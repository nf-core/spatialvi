#!/usr/bin/env python3
"""
Compute UMAP (Uniform Manifold Approximation and Projection) embedding.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import importlib.metadata
import yaml
import scanpy as sc
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"
min_dist = float("${min_dist}")
spread = float("${spread}")
key_added = "${key_added}"

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Computing UMAP for: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Parameters: min_dist={min_dist}, spread={spread}")

# Check if neighbors have been computed
if "neighbors" not in adata.uns:
    raise ValueError("Neighbor graph not found. Run scanpy.pp.neighbors first.")

# Compute UMAP
sc.tl.umap(
    adata,
    min_dist=min_dist,
    spread=spread,
    key_added=key_added
)

# Print summary
print(f"UMAP embedding shape: {adata.obsm['X_umap'].shape}")
print("UMAP coordinate ranges:")
print(f"  UMAP1: [{adata.obsm['X_umap'][:, 0].min():.2f}, {adata.obsm['X_umap'][:, 0].max():.2f}]")
print(f"  UMAP2: [{adata.obsm['X_umap'][:, 1].min():.2f}, {adata.obsm['X_umap'][:, 1].max():.2f}]")

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with UMAP to: {output_adata}")

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
