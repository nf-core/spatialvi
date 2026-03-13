#!/usr/bin/env python3
"""
Compute neighborhood enrichment by permutation test.
"""

import platform
import importlib.metadata
import yaml
import squidpy as sq
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"
cluster_key = "${cluster_key}"

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Computing neighborhood enrichment for: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Cluster key: {cluster_key}")

# Check if required data exists
if cluster_key not in adata.obs.columns:
    raise ValueError(f"Column '{cluster_key}' not found in adata.obs")
if "spatial_connectivities" not in adata.obsp:
    raise ValueError("Spatial connectivities not found. Run squidpy.gr.spatial_neighbors first.")

# Compute neighborhood enrichment
sq.gr.nhood_enrichment(
    adata,
    cluster_key=cluster_key
)

# Print summary
n_clusters = adata.obs[cluster_key].nunique()
print(f"Computed neighborhood enrichment for {n_clusters} clusters")
print(f"Results stored in adata.uns['{cluster_key}_nhood_enrichment']")

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with neighborhood enrichment to: {output_adata}")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "squidpy": importlib.metadata.version("squidpy"),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)