#!/usr/bin/env python3
"""
Compute spatial autocorrelation statistics to identify spatially variable genes.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import importlib.metadata
import yaml
import squidpy as sq
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"
output_csv = "${prefix}_svg.csv"
mode = "${mode}"
n_top_genes = int("${n_top_genes}")

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Computing spatial autocorrelation for: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Mode: {mode}")
print(f"Number of top genes to export: {n_top_genes}")

# Check if spatial neighbors exist
if "spatial_connectivities" not in adata.obsp:
    raise ValueError("Spatial connectivities not found. Run squidpy.gr.spatial_neighbors first.")

# Make var names unique
adata.var_names_make_unique()

# Compute spatial autocorrelation
sq.gr.spatial_autocorr(
    adata,
    mode=mode
)

# Get results key based on mode
if mode == "moran":
    results_key = "moranI"
elif mode == "geary":
    results_key = "gearyC"
else:
    raise ValueError(f"Unknown mode: {mode}. Use 'moran' or 'geary'.")

# Get results dataframe
svg_df = adata.uns[results_key]

# Print top SVGs
print(f"\\nTop {min(n_top_genes, len(svg_df))} spatially variable genes:")
print(svg_df.head(n_top_genes).to_string())

# Export to CSV
svg_df.to_csv(output_csv)
print(f"\\nExported SVG results to: {output_csv}")

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with spatial autocorrelation to: {output_adata}")

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
