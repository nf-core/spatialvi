#!/usr/bin/env python3
"""
Filter spots and genes from an AnnData object based on QC metrics.

Filtering steps (in order):
1. Remove spots outside tissue
2. Filter spots by minimum counts
3. Filter spots by minimum genes
4. Filter genes by minimum spots
5. Filter spots by mitochondrial content
6. Filter spots by ribosomal content
7. Filter spots by haemoglobin content
"""

import platform
import importlib.metadata
import json
import yaml
import scanpy as sc
import anndata as ad

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}.h5ad"
output_stats = "${prefix}_stats.json"

min_counts = int("${min_counts}")
min_genes = int("${min_genes}")
min_spots = int("${min_spots}")
mito_threshold = float("${mito_threshold}")
ribo_threshold = float("${ribo_threshold}")
hb_threshold = float("${hb_threshold}")

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Filtering AnnData from: {input_adata}")
print(f"Initial shape: {adata.shape}")

# Store initial counts
n_total_spots = adata.shape[0]
n_total_genes = adata.shape[1]

# Save a copy as restore-point if filtering results in 0 spots
adata_backup = adata.copy()

# Initialize statistics dictionary
stats = {
    "sample_id": "${meta.id}",
    "total_spots_before": n_total_spots,
    "total_genes_before": n_total_genes,
    "min_counts_threshold": min_counts,
    "min_genes_threshold": min_genes,
    "min_spots_threshold": min_spots,
    "mito_threshold": mito_threshold,
    "ribo_threshold": ribo_threshold,
    "hb_threshold": hb_threshold,
}

# Make var names unique before filtering
adata.var_names_make_unique()

#
# 1. Filter spots outside tissue
#
if "in_tissue" in adata.obs.columns:
    n_before = adata.shape[0]
    adata = adata[adata.obs["in_tissue"] == 1].copy()
    n_spots_outside_tissue = n_before - adata.shape[0]
    print(f"Removed {n_spots_outside_tissue} spots outside tissue")
else:
    n_spots_outside_tissue = 0
    print("Column 'in_tissue' not found, skipping tissue filtering")

stats["spots_filtered_outside_tissue"] = n_spots_outside_tissue

#
# 2. Filter spots by minimum counts
#
n_before = adata.shape[0]
sc.pp.filter_cells(adata, min_counts=min_counts)
n_spots_filtered_min_counts = n_before - adata.shape[0]
print(f"Removed {n_spots_filtered_min_counts} spots with < {min_counts} counts")

stats["spots_filtered_min_counts"] = n_spots_filtered_min_counts

#
# 3. Filter spots by minimum genes
#
n_before = adata.shape[0]
sc.pp.filter_cells(adata, min_genes=min_genes)
n_spots_filtered_min_genes = n_before - adata.shape[0]
print(f"Removed {n_spots_filtered_min_genes} spots with < {min_genes} genes")

stats["spots_filtered_min_genes"] = n_spots_filtered_min_genes

#
# 4. Filter genes by minimum spots
#
n_before = adata.shape[1]
sc.pp.filter_genes(adata, min_cells=min_spots)
n_genes_filtered_min_spots = n_before - adata.shape[1]
print(f"Removed {n_genes_filtered_min_spots} genes expressed in < {min_spots} spots")

stats["genes_filtered_min_spots"] = n_genes_filtered_min_spots

#
# 5. Filter spots by mitochondrial content
#
if "pct_counts_mt" in adata.obs.columns:
    n_before = adata.shape[0]
    adata = adata[adata.obs["pct_counts_mt"] <= mito_threshold].copy()
    n_spots_filtered_mito = n_before - adata.shape[0]
    print(f"Removed {n_spots_filtered_mito} spots with > {mito_threshold}% mito content")
else:
    n_spots_filtered_mito = 0
    print("Column 'pct_counts_mt' not found, skipping mito filtering")

stats["spots_filtered_mito"] = n_spots_filtered_mito

#
# 6. Filter spots by ribosomal content
#
if "pct_counts_ribo" in adata.obs.columns:
    n_before = adata.shape[0]
    adata = adata[adata.obs["pct_counts_ribo"] >= ribo_threshold].copy()
    n_spots_filtered_ribo = n_before - adata.shape[0]
    print(f"Removed {n_spots_filtered_ribo} spots with < {ribo_threshold}% ribo content")
else:
    n_spots_filtered_ribo = 0
    print("Column 'pct_counts_ribo' not found, skipping ribo filtering")

stats["spots_filtered_ribo"] = n_spots_filtered_ribo

#
# 7. Filter spots by haemoglobin content
#
if "pct_counts_hb" in adata.obs.columns:
    n_before = adata.shape[0]
    adata = adata[adata.obs["pct_counts_hb"] <= hb_threshold].copy()
    n_spots_filtered_hb = n_before - adata.shape[0]
    print(f"Removed {n_spots_filtered_hb} spots with > {hb_threshold}% hb content")
else:
    n_spots_filtered_hb = 0
    print("Column 'pct_counts_hb' not found, skipping hb filtering")

stats["spots_filtered_hb"] = n_spots_filtered_hb

#
# Check if filtering resulted in 0 spots or genes
#
filtering_failed = False
if adata.shape[0] == 0 or adata.shape[1] == 0:
    print("WARNING: Filtering resulted in 0 spots or genes remaining!")
    print("Restoring original unfiltered data.")
    adata = adata_backup
    filtering_failed = True

stats["filtering_failed"] = filtering_failed

#
# Calculate final statistics
#
n_remaining_spots = adata.shape[0]
n_remaining_genes = adata.shape[1]
n_spots_filtered = n_total_spots - n_remaining_spots
n_genes_filtered = n_total_genes - n_remaining_genes

stats["total_spots_after"] = n_remaining_spots
stats["total_genes_after"] = n_remaining_genes
stats["total_spots_filtered"] = n_spots_filtered
stats["total_genes_filtered"] = n_genes_filtered

# Store filtering stats in adata.uns for downstream use
adata.uns["filtering_stats"] = stats

print("\\nFiltering summary:")
print(f"  Spots: {n_total_spots} -> {n_remaining_spots} ({n_spots_filtered} removed)")
print(f"  Genes: {n_total_genes} -> {n_remaining_genes} ({n_genes_filtered} removed)")

# Write outputs
adata.write_h5ad(output_adata)
print(f"\\nWritten filtered AnnData to: {output_adata}")

with open(output_stats, "w") as f:
    json.dump(stats, f, indent=2)
print(f"Written filtering statistics to: {output_stats}")

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