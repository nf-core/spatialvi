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

import importlib.metadata
import json
import platform
import yaml

import anndata as ad
import scanpy as sc


def filter_outside_tissue(adata, stats):
    """Filter spots outside tissue based on 'in_tissue' column."""
    if "in_tissue" in adata.obs.columns:
        n_before = adata.shape[0]
        adata = adata[adata.obs["in_tissue"] == 1].copy()
        n_filtered = n_before - adata.shape[0]
        print(f"Removed {n_filtered} spots outside tissue")
    else:
        n_filtered = 0
        print("Column 'in_tissue' not found, skipping tissue filtering")

    stats["spots_filtered_outside_tissue"] = n_filtered
    return adata, stats


def filter_min_counts(adata, min_counts, stats):
    """Filter spots with fewer than min_counts total counts."""
    n_before = adata.shape[0]
    sc.pp.filter_cells(adata, min_counts=min_counts)
    n_filtered = n_before - adata.shape[0]
    print(f"Removed {n_filtered} spots with < {min_counts} counts")

    stats["spots_filtered_min_counts"] = n_filtered
    return adata, stats


def filter_min_genes(adata, min_genes, stats):
    """Filter spots with fewer than min_genes expressed genes."""
    n_before = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=min_genes)
    n_filtered = n_before - adata.shape[0]
    print(f"Removed {n_filtered} spots with < {min_genes} genes")

    stats["spots_filtered_min_genes"] = n_filtered
    return adata, stats


def filter_genes_min_spots(adata, min_spots, stats):
    """Filter genes expressed in fewer than min_spots spots."""
    n_before = adata.shape[1]
    sc.pp.filter_genes(adata, min_cells=min_spots)
    n_filtered = n_before - adata.shape[1]
    print(f"Removed {n_filtered} genes expressed in < {min_spots} spots")

    stats["genes_filtered_min_spots"] = n_filtered
    return adata, stats


def filter_mito(adata, mito_threshold, stats):
    """Filter spots with mitochondrial content above threshold."""
    if "pct_counts_mt" in adata.obs.columns:
        n_before = adata.shape[0]
        adata = adata[adata.obs["pct_counts_mt"] <= mito_threshold].copy()
        n_filtered = n_before - adata.shape[0]
        print(f"Removed {n_filtered} spots with > {mito_threshold}% mito content")
    else:
        n_filtered = 0
        print("Column 'pct_counts_mt' not found, skipping mito filtering")

    stats["spots_filtered_mito"] = n_filtered
    return adata, stats


def filter_ribo(adata, ribo_threshold, stats):
    """Filter spots with ribosomal content below threshold."""
    if "pct_counts_ribo" in adata.obs.columns:
        n_before = adata.shape[0]
        adata = adata[adata.obs["pct_counts_ribo"] >= ribo_threshold].copy()
        n_filtered = n_before - adata.shape[0]
        print(f"Removed {n_filtered} spots with < {ribo_threshold}% ribo content")
    else:
        n_filtered = 0
        print("Column 'pct_counts_ribo' not found, skipping ribo filtering")

    stats["spots_filtered_ribo"] = n_filtered
    return adata, stats


def filter_hb(adata, hb_threshold, stats):
    """Filter spots with haemoglobin content above threshold."""
    if "pct_counts_hb" in adata.obs.columns:
        n_before = adata.shape[0]
        adata = adata[adata.obs["pct_counts_hb"] <= hb_threshold].copy()
        n_filtered = n_before - adata.shape[0]
        print(f"Removed {n_filtered} spots with > {hb_threshold}% hb content")
    else:
        n_filtered = 0
        print("Column 'pct_counts_hb' not found, skipping hb filtering")

    stats["spots_filtered_hb"] = n_filtered
    return adata, stats


def filter_adata(adata, min_counts, min_genes, min_spots, mito_threshold, ribo_threshold, hb_threshold):
    """Apply all filtering steps to AnnData object."""
    print(f"Initial shape: {adata.shape}")

    # Store initial counts
    n_total_spots = adata.shape[0]
    n_total_genes = adata.shape[1]

    # Save a copy as restore-point if filtering results in 0 spots
    adata_backup = adata.copy()

    # Initialize statistics dictionary
    stats = {
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

    # Apply filtering steps in order
    adata, stats = filter_outside_tissue(adata, stats)
    adata, stats = filter_min_counts(adata, min_counts, stats)
    adata, stats = filter_min_genes(adata, min_genes, stats)
    adata, stats = filter_genes_min_spots(adata, min_spots, stats)
    adata, stats = filter_mito(adata, mito_threshold, stats)
    adata, stats = filter_ribo(adata, ribo_threshold, stats)
    adata, stats = filter_hb(adata, hb_threshold, stats)

    # Check if filtering resulted in 0 spots or genes
    filtering_failed = False
    if adata.shape[0] == 0 or adata.shape[1] == 0:
        print("WARNING: Filtering resulted in 0 spots or genes remaining!")
        print("Restoring original unfiltered data.")
        adata = adata_backup
        filtering_failed = True

    stats["filtering_failed"] = filtering_failed

    # Calculate final statistics
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

    print("Filtering summary:")
    print(f"  Spots: {n_total_spots} -> {n_remaining_spots} ({n_spots_filtered} removed)")
    print(f"  Genes: {n_total_genes} -> {n_remaining_genes} ({n_genes_filtered} removed)")

    return adata, stats


def write_stats(stats, output_path):
    """Write filtering statistics to JSON file."""
    with open(output_path, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Written filtering statistics to: {output_path}")


def write_versions(process_name):
    """Write software versions to YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "scanpy": importlib.metadata.version("scanpy"),
            "anndata": importlib.metadata.version("anndata"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Filter spots and genes from an AnnData object based on QC metrics."""

    # Template variables
    input_adata = "${adata}"
    output_adata = "${prefix}.h5ad"
    output_stats = "${prefix}_stats.json"
    process_name = "${task.process}"

    # Filtering parameters
    min_counts = int("${min_counts}")
    min_genes = int("${min_genes}")
    min_spots = int("${min_spots}")
    mito_threshold = float("${mito_threshold}")
    ribo_threshold = float("${ribo_threshold}")
    hb_threshold = float("${hb_threshold}")

    # Read AnnData
    print(f"Filtering AnnData from: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Apply filtering
    adata, stats = filter_adata(
        adata,
        min_counts=min_counts,
        min_genes=min_genes,
        min_spots=min_spots,
        mito_threshold=mito_threshold,
        ribo_threshold=ribo_threshold,
        hb_threshold=hb_threshold
    )

    # Add sample ID to stats
    stats["sample_id"] = "${meta.id}"

    # Write outputs
    adata.write_h5ad(output_adata)
    print(f"Written filtered AnnData to: {output_adata}")

    write_stats(stats, output_stats)
    write_versions(process_name)


if __name__ == "__main__":
    main()
