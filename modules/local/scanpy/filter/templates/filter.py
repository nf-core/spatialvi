#!/usr/bin/env python3
"""
Filter observations (cells for scRNA-seq data, spots for spatial data) and genes
from an AnnData object based on QC metrics.

Filtering steps (in order):
 1. Remove observations outside tissue (spatial data only)
 2. Filter observations by minimum counts
 3. Filter observations by minimum genes
 4. Filter genes by minimum observations
 5. Filter observations by maximum mitochondrial content
 6. Filter observations by minimum ribosomal content
 7. Filter observations by maximum haemoglobin content
"""

import importlib.metadata
import json
import platform

import anndata as ad
import scanpy as sc
import yaml


def filter_by_obs_column(adata, col, threshold, filter_below, stat_key, stats):
    """
    Filter observations based on a threshold for a given column.

    Parameters
    ----------
    adata : AnnData
        AnnData object to filter.
    col : str
        Column name in adata.obs to filter on.
    threshold : float
        Threshold value for filtering.
    filter_below : bool
        If True, filter obs with values smaller than the threshold, otherwise
        filter values larger than the threshold.
    stat_key : str
        Key for storing the filtered count in stats dictionary.
    stats : dict
        Statistics dictionary to update.

    Returns
    -------
    tuple
        Filtered AnnData and updated statistics.
    """
    if col not in adata.obs.columns:
        print(f"Column '{col}' not found, skipping filtering")
        stats[stat_key] = 0
        return adata, stats

    n_before = adata.shape[0]
    if filter_below:
        adata = adata[adata.obs[col] >= threshold].copy()
        symbol = "<"
    else:
        adata = adata[adata.obs[col] <= threshold].copy()
        symbol = ">"

    n_filtered = n_before - adata.shape[0]
    print(f"Removed {n_filtered} obs with {symbol} {threshold}% {col}")
    stats[stat_key] = n_filtered
    return adata, stats


def filter_outside_tissue(adata, stats):
    """Filter observations outside tissue based on 'in_tissue' column."""
    in_tissue_col = "in_tissue"
    if in_tissue_col not in adata.obs.columns:
        print(f"Column '{in_tissue_col}' not found, skipping tissue filtering")
        stats["obs_filtered_outside_tissue"] = 0
        return adata, stats

    n_before = adata.shape[0]
    adata = adata[adata.obs[in_tissue_col] == 1].copy()
    n_filtered = n_before - adata.shape[0]
    print(f"Removed {n_filtered} obs outside tissue")
    stats["obs_filtered_outside_tissue"] = n_filtered
    return adata, stats


def filter_min_counts(adata, min_counts, stats):
    """Filter observations with fewer than min_counts total counts."""
    n_before = adata.shape[0]
    sc.pp.filter_cells(adata, min_counts=min_counts)
    n_filtered = n_before - adata.shape[0]
    print(f"Removed {n_filtered} obs with < {min_counts} counts")
    stats["obs_filtered_min_counts"] = n_filtered
    return adata, stats


def filter_min_genes(adata, min_genes, stats):
    """Filter observations with fewer than min_genes expressed genes."""
    n_before = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=min_genes)
    n_filtered = n_before - adata.shape[0]
    print(f"Removed {n_filtered} obs with < {min_genes} genes")
    stats["obs_filtered_min_genes"] = n_filtered
    return adata, stats


def filter_genes_min_obs(adata, min_obs, stats):
    """Filter genes expressed in fewer than min_obs observations."""
    n_before = adata.shape[1]
    sc.pp.filter_genes(adata, min_cells=min_obs)
    n_filtered = n_before - adata.shape[1]
    print(f"Removed {n_filtered} genes expressed in < {min_obs} obs")
    stats["genes_filtered_min_obs"] = n_filtered
    return adata, stats


def filter_adata(
    adata,
    min_counts,
    min_genes,
    min_obs,
    mito_threshold,
    ribo_threshold,
    hb_threshold,
):
    """
    Apply all filtering steps to AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData object to be filtered.
    min_counts : int
        Minimum total counts per observation.
    min_genes : int
        Minimum genes expressed per observation.
    min_obs : int
        Minimum observations expressing each gene.
    mito_threshold : float
        Maximum mitochondrial content percentage.
    ribo_threshold : float
        Minimum ribosomal content percentage.
    hb_threshold : float
        Maximum haemoglobin content percentage.

    Returns
    -------
    tuple
        Filtered AnnData and filtering statistics dictionary.
    """
    print(f"Initial shape: {adata.shape}")

    n_total_obs = adata.shape[0]
    n_total_genes = adata.shape[1]

    stats = {
        "total_obs_before": n_total_obs,
        "total_genes_before": n_total_genes,
        "min_counts_threshold": min_counts,
        "min_genes_threshold": min_genes,
        "min_obs_threshold": min_obs,
        "mito_threshold": mito_threshold,
        "ribo_threshold": ribo_threshold,
        "hb_threshold": hb_threshold,
    }

    adata.var_names_make_unique()

    # Apply filtering steps
    adata, stats = filter_outside_tissue(adata, stats)
    adata, stats = filter_min_counts(adata, min_counts, stats)
    adata, stats = filter_min_genes(adata, min_genes, stats)
    adata, stats = filter_genes_min_obs(adata, min_obs, stats)
    adata, stats = filter_by_obs_column(
        adata,
        "pct_counts_mt",
        mito_threshold,
        filter_below=False,
        stat_key="obs_filtered_mito",
        stats=stats
    )
    adata, stats = filter_by_obs_column(
        adata,
        "pct_counts_ribo",
        ribo_threshold,
        filter_below=True,
        stat_key="obs_filtered_ribo",
        stats=stats
    )
    adata, stats = filter_by_obs_column(
        adata,
        "pct_counts_hb",
        hb_threshold,
        filter_below=False,
        stat_key="obs_filtered_hb",
        stats=stats
    )

    if adata.shape[0] == 0:
        raise ValueError("Filtering resulted in 0 observations remaining.")
    if adata.shape[1] == 0:
        raise ValueError("Filtering resulted in 0 genes remaining.")

    # Final statistics
    stats["total_obs_after"] = adata.shape[0]
    stats["total_genes_after"] = adata.shape[1]
    stats["total_obs_filtered"] = n_total_obs - adata.shape[0]
    stats["total_genes_filtered"] = n_total_genes - adata.shape[1]
    adata.uns["filtering_stats"] = stats

    print("Filtering summary:")
    print(f"  Obs: {n_total_obs} -> {stats['total_obs_after']} "
          f"({stats['total_obs_filtered']} removed)")
    print(f"  Genes: {n_total_genes} -> {stats['total_genes_after']} "
          f"({stats['total_genes_filtered']} removed)")

    return adata, stats


def write_stats(stats, output_path):
    """Write filtering statistics to a JSON file."""
    with open(output_path, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Written filtering statistics to: {output_path}")


def write_versions(process_name):
    """Write software versions to a YAML file."""
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
    """Filter observations and genes from an AnnData object."""
    # Template variables
    h5ad = "${adata}"
    sample_id = "${meta.id}"
    min_counts = int("${min_counts}")
    min_genes = int("${min_genes}")
    min_obs = int("${min_obs}")
    mito_threshold = float("${mito_threshold}")
    ribo_threshold = float("${ribo_threshold}")
    hb_threshold = float("${hb_threshold}")
    output_adata = "${prefix}.h5ad"
    output_stats = "${prefix}_stats.json"
    process_name = "${task.process}"

    print(f"Filtering AnnData from: {h5ad}")
    adata = ad.read_h5ad(h5ad)

    adata, stats = filter_adata(
        adata,
        min_counts=min_counts,
        min_genes=min_genes,
        min_obs=min_obs,
        mito_threshold=mito_threshold,
        ribo_threshold=ribo_threshold,
        hb_threshold=hb_threshold,
    )
    stats["sample_id"] = sample_id

    adata.write_h5ad(output_adata)
    print(f"Written filtered AnnData to: {output_adata}")
    write_stats(stats, output_stats)

    write_versions(process_name)


if __name__ == "__main__":
    main()
