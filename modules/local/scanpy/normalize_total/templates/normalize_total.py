 #!/usr/bin/env python3
"""
Normalize total counts per cell/spot using scanpy.

Normalizes each cell/spot to have the same total count after normalization.
By default, uses median total counts as the target sum.
"""

import importlib.metadata
import platform
import yaml

import anndata as ad
import scanpy as sc


def parse_target_sum(target_sum_str):
    """Parse target_sum parameter from string."""
    if target_sum_str.lower() in ["null", "none", ""]:
        return None
    return float(target_sum_str)


def normalize_adata(adata, target_sum):
    """Normalize total counts per cell/spot."""
    print(f"Shape: {adata.shape}")
    print(f"Target sum: {target_sum if target_sum else 'median of total counts'}")

    # Store pre-normalization statistics
    if "total_counts" in adata.obs:
        median_counts_before = adata.obs["total_counts"].median()
        print(f"Median total counts before normalization: {median_counts_before:.2f}")

    # Normalize total counts
    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)

    # Store post-normalization statistics
    if "total_counts" in adata.obs:
        # Recalculate after normalization
        if hasattr(adata.X, 'A1'):
            adata.obs["total_counts_normalized"] = adata.X.sum(axis=1).A1
        else:
            adata.obs["total_counts_normalized"] = adata.X.sum(axis=1)
        median_counts_after = adata.obs["total_counts_normalized"].median()
        print(f"Median total counts after normalization: {median_counts_after:.2f}")

    # Store normalization info in uns
    adata.uns["normalization"] = {
        "method": "normalize_total",
        "target_sum": target_sum if target_sum else "median",
    }

    return adata


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
    """Normalize total counts in an AnnData object."""

    # Template variables
    input_adata = "${adata}"
    output_adata = "${prefix}.h5ad"
    target_sum = "${target_sum}"
    process_name = "${task.process}"

    # Parse target_sum parameter
    target_sum = parse_target_sum(target_sum)

    # Read AnnData
    print(f"Normalizing AnnData from: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Normalize
    adata = normalize_adata(adata, target_sum)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written normalized AnnData to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
