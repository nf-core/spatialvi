#!/usr/bin/env python3
"""
Compute spatial autocorrelation statistics to identify spatially variable
genes.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform
import yaml

import anndata as ad
import squidpy as sq


def get_results_key(mode):
    """Get the uns key for results based on mode."""


def compute_spatial_autocorr(adata, mode):
    """Compute spatial autocorrelation statistics."""
    print(f"Shape: {adata.shape}")
    print(f"Mode: {mode}")

    # Validate input
    if "spatial_connectivities" not in adata.obsp:
        raise ValueError("Spatial connectivities not found. Run squidpy.gr.spatial_neighbors first.")

    # Make var names unique
    adata.var_names_make_unique()

    # Compute spatial autocorrelation
    sq.gr.spatial_autocorr(
        adata,
        mode=mode
    )

    return adata


def write_svg_to_csv(adata, mode, output_csv):
    """Export spatially variable genes results to CSV."""

    # Validate mode
    if mode == "moran":
        results_key = "moranI"
    elif mode == "geary":
        results_key = "gearyC"
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'moran' or 'geary'.")

    # Export genes to CSV
    svg_df = adata.uns[results_key]
    svg_df.to_csv(output_csv)
    print(f"Exported SVG results to: {output_csv}")


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "squidpy": importlib.metadata.version("squidpy"),
            "anndata": importlib.metadata.version("anndata"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Compute spatial autocorrelation for an AnnData object."""

    # Template variables
    input_adata = "${adata}"
    mode = "${mode}"
    output_adata = "${prefix}.h5ad"
    output_csv = "${prefix}_svg.csv"
    process_name = "${task.process}"

    # Read AnnData
    print(f"Computing spatial autocorrelation for: {input_adata}")
    adata = ad.read_h5ad(input_adata)

    # Compute spatial autocorrelation
    adata = compute_spatial_autocorr(adata, mode)

    # Write results to CSV
    write_svg_to_csv(adata, mode, output_csv)

    # Write output H5AD
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with spatial autocorrelation to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
