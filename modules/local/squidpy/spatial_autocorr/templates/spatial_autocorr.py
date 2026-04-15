#!/usr/bin/env python3
"""
Compute spatial autocorrelation statistics to identify spatially variable
genes.

Supports Moran's I and Geary's C statistics for identifying genes with
spatially variable expression patterns. Results are stored in adata.uns
and exported to CSV.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import logging
import platform

import anndata as ad
import squidpy as sq
import yaml

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def compute_spatial_autocorr(adata, mode):
    """Compute spatial autocorrelation statistics."""
    logger.info(f"Shape: {adata.shape}")
    logger.info(f"Mode: {mode}")

    # Validate input
    if "spatial_connectivities" not in adata.obsp:
        raise ValueError(
            "Spatial connectivities not found; "
            "run squidpy.gr.spatial_neighbors first."
        )

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
    if mode == "moran":
        results_key = "moranI"
    elif mode == "geary":
        results_key = "gearyC"
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'moran' or 'geary'.")

    svg_df = adata.uns[results_key]
    svg_df.to_csv(output_csv)
    logger.info(f"Exported SVG results to: {output_csv}")


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
    h5ad = "${adata}"
    mode = "${mode}"
    output_adata = "${prefix}.h5ad"
    output_csv = "${prefix}_svg.csv"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)
    logger.info(f"Computing spatial autocorrelation for: {h5ad}")

    adata = compute_spatial_autocorr(adata, mode)

    write_svg_to_csv(adata, mode, output_csv)

    adata.write_h5ad(output_adata)
    logger.info(f"Written AnnData with spatial autocorrelation to: {output_adata}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
