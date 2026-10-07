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

# Keep caches in the task's work directory, which is always writable and
# private to the task
os.environ["NUMBA_CACHE_DIR"] = os.path.join(os.getcwd(), ".cache", "numba")
os.environ["MPLCONFIGDIR"] = os.path.join(os.getcwd(), ".cache", "matplotlib")
os.environ["XDG_CACHE_HOME"] = os.path.join(os.getcwd(), ".cache")

import importlib.metadata
import logging
import pickle
import platform
from pathlib import Path

import anndata as ad
import squidpy as sq
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


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

    if not adata.var_names.is_unique:
        raise ValueError("Gene names in `adata.var_names` are not unique")

    # Compute spatial autocorrelation
    sq.gr.spatial_autocorr(
        adata,
        mode=mode
    )

    return adata


def get_results_key(mode):
    """Get the `uns` key that squidpy stores results under for a mode."""
    if mode == "moran":
        return "moranI"
    elif mode == "geary":
        return "gearyC"
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'moran' or 'geary'.")


def write_svg_to_csv(adata, results_key, output_csv):
    """Export spatially variable genes results to CSV."""
    svg_df = adata.uns[results_key]
    svg_df.to_csv(output_csv)
    logger.info(f"Exported SVG results to: {output_csv}")


def write_pickle(data, slot, name):
    """Write data to a `<slot>/<name>.pkl` pickle file."""
    Path(slot).mkdir(exist_ok=True)
    with open(f"{slot}/{name}.pkl", "wb") as f:
        pickle.dump(data, f, protocol=5)
    logger.info(f"Written slot data to: {slot}/{name}.pkl")


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
    write_adata = "${write_adata}" == "true"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)
    logger.info(f"Computing spatial autocorrelation for: {h5ad}")

    results_key = get_results_key(mode)
    adata = compute_spatial_autocorr(adata, mode)

    write_svg_to_csv(adata, results_key, output_csv)
    write_pickle(adata.uns[results_key], "uns", results_key)

    if write_adata:
        adata.write_h5ad(output_adata)
        logger.info(f"Written AnnData with spatial autocorrelation to: {output_adata}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
