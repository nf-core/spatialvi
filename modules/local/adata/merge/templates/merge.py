#!/usr/bin/env python3
"""
Merge multiple AnnData objects into one.

Gene annotations (`var` columns) are kept only if they are the same in every
object, so per-sample statistics are dropped. `preserve_spatial` keeps the
spatial data in `uns["spatial"]` of every object, and `layer` puts the merged
content of that layer into `X` (e.g. raw counts, to be normalised again),
keeping the layer itself.
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
import platform
from pathlib import Path

import anndata as ad
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def add_spatial(adata, adata_list):
    """
    Adds `.uns['spatial']` back into a merged AnnData objects from the original
    list of multiple AnnData objects.
    """
    merged_spatial = {}
    for adata_orig in adata_list:
        if "spatial" in adata_orig.uns:
            merged_spatial.update(adata_orig.uns["spatial"])
    if merged_spatial:
        adata.uns["spatial"] = merged_spatial
        logger.info("Preserved `.uns['spatial']` data")
    return adata


def validate_var_names(adata_list, keys):
    """Check that gene names are unique within each AnnData object."""
    for adata, key in zip(adata_list, keys):
        if not adata.var_names.is_unique:
            raise ValueError(
                f"Gene names in `adata.var_names` of '{key}' are not unique"
            )


def validate_layer(adata_list, keys, layer):
    """Check that every AnnData object has `layer`."""
    for adata, key in zip(adata_list, keys):
        if layer not in adata.layers:
            raise ValueError(
                f'`adata.layers["{layer}"]` not found in `{key}`'
            )


def merge_adata(adata_list, keys, join, label, preserve_spatial, layer):
    """
    Merge multiple AnnData objects into one, keeping the `var` columns that are
    the same in every object. Can optionally preserve `.uns['spatial']` and put
    a layer into `X` for the final merged object.
    """
    validate_var_names(adata_list, keys)
    if layer:
        validate_layer(adata_list, keys, layer)

    logger.info(f"Merging {len(adata_list)} AnnData objects using {join} join")
    adata = ad.concat(
        adata_list,
        join=join,
        merge="same",
        label=label,
        keys=keys,
        index_unique="-"
    )

    if preserve_spatial:
        adata = add_spatial(adata, adata_list)

    if layer:
        adata.X = adata.layers[layer]
        logger.info(f'Set `adata.X` to the merged `adata.layers["{layer}"]`')

    logger.info(f"Final merged AnnData {adata}")

    return adata


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Merge multiple AnnData objects into one."""

    # Template variables
    h5ads = "${h5ad}".split()
    join = "${join}"
    label = "${label}"
    preserve_spatial = "${preserve_spatial}" == "true"
    layer = "${layer}"
    output_file = "${prefix}.h5ad"
    process_name = "${task.process}"

    adata_list = []
    for h5ad in h5ads:
        adata = ad.read_h5ad(h5ad)
        adata_list.append(adata)
        logger.info(f"Read AnnData object {adata}")

    sample_names = [Path(h5ad).stem for h5ad in h5ads]
    adata_integrated = merge_adata(
        adata_list,
        sample_names,
        join=join,
        label=label,
        preserve_spatial=preserve_spatial,
        layer=layer
    )

    adata_integrated.write_h5ad(output_file)
    logger.info(f"Written merged AnnData to: {output_file}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
