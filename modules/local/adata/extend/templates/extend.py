#!/usr/bin/env python3
"""
Extend an AnnData object with pickled `obs`, `var`, `obsm`, `varm`, and `uns`
content produced by other modules; one directory per slot, with the file name as
the key.

Every slot index must match the input adata object exactly, and existing columns
or keys are not replaced. `allow_missing` re-indexes the slot index to the input
adata object instead, and `overwrite` allows existing columns and keys to be
replaced.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import logging
import pickle
import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(
    level=logging.INFO,
    format="%(name)s - %(levelname)s: %(message)s"
)
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def load_pickle(path):
    """Load a pickle file."""
    if path.suffix == ".pkl":
        with open(path, "rb") as f:
            loaded = pickle.load(f)
        return loaded
    else:
        raise ValueError(f"Unsupported file extension: `{path}`")


def handle_index_mismatches(adata, data, slot, name, allow_missing):
    """
    Check that the index of slot data matches the corresponding adata axis.

    `obs` and `obsm` are aligned to `obs_names`, `var` and `varm` to
    `var_names`. A mismatch (different values or order) raises an error, unless
    `allow_missing` is set, in which case the data is re-indexed to the axis and
    unmatched entries become missing values.
    """
    idx_name = slot[0:3]
    slot_idx = getattr(adata, idx_name).index
    if not data.index.equals(slot_idx):
        if allow_missing:
            data = data.reindex(slot_idx)
            logger.info(f"Re-indexed {name} with adata.{idx_name} index")
        else:
            raise ValueError(f"Index for {slot}.{name} differs from adata")
    return data


def handle_name_collisions(adata, slot, name, names, overwrite):
    """
    Handle overlapping names between data and `adata.{slot}`.

    Overwrites an already existing column/key in `adata.{slot}` if `overwrite`
    is set, otherwise raises an error.
    """
    if name in names:
        if overwrite:
            del getattr(adata, slot)[name]
            logger.warning(f"Overwrote existing `adata.{slot}[{name}]`")
        else:
            raise ValueError(f"`{name}` already exists in `adata.{slot}`")
    return adata


def extend_frame(adata, data, slot, name, allow_missing, overwrite):
    """
    Add the columns of a DataFrame to `adata.obs` or `adata.var`.
    """
    data = handle_index_mismatches(adata, data, slot, name, allow_missing)

    adata_cols = getattr(adata, slot).columns
    for data_col in data.columns:
        adata = handle_name_collisions(
            adata,
            slot,
            data_col,
            adata_cols,
            overwrite
        )

    setattr(adata, slot, pd.concat([getattr(adata, slot), data], axis=1))
    logger.info(f"Extended `adata.{slot}.{name}`")
    return adata


def extend_matrix(adata, data, slot, name, allow_missing, overwrite):
    """
    Add a matrix to `adata.obsm` or `adata.varm` under the key `name`.
    """
    data = handle_index_mismatches(adata, data, slot, name, allow_missing)
    names = getattr(adata, slot).keys()
    adata = handle_name_collisions(adata, slot, name, names, overwrite)
    getattr(adata, slot)[name] = data.to_numpy()
    logger.info(f"Extended `adata.{slot}.{name}`")
    return adata


def extend_uns(adata, data, name, overwrite):
    """
    Add unstructured data to `adata.uns` under the key `name`.
    """
    adata = handle_name_collisions(adata, "uns", name, adata.uns, overwrite)
    adata.uns[name] = data
    logger.info(f"Extended `adata.uns.{name}`")
    return adata


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
            "pandas": pd.__version__,
            "numpy": np.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Extend an adata object with slot data."""

    # Template variables
    h5ad = "${h5ad}"
    allow_missing = "${allow_missing}" == "true"
    overwrite = "${overwrite}" == "true"
    prefix = "${prefix}"
    process_name = "${task.process}"

    adata = ad.read_h5ad(h5ad)

    slot_dirs = {
        "obs": Path("obs/"),
        "var": Path("var/"),
        "obsm": Path("obsm/"),
        "varm": Path("varm/"),
        "uns": Path("uns/"),
    }

    for slot, directory in slot_dirs.items():
        for path in sorted(directory.glob("*")):
            data = load_pickle(path)
            name = path.stem
            if slot in ("obs", "var"):
                adata = extend_frame(
                    adata,
                    data,
                    slot,
                    name,
                    allow_missing,
                    overwrite
                )
            elif slot in ("obsm", "varm"):
                adata = extend_matrix(
                    adata,
                    data,
                    slot,
                    name,
                    allow_missing,
                    overwrite
                )
            elif slot == "uns":
                adata = extend_uns(
                    adata,
                    data,
                    name,
                    overwrite
                )

    adata.write_h5ad(f"{prefix}.h5ad")

    write_versions(process_name)

if __name__ == "__main__":
    main()
