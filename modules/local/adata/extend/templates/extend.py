#!/usr/bin/env python3
"""
Extend an AnnData object with pickled `obs`, `var`, `obsm`, `varm`, `obsp` and
`uns` content produced by other modules; one directory per slot, with the file
name as the key.

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
            # Re-indexing should fail with duplicated names
            if not (slot_idx.is_unique and data.index.is_unique):
                raise ValueError(
                    f"Can't re-index slot data `{slot}/{name}.pkl`: names in "
                    f"`adata.{idx_name}_names` or the slot data are duplicated"
                )
            # Log number of missing entries
            n_missing = (~slot_idx.isin(data.index)).sum()
            data = data.reindex(slot_idx)
            logger.warning(
                f"Re-indexed slot data `{slot}/{name}.pkl` to "
                f"`adata.{idx_name}_names`; {n_missing} of {len(slot_idx)} "
                "entries have no slot data"
            )
        else:
            raise ValueError(
                f"Index of slot data `{slot}/{name}.pkl` differs from "
                f"`adata.{idx_name}_names`"
            )
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
            logger.warning(f'Overwrote existing `adata.{slot}["{name}"]`')
        else:
            raise ValueError(f'`adata.{slot}["{name}"]` already exists')
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
    logger.info(
        f"Added columns {list(data.columns)} to `adata.{slot}` from slot data "
        f"`{slot}/{name}.pkl`"
    )
    return adata


def extend_matrix(adata, data, slot, name, allow_missing, overwrite):
    """Add a matrix to `adata.obsm` or `adata.varm` under the key `name`."""
    data = handle_index_mismatches(adata, data, slot, name, allow_missing)
    names = getattr(adata, slot).keys()
    adata = handle_name_collisions(adata, slot, name, names, overwrite)
    getattr(adata, slot)[name] = data.to_numpy()
    logger.info(f'Added `adata.{slot}["{name}"]`')
    return adata


def extend_pairwise(adata, data, name, overwrite):
    """
    Add a pairwise observation annotation to `adata.obsp` with `name` key.

    The `data` is required to be a dictionary with the "matrix" and "index"
    keys, storing the sparse matrix and the corresponding `obs_names`.
    Differences between the index and `adata.obs_names` are not allowed.
    """
    if not adata.obs_names.equals(data["index"]):
        raise ValueError(
            f"Index of slot data `obsp/{name}.pkl` differs from "
            "`adata.obs_names`; graphs can't be re-indexed"
        )
    names = adata.obsp.keys()
    adata = handle_name_collisions(adata, "obsp", name, names, overwrite)
    adata.obsp[name] = data["matrix"]
    logger.info(f'Added `adata.obsp["{name}"]`')
    return adata


def extend_uns(adata, data, name, overwrite):
    """Add unstructured data to `adata.uns` under the key `name`."""
    adata = handle_name_collisions(adata, "uns", name, adata.uns, overwrite)
    adata.uns[name] = data
    logger.info(f'Added `adata.uns["{name}"]`')
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

    slot_dirs = {
        "obs": Path("obs/"),
        "var": Path("var/"),
        "obsm": Path("obsm/"),
        "varm": Path("varm/"),
        "obsp": Path("obsp/"),
        "uns": Path("uns/"),
    }
    slot_paths = {
        slot: sorted(directory.glob("*"))
        for slot, directory in slot_dirs.items()
    }

    # Abort if no slot data is given
    if not any(slot_paths.values()):
        raise ValueError("No slot data given; there is nothing to extend")

    adata = ad.read_h5ad(h5ad)
    logger.info(f"Read base AnnData with shape {adata.shape}: {h5ad}")

    for slot, paths in slot_paths.items():
        for path in paths:
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
            elif slot == "obsp":
                adata = extend_pairwise(
                    adata,
                    data,
                    name,
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
    logger.info(f"Written extended AnnData to: {prefix}.h5ad")

    write_versions(process_name)

if __name__ == "__main__":
    main()
